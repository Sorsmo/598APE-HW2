#define STB_IMAGE_IMPLEMENTATION
#include "../external/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../external/stb_image_write.h"

#include "../src/he.h"
#include "../src/poly_utils.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  uint8_t *data;
  int width;
  int height;
  int channels;
} Image;

static Image load_image(const char *path) {
  Image img;
  img.data = stbi_load(path, &img.width, &img.height, &img.channels, 0);
  if (!img.data) {
    fprintf(stderr, "Failed to load image: %s\n", path);
    exit(1);
  }
  return img;
}

static void free_image(Image img) { stbi_image_free(img.data); }

static void save_image(const char *path, Image img) {
  stbi_write_png(path, img.width, img.height, img.channels, img.data,
                 img.width * img.channels);
}

static const int sobel_gx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
static const int sobel_gy[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

static void rgb_to_grayscale(uint8_t *input, uint8_t *output, int width,
                             int height, int channels) {
  for (int i = 0; i < width * height; i++) {
    if (channels >= 3) {
      uint8_t r = input[i * channels + 0];
      uint8_t g = input[i * channels + 1];
      uint8_t b = input[i * channels + 2];
      output[i] = (uint8_t)((r + g + b) / 3);
    } else {
      output[i] = input[i * channels];
    }
  }
}

static void sobel_plain(uint8_t *input, uint8_t *output, int width,
                        int height) {
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      int gx = 0, gy = 0;

      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          int pixel = input[(y + ky) * width + (x + kx)];
          gx += pixel * sobel_gx[ky + 1][kx + 1];
          gy += pixel * sobel_gy[ky + 1][kx + 1];
        }
      }

      int magnitude = sqrt(gx * gx + gy * gy);
      if (magnitude > 255)
        magnitude = 255;
      output[y * width + x] = (uint8_t)magnitude;
    }
  }

  for (int x = 0; x < width; x++) {
    output[x] = 0;
    output[(height - 1) * width + x] = 0;
  }
  for (int y = 0; y < height; y++) {
    output[y * width] = 0;
    output[y * width + (width - 1)] = 0;
  }
}

Ciphertext encode_zero(int64_t q, Poly poly_mod) {
  Ciphertext ct;
  ct.c0 = encode_plain_integer(q, 0);
  ct.c1 = encode_plain_integer(q, 0);
  return ct;
}

static void sobel_fhe(Ciphertext *input_enc, Ciphertext *output_enc, int width,
                      int height, int64_t q, int64_t t, Poly poly_mod) {
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      Ciphertext gx = encode_zero(q, poly_mod);
      Ciphertext gy = encode_zero(q, poly_mod);
  
      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          Ciphertext pixel = input_enc[(y + ky) * width + (x + kx)];

          int coeff_gx = sobel_gx[ky + 1][kx + 1];
          int coeff_gy = sobel_gy[ky + 1][kx + 1];

          if (coeff_gx != 0) {
            Ciphertext term = mul_plain(pixel, q, t, poly_mod, coeff_gx);
            gx = add_cipher(gx, term, q, poly_mod);
          }

          if (coeff_gy != 0) {
            Ciphertext term = mul_plain(pixel, q, t, poly_mod, coeff_gy);
            gy = add_cipher(gy, term, q, poly_mod);
          }
        }
      }

      output_enc[y * width + x] = add_cipher(gx, gy, q, poly_mod);
    }
  }

  Ciphertext zero = encode_zero(q, poly_mod);
  for (int x = 0; x < width; x++) {
    output_enc[x] = zero;
    output_enc[(height - 1) * width + x] = zero;
  }
  for (int y = 0; y < height; y++) {
    output_enc[y * width] = zero;
    output_enc[y * width + (width - 1)] = zero;
  }
}

int main(int argc, char **argv) {
  srand(42);
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input_image>\n", argv[0]);
    return 1;
  }

  const char *input_path = argv[1];

  size_t n = 1u << 4;
  int64_t q = 1ll << 30;
  int64_t t = 1ll << 10;

  printf("Loading image: %s\n", input_path);
  Image img = load_image(input_path);
  printf("Image size: %dx%d, channels: %d\n", img.width, img.height, img.channels);

  uint8_t *gray = (uint8_t *)malloc(img.width * img.height * sizeof(uint8_t));
  rgb_to_grayscale(img.data, gray, img.width, img.height, img.channels);

  int total_pixels = img.width * img.height;

  Poly poly_mod = create_poly();
  set_coeff(&poly_mod, 0, 1);
  set_coeff(&poly_mod, n, 1);

  printf("Generating keys...\n");
  KeyPair keys = keygen(n, q, poly_mod);
  PublicKey pk = keys.pk;
  SecretKey sk = keys.sk;

  const int chunk_size = 32; 
  uint8_t *fhe_sobel = (uint8_t *)malloc(total_pixels * sizeof(uint8_t));
  clock_t fhe_start = clock();
  for (int c = 0; c < img.height; c += chunk_size) {
    int start_y;
    if (c == 0) {
      start_y = 0;
    } else {
      start_y = c - 1;
    }

    int end_y; 
    if (c + chunk_size >= img.height) {
      end_y = img.height;
    } else {
      end_y = c + chunk_size + 1;
    }

    int rows = end_y - start_y;

    Ciphertext *gray_chunk_enc = (Ciphertext *)malloc(rows * img.width * sizeof(Ciphertext));
    Ciphertext *sobel_chunk_enc = (Ciphertext *)malloc(rows * img.width * sizeof(Ciphertext));

    for (int i = 0; i < rows; i++) {
      int global_y = start_y + i;
      for (int j = 0; j < img.width; j++) {
        int idx = global_y * img.width + j;
        gray_chunk_enc[i * img.width + j] = encrypt(pk, n, q, poly_mod, t, gray[idx]);
      }
    }

    sobel_fhe(gray_chunk_enc, sobel_chunk_enc, img.width, rows, q, t, poly_mod);

    int copy_start;
    if (c == 0) {
      copy_start = rows;
    } else {
      copy_start = rows - 1;
    }

    int copy_end;
    if (c + chunk_size >= img.height) {
      copy_end = rows;
    } else {
      copy_end = rows - 1;
    }

    for (int i = copy_start; i < copy_end; i++) {
      int global_y = start_y + i;
      for (int j = 0; j < img.width; j++) {
        int idx = global_y * img.width + j;
        int64_t val = decrypt(sk, n, q, poly_mod, t, sobel_chunk_enc[i * img.width + j]);
        if (val > t / 2)
          val = t - val;
        if (val > 255)
          val = 255;
        fhe_sobel[idx] = (uint8_t)val;
      }
    }

    free(gray_chunk_enc);
    free(sobel_chunk_enc);
  }

  clock_t fhe_end = clock();
  double fhe_time = ((double)(fhe_end - fhe_start)) / CLOCKS_PER_SEC;

  // Compute plaintext Sobel
  uint8_t *plain_sobel = (uint8_t *)malloc(total_pixels * sizeof(uint8_t));
  clock_t plain_start = clock();
  sobel_plain(gray, plain_sobel, img.width, img.height);
  clock_t plain_end = clock();
  double plain_time = ((double)(plain_end - plain_start)) / CLOCKS_PER_SEC;

  printf("\n=== Results ===\n");
  printf("FHE Sobel time: %.4f s\n", fhe_time);
  printf("Plaintext Sobel time: %.4f s\n", plain_time);

  save_image("output/sobel_fhe.png", (Image){fhe_sobel, img.width, img.height, 1});
  save_image("output/sobel_plain.png", (Image){plain_sobel, img.width, img.height, 1});
  save_image("output/grayscale.png", (Image){gray, img.width, img.height, 1});

  printf("\nSaved outputs:\n");
  printf("  output/grayscale.png     (grayscale input)\n");
  printf("  output/sobel_fhe.png     (FHE Sobel edges)\n");
  printf("  output/sobel_plain.png   (plaintext Sobel edges)\n");

  free(fhe_sobel);
  free(plain_sobel);
  free(gray);
  free_image(img);

  return 0;
}