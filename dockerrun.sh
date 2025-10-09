# docker run -it --security-opt seccomp=unconfined -v `pwd`:/host wsmoses/598ape /bin/bash
docker run -it --security-opt seccomp=unconfined -v "$(pwd)":/host gis3/598ape /bin/bash
