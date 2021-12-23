#!/bin/bash
project_path=$(realpath $(pwd)/../)

vuid=$(id -u)
vgid=$(id -g)

docker run -it --rm --net host --hostname tpolylower -e LOCAL_USER_ID=$vuid -e LOCAL_USER_GID=$gvid -v $(pwd):/home/pl \
       tpolylower /bin/bash
