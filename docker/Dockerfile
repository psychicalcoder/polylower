FROM ubuntu:20.04

ENV TZ=Asia/Taipei
RUN ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo $(readlink /etc/localtime)

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y python3 python3-pip
RUN apt-get install -y cmake
RUN pip install numpy pytest
RUN apt-get install -y gosu

ARG CONTAINER_USERNAME=pl
ARG CONTAINER_GROUPNAME=plg
ARG CONTAINER_HOMEDIR=/home/pl

ENV DOCKER_USERNAME_PASSIN ${CONTAINER_USERNAME}
ENV DOCKER_GROUPNAME_PASSIN ${CONTAINER_GROUPNAME}
ENV DOCKER_HOMEDIR_PASSIN ${CONTAINER_HOMEDIR}

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

CMD ["/bin/bash"]
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
