#!/bin/bash

set -x

USER_ID=${LOCAL_USER_ID:-9001}
USER_GID=${LOCAL_USER_GID:-9001}

USER_NAME=`echo $DOCKER_USERNAME_PASSIN`
GROUP_NAME=`echo $DOCKER_GROUPNAME_PASSIN`
HOMEDIR=`echo $DOCKER_HOMEDIR_PASSIN`

groupadd -g ${USER_GID} ${GROUP_NAME}

if [[ -d ${HOMEDIR} ]]
then
    # supress home already exist warning
    useradd --shell /bin/bash -u $USER_ID -o -c "" -m ${USER_NAME}  2>/dev/null
else
    useradd --shell /bin/bash -u $USER_ID -o -c "" -m ${USER_NAME}
fi

usermod -g ${GROUP_NAME} ${USER_NAME}

cd ${HOMEDIR}
exec /usr/sbin/gosu ${USER_NAME} "$@"
