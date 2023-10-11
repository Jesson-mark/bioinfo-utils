#!/usr/bin/bash

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

# arguments
if [[ $# -eq 0 ]]; then
    echo "Usage: "
    exit 0
fi

# parse arguments
until [ $# -eq 0 ]
do
    name=${1:1}; shift;
    if [[ -z "$1" || $1 == -* ]] ; then 
        eval "$name=true"; 
    else 
        eval "$name=$1"; 
        shift; 
    fi  
done
