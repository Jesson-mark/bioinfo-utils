#!/usr/bin/bash 

set -ue
source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"

if [ $# -eq 0 ]; then
    echo "Usage: print-var-info.sh input_file"
    echo "Example: "
    echo '    echo "genehancer_ids=tmp/genehancer_ids.txt" > tmp/a'
    echo '    print-var-info.sh tmp/a'
    exit
fi

########################## Parse arguments ##########################
input_file=$1


########################## Begin ##########################

sed 's/=.*//' "$input_file" | awk 'NF>0{print "echo \""$1" is $"$1"\""}'

sed 's/=.*//' "$input_file" | awk 'NF>0{printf "\"$"$1"\" "} END{print ""}'

########################## Done ##########################
