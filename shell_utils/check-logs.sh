#!/usr/bin/bash 

if [ $# -eq 0 ]; then
    echo "Usage: check_logs.sh logfile1 logfil2 ... "
    echo "Example1: get_log_time.sh logs/A.log logs/B.log"
    echo 'Example2: get_log_time.sh logs/AAA*'
    exit 0
fi

#usage: check_logs.sh 'logs/convert_chr*'
exp=$1
files=$(ls $exp)
#echo $files
#echo ${#files[*]}

function check_log(){
	#files=$(echo $1 | tr " " "\n")
	#echo ${files[0]}
	#nums=${#files[*]}
	#echo $nums
	file=$1
	a=`grep "Program is done" $file`
	if [ $? -eq 0 ]; then
		echo True
	else
		echo False
	fi
}

num_true=0
num_false=0

for logfile in "$@"; 
do
	status=`check_log $logfile`

	if [ $status == "True" ]; then
		num_true=$((num_true+1))
	elif [ $status == "False" ]; then
		num_false=$((num_false+1))
	else
		echo "Wrong status"
		exit
	fi
	echo -e $logfile"\t"$status
done
echo "Total successes: $num_true, total fails: $num_false"


