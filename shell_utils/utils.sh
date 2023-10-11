#util1: get time
# time=`date '+%Y-%m-%d %H:%M:%S'`

#util2: log file format
# log="name_"`date '+%Y-%m%d-%H%M'`".log"

# args specification of PBS job
#args=$(echo work_dir=$work_dir,bam=$bam,ref=$ref, \
#            var_catalog=$var_catalog,out_prefix=$out_prefix, \
#			            threads=$threads | tr -d ' ')

# pbs utils
# log="$log_dir/call_hpgp_ngs_vntr_"`date "+%Y-%m%d-%H%M"`".log"
# args=$(echo work_dir=$work_dir,var_catalog=$var_catalog, \
#         sample_info=$sample_info,ref_genome=$ref_genome,threads=$threads, \
#         out_dir=$out_dir | tr -d ' ' )
# # echo qsub -N ${sample_id} -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script
# echo qsub -N call_vntr -q $job_queue -l nodes=1:ppn=$threads -l mem=$applied_memory -o $log -v $args $script

function source_straglr_env(){
    echo "Setting following path to PATH"
    echo 'export PATH=/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin:$PATH'
    export PATH=/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin:$PATH
}

function python_utils(){
    echo 'import sys'
    echo 'sys.path.append("/public/home/fan_lab/wangjie//utils/py_utils/")'
}

function calcul(){
    expr=$1
	echo "expr is $expr"
    python -c 'print('"$expr"')'
}

#@ text_util
function mgrep(){
    if [ $# -eq 0 ]; then
        echo "Usage: mgrep text1 text2 ... target_file "
        return
    fi
    
    num_args=$#
    
    # return 0
    search_text=""
    num=0
    for text in $@
    do
        num=$((num+1))
        if [[ $num -ne $num_args ]]; then
            search_text=$search_text"|"$text
        else
            target_file=$text
        fi
    done

    search_text_final=$(echo "$search_text" | sed 's/|//' ) # 将第一个|去掉

    # echo "text is $search_text_final, target_file=$target_file"
    grep -E "$search_text_final" "$target_file"
}

#@ fileop
function get_filename(){
    if [ $# -eq 0 ]; then
        echo "Usage: get_filename filepath "
        return
    fi

    filepath=$1
    filename=${filepath##*/}
    echo "$filename"
}

#@ fileop
function get_filedir(){
    if [ $# -eq 0 ]; then
        echo "Usage: get_filedir filepath "
        echo "get dir of a filepath"
        return
    fi

    filepath=$1
    filedir=${filepath%/*}
    echo "$filedir"
}

function find_log_file(){
    sample_id=$1
    log_dir=$2

    ls_tmp_err=/tmp/ls.error
    log=$(ls "$log_dir"/*"${sample_id}"* 2>$ls_tmp_err )
    
    # 判断log文件是否存在
    grep -q "No such file or directory" $ls_tmp_err
    if [[ $? -eq 0 ]]; then
        echo "Log file of $sample_id does not exist!"
        rm -f $ls_tmp_err
        return 1
    fi

    # 判断该样本是否存在多个log文件
    log_arr=($log)    
    if [[ ${#log_arr[*]} -ne 1 ]]; then 
        echo "This sample has many log files. Please check it. $sample_id"
        rm -f $ls_tmp_err
        return 1
    fi

    log=${log_arr[0]}
    echo "$log"
    rm -f $ls_tmp_err

    return 0
}

function log_is_done(){
    log=$1
    is_done=$(check_logs.sh "$log" | head -1 | cut -f2)
    if [[ $is_done == "True" ]]; then
        echo "done"
    else
        echo "not done"
    fi
}

function get_log_time() {
    # usage: $(get_log_time $sample_id $log_dir )
    sample_id=$1
    log_dir=$2

    log=$(find_log_file "$sample_id" "$log_dir")

    if [[ $? -eq 0 ]];then
        # 判断是否跑完
        grep -q "Program is done" "$log"
        if [ $? -eq 0 ]; then
            total_time=$(grep -v '+' "$log" | grep "Total time" | sed 's/Total time: //' )
            echo "$total_time"
            return 0
        else
            echo "file $log isn't completed"
            return 1
        fi
    
    else
        echo "$log"
    fi

}

function icd(){
    # my_cd: 
    if [[ $# -eq 0 ]]; then
        echo "Usage: icd [ dir | file ] "
        echo "icd means integrated cd"
        return 1
    fi

    dir_or_file=$1

    if [[ $dir_or_file == "-" ]]; then
        cd - || exit
    elif [[ -d $dir_or_file ]]; then
        # echo "dir"
        cd $dir_or_file || exit
    elif [[ -f $dir_or_file ]]; then
        filedir=$(dirname "$dir_or_file")
        cd $filedir || exit
    else
        echo "Error! Wrong input! Must be a dir or a file"
    fi

}

function jobnode() {
    # 打印每个任务地详细信息并将其计算节点也打印出来
    qstat -a | grep wangjie > /tmp/job_status
    cat /tmp/job_status | while read line
    do
        jobid=$(echo "$line" | awk '{split($1, res,"\t"); split(res[1], res2, "."); print res2[1]}')

        nodeid=$(pestat | grep "$jobid" | awk '{print $1}')
        echo -e "$line\t$nodeid"
    done
	rm -f /tmp/job_status
}

function tar_dir() {
    if [ $# -eq 0 ]; then
        echo "Usage: tar_dir dir"
        echo 'tar all sub_dirs, use: ls | while read dir; do if [[ -d $dir ]]; then tar_dir $dir & fi; done'
        return
    fi

    dir=$1
    adir=${dir%%/} # 将路径末尾的/去掉
    echo "Dir is $adir"

    if [ -d "$dir" ]; then
        echo "taring $adir to ${adir}.tar"
        tar -cf "${adir}".tar "${adir}"
        echo "removing $adir"
        rm -rf "${adir}"
    else
        echo "$adir is not a directory!"
    fi
}

function chmod_dirs(){
	# main: 可以统一修改当前目录下面所有的目录的权限
	# usage: chmod_dirs . 755
	if [ $# == 0 ]; then
		echo "Usage: chmod_dirs dir mode, such as chmod_dirs . 755"
		return
	fi

	work_dir=$1
	mode=$2
    sub_dirs=$(find $work_dir)
    ssub_dirs=($sub_dirs) # 切分成数组
    for dir in ${ssub_dirs[*]}
    do
        if [ -d "$dir" ]; then
            echo "changing mode of ${dir} to ${mode}"
            chmod "${mode}" "$dir"
        fi
    done

}

function time_template() {
	starttime=$(date +'%Y-%m-%d %H:%M:%S')
    #执行程序
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    
    start_seconds=$(date --date="$starttime" +%s); # 将时间转换成秒数
    end_seconds=$(date --date="$endtime" +%s);

    echo "$start_seconds" "$end_seconds" | awk '{total_sec=$2-$1; min=int(total_sec/60); lef_sec=total_sec%60; printf("Total time: %d minutes, %d seconds.\n", min, lef_sec); }'
}

function cmd_templates() {
	echo "1. get time"
	echo 'time=`date "+%Y-%m-%d %H:%M:%S"`'
	echo "2. log file format"
	echo 'log="name_"`date "+%Y-%m%d-%H%M"`".log"'

	echo "3. qsub template"
	echo 'qsub -N name -o $log -q mini -l nodes=1:ppn=1 -v work_dir=$work_dir,name=$name,genome=$genome,threads=$threads pbs_script.sh'
}

function sam2fastq() {
	awk '{print "@"$1"\n"$10"\n+\n"$11}' "$1"
}

function dir_exists() {
    if [ ! -d "$1" ]; then
        echo "dir $1 not exists, creating it..."
        mkdir -p "$1"
    fi
}

function count_files() {
	# usage: count_files dir
    if [[ $# -eq 0 ]]; then
        echo -e "Usage: count_files dir1 dir2 ... "
        echo "for quickly checking dirs that have most files, use: count_files * |sort -k2,2n"
        return
    fi

    # 循环对每个目录计数
    for dir in "$@"
    do
        file_num=$(find "$dir" | wc -l)
        echo "$dir"" ""$file_num"
    done

}

function add_path(){
    new_path=$1
	echo "adding path: ""$new_path"" to ~/utils/shell_utils/envs.sh"
	echo 'added cmd is: export PATH='"$new_path"':$PATH'
	echo 'export PATH='"$new_path"':$PATH' >> ~/utils/shell_utils/envs.sh
}

function myprint() {
  time=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[ INFO $time ] [ $1 ]"
}

function loginfo() {
  time=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[ INFO $time ] [ $1 ]"
}

function logerror() {
  time=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[ ERROR $time ] [ $1 ]"
}

function file_exists(){
	if [ -f "$1" ]; then
		myprint "Error! File ""$1"" already exists!"
		exit
	fi
}
# logfile format
#log="logs/anno_${name}_"`date '+%Y%m%d-%H-%M'`".log"

### 集合操作

function setdiff() {
    # 差集
    if [[ $# -eq 0 ]]; then
        echo "Usage: setdiff a.txt b.txt [col1 col2]"
        return
    fi

    local curr_time=$(date "+%Y-%m%d-%H%M%S-%N")

    if [[ $# -eq 2 ]]; then
        file1=$1
        file2=$2
        col1=1
        col2=1

    elif [[ $# -eq 4 ]]; then # 指定根据文件的哪些列进行计算
        file1=$1
        file2=$2
        col1=$3
        col2=$4

    else
        echo "Wrong inputs"
        echo "Usage: setdiff a.txt b.txt [col1 col2]"
        return 1
    fi

    tmp_diff_res=/tmp/tmp_diff_res_${curr_time}.txt
    tmp_file1_col=/tmp/tmp_file1.col.tmptmptmp_${curr_time}
    tmp_file2_col=/tmp/tmp_file2.col.tmptmptmp_${curr_time}
    tmp_file1_sorted_col=/tmp/tmp_file1_sorted.tmptmptmp_${curr_time}
    tmp_file2_sorted_col=/tmp/tmp_file2_sorted.tmptmptmp_${curr_time}

    # fetch cols
    cut -f "$col1" "$file1" > "${tmp_file1_col}"
    cut -f "$col1" "$file2" > "${tmp_file2_col}"

    # sort
    sort -k1,1 "${tmp_file1_col}" > "${tmp_file1_sorted_col}"
    sort -k1,1 "${tmp_file2_col}" > "${tmp_file2_sorted_col}"

    comm -23 "${tmp_file1_sorted_col}" "${tmp_file2_sorted_col}" > "$tmp_diff_res"

    # 按照原始文件的顺序打印结果
    # grep -f $tmp_diff_res ${tmp_file1_col} | cat
    cat "$tmp_diff_res"

    rm -f "$tmp_diff_res"
    rm -f "${tmp_file1_col}" "${tmp_file2_col}" 
    rm -f "${tmp_file1_sorted_col}" "${tmp_file2_sorted_col}"

    return 0
}

function intersect() {
    # 交集
    if [ $# -eq 0 ]; then
        echo "Usage: intersect a.txt b.txt"
        return
    fi

    local curr_time=$(date "+%Y-%m%d-%H%M")

    file1=$1
    file2=$2

    sort -k1,1 "$file1" > /tmp/file1_sorted.tmptmptmp_"${curr_time}"
    sort -k1,1 "$file2" > /tmp/file2_sorted.tmptmptmp_"${curr_time}"

    comm -12 /tmp/file1_sorted.tmptmptmp_"${curr_time}" /tmp/file2_sorted.tmptmptmp_"${curr_time}"

    rm /tmp/file1_sorted.tmptmptmp_"${curr_time}" /tmp/file2_sorted.tmptmptmp_"${curr_time}"
}

function union() {
    # 并集
    if [ $# -eq 0 ]; then
        echo "Usage: union a.txt b.txt"
        return
    fi

    file1=$1
    file2=$2

    sort "$file1" "$file2" | uniq
}

function qmefull(){
    qstat -a |grep "${USER}" |awk '{print $1}' | while read jobid
    do
        tsk_info=$(qstat -a "$jobid" | tail -n1)
        jobname=$(qstat -f "$jobid" | grep Job_Name | sed 's/.*= //')
        echo -e "$tsk_info\t$jobname"
    done
}


