#!/usr/bin/bash

# function left_join() {}

if [[ $# -eq 0 ]]; then
    echo "Usage: left_join.sh query_file target_file [qry_idx, tgt_idx]"
    echo "  If qry_idx and tgt_idx not provided, this script will use the first column of both file to left join. "
    exit
fi

if [[ $# -eq 2 ]]; then
    qry_file=$1
    tgt_file=$2
    qry_col_idx=1
    tgt_col_idx=1

elif [[ $# -eq 4 ]]; then
    qry_file=$1
    tgt_file=$2
    qry_col_idx=$3
    tgt_col_idx=$4

else
    echo "Wrong number of arguments! Please see the usage carefully! "
    exit
fi

# temporary files
curr_time=`date "+%Y-%m%d-%H%M"`

query_col=/tmp/query_col_${curr_time}
left_cols=/tmp/left_cols_${curr_time}
target_col=/tmp/target_col_${curr_time}
tgt_res=/tmp/tmp_res_tgt_${curr_time}
qry_res=/tmp/tmp_res_qry_${curr_time}

# process query file
# cut -f 1 $qry_file > query_col
cut -f $qry_col_idx $qry_file > $query_col
cat $qry_file | awk -F '\t' '{
    for(i=1;i<=NF;i++){
        if(i != "'$qry_col_idx'") {
            if(i == NF){
                printf("%s",$i); 
            } else {
                printf("%s\t",$i);
            }
        }
    }
    print "" ;
    }' > $left_cols
    # awk -F '\t' '{for(i=1;i<=NF;i++){if(i != "'$qry_col_idx'") {printf("%s\t",$i);} } print "" }' $qry_file > $left_cols

# process target file
cut -f $tgt_col_idx $tgt_file > $target_col

# 开始匹配

# 搜索target中匹配的行
lines_in_tgt=$(grep -w -n -f $query_col $target_col | awk '{split($1,res,":"); printf("%sp;",res[1])}END{print ""}')
sed -n "$lines_in_tgt" $tgt_file > $tgt_res

# 搜索匹配target中的query中的行
lines_in_qry=$(cut -f $tgt_col_idx $tgt_res | grep -w -n -f - $qry_file | awk '{split($1,res,":"); printf("%sp;",res[1])}END{print ""}')
sed -n "$lines_in_qry" $left_cols > $qry_res

# print result
paste $tgt_res $qry_res

# remove temporary files
# remove=true
remove=false
if [[ $remove == "true" ]]; then
    rm $query_col
    rm $left_cols
    rm $target_col
    rm $tgt_res
    rm $qry_res
fi

