awk 'NR==FNR{arr[$1]=$0; next } NR>FNR{if($1 in arr){print $0}}' qry_file.txt tgt2_file.txt

# 这一句命令可以直接实现left_join，只是多输出了一列，再存储的时候把那一列去掉就行了
awk -v OFS="\t" 'NR==FNR{arr[$1]=$0; next } NR>FNR{if($1 in arr){print $0,arr[$1]}}' tgt2_file.txt qry_file.txt


