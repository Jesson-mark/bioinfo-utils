#task3: run annotate_variation.pl for all chrs avinput file
work_dir=$HOME"/han_data"
data_dir="avinput"
out_dir="geneanno"

for name in `cat chrs_name.txt`
do
	avinput=$data_dir/${name}.avinput
	job_name=annotate_$name
	log_file="logs/${job_name}_"`date '+%Y-%m%d-%H-%M'`".log"
	echo "qsub -N $job_name -o $log_file -v work_dir=$work_dir,out_dir=$out_dir,avinput=$avinput,name=$name annotate_variation.sh"
done

