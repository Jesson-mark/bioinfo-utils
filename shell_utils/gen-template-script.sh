#!/usr/bin/bash

if [[ $# -ne 2 ]]; then
    echo "Usage: gen-script-template.sh script_type out_script"
    echo "Examples: "
    echo "    gen-script-template.sh R test.R"
    echo "    gen-script-template.sh r test.R"
    echo "    gen-script-template.sh py test.py"
    echo "    gen-script-template.sh shell test.sh"
    echo "    gen-script-template.sh pbs pbs.sh"
    exit 
fi

file_exists(){
    if [[ -f $1 ]]; then
        echo "Error! File $1 already exist! Please specify a new filename"
        exit 1
    fi
}

newline(){
    echo ""
}

generate_r_template(){
    echo "library(tidyverse)"
    echo "library(data.table)"
    echo "library(argparse)"
    echo ""
    echo "source('/public/home/fan_lab/wangjie/utils/r_utils/set_gg_themes.R')"; echo ""
    work_dir=$PWD
    printf 'work_dir <- "%s"\n' "$work_dir"
    echo "setwd(work_dir)"
    echo "getwd()"
    echo ""

    echo '### get arguments'
    echo 'parser <- ArgumentParser()'
    echo 'parser$add_argument("--input_file", type = "character", required = TRUE)'
    echo 'parser$add_argument("--out_file", type = "character", required = TRUE)'
    echo 'args <- parser$parse_args()'; echo ""

    echo 'input_file <- args$input_file'
    echo 'out_file <- args$out_file'
    echo 'cat("input_file is", input_file, "\n")'
    echo 'cat("out_file is", out_file, "\n")'; echo ""

    echo '# tmp'
    echo 'tmp <- function(){'; newline; newline; echo '}'

}

generate_sh_template(){
    echo "#!/usr/bin/bash "; echo ""
    echo "set -ue"
    echo 'source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"'
    newline
    echo 'if [ $# -eq 0 ]; then'
    echo '    echo "Usage: "'
    echo '    exit'
    echo 'fi'; newline

    echo "########################## Parse arguments ##########################"
    echo 'input_file=$1'
    echo 'out_prefix=$2'
    echo 'threads=$3'
    newline; newline

    echo "########################## Begin ##########################"
    echo 'starttime=$(date +"%Y-%m-%d %H:%M:%S")'
    echo 'myprint "Begin to analysis"'; echo -e "\n\n"
    echo "########################## Done ##########################"
    echo "# 计时结束"
    echo 'endtime=$(date +"%Y-%m-%d %H:%M:%S")'
    echo 'cal_used_time.sh "$starttime" "$endtime"'
    echo 'myprint "Program is done"'

}

generate_py_template(){
    echo "#!/usr/bin/env python"; newline
    echo -e "import os\nimport sys\nimport gzip\nimport argparse\nimport subprocess"; newline

    echo 'class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):'
    echo '    """'
    echo '        可以将参数的默认选项打印在帮助文档中'
    echo '    """'
    echo '    pass'; newline

    echo 'def run_shell_cmd(cmd):'
    echo '    with subprocess.Popen('
    echo '        cmd, '
    echo '        shell=True, '
    echo '        bufsize=1, '
    echo '        executable = "/usr/bin/bash", '
    echo '        stdout=subprocess.PIPE,'
    echo '        stderr=subprocess.STDOUT, '
    echo '        text=True'
    echo '    ) as proc:'
    echo '        for line in proc.stdout:'
    echo '            print(line, end="")'
    newline
    echo '    if proc.returncode != 0:'
    echo '        sys.exit(proc.returncode)'; newline

    echo 'def check_files(*args):'
    echo '    for afile in args:'
    echo '        if not os.path.exists(afile):'
    echo '            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))'
    newline

    echo 'def open_file(afile):'
    echo '    if afile.endswith("gz"):'
    echo '        file_reader = gzip.open(afile, "rt")'
    echo '    else:'
    echo '        file_reader = open(afile)'
    echo '    return file_reader'; newline

    echo 'def myfunc(args):'
    echo '    input_file = args.input_file'
    echo '    out_file = args.out_file'; newline
    echo '    pass'; newline

    echo 'def get_args():'
    echo '    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)'
    echo '    parser.add_argument("--input_file", help="Project dir", required=True)'; newline
    echo '    parser.add_argument("--out_file", help="Project dir", required=True)'; newline
    echo '    if len(sys.argv)==1:'
    echo '        parser.print_help(sys.stderr)'
    echo '        sys.exit(1)'; newline
    echo '    args = parser.parse_args()'; newline
    echo '    myfunc(args)'; newline
    echo '    return args'; newline

    echo 'def main():'
    echo '    args = get_args()'; newline
    echo 'if __name__ == "__main__":'
    echo '    main()'

}

generate_pbs_template(){
    echo -e "#PBS -j oe\n"
    echo -e "#PBS -N test\n"
    echo -e "#PBS -l nodes=1:ppn=1\n"
    echo "########################## Presets ##########################"
    echo "set -ue"
    echo 'source "/public/home/fan_lab/wangjie/utils/shell_utils/utils.sh"'; echo "";
    echo "########################## Arguments ##########################"
    echo "# work_dir "; 
    echo "threads=1"; echo ""; 
    echo "########################## softwares ##########################"
    echo 'export PATH=/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin:$PATH '; echo "";
    echo '########################## Change workdir ##########################'
    echo 'myprint "Working dir is: ""${work_dir}"'
    echo 'cd "${work_dir}"'; echo ""; 
    echo "########################## Begin ##########################"
    echo 'starttime=$(date +"%Y-%m-%d %H:%M:%S")'
    echo 'myprint "Begin to analysis"'; echo -e "\n\n"
    echo "########################## Done ##########################"
    echo "# 计时结束"
    echo 'endtime=$(date +"%Y-%m-%d %H:%M:%S")'
    echo 'cal_used_time.sh "$starttime" "$endtime"'
    echo 'myprint "Threads is $threads"'
    echo 'myprint "Program is done"'
}

script_type=$1
out_script=$2
echo "script_type is $script_type"
echo "out_script is $out_script"
file_exists "${out_script}"

case $script_type in 
    "R" | "r" | "rscript")
        echo "Generating a R template script"
        generate_r_template > "${out_script}"
        echo "Done"
        ;;
    
    "py" | "Py" | "python" | "Python")
        echo "Generating a Python template script"
        generate_py_template > "${out_script}"
        chmod +x "${out_script}"

        ;;

    "sh" | "shell" | "Shell")
        echo "Generating a Shell template script"
        generate_sh_template > "${out_script}"
        chmod +x "${out_script}"
        ;;

    "pbs" | "PBS")
        echo "Generating a PBS template script"
        generate_pbs_template > "${out_script}"
        echo "Done"
        ;;

    *)
        echo "Error! $script_type not supported"
        exit 1
        ;;

esac


