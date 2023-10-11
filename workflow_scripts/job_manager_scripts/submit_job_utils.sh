function if_submit(){
    # 询问是否提交任务
    read -p "Do you want to run program for those samples? (y/n): " proceed
    if [[ $proceed != 'y' ]]; then
        echo "Bye."
        exit 1
    fi
}

function read_arguments(){
    # batch_size queue threads dry_run
    read -p "Please input batch_size(n or start-end), jobqueue(dft, fat, fsh, mini), threads, (dry_run:y|n): " batch_size queue threads dry_run
    if [[ $queue == "dft" ]]; then
        job_queue="default"
    elif [[ $queue == "fat" ]]; then
        job_queue="fat"
    elif [[ $queue == "fsh" ]]; then
        job_queue="fsh_team"
    elif [[ $queue == "mini" ]]; then
        job_queue="mini"
    else
        echo "Wrong jobqueue! Must be one of dft, fat, fsh!"
        exit 1
    fi

    if [[ $batch_size =~ "-" ]]; then # 判断用户是指定了范围还是指定了提交任务数目
        bstart=$(echo $batch_size | cut -d- -f1) # bstart: batch_start
        bend=$(echo $batch_size | cut -d- -f2)
    else
        bstart=1
        bend=$batch_size
    fi

    if [[ $bend -gt $display_batch_size ]]; then
        echo "batch_size is out of range!"
        exit 1
    fi

    # 是否dry_run? 若dry_run则只输出运行的命令，但并不提交任务
    if [[ $dry_run ]]; then
        if [[ $dry_run == "y" ]]; then
            dry_run="y"
            # echo "dry_run mode is on. Will print cmds of jobs, but not submit them"
        elif [[ $dry_run == "n" ]]; then
            dry_run="n"
        else
            echo "wrong dry_run argument! Must be y or n"
            exit 1
        fi
    else
        dry_run="n"
    fi

    echo -e "$job_queue\t$threads\t$bstart\t$bend\t$dry_run"

}

