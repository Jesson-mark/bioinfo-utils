#!/usr/bin/env python

#%%
import os
import sys
import time
import json
import argparse
import logging
import subprocess
import pathlib

# logging.basicConfig(level = logging.INFO, 
#                     format = '%(levelname)s %(asctime)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s', 
#                     datefmt = '%Y-%m-%d %H:%M:%S'
#                     )

logger = ''

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    '''
        可以将参数的默认选项打印在帮助文档中
    '''
    pass

def timmer(func):
    '''
        Run function and report used time
    '''
    def wrapper():
        start_time = time.time()

        func()

        end_time = time.time()
        used_time = end_time - start_time
        used_mins = used_time / 60
        print('Total time: %.3s minutes'%(used_mins))
    return wrapper

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception('Error! File `%s` not exists! Please check it!'%(afile))

class Project:
    def __init__(self, project_dir):
        self.project_dir = project_dir
        self.task1_dir = os.path.join(project_dir, 'task1')
        self.task2_dir = os.path.join(project_dir, 'task2')
        self.logs_dir = os.path.join(project_dir, 'logs')
        self.task_status_dir = os.path.join(project_dir, 'task_status')
        self._create_dirs()

    def _create_dirs(self):
        for dirname, dir in self.__dict__.items():
            if not os.path.exists(dir):
                print('Making dir %s for %s'%(dir, dirname))
                os.makedirs(dir)

    def __repr__(self):
        return f'Project({self.project_dir})'

#%%
def task1(project, args):
    '''
        Return commands associated with task1
        Multiple commands are joined by \n.
    '''
    cmd = 'echo task1\n'
    cmd += 'pwd\n'
    cmd += 'cd && ls \n'
    return cmd

def task2(project, args):
    cmd = 'echo task2'
    return cmd

def task3(project, args):
    cmd = 'echo task3'
    return cmd

def task4(project, args):
    # print('hello')
    logger.info('hello')
    cmd = 'echo task4'
    return cmd

def get_task_done_file(project, task_name):
    done_file = os.path.join(project.task_status_dir, f'{task_name}.done')
    return done_file

def check_tasks(project, task_checklist):
    logger.info('************ Checking tasks ************\n')
    num_done_tasks = 0
    done_tasks = []
    not_done_tasks = []
    for task_name in task_checklist:
        task_done_file = get_task_done_file(project, task_name)
        if os.path.exists(task_done_file):
            # logger.info('Task %s done'%(task_name))
            num_done_tasks += 1
            done_tasks.append(task_name)
        else:
            # logger.info('Task %s not done'%(task_name))
            not_done_tasks.append(task_name)

    task_status = 'Not done'
    if num_done_tasks == len(task_checklist):
        logger.info('************ All tasks are done. Bye! ************\n')
        task_status = 'Done'

    else:
        logger.info('************ %s / %s tasks are done. ************\n'%(num_done_tasks, len(task_checklist)))
        logger.info('Done tasks: %s'%(', '.join(done_tasks)))
        logger.info('Not done tasks: %s'%(', '.join(not_done_tasks)))

    return task_status, done_tasks

def check_tasks_light(project, task_checklist):
    done_tasks = []
    not_done_tasks = []
    for task_name in task_checklist:
        task_done_file = get_task_done_file(project, task_name)
        if os.path.exists(task_done_file):
            done_tasks.append(task_name)
        else:
            not_done_tasks.append(task_name)
    print('Project_dir\tDone_tasks\tNot_done_tasks')
    print('%s\t%s\t%s'%(project.project_dir, ','.join(done_tasks), ','.join(not_done_tasks)))

def run_task(task_name, cmd, done_file, force=False, dry_run=False):
    if not force and os.path.exists(done_file):
        logger.info(' ************ Task %s has done ************ \n'%(task_name))
    else:
        # run task
        logger.info(' ************ Running task %s with command: \n%s\n'%(task_name, cmd))

        # fp = open(logfile, 'a+') if logfile is not None else None
        if not dry_run:
            with subprocess.Popen(
                cmd, 
                shell=True, 
                bufsize=1, 
                executable = '/usr/bin/bash', 
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, 
                text=True
            ) as proc:
                for line in proc.stdout:
                    # print(line, end='')
                    logger.info(line.rstrip())
            # if fp is not None:
            #     fp.close()

            if proc.returncode != 0:
                sys.exit(proc.returncode)
            else:
                # if task is done, touch done file
                pathlib.Path(done_file).touch()

def get_task_checklist():
    '''
        Users need to specify their tasks in this function and implement
        the corresponding function named by task name. These functions will
        generate (shell) commands for each task, when supplied required 
        arguments (such as project, args). 
        If not all tasks have their functions, an error will be raised. 
    '''
    task_checklist = ['task1', 'task2', 'task3', 'task4']

    ### check if specified tasks have the corresponding function to
    # generate cmds for them
    
    return task_checklist

def get_task_dependency():

    task_dependency = []
    task_dependency.append(['start', 'task1'])
    task_dependency.append(['start', 'task3'])
    task_dependency.append(['task1', 'task2'])
    task_dependency.append(['task2', 'task4'])
    task_dependency.append(['task3', 'task4'])
    task_dependency.append(['task4', 'end'])

    return task_dependency

def get_next_tasks(task_dependency, done_tasks):
    # next_tasks = [x[1] for x in task_dependency if x[0] == curr_task]
    next_tasks = []
    for done_task in done_tasks:
        _next_tasks = [x[1] for x in task_dependency if x[0] == done_task and x[1] not in done_tasks]
        next_tasks.extend(_next_tasks)

    return next_tasks

def get_logger(logfile = None):
    logger = logging.getLogger()
    logger.setLevel(level=logging.INFO)

    formatter = logging.Formatter('%(levelname)s %(asctime)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s')

    # write into file
    if logfile is not None:
        file_handler = logging.FileHandler(logfile, mode='w')
        file_handler.setLevel(level=logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # write to stdout
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument('-p', '--project-dir', help='Project dir', required=True)
    parser.add_argument('-r', '--run', help='Can be `all`, `check`, or some tasks. Use comma to run multiple tasks such as `task1,task2`', required=False, default='all')
    parser.add_argument('-d', '--dry-run', help='If specified, only print the commands but not execute commands', action='store_true', required=False)
    parser.add_argument('-f', '--force', help='Run tasks even if they are already done.', action='store_true', required=False)

    if len(sys.argv)==1:   
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

#%%
@timmer
def main():
    ### get args
    args = get_args()

    project_dir = args.project_dir
    run = args.run
    dry_run = args.dry_run
    force = args.force

    ### Create dirs and get logfile
    project = Project(project_dir)
    project_done_file = os.path.abspath(project_dir) + '.done'

    ### light check (print a more user-friendly format)
    if run == 'check_light':
        task_checklist = get_task_checklist()
        check_tasks_light(project, task_checklist)
        sys.exit(0)
    
    ### get logfile
    if dry_run or run == 'check':
        logfile = None
    
    else:
        tasks = run.replace(',', '-')
        curr_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logfile = os.path.join(project.logs_dir, 'run.%s.%s.log'%(curr_time, tasks))
    
    ### get logger
    global logger
    logger = get_logger(logfile)

    ### print dry-run mode    
    if dry_run:
        logger.info('************ Dry-run mode is on. Only list commands but do not execute them. Also not record logs. ************\n')

    ### print arguments
    logger.info('Arguments:')
    args_str = json.dumps(vars(args), indent=4)
    logger.info(args_str)

    logger.info('project_dir is %s'%(project_dir))
    logger.info('run is %s'%(run))
    logger.info('dry_run is %s'%(dry_run))

    ### check if this project is done
    if not force and os.path.exists(project_done_file):
        logger.info('************ All tasks are done. Bye! ************\n')
        sys.exit(0)

    ### Get tasks list
    task_checklist = get_task_checklist()
    task_dependency = get_task_dependency()

    ### Check run command
    if run == 'check':
        task_status, done_tasks = check_tasks(project, task_checklist)
        sys.exit(0)

    else:
        # get tasks needed to run
        tasks_to_run = []
        if run == 'all':
            tasks_to_run = task_checklist
        else:
            task_names_check = run.split(',')
            for task_name in task_names_check:
                if task_name in task_checklist:
                    tasks_to_run.append(task_name)
                else:
                    raise Exception('Error! Task %s is not in checklist'%(task_name))
        logger.info('Will running the following tasks: %s'%(', '.join(tasks_to_run)))

        # run tasks
        for task_name in tasks_to_run:
            task_cmd = eval('%s(project, args)'%(task_name))
            task_done_file = get_task_done_file(project, task_name)

            # run task
            run_task(task_name, task_cmd, task_done_file, force, dry_run)

        # check tasks
        task_status, done_tasks = check_tasks(project, task_checklist)
        if task_status == 'Done':
            pathlib.Path(project_done_file).touch()

        # print next tasks
        next_tasks = get_next_tasks(task_dependency, done_tasks)
        logger.info('Following tasks that can be run are: \n%s'%(', '.join(next_tasks)))


if __name__ == '__main__':
    main()

