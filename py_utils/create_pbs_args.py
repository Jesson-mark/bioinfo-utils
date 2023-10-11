#!/public/home/fan_lab/wangjie/miniconda3/envs/straglr/bin/python

import sys
if len(sys.argv)==1:
    print('Usage: create_pbs_args.py "args"')
    print('Example: create_pbs_args.py "workdir input_file out_file"')
    sys.exit(1)

args = sys.argv[1]
args = args.split()
# print(args)

p_args = []
for arg in args:
    p_args.append('%s=$%s'%(arg, arg))
print('args=$(echo %s)'%(','.join(p_args)))

# args_file = sys.argv[1]
# args_file = 'scripts/a'
# args = [l.strip().split() for l in open(args_file)]


# for arg in args:
#     p_args = []
#     for a in arg:
#         p_args.append('%s=$%s'%(a, a))
#     print('args=$(echo %s)'%(','.join(p_args)))

