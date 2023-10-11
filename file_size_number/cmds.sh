# 提交任务，该任务会每天计算一次服务器存储占用
nohup bash cal_storage_per_day.sh > logs/cal_storage_per_day.log 2>&1 &


