alias pwdp='pwd -P'

# pbs
alias qme='qstat -a |grep wangjie'
alias qmer='qstat -r -u wangjie |grep wangjie'
alias peme='pestat |grep "wangjie"'
alias qfan='qstat -a | grep shfan'
alias qpan='qstat -a | grep mppan'

# simpler commands
alias ds='du -sh'
alias wl='wc -l'
alias le='less'
alias len='less -N'
alias lt='ls -lhtr'
alias ll='ls -lh --time-style long-iso'
alias ca='cat'
alias he='head'
alias grepn='grep -n'
alias tl='tail'
alias catt='cat -T'
alias mv="mv -i"
alias cp="cp -i"
alias rm="rm -i"

#cd
alias ..='cd .. && ls -lh'
alias ...='cd ../..'

alias condaac='conda activate'
alias condade='conda deactivate'
alias condain='conda info -e'

alias lln="ls -lhtr  --time-style long-iso | tac | cat -n | tac | sed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
function lf() {
    if [ "x${1}" == "x" ]
    then
        n=1 
    else
        n="${1}"
    fi  
    ls -rt1 | tail -n ${n} | head -n 1
}

alias topu='top -u wangjie'

# software aliases
#alias R=/public/home/fan_lab/wangjie/miniconda3/envs/r4/bin/R

alias qfan='qstat -a | grep shfan'
alias grep='LC_ALL=C grep --color=auto'

alias psme='ps aux |grep wangjie'

alias df_group_share='df -h /public/group_share_data/fan_lab/'
alias mambaac='mamba activate'
alias he3='head -n3'
alias tl3='tail -n3'
alias l='ls -lh --time-style long-iso'

