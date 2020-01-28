#!/bin/bash
shopt -s extglob
shopt -s expand_aliases
# Remarks:
# pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname $depend
# note PBS scripts won't expand aliases
# to enforce aliases, use:
# shopt -s expand_aliases

# customize this to suit the computing environment (queues etc)
alias pbsgen=${BASH_SOURCE%/*}/../pbsgen/pbsgen_respublica.py

#########################
## Usage
#########################
function pipeline_template {
  cat <<'EOF'
source $WZSEQ_ENTRY
wzref_hg19
pipeline_prepare

while read sname; do
  jump_comments

done <<EOM
sample1
sample2
EOM
EOF
}

###########################################
## Helper functions
###########################################

_pipeline_prepare=$(cat <<'EOF'
callargarray=("$@"); base=$(pwd); [[ -d pbs ]] || mkdir pbs; leveljobids=();
# echo ${#callargarray[@]},0
if [[ ${#callargarray[@]} > 0 && ${callargarray[-1]} == "do" ]]; then
  pipeline_submit=true
  unset 'callargarray[${#callargarray[@]}-1]'
else
  pipeline_submit=false
fi
# echo ${#callargarray[@]},1

if [[ ${#callargarray[@]} > 0 && ${callargarray[-1]} =~ ^[0-9]+$ ]]; then
  pipeline_select=${callargarray[-1]}
elif [[ ${#callargarray[@]} > 0 && "${callargarray[-1]}" =~ ([0-9]*)-([0-9]*) ]]; then
  pipeline_select=${callargarray[-1]}
  pipeline_range=true
  pipeline_range_start=${BASH_REMATCH[1]}
  pipeline_range_end=${BASH_REMATCH[2]}
else
  pipeline_select="all"
fi
# echo ${#callargarray[@]},2
submitted_jobids=()
level_jobids=()

# echo ${#callargarray[@]},3
# echo "submit:" $pipeline_submit
# echo "select:" $pipeline_select

# pipeline_submit=false
hour=24; memG=10; ppn=1; queue=default
EOF
)

alias pipeline_prepare='eval "$_pipeline_prepare"'
alias jump_comments='sname_re="^#"; [[ "$sname" =~ $sname_re ]] && continue; depend="";'

## define pipeline_select and pipeline_submit
function pipeline_eval {
  pipeline_component=$1
  if [[ "$pipeline_select" == "all" || "$pipeline_component" == "$pipeline_select" || ( "${pipeline_select: -1}" == "+" && "$pipeline_component" -ge "${pipeline_select::-1}" ) || ($pipeline_range && $pipeline_component -ge $pipeline_range_start && $pipeline_component -le $pipeline_range_end) ]]; then
    echo "submit:" $pipeline_submit
    echo "select:" $pipeline_select
    echo "component:" $pipeline_component

    $2
    pbsfn=$base/pbs/${jobname}.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour $hour -memG $memG -ppn $ppn -queue $queue -workd $(pwd)

    ## whether to submit
    if $pipeline_submit; then
      jobid=$(qsub -terse $depend $pbsfn)
      # jobid=`date "+%Y-%m-%d_%H-%M-%S"`
      echo "depend: " $depend
      echo "submitted jobid: " $jobid
      echo

      # by default, let next job depend on this one.
      # depend="-W depend=afterok:$jobid"
      depend="-hold_jid $jobid"

      # collect job id
      submitted_jobids[$pipeline_component]=$jobid
      
      # collect all jobids on the level
      if [[ -z ${leveljobids[$pipeline_component]} ]]; then
        leveljobids[$pipeline_component]=$jobid
      else
        leveljobids[$pipeline_component]="${leveljobids[$pipeline_component]}:$jobid"
      fi
    fi
  fi
}

function pipeline_depend {
  depended_component=$1

  if [[ $depended_component == 'none' ]]; then
    depend=''
  fi

  # if dependency is not submitted, assume it's done
  if [[ -z ${submitted_jobids[$depended_component]} ]]; then
    depend=""
  else
    depend="-hold_jid $jobid"
    # "-W depend=afterok:"${submitted_jobids[$depended_component]}
  fi
}

function pipeline_dependlevel {
  if [[ $# -lt 1 ]]; then
    depend=""
  elif [[ -z "${leveljobids[$1]}" ]]; then
    depend=""
  else
    depend="-W depend=afterok:${leveljobids[$1]}"
  fi
}

##############################################
# source the files under references/ and src/
##############################################
source ${BASH_SOURCE%/*}/../references/CHOP_HPC.sh
# for fn in ${BASH_SOURCE%/*}/../references/*.sh; do
#   source $fn;
# done

for fn in ${BASH_SOURCE%/*}/../src/**/*.sh; do
  source $fn;
done

