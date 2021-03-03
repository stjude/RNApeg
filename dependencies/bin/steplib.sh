#!/bin/bash
# Provides functions for use by setup scripts for step-based pipelines.
#
# IMPORTANT: This script is intended to be sourced.
#
# The general use pattern is:
# 1. Source this script
# 2. call set_step_script_dir DIR RUN_TYPE
# 3. For each step, do:
# 3a. init_step name, where name is the name of the step
# 3b. Use the "getters" to get the names of files and directories for the
#     step, and generate the scripts, writing to those locations.
# 3c. For the submission script, you can use write_step_submit_script, if you
#     are following that method's conventions.
#
# Then you are done!
#
# Example:
#   set_step_script_dir /path/to/run
#   init_step finishalign
#   cat > `get_step_qc_script` <<EOF
#   #!/bin/bash
#   # QC:
#   anyfail=no
#   for sample in $SAMPLES
#   do
#     if ! qcquiet.sh $failed_qc_dir/\$sample qc_...
#     then anyfail=yes
#     fi
#   done
#   EOF
#   cat > `get_step_make_cmds_script` <<EOF
#   #!/bin/bash
#   make_script_... > `get_step_cmds_file`
#   EOF
#   write_step_submit_script

_STEPLIB_RUN_DIR=`pwd`
_STEPLIB_CUR_STEP_NUM=0
_STEPLIB_CUR_STEP_NAME=
_STEPLIB_RUN_TYPE=

# Sets the run directory into which scripts are written
# $1 = Run directory 
# $2 = Run type (Ex: "mapping-standard", "mapping-rnaseq", "std-analysis") 
function set_step_script_dir {
  _STEPLIB_RUN_DIR=$1
  _STEPLIB_RUN_TYPE=$2
  mkdir -p $_STEPLIB_RUN_DIR
}

# Initializes a step and returns the script name to use for the step.
# This function:
# 1. Increments the step number, making the appropriate step information
#    available from the getters.
# 2. Creates the files directory with subdirs
# 3. Writes the last_step.txt file
#
# $1 = step name
function init_step {
  let ++_STEPLIB_CUR_STEP_NUM
  _STEPLIB_CUR_STEP_NAME=$1
  mkdir -p `get_step_log_dir` `get_step_failed_qc_dir`
  echo $_STEPLIB_CUR_STEP_NAME > $_STEPLIB_RUN_DIR/last_step.txt
}

# Returns the full path to the script that should be written to perform QC
function get_step_qc_script {
  _get_step_subscript "c"
}

# Returns the full path to the script that can be used by external code for
# custom hooks.  Note that the hook runs after QC of the preceeding step's
# results, but before any other work on this step.
function get_step_hook_script {
  _get_step_subscript "h"
}

# Returns the full path to the script that should be written to perform local
# work
function get_step_local_work_script {
  _get_step_subscript "l"
}

# Returns the full path to the script that should be written to generate
# commands for submission
function get_step_make_cmds_script {
  _get_step_subscript "q"
}

# Returns the full path to the script that should be written to submit commands
function get_step_submit_script {
  _get_step_subscript "y"
}


# Automatically writes the step submit script.
# This uses:
# - AFC_ARGS_## for any step-specific arguments to submit
# - the step log dir for the log template
# - the commands file
# All of these are used automatically, and there is no way to override.
function write_step_submit_script {
  args_var=`printf "AFC_ARGS_%02d" $_STEPLIB_CUR_STEP_NUM`
  args_var=\$$args_var
  cat > `get_step_submit_script` <<EOF
#!/bin/bash
sub_array_for_cmdfile.sh `get_step_cmds_file` --log-tpl `get_step_log_dir`/%W.%J.%I.txt `eval echo $args_var`
EOF
}

# Returns the full path to the log dir
function get_step_log_dir {
  _get_step_subdir "logs"
}

# Returns the full path to the failed_qc dir
function get_step_failed_qc_dir {
  _get_step_subdir "failed_qc"
}

# Returns the full path to the commands file. Accepts an optional 
# argument that can be included in the commands file name.
# $1 = Name suffix (optional) (Used in RNA-Seq mapping to include DB name in the file name)
function get_step_cmds_file {
  if [ "$1" != "" ]
  then
    printf "%s/cmds-%02d-%s.sh" $_STEPLIB_RUN_DIR $_STEPLIB_CUR_STEP_NUM $1
  else
    printf "%s/cmds-%02d.sh" $_STEPLIB_RUN_DIR $_STEPLIB_CUR_STEP_NUM
  fi
}

# Returns the full path to one of the scripts, based on letter
# Parameters:
# $1 = letter to append
function _get_step_subscript {
  printf "%s/%02d%s.sh" $_STEPLIB_RUN_DIR "${_STEPLIB_CUR_STEP_NUM#0}" $1
}

# Returns the full path to one of the directories, based on subdir
# Parameters:
# $1 = subdir
function _get_step_subdir {
  printf "%s/files-%02d-%s/%s" $_STEPLIB_RUN_DIR $_STEPLIB_CUR_STEP_NUM $_STEPLIB_CUR_STEP_NAME $1
}

# To be called after initial setup, returns the final postprocess script
# DEPRECATED: use get_end_hook_script
function get_final_postprocess_script {
  get_end_hook_script
}

# Gets the hook script for after all other processing is complete.
# Only call this after the setup script has been run.
function get_end_hook_script {
  _STEPLIB_CUR_STEP_NUM=$(ls $_STEPLIB_RUN_DIR/*.sh | sed -re "s#$_STEPLIB_RUN_DIR/##" | grep -e "^[0-9]\+.\.sh" | sed -re "s/.\.sh//" | tail -n 1 | sed 's/^0*//')
  init_step end_hook
  get_step_hook_script
}

function get_failure_hook_script {
  printf "%s/.failure.sh" $_STEPLIB_RUN_DIR 
}

function get_preprocess_hook_script {
  printf "%s/.preprocess.sh" $_STEPLIB_RUN_DIR
}

function run_step {
  STEP=$1
  FORCE=$2
  _STEPLIB_CUR_STEP_NUM=$STEP
  for script in $_STEPLIB_RUN_DIR/$STEP*.sh
  do
    # Run the script, falling back to bash if eval fails (e.g. not executable)
    eval $script
    exitcode=$?
    if [ $exitcode == 126 ];
    then
      echo "Trying to recover from permission denied by running as a bash script" >&2
      bash $script
      exitcode=$?
    fi
    
    # Exit if appropriate
    if [ $exitcode != 0 ]
    then
      if [ "$FORCE" == "force" -a "$(get_step_qc_script)" == "$script" ]
      then :
      else exit $? 
      fi
    fi
  done
}
