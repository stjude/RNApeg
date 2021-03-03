#!/bin/bash
# Imports one or more values from a config file; this must be sourced to have
# any effect.
#
# Example:
# . import_config.sh genome GRCh37-lite
#
# Note that this will exit on error; since you are sourcing it that means that
# it has the ability to abort the calling script.
#
# $1 = category name
# $2 = config name 
# $3,... = optional list of variables to import; if none are specified, then all
#          are imported by sourcing the config file

# Show usage information if no parameters were sent
if [[ "$#" == 0 ]]; then about.sh $0; exit 1; fi

# Find the config file
_CONFIG_FILE=`which_config.sh $1 $2`
if [[ ! -f "$_CONFIG_FILE" ]]; then echo "No config found for $1 $2; exiting"; exit 65; fi

# Now, process variable names
shift 2
IFS=$'\t'
# If there are no variable names, then load everything, otherwise, load
# selectively
if [[ "$#" == 0 ]]
then
  # Load everything by sourcing the config file
  #. $_CONFIG_FILE
  while read a b 
  do 
   #ignore comment lines and blank line
   read ${a} <<IN
$b
IN
  done <<< "`cat $_CONFIG_FILE | grep -v "^#" | grep -v "^$"`"
else
  # Load selectively
  _var=$1
  while [[ -n "$_var" ]]
  do
    _search="^$_var\t"
    _assign=`grep -Pe "$_search" $_CONFIG_FILE`
    if [[ -z "$_assign" ]]
    then
      echo "Could not find assignment for variable $_var in config file: $_CONFIG_FILE" >&2
      exit 64
    fi
    a=$_var 
    b=$( echo "$_assign" | sed -re "s/$_search//" )
    read ${a} <<EOF
$b
EOF
    shift
    _var=$1
  done
fi
unset IFS
