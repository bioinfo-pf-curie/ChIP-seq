## utils.inc.sh
##
## Copyright (c) 2017 Institut Curie                               
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

set -o pipefail  # trace ERR through pipes                                                                                                                                                                 
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

###########################
## Subroutine for pipelines
###########################

die() 
{ 
    echo "Exit: $@" 1>&2 
    exit 1
}


exec_cmd()
{
    echo $*
    if [ -z "$DRY_RUN" ]; then
	eval "$@" || die 'Error'
    fi
}

exec_ret()
{
    if [ -z "$DRY_RUN" ]; then
	eval "$@" || die 'Error'
    fi
}

abspath() 
{
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

read_config()
{
    eval "$(sed -e '/^$/d' -e '/^#/d' -e 's/ =/=/' -e 's/= /=/' $1 | awk -F"=" '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)}')"
}

is_in_path()
{
    type -P $1 > /dev/null && echo 1 || echo 0
}


