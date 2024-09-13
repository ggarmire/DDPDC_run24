#!/usr/bin/env bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL="$HOME/install"

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL



filename=$1
outname=$2
root -b -l -q "/direct/sphenix+tg+tg01/jets/ajwood3/pi0DataProduction/macro/pi0_position_analysis.c(\"$filename\",\"$outname\")"

