#!/usr/bin/env bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL="$HOME/install"

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

events=$1
filename=$2
outname1=$3
root -b -l -q "macro/Fun4All_CaloTreeGen.C($events,\"$filename\",\"$outname1\")"

