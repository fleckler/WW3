#!/bin/bash

#set -e
#set -x

# get arguments
if [ $# == 2 ]; then
  switch=$1
  debug=$2
else
  echo "ERROR: 2 argument [SWITCH, DEBUG] are requiere!!!"
  exit 10
fi
  
# switch file
switch_file="switch_${switch}"
if [ ! -f ${switch_file} ]; then
  echo "ERROR: Switch file ${switch_file} NOT FOUND !!!"
  exit 10
fi

# debug ?
if [ $debug == 0 ]; then
  comp='datarmor_mpt'
  suffix="${switch}_mpt"
elif [ $debug == 1 ]; then
  comp='datarmor_mpt_debug'
  suffix="${switch}_mpt_debug"
else
  echo "ERROR: DEBUG argument must be 0 or 1 !!!"
  exit 10
fi

# check env. var. $PATH
if [ -z "$( echo ${PATH} | grep 'NETCDF2019/MPT')" ]; then
  echo "ERROR: NetCDF (MPT) not found in path"
  exit 250
elif  [ ! -z "$( echo ${PATH} | grep 'NETCDF2019/INTEL')" ]; then
  echo "ERROR: NetCDF for both INTEL and MPT is found in path"
  exit 250
fi

# remove old files
rm -fr ../exe_${suffix}

# compilation
if [ -f "comp_${comp}" -a -f "link_${comp}" ]; then
  ./w3_clean
  cp comp_${comp}      comp
  cp link_${comp}      link
  cp switch_${switch}  switch
else
  ./w3_clean -c
  ./w3_setup -c "${comp}" -s ${switch} ../
  cp comp comp_${comp}
  cp link link_${comp} 
fi

./w3_automake
if [ $? != 0 ]; then
  exit 100
fi

# rename directories
cp -vr ../exe ../exe_${suffix}
cp -vr ../tmp_SEQ ../exe_${suffix}/tmp_SEQ
if [ "${switch}" != 'STX' -a  ${switch} != 'STX_TPROF' ]; then
  cp -vr ../tmp_MPI ../exe_${suffix}/tmp_MPI
fi



