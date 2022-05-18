#!/bin/bash

# first arg is the key
key=$1

# second arg is the file
file=$2

# #idef depth counter
num_in=0
in_key=0

echo KEY  : $key
echo FILE : $file

# loop over file
iline=0
while read line; do

  # line counter
  iline=$(($iline+1))
  #echo "LINE ($iline) : $line"

  # if not currently in "part", ...
  if [ $in_key == 0 ]; then

    # check if a part is starting
    if [ "${line:0:3}" == '#if' ]; then
      if [ ! -z "$(echo -e "$line" |grep -E "$key")" ]; then
        in_key=1
        echo "------- line $iline (START) ------------"
        echo -e "${line}"
        continue
      fi
    fi

  # if in "part", ...
  else

    # display current line
    echo -e "$line"
  
    # check if a sub-part is starting in the main part
    if [ "${line:0:3}" == '#if' ]; then

      # open new sub-part
      num_in=$(($num_in + 1))

    # check for "#endif"
    elif [ "${line:0:6}" == '#endif' ]; then

      # if all sub-parts are closed ?
      if [ $num_in == 0 ]; then
        
        # close current part
        in_key=0
        echo "------- line $iline (END) ------------"

      # esle close a sub-part
      else
        num_in=$(($num_in-1))
      fi
    fi
  fi

done < $file

