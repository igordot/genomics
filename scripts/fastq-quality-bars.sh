#!/bin/bash
#
# Visualize FASTQ quality using bars and "animate" them by looping through the individual reads.
# Demo: https://asciinema.org/a/194133
# Usage: cat <file.fastq> | ./fastq-quality-bars.sh
#


# check if a stdin pipe exists
if [ -p /dev/stdin ]; then

  # blank line for cleaner output
  echo ""

  # initiate line counter
  n=0

  while read -r line; do

    # update line counter
    n=$(($n+1))

    # use only every fourth line (quality scores)
    if [ "$n" -eq 4 ]; then

      # clear screen
      printf "\033c"

      # bin quality scores
      qual8=$(echo "$line" | sed -e 'y/!"#$%&'\''()*+,-.\/0123456789:;<=>?@ABCDEFGHIJKL/00001111122222333334444455555666667777788888/')

      # convert binned quality scores to vertical bars
      awk -v q="$qual8" 'BEGIN {
        # number of bases
        len = split(q, qarr, "");
        # height of bars (plus an extra row on top and bottom)
        h = 8 + 2;
        # matrix for output characters
        for (i = 1; i <= len; i++) {
          # top border
          a[i,h] = "▁";
          # bottom border
          a[i,1] = "▔";
          # add quality bars
          for (j = 1; j <= qarr[i]; j++) {
            a[i,j+1] = "▊";
          }
        }
        # transpose matrix and print from top
        for (j = h; j >= 1; j--) {
          out = "";
          for (i = 1; i <= len; i++) {
            if (a[i,j] == "") a[i,j] = " ";
            out = out a[i,j];
          }
          print out;
        }
      }'

      # pause
      sleep 0.2

      # reset line counter
      n=0

    fi

  done

  echo "fastq ended"

else

  # show usage if nothing was piped in
  echo -e "\n Usage: cat <file.fastq> | ./fastq-quality-bars.sh \n"

fi



# end
