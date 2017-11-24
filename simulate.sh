#!/bin/bash
py main.py $1 $2
./test $2
yes | ffmpeg -i plane_%d.png simul.gif

