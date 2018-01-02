#!/bin/bash
py main.py $1 $2
./draw 0 $2
yes | ffmpeg -i plane_%d.png simul.gif

