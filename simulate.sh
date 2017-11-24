#!/bin/bash
py main.py INPUT.txt $1
./test $1
yes | ffmpeg -i plane_%d.png simul.gif

