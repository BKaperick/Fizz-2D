#!/bin/bash
rm *png
rm *txt
py main.py $1 $2
./draw $2
yes | ffmpeg -i plane_%d.png simul.gif

