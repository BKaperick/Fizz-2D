#!/bin/bash
py ../../extract_skeleton.py -d docs.md
./simulate.sh INPUT.in 125 
git add .
git commit -m "$1"
git push

