#!/bin/bash
py ../../extract_skeleton.py -d docs.md
./simulate.sh INPUT.in $1
git add .
git commit -m "$2"
git push

