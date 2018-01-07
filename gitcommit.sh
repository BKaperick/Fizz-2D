#!/bin/bash
py ../../extract_skeleton.py -d docs.md
./simulate.sh $1 $2
git add .
git commit -m "$3"
git push

