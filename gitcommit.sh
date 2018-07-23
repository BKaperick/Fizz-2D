#!/bin/bash
py ../../extract_skeleton.py -d docs.md
./simulate.sh $1 $2
sed -i -e 's/py ./python ./g' simulate.sh
git add .
git commit -m "$3"
git push
sed -i -e 's/python ./py ./g' simulate.sh
