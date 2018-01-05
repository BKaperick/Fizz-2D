#!/bin/bash
rm *png
rm *txt
py main.py $1 $2
start=0
end=126
while [ $end -le $2 ]; do
    ./draw $start $end
    (( start+=126 ))
    (( end+=126 ))
done
echo $end
./draw $start $2
yes | ffmpeg -i plane_%03d.png simul.gif

