#! /bin/bash
if pdal --drivers | grep sritrajectory > /dev/null; then
    :
else
    echo "sritrajectory filter not available" 1>&2
    exit 1
fi

while read dataset prefix params extrafilters; do
    let i=0
    for f in $dataset/$prefix*.las; do
	pdal translate $extrafilters \
	     -f filters.sritrajectory \
	     --filters.sritrajectory.params="$params" \
	     --writers.text.precision=6 \
	     --writers.text.order=GpsTime,X,Y,Z \
	     --writers.text.keep_unspecified=false \
	     --writers.text.write_header=false --writers.text.delimiter=' ' \
	     "$f" $dataset/traj$i.txt
	let ++i
    done
done <<EOF
ak_sitka  auto dt=0.005,minsep=0.5 -f filters.sort --filters.sort.dimension=GpsTime
uh_campus C2_L dt=0.001,minsep=0.1
EOF

exit
