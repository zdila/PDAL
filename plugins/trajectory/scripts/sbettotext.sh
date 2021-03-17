#! /bin/bash
while read dataset zone; do
    pdal translate -r sbet -f filters.reprojection \
	 --filters.reprojection.in_srs=EPSG:4326 \
	 --filters.reprojection.out_srs=EPSG:326$zone \
	 --writers.text.precision=6 \
	 --writers.text.order=GpsTime,X,Y,Z \
	 --writers.text.write_header=false \
	 --writers.text.keep_unspecified=false \
	 --writers.text.delimiter=' ' \
	 $dataset/*.out $dataset/sbet.txt
done <<EOF
uh_campus 15
ak_sitka  08
EOF
