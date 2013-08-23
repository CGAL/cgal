#!/bin/bash
dir="/home/kan/cgal-gsoc/Visibility_2/data/example"
nameend=".environment"
type1=T
type2=T
report_name=benchmark_report.txt
touch $report_name
for (( i=1; i<6; i=i+1)); do
	file="${dir}${i}${nameend}"
	for query_type in face edge vertex; do
		echo "file_name: ${file}" >>$report_name
		echo "query_type: ${query_type}" >>$report_name
		./benchmark $file $type1 $type2 $query_type >>$report_name
	done	
done
