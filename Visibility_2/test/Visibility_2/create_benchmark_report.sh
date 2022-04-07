#!/bin/bash
dir="../../data/example"
nameend=".environment"
type1=S
type2=T
report_name="${type1}${type2}_report.txt"
echo "" >$report_name
for (( i=1; i<9; i=i+1)); do
	file="${dir}${i}${nameend}"
	for query_type in face edge vertex; do
		echo "file_name: ${file}" >>$report_name
		echo "query_type: ${query_type}" >>$report_name
		echo "regularize: yes" >>$report_name
		./simple_benchmark $file $type1 $type2 $query_type true >>$report_name
		echo "regularize: no" >>$report_name
		./simple_benchmark $file $type1 $type2 $query_type false >>$report_name
	done
done