#!/bin/sh

log="deimos.test_min.sqr.log"

for n in 100 200 300 400 500 600 700 800 900 \
         1000 2000 3000 4000 5000 6000 7000 8000 9000 \
         10000
do
  for i in 1 2 3
  do
    echo "${i}: time ./test_min_sqr $n" 
    time ./test_min_sqr $n 
  done
done
   
   
