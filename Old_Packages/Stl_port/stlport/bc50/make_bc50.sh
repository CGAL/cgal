#!/bin/sh
for file in `cat ../export_names`
do
rm -fr $file.h
cat stl_tmpl.h | sed -e "s/REPLACEME/$file/g" > $file.h
done

mv algorithm.h algorith.h
mv functional.h function.h
mv stdexcept.h stdexcep.h
mv streambuf.h streambu.h
mv strstream.h strstrea.h
