#!/bin/bash

build_dir=build
for f in `ls *.F90 *.h *.cpp *.c *.h *.hpp 2>/dev/null`
do
    if [ ! -d "$build_dir" ]; then
	mkdir ${build_dir}
    fi
    filename="${f%.*}"
    extension="${f##*.}"

    #m4 ${filename}.${extension} > ${build_dir}/$filename.${extension}
    ./fypp ${filename}.${extension} > ${build_dir}/$filename.${extension} -m re
    cp makefile  ${build_dir}/makefile

done
