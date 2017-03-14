#!/bin/bash

build_dir=./build

if [ ! -d "$build_dir" ]; then
    mkdir ${build_dir}
fi

for d in `find . -type d -maxdepth 1`
do
    if [ "$d" == "$build_dir" ]
    then
	continue
    fi

    if [ ! -d "${build_dir}/${d}" ]; then
    	mkdir ${build_dir}/${d}
    fi
    
    for f in `find $d -maxdepth 1 -type f \
              \( -iname \*.F90 -o -iname \*.h -o -iname \*.cpp \) \
              2>/dev/null`
    do
    	filename="${f%.*}"
    	extension="${f##*.}"

    	#m4 ${filename}.${extension} > \
	#   ${build_dir}/$filename.${extension}
    	#./fypp ${filename}.${extension} > \
	#  ${build_dir}/$filename.${extension} -m re
	
    	cp makefile  ${build_dir}/makefile
    	#echo ./${build_dir}/$f
    	#echo $f | cut -d '\.\/' -f 2
	./fypp -n  -p -m re --create-parents $f > \
	       ${build_dir}/`echo $f | cut -c3-`
	
	echo " >>>processing file " \
	     ${build_dir}/`echo $f | cut -c3-`
    done
done
