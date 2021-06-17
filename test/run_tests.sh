#!/bin/bash


for t in ${SRT_PRIVATE_CONTEXT}/bin/${SRT_SUBDIR}/test_*; do
    $t
    if [ $? -ne 1 ]; then
	echo "*** $t failed ***"
    else
	echo "$t...passed"
    fi
done

