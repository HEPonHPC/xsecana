#!/bin/bash



for t in $1/test_*; do
    $t
    if [ $? -ne 1 ]; then
	echo "*** $t failed ***"
    else
	echo "$t...passed"
    fi
done

