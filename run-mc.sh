#!/bin/bash
cd ./output
#cp timestamp_cpp.root timestamp_cpp_backup.root
./run.out --input timestamp_cpp.root 2>/dev/null
#./run.out --input timestamp_cpp.root 2>/dev/null
cd ..
