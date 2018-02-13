#!/bin/bash
cd ./output
mkdir -p ./anode_histo_good
mkdir -p ./anode_histo_fail
mkdir -p ./differential_histo_good
mkdir -p ./differential_histo_fail
#cp timestamp_cpp.root timestamp_cpp_backup.root
./run.out --input timestamp_cpp.root 2>/dev/null
#./run.out --input timestamp_cpp.root 2>/dev/null
cd ..
