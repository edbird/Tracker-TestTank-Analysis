#!/bin/bash
cd ./output
mkdir -p ./anode_histo_good
mkdir -p ./anode_histo_fail
mkdir -p ./differential_histo_good
mkdir -p ./differential_histo_fail
#cp cell8.root cell8_backup.root
./run.out 2>/dev/null
cd ..
