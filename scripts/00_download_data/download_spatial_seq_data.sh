#!/bin/bash

# download dataset GSE317410 from GEO

mkdir -p ../../data/spatial_seq_data/
cd ../../data/spatial_seq_data/
wget -nc -O "GSE317410_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE317410&format=file"
tar -xvf GSE317410_RAW.tar

# Extract all .gz files

gunzip *.gz
