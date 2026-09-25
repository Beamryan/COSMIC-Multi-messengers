#!/bin/bash

head -n 101 gx_dat_full_monochrome.csv > gx_dat_full_300.csv

head -n 101 gx_dat_full_chirping.csv | tail -100 >> gx_dat_full_300.csv

head -n 101 gx_dat_full_confused.csv | tail -100 >> gx_dat_full_300.csv

# sll versions of the files
head -n 101 sll_gx_dat_full_monochrome.csv > sll_gx_dat_full_300.csv

head -n 101 sll_gx_dat_full_chirping.csv | tail -100 >> gx_dat_full_300.csv

head -n 101 sll_gx_dat_full_confused.csv | tail -100 >> gx_dat_full_300.csv

