#!/bin/bash
set -x
cp -f ../code/weno7.out ./
bsub -o run.log -b -q q_sw_expr -n 4 -np 4 -cgsp 64 -share_size 6000 -host_stack 512 -priv_size 16 -pe_stack 3 ./weno7.out
