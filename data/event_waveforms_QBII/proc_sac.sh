#!/bin/bash
for sacfile in *SAC
do
sac<<!
r $sacfile
decimate 5
decimate 2
w over 
q
!
done
