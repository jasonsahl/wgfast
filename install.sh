#!/bin/bash -e -x
cd wgfast/bin/standard-raxml
# make raxmlHPC executable
make -f Makefile.SSE3.gcc
make -f Makefile.SSE3.PTHREADS.gcc
rm *.o
chmod +x raxmlHPC-PTHREADS-SSE3
chmod +x raxmlHPC-SSE3

mv raxmlHPC-PTHREADS-SSE3 ../
mv raxmlHPC-SSE3 ../





