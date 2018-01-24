#1/bin/sh
rm -f dvbs2_V1.3

gcc -g -Wall -Wconversion -Wshadow -O3 -march=native -mfpmath=sse -funroll-loops  -fomit-frame-pointer -o dvbs2_V1.3  DVBS2_Xu_V1.3.c -lm

./dvbs2_V1.3
