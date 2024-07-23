
CC99='mpicc -std=c99' qcc -disable-dimensions -grid=octree -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DMTRACE=3 -DTREE -o field_output field_output.c -L$BASILISK/gl -lfb_tiny -lm
