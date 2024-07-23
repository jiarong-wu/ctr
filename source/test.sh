# CC99='mpicc -std=c99' qcc -grid=octree -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DTREE -o NWP NWP.c -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
CC99='mpicc -std=c99' qcc -disable-dimensions -grid=octree -autolink -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -DTREE -o NWP NWP.c -L$BASILISK/gl -lfb_tiny -lm
