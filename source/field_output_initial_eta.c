#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "output_def.h"  // some output functions

/**
    A few parameters controlling the output size. */
double snapshot_time = 0;
int OUTLEVEL = 7;
int Nslice = 512;
int MAXLEVEL = 10;

scalar pair[];
vector accel[];

int main(int argc, char *argv[]) {
  if (argc > 1)
    snapshot_time = atof(argv[1]);
  if (argc > 2)
    OUTLEVEL = atoi(argv[2]);
  if (argc > 3)
    Nslice = atoi(argv[3]);

  u.r[top] = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top] = dirichlet(0); // This is supposed to be neumann 
  u.n[bottom] = dirichlet(0);
  u.t[top] = neumann(0);
  u.t[bottom] = neumann(0);
  // TO-DO: Test if setting to neumann change 
  periodic (right);
  periodic (front);
  run();
}

event init(i=0)
{
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
    fprintf(ferr, "lsNot restored!\n");
    return 1;
  }
  if (pid() == 0) {
    fprintf(stderr, "Field restored!\n");
  }

  // Do we need to call boundary to make linear interpolation ok?
  boundary ((scalar *){u});

  char filename[100];
  bool do_linear = true;
  bool print_bin = true;
  double L0 = 2*pi;

  // Surface
  output_twophase (t, MAXLEVEL);
}




