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
    fprintf(ferr, "Not restored!\n");
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

  // scalar pair[];
  // foreach () {
  //   pair[] = p[]*(1-f[]);
  // }

  // Spanwise averaged
  int res = 9;
  sprintf (filename, "./field/ux_2d_avg_t%g.bin", snapshot_time); // x-vel
  output_2d_span_avg (filename,u.x  ,res, do_linear, print_bin);
  sprintf (filename, "./field/uy_2d_avg_t%g.bin", snapshot_time); // y-vel
  output_2d_span_avg (filename,u.y  ,res, do_linear, print_bin);
  sprintf (filename, "./field/uz_2d_avg_t%g.bin", snapshot_time); // z-vel
  output_2d_span_avg (filename,u.z  ,res, do_linear, print_bin);
  sprintf (filename, "./field/f_2d_avg_t%g.bin", snapshot_time); // VoF
  output_2d_span_avg (filename,f    ,res, do_linear, print_bin);
  // sprintf (filename, "./field/p_2d_avg_t%g.bin", snapshot_time); // pressure
  // output_2d_span_avg (filename,pair ,res, do_linear, print_bin);
  fprintf (stderr, "Spanwise average output finished!\n");

  // Slices
  // double zslice = -L0/2.+L0/2./Nslice;
  // for (int i=0; i<Nslice; i++) {
  //   sprintf (filename, "./field/ux_t%g_slice%d", snapshot_time, i);
  //   sliceXY (filename,u.x,zslice,OUTLEVEL,do_linear);
  //   sprintf (filename, "./field/uy_t%g_slice%d", snapshot_time, i);
  //   sliceXY (filename,u.y,zslice,OUTLEVEL,do_linear);
  //   sprintf (filename, "./field/uz_t%g_slice%d", snapshot_time, i);
  //   sliceXY (filename,u.z,zslice,OUTLEVEL,do_linear);
  //   sprintf (filename, "./field/f_t%g_slice%d", snapshot_time, i);
  //   sliceXY (filename,f,zslice,OUTLEVEL,do_linear);
  //   // sprintf (filename, "./field/pair_t%g_slice%d", t, i);
  //   // sliceXY (filename,pair,zslice,OUTLEVEL,do_linear);
  //   zslice += L0/Nslice;
  // }
  // fprintf (stderr, "Slices output finished!\n");

  // Surface
  output_twophase (t, MAXLEVEL);
}




