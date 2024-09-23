//#include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"
#include "maxruntime.h"


//#include "adapt_wavelet_limited.h"
//#include "sandbox/frac-dist.h" // extra headerfiles used in profiling function
//#include "profile6.h"    // from Antoon

#define POPEN(name, mode) fopen (name ".ppm", mode)

double RELEASETIME = 200; 
int MAXLEVEL = dimension == 2 ? 10 : 5; // max level if not use limited refinement

double RE_tau = 180.; // Friction Reynolds number that we don't easily change 
double BO = 200; // Default Bond number
double RE = 100.; // Wave Reynolds number that is dependent on c
double k_ = 1;
double g_ = 1;
double h_ = 1;
double ak = 0.1;
double c_;
double UstarRATIO = 0.5;
double Ustar;

double uemax = 0.01;
double femax = 0.0001;
double uwemax = 0.001;
double uemaxRATIO = 0.01;

#define RATIO 1.225/1000. //density ratio, air to water
#define MURATIO 18.31e-6/10.0e-4 //dynamic viscosity ratio, air to water
// kinematic viscosity air = 16*water
double alter_MU = 1.; //not matching MURATIO extra factor. Default is 1, which 
                      //means MURATIO is used.

vector u_water[];

#include "wave_def.h"    // defines initial wave motion
#include "output_def.h"  // some output functions

/**
    We need to store the variable forcing term. */
double amp_force = 0.1; //amplitude of the forcing
face vector av[];

int main(int argc, char *argv[]) {
  maxruntime (&argc, argv); // from Nicolo
  if (argc > 1)
    RE_tau  = atof (argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    g_ = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    RELEASETIME = atof(argv[6]);
  if (argc > 7)
    uemaxRATIO = atof(argv[7]);
  if (argc > 8)
    alter_MU = atof(argv[8]);

  L0 = 2*pi;
  h_ = 1; // Water depth
  k_ = 4; // Four waves per box
  origin (-L0/2., 0, -L0/2.);
  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top] = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top] = dirichlet(0); // This is supposed to be neumann 
  u.n[bottom] = dirichlet(0);
  u.t[top] = neumann(0);
  u.t[bottom] = neumann(0);
  // TO-DO: Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  f.sigma = g_*rho1/(BO*sq(k_));
  G.y = -g_;
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  // Ustar = c_*UstarRATIO; // Depleted
  Ustar = 0.25; // Pick a fixed value
  mu2 = Ustar*rho2*(L0-h_)/RE_tau;
  // mu1 = (2.*pi/k_)*c_*rho1/RE; // Depleted: using wavelength as length scale
  mu1 = mu2/(MURATIO)/alter_MU;
  RE = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c 
  fprintf (stderr, "g = %g, c = %g, Ustar = %g, MURATIO = %g, mu_w = %g, rho_w = %g, mu_a = %g, rho_a = %g, sigma = %g, Bo = %g, RE = %g, Re_tau = %g\n", g_, c_, Ustar, MURATIO*alter_MU, mu1, rho1, mu2, rho2, f.sigma, BO, RE, RE_tau);
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << 7);
  // Refine according to 
  uemax = uemaxRATIO*Ustar;
  uwemax = 0.001*c_;
  fprintf (stderr, "RELEASETIME = %g, uemax = %g \n", RELEASETIME, uemax);
  run();
}



event init (i = 0) {
  if (!restore("restart")){
    // double rand = 0;
    double ytau = (mu2/rho2)/Ustar;
    do {
      fraction (f, WaveProfile(x,z)-y);
      foreach() {
        // rand = randInRange (0,0.1);
        // Initialize with a more accurate profile
        if ((y-WaveProfile(x,z)) > 0.05)
          u.x[] = (1-f[])*(log((y-WaveProfile(x,z))/ytau)*Ustar/0.41);
        else
          u.x[] = 0.;
        // Random noise gets killed by adaptation anyway. We wait for the instability to naturally develop instead.
        // u.x[] = (2+rand)*(1.-f[]); 
        u.y[] = 0.;
        u.z[] = 0.;
      }
      boundary ((scalar *){u});
    }
#if TREE
    // No need for adaptation when starting 
    while (0);
#else
    while (0);
#endif
  }
}

/**
  Set the wave velocity 0. */
event set_wave(i=0; i++; t<RELEASETIME) {
  fraction (f, WaveProfile(x,z)-y);
  foreach(){
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

event start(t = RELEASETIME) {
  // A slightly changed version of stokes wave as y = 0 at the bottom now so y+h -> y
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  boundary ((scalar *){u});
}

/**
   Forcing term equivalent to pressure gradient in x. */
event acceleration (i++) {
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);
}

/** 
    Output video and field. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 0.1) {

  /**
     We first do simple movies of the volume fraction, level of
     refinement fields. In 3D, these are in a $z=0$ cross-section. */

  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 40, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -3.1415);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -3.1415);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  char s[80];
  sprintf (s, "t = %0.1f", t);
  draw_string (s, size = 30);
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
  }
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 40, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.y", linear = true, n = {0,0,1}, alpha = -3.1415);
  squares ("u.x", linear = true, n = {1,0,0}, alpha = -3.1415);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  isosurface ("l2", -1);
  draw_string (s, size = 30);
  {
    static FILE * fp = POPEN ("vortex", "a");
    save (fp = fp);
  }
}

/**
   Generate averaged profile in y direction on the fly. */
// event profile_output (t += 0.1) { 
//    char file[99]; 
//    sprintf (file, "prof_%g", t); 
//    scalar uxuy[],uxux[],uyuy[],uzuz[]; 
//    foreach () { 
//      uxuy[] = u.x[]*u.y[]; 
//      uxux[] = u.x[]*u.x[]; 
//      uyuy[] = u.y[]*u.y[]; 
//      uzuz[] = u.z[]*u.z[]; 
//    } 
//    vertex scalar phi[]; 
//    foreach_vertex () 
//      phi[] = y; 
//    // default phi is y 
//    profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5, fname = file, min = 0.8, max = 2.*pi); 
//  } 

event eta_output (t += 0.1) {
  if (t > RELEASETIME) // Event does not take variable time condition 
    output_twophase (t, MAXLEVEL);
}

scalar pair[];
event turbulence_stat (t += 0.1) {
  char filename[100];
  bool do_linear = true;
  bool print_bin = true;
  int Nslice = 256;
  int OUTLEVEL = 9;
  double L0 = 2.*pi;

  foreach()
    pair[] = p[]*(1-f[]);

  // Spanwise averaged
  int res = 9;
  sprintf (filename, "./field/ux_2d_avg_t%g.bin", t); // x-vel
  output_2d_span_avg (filename,u.x  ,res, do_linear, print_bin);
  sprintf (filename, "./field/uy_2d_avg_t%g.bin", t); // y-vel
  output_2d_span_avg (filename,u.y  ,res, do_linear, print_bin);
  sprintf (filename, "./field/uz_2d_avg_t%g.bin", t); // z-vel
  output_2d_span_avg (filename,u.z  ,res, do_linear, print_bin);
  sprintf (filename, "./field/f_2d_avg_t%g.bin", t); // VoF
  output_2d_span_avg (filename,f    ,res, do_linear, print_bin);
  sprintf (filename, "./field/p_2d_avg_t%g.bin", t); // pressure
  output_2d_span_avg (filename,pair ,res, do_linear, print_bin);

  // Slices: too expensive to output for now
  if (fmod(t, 1.0) == 0.0) {
    double zslice = -L0/2.+L0/2./Nslice;
    for (int i=0; i<Nslice; i++) {
      sprintf (filename, "./field/ux_t%g_slice%d", t, i);
      sliceXY (filename,u.x,zslice,OUTLEVEL,do_linear);
      sprintf (filename, "./field/uy_t%g_slice%d", t, i);
      sliceXY (filename,u.y,zslice,OUTLEVEL,do_linear);
      sprintf (filename, "./field/uz_t%g_slice%d", t, i);
      sliceXY (filename,u.z,zslice,OUTLEVEL,do_linear);
      sprintf (filename, "./field/f_t%g_slice%d", t, i);
      sliceXY (filename,f,zslice,OUTLEVEL,do_linear);
      sprintf (filename, "./field/pair_t%g_slice%d", t, i);
      sliceXY (filename,pair,zslice,OUTLEVEL,do_linear);
      zslice += L0/Nslice;
    }
  }
}



/**
   ## End or dump regularly. */

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

/* event dumpforrestart (t += 0.1) { */
/*   char dname[100]; */
/*   u_water.x.nodump = true; */
/*   u_water.y.nodump = true; */
/*   u_water.z.nodump = true; */
/*   // p.nodump = false; */
/*   // sprintf (dname, "restart"); */
/*   dump (); */
/* } */

// event dumpforrestart (t += 1) {
//   char dname[100];
//   u_water.x.nodump = true;
//   u_water.y.nodump = true;
//   u_water.z.nodump = true;
//   // p.nodump = false;
//   sprintf (dname, "restart");
//   dump (dname);
// }

event dumpstep (t += 0.1) {
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  // pair.nodump = true;
  // p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

/** 
    Adaptive function. uemax is tuned. We need a more strict criteria for water speed once the waves starts moving. */ 
#if TREE
event adapt (i++) {
  if (i == 5)
    fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO);
  if (t < RELEASETIME)
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL);
  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,uwemax,uwemax,uwemax}, MAXLEVEL);
  }
}
#endif

