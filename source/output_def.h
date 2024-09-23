
/** 
    Outputting slices on the fly. */
void sliceXY(char * fname, scalar s, double zp, int maxlevel, bool do_linear) {

  FILE * fpver = fopen (fname,"w");
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;

  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      field[i][j] = interpolate(s,xp,yp,zp,do_linear);
    }
  }

  if (pid() == 0) { // only the master prints!
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	      fwrite ( &field[i][j], sizeof(double), 1, fpver);
      }
    }
    fflush (fpver);
  }
  matrix_free (field);
  fclose (fpver); // we close at the end
}

void sliceXZ(char * fname, scalar s, double yp, int maxlevel, bool do_linear) {

  FILE * fpver = fopen (fname,"w");
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;

  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double zp = stp*j + Z0 + stp/2.;
      field[i][j] = interpolate(s,xp,yp,zp,do_linear);
    }
  }

  if (pid() == 0) { // only the master prints!
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	      fwrite ( &field[i][j], sizeof(double), 1, fpver);
      }
    }
    fflush (fpver);
  }
  matrix_free (field);
  fclose (fpver); // we close at the end
}

/** 
   Perform the spanwise averaging in a memory efficient way. */

void output_2d_span_avg (char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary

  //fprintf(stderr, "I am here 1\n"), fflush (stderr);
  double ** field = (double **) matrix_new (nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      field[i][j] = 0.0;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
        double val = nodata;
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
	}
	field[i][j] += val != nodata ? val/(double)nn : 0.0; 
	/*
	if (val != nodata) {
          field[i][j] += val/(double)nn;
	}
	*/
      }
    }
  }

  //fprintf(stderr, "I am here 2\n"), fflush (stderr);
  FILE * fpver = fopen (fname,"w");
  if (pid() == 0) { // master

    // reduce it to the first pid()
  #if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
  #endif

    // print
    if(print_bin) { // print in binary format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
	  fwrite ( &field[i][j], sizeof(double), 1, fpver );
        }
      }
    }
    else { // print in ascii format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          fprintf(fpver, "%8E", field[i][j]);
          fputc('\n', fpver);  // not double quotation""
        }
      }
    }
    fflush(fpver);
    //fprintf(stderr, "I am here 4\n"), fflush (stderr);

  }
  #if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, sq(nn), MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
  #endif
  //fprintf(stderr, "I am here 5\n"), fflush (stderr);

  matrix_free (field);
  fclose (fpver); // we close at the end
}



/**
   Output eta on the fly. */

void output_twophase (double snapshot_time, int MAXLEVEL) {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,epsilon\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  foreach(){
    if (interfacial (point, f)){
      if (point.level == MAXLEVEL) {
	      coord n = mycs (point, f);
	      double norm_2 = sqrt(sq(n.x) + sq(n.y));
      	double eta = pos[];
        fprintf (feta, "%g,%g,%g,%g\n", x, z, eta, -n.x/n.y);
      }
    }
  }
  fclose (feta);
} 

/** 
   Creation of a 3D matrix. */

void * matrix_new_3d (int nx, int ny, int nz, size_t size) {
  void ** m = qmalloc(nx, void *);  //Define a pointer that points to every x coordinate.
  char * a  = qmalloc(nx*ny*nz*size, char);
  for (int i=0; i<nx; i++) {
    m[i] = a+i*ny*nz*size;
  }
  if (pid() == 0) {
    fprintf(stderr, "Matrix created!\n");
  }
  return m;
}

/** 
   Output a 3d field in a linearly interpolated uniform grid. */

void output_3d(char * fname, scalar s, int maxlevel, bool do_linear, bool print_bin) {

  FILE * fpver = fopen (fname,"w");
  int nn = (1<<maxlevel);
  double stp = 0.999999*(L0+X0-X0)/(double)nn; // to avoid interpolated point coincides with fine grid boundary
  double ** field = (double **) matrix_new_3d (nn, nn, nn, sizeof(double));
  for (int i = 0; i < nn; i++) {
    double xp = stp*i + X0 + stp/2.;
    for (int j = 0; j < nn; j++) {
      double yp = stp*j + Y0 + stp/2.;
      for (int k = 0; k < nn; k++) {
        double zp = stp*k + Z0 + stp/2.;
	      double val = nodata;
        foreach_point (xp, yp, zp, serial) {
          val = do_linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
	}
        field[i][j*nn+k] = val;
      }
    }
  }
  if (pid() == 0) { // master
  #if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  #endif
    if(print_bin) { // print in binary format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          for (int k = 0; k < nn; k++) {
            fwrite ( &field[i][j*nn+k], sizeof(double), 1, fpver );
          }
        }
      }
    }
    else { // print in ascii format
      for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
          for(int k = 0; k < nn; k++) {
            fprintf(fpver, "%.10e", field[i][j*nn+k]);
            fputc('\n', fpver);  // not double quotation""
          }
        }
      }
    }
    fflush(fpver);
  }
  #if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, cube(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
  #endif
  matrix_free (field);
  fclose (fpver); // we close at the end
}


// void output_twophase_locate (double snapshot_time) {
//   scalar pos[];
//   coord G = {0.,1.,0.}, Z = {0.,0.,0.};
//   position (f, pos, G, Z);
//   char etaname[100];
//   sprintf (etaname, "./eta/eta_loc_t%g_%d", snapshot_time, pid());
//   FILE * feta = fopen (etaname, "w");
//   fprintf(feta, "x,z,pos,epsilon,p,dudy,dvdx,dudx,dvdy\n");
//   // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
//   double stp = L0/256.;
//   foreach(){
//     if (interfacial (point, f)){
//       if (point.level == MAXLEVEL) {
// 	coord n = mycs (point, f);
// 	double norm_2 = sqrt(sq(n.x) + sq(n.y));
// 	double eta = pos[];
// 	double yp = y + stp;
// 	point = locate (x, yp, z);
// 	// Getting the local normal vector
// 	// n is norm 1 and has to be normalized to norm 2
// 	// double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
// 	// fprintf (stderr, "Inline 1! \n");
// 	if (point.level > 0) {
// 	  POINT_VARIABLES;
// 	  /* fprintf (stderr, "Inline 2! \n"); */
// 	  /* fprintf (stderr, "Above cell x = %g, y = %g \n", x, y); */
// 	  /* fprintf (stderr, "Above cell Delta = %g \n", Delta); */
// 	  /* fprintf (stderr, "Above cell ux = %g \n", u.x[]); */
// 	  double dudy1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
// 	  //double dudy2 = (u.x[0,2] - u.x[0,0])/(2.*Delta);
// 	  // double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
// 	  double dvdx1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
// 	  //double dvdx2 = (u.y[1,1] - u.y[1,-1])/(2.*Delta);
// 	  double dudx1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
// 	  //double dudx2 = (u.x[1,1] - u.x[1,-1])/(2.*Delta);
// 	  double dvdy1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
// 	  //double dvdy2 = (u.y[0,2] - u.y[0,0])/(2.*Delta);
// 	  /*  tau.x[] = 2*mu2*(SDeform.x.x[]*n.x + SDeform.y.x[]*n.y)/norm_2; */
// 	  /*  tau.y[] = 2*mu2*(SDeform.x.y[]*n.x + SDeform.y.y[]*n.y)/norm_2;  */
// 	  /* double tau1 = 2*mu2*(dudy1 + dvdx); */
// 	  /* double tau2 = 2*mu2*(dudy2 + dvdx); */
// 	  fprintf (feta, "%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
// 		   x, z, eta, -n.x/n.y, p[0,0], dudy1, dvdx1, dudx1, dvdy1);
// 	/* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
// 	}
//       }
//     }
//   }
//   fflush (feta);
//   fclose (feta);
// }
