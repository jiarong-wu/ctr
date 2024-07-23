/** 
    Specify the interface shape. Since we are considering ak=0.1/0.2/0.3, 
    use Stokes instead of linear profile. */

double WaveProfile_linear (double x, double z) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  return eta1 + h_;
}

double WaveProfile (double x, double z)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;
}

/**
   Start the wave at RELEASETIME. We don't do any adaptation at this step. 
   And we use linear wave instead of stokes. (UPDATE: use stokes for 0.1/0.2/0.3) */
// #include "test/stokes.h"
double u_x_linear (double x, double y) {
  return ak*c_*cos(x*k_)*exp(y*k_);
}
double u_y_linear (double x, double y) {
  return ak*c_*sin(x*k_)*exp(y*k_);
}

// TO-DO: maybe include g correction by Bond number. High Bond should be fine.
double u_x (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*3.*k_*cos(3.*k_*x);
}

double u_y (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*sin(3.*k_*x);
}

/**
   A function to generate a random field initially. */
double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}