
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

const int NTHETA = 600; //350; //200;
const int NPHI = 200; //350; //250;

inline double sph_to_cart_x(double k, double theta_k, double phi_k);
inline double sph_to_cart_y(double k, double theta_k, double phi_k);
inline double sph_to_cart_z(double k, double theta_k);

inline double window_tophat(
                     double kx, double ky, double kr, 
                     double x0, double y0, double r0,
                     double sgn, double w, double dr);

double angle_avg_window(double q, double kx, double ky, double kr,
                        double kxp, double kyp, double krp,
                        double w, double dr);
