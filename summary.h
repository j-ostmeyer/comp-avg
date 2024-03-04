#ifndef COMP_AVG
#define COMP_AVG

double arith_mittel(double *x, unsigned n);
double stabw(double *x, unsigned n, double mu);
void transpose(double *x, double *y, unsigned n, unsigned length);
void transpose_sym(double *x, double *y, unsigned n, unsigned length);

#ifndef STANDALONE
#ifndef GAUGE_FFTW3
#define GAUGE_FFTW3
#include <fftw3.h>
#endif
#else
unsigned hoch2(unsigned n);
double complex *vfft(double *f, double complex *h, double complex *g, unsigned n, int hin);
double complex *cvfft(double complex *f, double complex *h, double complex *g, unsigned n, int hin);
double complex *chirp_z_vfft(double *f, double complex *h, double complex *a, double complex *b, double complex *bft, double complex *dummy, unsigned n, unsigned m, int hin);
#endif

double *fast_auto_cov(double *x, double mu, unsigned n);
double local_auto_cov(double *x, double mu, unsigned n, int t);
double tau_int_error(double t_int, unsigned n, unsigned t_max);

double error_on_the_error(double error, double t_int, double t_int_err);
double error_auto_weight_simple_naive(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max);
double error_auto_weight_naive(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max);
double error_auto_weight_simple(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max);
double error_auto_weight(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max);

#endif
