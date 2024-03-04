#include <stdio.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "summary.h"

double arith_mittel(double *x, unsigned n){
	unsigned i;
	double sum=0;
	for(i = 0; i < n; i++) sum += x[i];
	return sum/n;
}

double stabw(double *x, unsigned n, double mu){
	unsigned i;
	double sum=0, summand;
	for(i = 0; i < n; i++){
		summand = x[i]-mu;
		sum += summand*summand;
	}
	return sqrt(sum/(n-1));
}

void transpose(double *x, double *y, unsigned n, unsigned length){
	for(unsigned k = 0; k < length; k++){
		const unsigned shift = k*n;
		for(unsigned i = 0; i < n; i++) y[shift + i] = x[i*length + k];
	}
}

void transpose_sym(double *x, double *y, unsigned n, unsigned length){
	for(unsigned i = 0; i < n; i++) y[i] = x[i*length];

	unsigned k;
	for(k = 1; k < (length+1)/2; k++){
		const unsigned shift = k*n;
		for(unsigned i = 0; i < n; i++){
			const unsigned shiftI = i*length;
			y[shift + i] = .5*(x[shiftI + k] + x[shiftI + length-k]);
		}
	}

	if(k*2 == length){ // only need this in case length is even
		const unsigned shift = k*n;
		for(unsigned i = 0; i < n; i++) y[shift + i] = x[i*length + k];
	}
}

////////////////////////////////////////////////////////////////////////////////////
// Several functions that can be treated as a black box.
// They are tested and correctly perform needed Fourier trafos.
////////////////////////////////////////////////////////////////////////////////////
#ifdef STANDALONE

unsigned hoch2(unsigned n){
	unsigned potenz=1, x=2;
	while(n){
		if(n%2)	potenz *= x;
		x *= x;
		n /= 2;
	}
	return potenz;
}

	/*very fast transform: no time waisted with memory allocation.
	needs input, output and dummy arrays of lengths 2^r*/
double complex *vfft(double *f, double complex *h, double complex *g, unsigned n, int hin){
	unsigned r = (unsigned)(log2(n)*(1+DBL_EPSILON)), m=n/2, l=1;
	unsigned i, k, j;
	unsigned a, b;
	double complex w, wSchritt;
	const double wurzelN = 1/sqrt(n);
	const double complex phase=-hin*M_PI*I;

	for(j = 0; j < n; j++) g[j] = f[j]*wurzelN;
	
	for(i = 0; i < r; i++){
		w = cexp(phase/m);
		for(k = 0; k < l; k++){
			wSchritt = 1;
			for(j = 0; j < m; j++){
				a = 2*k*m + j;
				b = a+m;
				g[a] += g[b];
				g[b] = wSchritt*(g[a] - 2*g[b]);
				wSchritt *= w;
			}
		}
		m /= 2;
		l *= 2; 
	}

		/*sortieren*/
	j = 0;
	for(i = 0; i < n-1; i++){
		h[i] = g[j];
		for(m = n/2; m <= j; m /= 2){
			j -= m;
		}
		j += m;
	}
	h[n-1] = g[n-1];

	return h;
}

	/*very fast transform: no time waisted with memory allocation.
	needs input, output and dummy arrays of lengths 2^r
	input and output may be identical, dummy array may not*/
double complex *cvfft(double complex *f, double complex *h, double complex *g, unsigned n, int hin){
	unsigned r = (unsigned)(log2(n)*(1+DBL_EPSILON)), m=n/2, l=1;
	unsigned i, k, j;
	unsigned a, b;
	double complex w, wSchritt;
	const double wurzelN = 1/sqrt(n);
	const double complex phase=-hin*M_PI*I;

	for(j = 0; j < n; j++) g[j] = f[j]*wurzelN;
	
	for(i = 0; i < r; i++){
		w = cexp(phase/m);
		for(k = 0; k < l; k++){
			wSchritt = 1;
			for(j = 0; j < m; j++){
				a = 2*k*m + j;
				b = a+m;
				g[a] += g[b];
				g[b] = wSchritt*(g[a] - 2*g[b]);
				wSchritt *= w;
			}
		}
		m /= 2;
		l *= 2; 
	}

		/*sortieren*/
	j = 0;
	for(i = 0; i < n-1; i++){
		h[i] = g[j];
		for(m = n/2; m <= j; m /= 2){
			j -= m;
		}
		j += m;
	}
	h[n-1] = g[n-1];

	return h;
}

	/*input- and output-arrays of lengths n
	and four dummy arrays of lengths m=2^log(2n-1)*/
double complex *chirp_z_vfft(double *f, double complex *h, double complex *a, double complex *b, double complex *bft, double complex *dummy, unsigned n, unsigned m, int hin){
	unsigned k, test=hoch2((unsigned)(log2(n)*(1+DBL_EPSILON)));
	const double norm = sqrt((double)m/n);
	const double complex phase = hin*M_PI*I/n;

	if(n == test) return vfft(f, h, dummy, n, hin);

	for(k = 0; k < n; k++){
		b[k] = cexp(phase*k*k);
		a[k] = f[k]*conj(b[k]);
	}
	for(; k <= m-n; k++){
		a[k] = 0;
		b[k] = 0;
	}
	for(; k < m; k++){
		a[k] = 0;
		b[k] = b[m-k];
	}

	cvfft(a, a, dummy, m, 1);
	cvfft(b, bft, dummy, m, 1);
	for(k = 0; k < m; k++) a[k] *= bft[k];
	cvfft(a, a, dummy, m, -1);
	for(k = 0; k < n; k++) h[k] = a[k]*conj(b[k])*norm;

	return h;
}

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relevant functions start here
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double *fast_auto_cov(double *x, double mu, unsigned n){
	double complex *cov;
	double complex *ab=NULL, *bft=NULL, *dummy=NULL;
	double *abs_sq;
	double a, b, norm;
	unsigned i, m;

	abs_sq = malloc(n*sizeof(double));

#ifndef STANDALONE
	const double oneoverN2=1./n/n;

	m = n/2+1;
	cov = malloc(m*sizeof(double complex));

	const fftw_plan fft = fftw_plan_dft_r2c_1d(n, abs_sq, cov, FFTW_ESTIMATE);
#else
	const double wurzelN=1/sqrt(n);

	m = hoch2((unsigned)(log2(2*n-1)*(1-DBL_EPSILON)+1));

	cov = malloc(m*sizeof(double complex));
	ab = malloc(m*sizeof(double complex));
	bft = malloc(m*sizeof(double complex));
	dummy = malloc(m*sizeof(double complex));
#endif

	for(i = 0; i < n; i++) abs_sq[i] = x[i]-mu;

#ifndef STANDALONE
	fftw_execute(fft);

	abs_sq[0] = pow(creal(cov[0]), 2);
	for(i = 1; i < m; i++){
		a = creal(cov[i]);
		b = cimag(cov[i]);
		norm = a*a+b*b;
		abs_sq[i] = norm;
		abs_sq[n-i] = norm;
	}

	fftw_execute(fft);
	fftw_destroy_plan(fft);

	abs_sq[0] = creal(cov[0])*oneoverN2;
	for(i = 1; i < m; i++){
		norm = creal(cov[i])*oneoverN2;
		abs_sq[i] = norm;
		abs_sq[n-i] = norm;
	}
#else
	chirp_z_vfft(abs_sq, cov, cov, ab, bft, dummy, n, m, -1);

	for(i = 0; i < n; i++){
		a = creal(cov[i]);
		b = cimag(cov[i]);
		abs_sq[i] = a*a+b*b;
	}

	chirp_z_vfft(abs_sq, cov, cov, ab, bft, dummy, n, m, 1);

	for(i = 0; i < n; i++) abs_sq[i] = creal(cov[i])*wurzelN;
#endif

	if(cov) free(cov);
	if(ab) free(ab);
	if(bft) free(bft);
	if(dummy) free(dummy);

	return abs_sq;
}

double local_auto_cov(double *x, double mu, unsigned n, int t){
	unsigned i;
	double cov=0;

	for(i = 0; i < n-t; i++) cov += (x[i]-mu)*(x[i+t]-mu);
	// Don't use periodic boundaries.
	//t -= n;
	//for(; i < n; i++) cov += (x[i]-mu)*(x[i+t]-mu);

	return cov/(n-t);
}

double tau_int_error(double t_int, unsigned n, unsigned t_max){
	double tau = 0.5*t_int*(1+exp(-2.*t_max/t_int));
	double error_stat = 2*sqrt((t_max+0.5-t_int)/n);
	double error_syst = exp(-1.*t_max/tau);

	return t_int*(error_stat + error_syst);
}

double error_on_the_error(double error, double t_int, double t_int_err){
	return error*t_int_err/(2*t_int);
}

double error_auto_weight_simple_naive(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max){
	unsigned i;
	double cov, cov_int;
	double var=local_auto_cov(x, mu, n, 0);

	if(std_dev) *std_dev = sqrt(var*n/(n-1));

	cov_int = var*0.5;
	for(i = 1; i < n/2 && (cov=local_auto_cov(x, mu, n, i)) > 0; i++){
		cov_int += cov;
	}
	cov_int *= 1./(1-(2*i+1.)/n); // Last factor is a correction for 1/n bias.
	if(t_corr) *t_corr = cov_int/var;
	if(t_max) *t_max = i-1;

	return sqrt(2*cov_int/n);
}

double error_auto_weight_naive(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max){
	unsigned i;
	double t_int=0.5, t, g=1;
	double cov_0_inv;

	cov_0_inv = 1/local_auto_cov(x, mu, n, 0);
	if(std_dev) *std_dev = sqrt(n/(cov_0_inv*(n-1)));

	for(i = 1; i < n/2 && g > 0; i++){
		t_int += local_auto_cov(x, mu, n, i)*cov_0_inv;
		t = 1/log((2*t_int+1)/(2*t_int-1));
		g = exp(-(i/t))-t/sqrt(i*n);
	}
	t_int *= 1./(1-(2*i+1.)/n); // Last factor is a correction for 1/n bias.
	if(t_int < 0) t_int = 0.5; // Data is extremely anti-correlated. This is probably not the best solution.
	if(t_corr) *t_corr = t_int;
	if(t_max) *t_max = i-1;

	return sqrt(2*t_int/(cov_0_inv*n));
}

double error_auto_weight_simple(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max){
	unsigned i;
	double cov_int;
	double *cov=fast_auto_cov(x, mu, n);

	if(std_dev) *std_dev = sqrt(cov[0]*n/(n-1));

	cov_int = cov[0]*0.5;
	for(i = 1; i < n/2 && cov[i] > 0; i++){
		cov_int += cov[i];
	}
	cov_int *= 1./(1-(2*i+1.)/n); // Last factor is a correction for 1/n bias.
	if(t_corr) *t_corr = cov_int/cov[0];
	if(t_max) *t_max = i-1;
	free(cov);

	return sqrt(2*cov_int/n);
}

double error_auto_weight(double *x, double mu, unsigned n, double *std_dev, double *t_corr, unsigned *t_max){
	unsigned i;
	double t_int=0.5, t, g=1;
	double cov_0_inv;
	double *cov=fast_auto_cov(x, mu, n);

	cov_0_inv = 1/cov[0];
	if(std_dev) *std_dev = sqrt(n/(cov_0_inv*(n-1)));

	for(i = 1; i < n/2 && g > 0; i++){
		t_int += cov[i]*cov_0_inv;
		t = 1/log((2*t_int+1)/(2*t_int-1));
		g = exp(-(i/t))-t/sqrt(i*n);
	}
	t_int *= 1./(1-(2*i+1.)/n); // Last factor is a correction for 1/n bias.
	if(t_int < 0) t_int = 0.5; // Data is extremely anti-correlated. This is probably not the best solution.
	if(t_corr) *t_corr = t_int;
	if(t_max) *t_max = i-1;
	//for(i = 0; i < n; i++) printf("%g\n", cov[i]*cov_0_inv);
	free(cov);

	return sqrt(2*t_int/(cov_0_inv*n));
}

int main(int argc, char **argv){
	unsigned i, n;
	unsigned t_max;
	double *x, *x0;
	double mu, err, std_dev, t_corr;
	double t_int_err, err_err;
	int scheme=0, precise_out=0, length=1, sym=0;

	if(argc >= 2) scheme = atoi(argv[1]);
	if(argc >= 3) precise_out = atoi(argv[2]);
	if(argc >= 4) length = atoi(argv[3]);
	if(argc >= 5) sym = atoi(argv[4]);
	if(argc > 5){
		printf("Error: Too many parameters!\n");
		return 0;
	}

	/* Get data from STDIN */
	n = 1024;
	x = malloc(n*length * sizeof(double));
	for(i = 0; scanf("%lf ", x+i*length) > 0;){
		for(unsigned k = 1; k < length; k++) scanf("%lf ", x+i*length+k);
		if(++i == n){
			n *= 2;
			x = realloc(x, n*length*sizeof(double));
		}
	}
	n = i;

	if(n < 2){
		printf("Error: An argument is needed!\n");
		return 0;
	}

	x0 = x;
	if(length > 1){
		x0 = realloc(x, 2*n*length*sizeof(double));
		x = x0 + n*length;
		if(sym){
			transpose_sym(x0, x, n, length);
			length = length/2 + 1;
		}else transpose(x0, x, n, length);
	}

	for(unsigned k = 0; k < length; k++, x += n){
		mu = arith_mittel(x, n);
		// The autocorrelation can be calculated directly up to the needed point (naive) or completely using Fourier trafos (default).
		// The integrated autocorrelation can be summed to first zero crossing (simple) or using the Ulli Wolff method (default).
		switch(scheme){
			case 1:
				err = error_auto_weight_naive(x, mu, n, &std_dev, &t_corr, &t_max);
				break;
			case 2:
				err = error_auto_weight_simple(x, mu, n, &std_dev, &t_corr, &t_max);
				break;
			case 3:
				err = error_auto_weight_simple_naive(x, mu, n, &std_dev, &t_corr, &t_max);
				break;
			default:
				err = error_auto_weight(x, mu, n, &std_dev, &t_corr, &t_max);
		}

		t_int_err = tau_int_error(t_corr, n, t_max);
		err_err = error_on_the_error(err, t_corr, t_int_err);

		if(precise_out)
			printf("%d\t%.15g\t%.15g\t%.15g\t%.15g\t%d\t%.15g\t%.15g\n", length>1?k:n, mu, std_dev, err, t_corr, t_max, t_int_err, err_err);
		else
			printf("%d\t%g\t%g\t%g\t%g\t%d\t%g\t%g\n", length>1?k:n, mu, std_dev, err, t_corr, t_max, t_int_err, err_err);
	}

	free(x0);

	return 0;
}
