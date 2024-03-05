# Compare Averages of time series and more

by Johann Ostmeyer

A light-weight tool to get the most important info from a time series and to COMPare the AVeraGe of different time series.

Get the average, standard deviation, standard error, integrated autocorrelation time and corresponding errors quickly and reliably with a simple bash command.

For questions concerning the code contact [ostmeyer@hiskp.uni-bonn.de](mailto:ostmeyer@hiskp.uni-bonn.de).

## Installation

The entire analysis is implemented in `C` with a `bash` front-end.

After downloading the code, enter the newly created directory and type
```
make
```
to create the executable. The tool is now ready to be used. It is recommended, however, to copy the two files `comp-avg` and `summary` into a directory that is part of your environment `$PATH`, e.g. `~/bin/`.

The default code depends on `FFTW3`. If you do not have `FFTW3` installed and cannot install it easily, compile with
```
make MYFFT=-DSTANDALONE
```
instead. This will use an integrated implementation for the FFTs. They are slower than those in `FFTW3`, but typically the difference is not important.

## Usage

A detailed description of the options and results is displayed on (see last section for output):
```
comp-avg -h
```

There are three ways in which `comp-avg` can be used. In all cases it accepts any number of numeric values (scientific notation is supported) separated by arbitrary white space characters (unless otherwise specified with the `-d` option). In a file with several columns you can select one of them using `-f`. See included `foo.txt` and `bar.txt` files containing random numbers for examples.

1. Pipe the time series into `comp-avg`, e.g.
```
echo {1..1000} | comp-avg -v
# lines   mean    std. dev    error      t_int      t_max    t_int_err    err_err
# 1000    500.5   288.819     130.657    102.427    115      36.8677      23.5144
```
2. Analyse a single file
```
comp-avg foo.txt
# 300     -0.0349672      0.98703 0.0624074       0.601662        2       0.0965079       0.00500514
```
3. Compare two files (e.g. the second columns of both)
```
comp-avg -vf2 foo.txt bar.txt
# lines   mean         std. dev    error        t_int       t_max    t_int_err    err_err       file
# 100     -0.104779    1.01967     0.10424      0.527818    1        0.117064     0.0115596     foo.txt
# 200     0.0550419    0.946643    0.0669586    0.502825    1        0.0811432    0.00540271    bar.txt
# 1.28 sigma deviation.
```

### Why an FFT?
The autocorrelation function required for a reliable estimate of the error on the average via the integrated autocorrelation time can be estimated time-slice by time-slice (done by `comp-avg -n`). If autocorrelations are long and $O(n)$ time slices are needed, this leads to a runtime in $O(n^2)$. This is where the Fast Fourier Transformation (FFT) comes into play. It allows to calculate the entire autocorrelation function in $O(n \log(n))$ runtime.

The FFT-based autocorrelation calculation is the default because it guarantees runtime in $O(n \log(n))$. For time series with short autocorrelation this might not be the best option and `comp-avg -n` could be faster, ideally in $O(n)$.

Note also that the FFT assumes periodic boundary conditions. In most cases this introduces a very negligible error, but technically it is not correct. Of course, `comp-avg -n` does not have this problem.

## Reference

All the algorithms used here are explained in:

U. Wolff, “Monte Carlo errors with less errors”, [Computer Physics Communications 156, 143–153 (2004)](https://www.sciencedirect.com/science/article/pii/S0010465503004673).

## Stable Releases

`v2.1.1` first version made publicly available with reasonably comprehensive documentation.\

## Help Message

```
Calculate most important statistical information from one or more time-series:
 n, mu, sd, err, t_int, t_max, t_int_err, err_err

 n:             No. of lines
 mu:            mean
 sd:            standard deviation
 err:           standard error
 t_int:         integrated autocorrelation time
 t_max:         stopping time for summation of autocorrelation
 t_int_err:     error on tau_int
 err_err:       error on the error

 The time series can be provided as a file or piped from STDIN.
 If two files are provided, results are calculated for both and the relative deviation is given.

 -e: using only every ...-th line
 -s: summing autocovariance only up to first zero-crossing, default is Ulli Wolff method
 -n: no FFT, calculate covariance time slice by time slice instead
 -f: use ...-th field/column only
 -d: separate fields by delimiter ..., default is tab
 -p: output with machine precision, default is 5 digits
 -t: number of observables, input has to contain one column of each observable, n replaced by index (0..) of observable
 -y: symmetrise, only if '-t' option is used for correlator type data
 -v: verbose, display header line
```
