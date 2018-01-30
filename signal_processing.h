#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <math.h>
#include <float.h>
#include <vector>
#include <dirent.h>
#include <cstring>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <cvplot.h>

#define DSP_FILTER_LOWPASS      0
#define DSP_FILTER_HIGHPASS     1
#define DSP_FILTER_BANDPASS     2
#define DSP_FILTER_BANDSTOP     3

#define GSL_FFT_FORWARD 1
#define GSL_FFT_REVERSE 2

class Signal_processing
{
public:
    Signal_processing( int MaxSP, int size_signal );
    ~Signal_processing();

    void add_frame( int* parent, float* c_1, float* c_2, float* c_3, int nb_SP );
    float* get_SNR() { return SNR; }
    void process();

private:
    int size_signals;
    int max_SP;
    int nb_sp;
    int index_circular;
    float* buff_signals;
    float* circular_data;
    int* circular_parent;
    float* SNR;
    bool ready;

    void construct_signal();

    void fft(double*, int, int);
    void filter(float *signal, int length, float fs, float *output, int order, float f1, float f2=0);// The order of the filter must be a multiple of 4.
    void _CHROM(double* R_trace, double* G_trace, double* B_trace, double* output, int fft_size);
    float compute_SNR(double* input, int n_input, int dirac_width, int nb_harmonic, int visu);

};

#endif // SIGNAL_PROCESSING_H
