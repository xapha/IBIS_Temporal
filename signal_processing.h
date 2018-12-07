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
#include <utils.h>

#define DSP_FILTER_LOWPASS      0
#define DSP_FILTER_HIGHPASS     1
#define DSP_FILTER_BANDPASS     2
#define DSP_FILTER_BANDSTOP     3

#define GSL_FFT_FORWARD 1
#define GSL_FFT_REVERSE 2

#define FS 30

class Signal_processing
{
public:
    Signal_processing( int MaxSP, int size_signal );
    ~Signal_processing();

    void add_frame( int* parent, float* c_1, float* c_2, float* c_3, int nb_SP );
    const float* get_SNR();
    int get_HR() { return HR; }
    void process();

    float* get_C1() { return buff_signals_c1; }
    float* get_C2() { return buff_signals_c2; }
    float* get_C3() { return buff_signals_c3; }
    int* get_color() { return color; }

private:
    int index_processed;
    int size_signals;
    int max_SP;
    int nb_sp;
    int index_circular;
    int HR;
    float* buff_signals;
    float* buff_signals_c1;
    float* buff_signals_c2;
    float* buff_signals_c3;
    float* circular_data_c1;
    float* circular_data_c2;
    float* circular_data_c3;
    int* circular_parent;
    float* circular_SNR;
    float* SNR;
    float* mean_HR;
    double* buff_HR;
    int* HR_final;
    double* circular_fundamental;
    float* fundamental;
    float* signal;
    float* visu_signal;
    float* output;
    double* fft_data;
    double* fft_buff;
    double* variance;
    int* color;

    //iPBV
    double* PBV;

    bool ready;
    bool complete_SNR;
    int count_SNR;

    void construct_signal();

    void fft(double*, int, int);
    void filter(float *signal, int length, float fs, float *output, int order, float f1, float f2=0);// The order of the filter must be a multiple of 4.
    void filter(double *signal, int length, float fs, float *output, int order, float f1, float f2=0);// The order of the filter must be a multiple of 4.
    float compute_SNR(double* input, int n_input, int dirac_width, int nb_harmonic, int visu);

};

#endif // SIGNAL_PROCESSING_H
