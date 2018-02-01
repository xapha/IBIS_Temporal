#include "signal_processing.h"

Signal_processing::Signal_processing( int MaxSP, int size_signal)
{
    size_signals = size_signal;
    max_SP = MaxSP;
    index_circular = 0;
    count_SNR = 0;
    ready=false;
    complete_SNR=false;

    buff_signals = new float[ max_SP * size_signals ];
    circular_data = new float[ max_SP * size_signals ];
    circular_parent = new int[ max_SP * size_signals ];
    circular_SNR = new float[ max_SP * size_signals ];
    circular_fundamental = new double[ max_SP * size_signals ];
    fundamental = new float[ max_SP ];
    SNR = new float[ max_SP ];

    memset(buff_signals, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_data, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_SNR, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_parent, 0, sizeof(int) * max_SP * size_signals );

}

Signal_processing::~Signal_processing() {
    delete[] buff_signals;
    delete[] circular_data;
    delete[] circular_parent;
    delete[] circular_SNR;
    delete[] circular_fundamental;
    delete[] fundamental;
    delete[] SNR;

}

void Signal_processing::process() {
    // construct signals based on inheritance table
    if( ready ) {
        construct_signal();
        int temp_index = index_circular - 1;
        if( temp_index < 0 )
            temp_index += size_signals;

        float* signal = new float[ size_signals ];
        float* output = new float[ size_signals ];
        double* fft_data = new double[ size_signals ];

        for( int i=0; i<nb_sp; i++ ) {
            float mean_val = 0;

            for( int j=0; j<size_signals; j++ ) {
                signal[ j ] = buff_signals[ max_SP * j + i ];
                mean_val += signal[ j ];

            }

            for( int j=0; j<size_signals; j++ )
                signal[ j ] -= mean_val/size_signals;

            signal[ 0 ] = 0.f;

            // filter signal
            // filter( signal, size_signals, 30.f, output, 4, 0.4f, 4.f);

            // compute fft
            for( int j=0; j<size_signals; j++ )
                fft_data[j] = double( signal[j] );

            fft(fft_data, size_signals, GSL_FFT_FORWARD);
            for( int j=0; j<size_signals; j++ )
                fft_data[j] *= fft_data[j];

            // get fundamental freq
            circular_fundamental[ max_SP * temp_index + i ] = gsl_stats_max_index( fft_data, 1, size_signals ) / 30 * 100;
            // get SNR value
            circular_SNR[ max_SP * temp_index + i ] = compute_SNR(fft_data, size_signals, 5, 2, 0);
            SNR[i] += circular_SNR[ max_SP * temp_index + i ];


        }

        if( count_SNR < size_signals )
            count_SNR++;

        for( int i=0; i<nb_sp; i++ ) {
            SNR[i] /= count_SNR;

        }

        int repart[ 200 * 150 ] = {0};
        for( int i=0; i<nb_sp; i++ ) {
            if( int(round(circular_fundamental[ max_SP * temp_index + i ])) > 0 && int(round(circular_fundamental[ max_SP * temp_index + i ])) < 2 && int( round( double( 10*SNR[i] ) ) ) > 0 ) {
                int temp = int(round(circular_fundamental[ max_SP * temp_index + i ]));
                repart[ int(round(circular_fundamental[ max_SP * temp_index + i ])) + 200 * ( 150-int( round( double( 10*SNR[i] ) ) ) ) ] = i;

            }
            else if( int(round(circular_fundamental[ max_SP * temp_index + i ])) < 0 )
                printf("invalid freq measure : %i\n", int(round(circular_fundamental[ max_SP * temp_index + i ] * 30)));

        }

        imagesc( "repart freq", repart, 100, 100 );

        CvPlot::clear("signal");
        CvPlot::plot( "signal", SNR, nb_sp );
        cv::waitKey( 1 );

        delete[] signal;
        delete[] output;
        delete[] fft_data;

    }

}

float* Signal_processing::get_SNR() {
    return SNR;

}


void Signal_processing::add_frame( int* parent, float* c_1, float* c_2, float* c_3, int nb_SP ) {
    nb_sp = nb_SP;

    for( int i=0; i<nb_sp; i++ ) {
        circular_data[ max_SP * index_circular + i ] = c_2[ i ];
        circular_parent[ max_SP * index_circular + i ] = parent[ i ];

    }

    index_circular++;
    if( index_circular >= size_signals ) {
        index_circular = 0;
        ready = true;

    }

}

void Signal_processing::construct_signal() {
    int parent_index;
    int temp_index;
    memset( SNR, 0, sizeof(float) * max_SP );

    for( int j=0; j<nb_sp; j++ ) {
         parent_index = j;

        for( int i=0; i<size_signals-1; i++ ) {
            temp_index = index_circular - 1 - i;
            if( temp_index < 0 )
                temp_index += size_signals;

            // construct signals
            buff_signals[ max_SP * ( size_signals - 1 - i ) + j ] = circular_data[ max_SP * temp_index + parent_index ];
            parent_index = circular_parent[ max_SP * temp_index + parent_index ];

            // construct SNR_history
            if( circular_SNR[ max_SP * temp_index + parent_index ] != 0 ) {
                SNR[ j ] += circular_SNR[ max_SP * temp_index + parent_index ];

            }

        }

    }

}

void Signal_processing::fft(double* data, int n, int type)
{
    gsl_fft_real_wavetable * real;
    gsl_fft_halfcomplex_wavetable * hc;
    gsl_fft_real_workspace * work;

    if(type == 1) {
        work = gsl_fft_real_workspace_alloc(n);
        real = gsl_fft_real_wavetable_alloc(n);

        gsl_fft_real_transform(data, 1, n, real, work);

        gsl_fft_real_wavetable_free(real);
        gsl_fft_real_workspace_free (work);

        //gsl_fft_real_radix2_transform(data, 1, n);

        /*for(int i=0; i<n/2; i++) {
            data[i] = pow(data[i], 2) + pow(data[n-i], 2);
        }*/
    }
    else if(type == 2) {
        work = gsl_fft_real_workspace_alloc(n);
        hc = gsl_fft_halfcomplex_wavetable_alloc (n);

        gsl_fft_halfcomplex_inverse (data, 1, n, hc, work);

        gsl_fft_halfcomplex_wavetable_free (hc);
        gsl_fft_real_workspace_free (work);
    }
    else
        printf("Wrong fft flag: < GSL_FFT_FORWARD || GSL_FFT_REVERSE >");
}

void Signal_processing::filter(float* signal, int length, float s, float* output, int n, float f1, float f2)
{
    f1 *= 2;
    f2 *= 2;
    float a = cos(M_PI*(f1+f2)/s)/cos(M_PI*(f2-f1)/s);
    float a2 = a*a;
    float b = tan(M_PI*(f2-f1)/s);
    float b2 = b*b;
    float r;

    int i,j;
    n = n/4;
    float *A = (float *)malloc(n*sizeof(float));
    float *d1 = (float *)malloc(n*sizeof(float));
    float *d2 = (float *)malloc(n*sizeof(float));
    float *d3 = (float *)malloc(n*sizeof(float));
    float *d4 = (float *)malloc(n*sizeof(float));
    float *w0 = (float *)calloc(n, sizeof(float));
    float *w1 = (float *)calloc(n, sizeof(float));
    float *w2 = (float *)calloc(n, sizeof(float));
    float *w3 = (float *)calloc(n, sizeof(float));
    float *w4 = (float *)calloc(n, sizeof(float));

    float C1 = 1.0;
    float C2 = 2.0;
    float C4 = 4.0;

    for(i=0; i<n; ++i) {
        r = sin(M_PI*(C2*i+C1)/(C4*n));
        s = b2 + C2*b*r + C1;

        A[i] = b2/s;
        d1[i] = C4*a*(C1+b*r)/s;
        d2[i] = C2*(b2-C2*a2-C1)/s;
        d3[i] = C4*a*(C1-b*r)/s;
        d4[i] = -(b2 - C2*b*r + C1)/s;
    }
    for(j=0; j<length; j++) {
        for(i=0; i<n; ++i) {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + signal[j];
            output[j] = A[i]*(w0[i] - C2*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
    }

    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
}

float Signal_processing::compute_SNR(double* input, int n_input, int dirac_width, int nb_harmonic, int visu) {
    int i, j;
    int fft_plot[n_input];
    double noise_signal[n_input] = {0};
    double pure_signal[n_input] = {0};

    double int_signal = 0;
    double int_noise = 0;

    float SNR;

    //define a model
    double fft_model[n_input] = {0};
    int demi_width_dirac = dirac_width/2;

    input[0] = 0;
    int SNR_max_index = gsl_stats_max_index(input, 1, n_input);

    for(i=0; i<n_input; i++) {
        for(j=0; j<nb_harmonic; j++) {
            if(  i > ( ( (j+1) * SNR_max_index ) - demi_width_dirac ) && i < ( ( (j+1) * SNR_max_index ) + demi_width_dirac ) ) {
                fft_model[i] = 1;

                if(j>0 && SNR_max_index > 0)
                    SNR = 0;
            }
        }
    }

    for(i=0; i<n_input; i++) {
        if( fft_model[i] > 0 ) {
            int_signal += input[i];
            pure_signal[i] = input[i];
        }
        else {
            int_noise += input[i];
            noise_signal[i] = input[i];
        }
    }

    SNR = float(10*log10(int_signal / int_noise));

    if(visu > 0) {
        CvPlot::clear("Visu fft");

        for(i=0; i<n_input; i++)
            fft_plot[i] = int(pure_signal[i] * 1000);

        CvPlot::plot("Visu fft", fft_plot, n_input);

        for(i=0; i<n_input; i++)
            fft_plot[i] = int(noise_signal[i] * 1000);

        CvPlot::plot("Visu fft", fft_plot, n_input);
    }

    return SNR;
}

void Signal_processing::_CHROM(double* R_trace, double* G_trace, double* B_trace, double* output, int fft_size)
{
    double Xs[fft_size];
    double Ys[fft_size];
    double Xf[fft_size];
    double Yf[fft_size];
    double alpha;

    //Xs && Ys
    for(int i=0; i<fft_size; i++) {
        Xf[i] = 3*R_trace[i]-2*G_trace[i];
        Yf[i] = 1.5*R_trace[i]+G_trace[i]-1.5*B_trace[i];
    }

    double std_xf = gsl_stats_sd(Xf, 1, fft_size);
    double std_yf = gsl_stats_sd(Yf, 1, fft_size);

    alpha = std_xf / std_yf;
    for(int i=0; i<fft_size; i++)
        output[i] = Xf[i] - alpha*Yf[i];
}
