#include "signal_processing.h"

Signal_processing::Signal_processing( int MaxSP, int size_signal)
{
    size_signals = size_signal;
    max_SP = MaxSP;
    index_circular = 0;
    HR = 0;
    count_SNR = 0;
    ready=false;
    complete_SNR=false;
    index_processed=0;

    buff_signals = new float[ max_SP * size_signals ];
    buff_signals_c1 = new float[ max_SP * size_signals ];
    buff_signals_c2 = new float[ max_SP * size_signals ];
    buff_signals_c3 = new float[ max_SP * size_signals ];
    circular_data_c1 = new float[ max_SP * size_signals ];
    circular_data_c2 = new float[ max_SP * size_signals ];
    circular_data_c3 = new float[ max_SP * size_signals ];
    circular_parent = new int[ max_SP * size_signals ];
    circular_SNR = new float[ max_SP * size_signals ];
    circular_fundamental = new double[ max_SP * size_signals ];
    buff_HR = new double[ max_SP * size_signals ];
    fundamental = new float[ max_SP ];
    SNR = new float[ max_SP ];
    mean_HR = new float[ max_SP ];
    HR_final = new int[ size_signals ];
    signal = new float[ size_signals ];
    visu_signal = new float[ size_signals ];
    output = new float[ size_signals ];
    fft_data = new double[ size_signals ];
    fft_buff = new double[ size_signals * max_SP ];
    variance = new double[ max_SP ];

    memset(mean_HR, 0, sizeof(float) * max_SP);
    memset(buff_signals, 0, sizeof(float) * max_SP * size_signals );
    memset(buff_signals_c1, 0, sizeof(float) * max_SP * size_signals );
    memset(buff_signals_c2, 0, sizeof(float) * max_SP * size_signals );
    memset(buff_signals_c3, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_data_c1, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_data_c2, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_data_c3, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_SNR, 0, sizeof(float) * max_SP * size_signals );
    memset(circular_parent, 0, sizeof(int) * max_SP * size_signals );
    memset(HR_final, 0, sizeof(int) * size_signals );
    memset(buff_HR, 0, sizeof(int) * max_SP * size_signals );
    memset(circular_fundamental, 0, sizeof(double) * max_SP * size_signals);

    PBV = new double[3*max_SP];
    color = new int[ max_SP ];

}

Signal_processing::~Signal_processing() {
    delete[] buff_signals;
    delete[] buff_signals_c1;
    delete[] buff_signals_c2;
    delete[] buff_signals_c3;
    delete[] circular_data_c1;
    delete[] circular_data_c2;
    delete[] circular_data_c3;
    delete[] circular_parent;
    delete[] circular_SNR;
    delete[] circular_fundamental;
    delete[] fundamental;
    delete[] SNR;
    delete[] buff_HR;
    delete[] HR_final;
    delete[] mean_HR;

    delete[] signal;
    delete[] visu_signal;
    delete[] output;
    delete[] fft_data;
    delete[] fft_buff;
    delete[] variance;
    delete[] PBV;
    delete[] color;

}

void Signal_processing::process() {
    // construct signals based on inheritance table
    if( ready ) {
        index_processed++;

        construct_signal();

        int temp_index = index_circular - 1;
        if( temp_index < 0 )
            temp_index += size_signals;

        double max_variance=0;

        /*for( int i=0; i<nb_sp; i++ )
            fft_data[ i ] = double(SNR[ i ]);*/
        int visu_snr = -1; //gsl_stats_max_index( fft_data, 1, nb_sp );

        for( int i=0; i<nb_sp; i++ ) {

            // succesive HR estimation variance
            for( int j=0; j<size_signals; j++ ) {
                if( buff_HR[ max_SP * j + i ] > 40 && buff_HR[ max_SP * j + i ] < 140 )
                    fft_data[ j ] = buff_HR[ max_SP * j + i ];
                else
                    fft_data[ j ] = 0;

            }

            bool stat = false;
            for( int j=0; j<size_signals; j++ ) {
                if( fft_data[j] > 0  ) {
                    stat = true;
                    break;

                }

            }

            if( stat )
                variance[ i ] = gsl_stats_variance( fft_data, 1, size_signals );
            else
                variance[ i ] = 0;

            if( variance[ i ] > max_variance )
                max_variance = variance[ i ];

            for( int j=0; j<size_signals; j++ ) {
                signal[ j ] = buff_signals[ max_SP * j + i ];

            }

            signal[ 0 ] = 0.f;
            signal[ size_signals-1 ] = 0.f;

            // filter signal
            //filter( signal, size_signals, FS, output, 4, 0.4f, 4.f);

            // compute fft
            for( int j=0; j<size_signals; j++ )
                fft_data[j] = double( signal[j] );

            fft(fft_data, size_signals, GSL_FFT_FORWARD);
            for( int j=0; j<size_signals; j++ ) {
                if( j > 0.8 * 2 / float( FS / float( size_signals ) ) && j < 4 * 2 / float( FS / float( size_signals ) ) )
                    fft_data[j] *= fft_data[j];
                else
                    fft_data[j] = 0;

            }

            for( int j=0; j<size_signals; j++ )
                fft_buff[ nb_sp * j + i ] = fft_data[ j ];

            // get fundamental freq
            circular_fundamental[ max_SP * temp_index + i ] = double(gsl_stats_max_index( fft_data, 1, size_signals ));

            // get SNR value
            float d_width = 0.5 * 2 / float( FS / float( size_signals ) );

            if( i == visu_snr )
                circular_SNR[ max_SP * temp_index + i ] = compute_SNR( fft_data, size_signals, int( round( d_width ) ), 2, 1 );
            else
                circular_SNR[ max_SP * temp_index + i ] = compute_SNR( fft_data, size_signals, int( round( d_width ) ), 2, 0 );

            SNR[ i ] += circular_SNR[ max_SP * temp_index + i ];
            mean_HR[ i ] += circular_fundamental[ max_SP * temp_index + i ];

        }

        if( count_SNR < size_signals )
            count_SNR++;

        for( int i=0; i<nb_sp; i++ ) {
            if( variance[ i ] > 0 )
                variance[ i ] = 1 - variance[ i ]/max_variance;
            else
                variance[ i ] = 0;

        }

        /*float max_SNR=0;
        float weight[nb_sp] = {0.f};
        float sum_weight=0.f;*/
        for( int i=0; i<nb_sp; i++ ) {
            SNR[i] = SNR[i] / count_SNR;
            mean_HR[ i ] /= count_SNR;

            //weight[i] = pow( 10.0, double(SNR[i])/6 );
            //sum_weight += weight[i];

        }

        /*memset( fft_data, 0, sizeof(double)*size_signals );

        for( int i=0; i<nb_sp; i++ ) {
            for( int j=0; j<size_signals; j++ ) {
                if( SNR[i] > 0 ) {
                    fft_data[j] += SNR[i] * fft_buff[ nb_sp * j + i ];

                }

            }

        }

        CvPlot::clear("fft");
        CvPlot::plot( "fft", fft_data, size_signals );*/

        //HR = round( float(gsl_stats_max_index( fft_data, 1, size_signals )) * 30.f * float( FS / float( size_signals ) ) );

        for( int j=0; j<size_signals; j++ ) {
            visu_signal[ j ] = buff_signals[ max_SP * j + 75 ];

        }
        visu_signal[0] = 0;

        memset( HR_final, 0, sizeof( int ) * size_signals );
        int count_sp=0;
        for( int i=0; i<nb_sp; i++ ) {
            if( SNR[i] > 0 ) {
                HR_final[ int(round(mean_HR[ i ])) ] += 6*int( round( double( SNR[i] ) ) );
                count_sp++;

            }

        }

        // PCA
        /*cv::Mat A( count_sp, 4, CV_32F );
        cv::Mat C( count_sp, count_sp, CV_32F );
        cv::Mat u( count_sp, 1, CV_32F );
        cv::Mat v( count_sp, count_sp, CV_32F );
        cv::Mat w( count_sp, count_sp, CV_32F );

        int index_sp=0;
        int corresp[ count_sp ];
        for( int i=0; i< nb_sp; i++ ) {
            if( SNR[i] > 0 ) {
                A.at<float>( index_sp, 0 ) = SNR[i] * ( circular_data_c1[ max_SP * temp_index + i ] / 255 );
                A.at<float>( index_sp, 1 ) = SNR[i] * ( circular_data_c2[ max_SP * temp_index + i ] / 255 );
                A.at<float>( index_sp, 2 ) = SNR[i] * ( circular_data_c3[ max_SP * temp_index + i ] / 255 );
                A.at<float>( index_sp, 3 ) = SNR[i] * ( mean_HR[ i ] / size_signals );

                corresp[ index_sp ] = i;

                if( index_sp >= count_sp )
                    printf("!!!!!!!!!\n");

                index_sp++;

            }

        }

        C = A * A.t();
        cv::SVD::compute( C, u, v, w );
        //cv::eigen( C, u, v );

        CvPlot::clear("SVD value");
        int ii = 0;
        memset(color, 0, sizeof(int) * max_SP);

        int count_signature = 0;
        int signatures[ count_sp ];
        for( int j=1; j<count_sp; j++ ) {

            bool pos=true;
            for( int k=1; k<count_sp; k++ ) {
                if( v.at<float>( k,j ) < 0 ) {
                    pos = false;
                    break;
                }
            }

            if( pos ) {
                signatures[count_signature] = j;
                count_signature++;
            }
            else
                break;

        }

        cv::Mat ref_S( count_signature, 4, CV_32F );

        float color_dist[count_sp] = {0.f};
        for( int k=0; k<count_signature; k++ ) {
            cv::Mat S( 1, 4, CV_32F );
            cv::Mat u_0( count_sp, 1, CV_32F );
            cv::Mat out_0( 1, count_sp, CV_32F );
            for( int i=0; i<count_sp; i++ )
                u_0.at<float>( i, 0 ) = v.at<float>(i,signatures[k]);

            S = u_0.t() * A;
            float norm_s = float( sqrt( double( S.at<float>(0, 0) * S.at<float>(0, 0) +
                                                S.at<float>(0, 1) * S.at<float>(0, 1) +
                                                S.at<float>(0, 2) * S.at<float>(0, 2) +
                                                S.at<float>(0, 3) * S.at<float>(0, 3) ) ) );
            S = 1/norm_s * S;

            ref_S.at<float>(ii, 0) = S.at<float>(0, 0);
            ref_S.at<float>(ii, 1) = S.at<float>(0, 1);
            ref_S.at<float>(ii, 2) = S.at<float>(0, 2);
            ref_S.at<float>(ii, 3) = S.at<float>(0, 3);

            out_0 = S * A.t();
            norm_s = 0.f;
            for( int i=0; i<count_sp; i++ )
                norm_s += out_0.at<float>(0, i) * out_0.at<float>(0, i);

            norm_s = float( sqrt( double( norm_s ) ) );
            out_0 /= norm_s;

            float out_hr = 0.f;
            float sum_hr = 0.f;
            for( int i=0; i<count_sp; i++ ) {
                signal[ i ] = out_0.at<float>(0,i);
                if( signal[ i ] < 0 )
                    signal[ i ] = -signal[ i ];

                float weight = out_0.at<float>(0,i);

                if( weight < 0 )
                    weight = -weight;

                out_hr += weight * mean_HR[ corresp[ i ] ];
                sum_hr += weight;

                if( color_dist[ i ] < weight ) {
                    color_dist[ i ] = weight;
                    color[ corresp[i] ] == ii;

                }

            }

            CvPlot::plot( "SVD value", signal, count_sp );

            printf(" ---- > Signature %i : (%f, %f, %f, %f) \n",
                   ii,
                   S.at<float>(0,0),
                   S.at<float>(0,1),
                   S.at<float>(0,2),
                   S.at<float>(0,3) );
            printf(" ---- -> HR : %f bpm\n", ( out_hr / sum_hr ) * 30.f * float( FS / float( size_signals ) ) );

            ii++;

        }*/

        // check correlation between signature and initial vectors
        /*int nb_s[ count_signature ] = {0};
        for( int j=0; j<count_sp; j++ ) {
            float c_dist = 0.f;
            int best_color;

            for( int i=0; i<count_signature; i++ ) {*/
                /*float correlation = ( ref_S.at<float>( i,0 ) - A.at<float>( j,0 ) ) * ( ref_S.at<float>( i,0 ) - A.at<float>( j,0 ) ) +
                                    ( ref_S.at<float>( i,1 ) - A.at<float>( j,1 ) ) * ( ref_S.at<float>( i,1 ) - A.at<float>( j,1 ) ) +
                                    ( ref_S.at<float>( i,2 ) - A.at<float>( j,2 ) ) * ( ref_S.at<float>( i,2 ) - A.at<float>( j,2 ) ) +
                                    ( ref_S.at<float>( i,3 ) - A.at<float>( j,3 ) ) * ( ref_S.at<float>( i,3 ) - A.at<float>( j,3 ) );
*/
                /*cv::Mat ref( 1, 4, CV_32F );
                cv::Mat te( 1, 4, CV_32F );

                ref.at<float>( 0,0 ) = ref_S.at<float>( i,0 );
                ref.at<float>( 0,1 ) = ref_S.at<float>( i,1 );
                ref.at<float>( 0,2 ) = ref_S.at<float>( i,2 );
                ref.at<float>( 0,3 ) = ref_S.at<float>( i,3 );

                te.at<float>( 0,0 ) = A.at<float>( j,0 );
                te.at<float>( 0,1 ) = A.at<float>( j,1 );
                te.at<float>( 0,2 ) = A.at<float>( j,2 );
                te.at<float>( 0,3 ) = A.at<float>( j,3 );

                float correlation = float( ref.dot( te ) );

                if( i == 0 || correlation < c_dist ) {
                    best_color = i;
                    c_dist = correlation;

                }

            }

            color[ corresp[ j ] ] = best_color;
            nb_s[ best_color ]++;

        }

        cv::waitKey(1);*/

        // detect independent HR
        /*int tentative_HR[20] = {0};
        int nb_tentative = 0;
        int energic=0;
        int energic_id=0;
        for( int i=0; i<size_signals; i++ ) {
            if( HR_final[i] > 0 ) {
                tentative_HR[ nb_tentative ] = HR_final[i];
                nb_tentative++;

                if( HR_final[i] > energic ) {
                    energic = HR_final[i];
                    energic_id = i;
                }

            }

        }
        HR = int( round( energic_id * 30.f * float( FS / float( size_signals ) ) ) );
        printf(" tentative HR %i, most energic: %i bpm\n", nb_tentative, HR );

        int repart[ 200 * 200 ] = {0};
        for( int i=0; i<nb_sp; i++ ) {
            if( int(round(mean_HR[ i ])) > 0 && int(round(mean_HR[ i ])) < 200 && int( round( double( SNR[i] ) ) ) > 0 ) {
                int temp = int(round(mean_HR[ i ]));
                //repart[ int(round(mean_HR[ i ])) + 200 * ( 200 - int( round( double( 10*SNR[i] ) ) ) ) ] = 1;
                repart[ int(round(mean_HR[ i ])) + 200 * ( 200 - int( round( double( min<float>(200, 10*SNR[i]) ) ) ) ) ] = 1;

            }
            else if( int(round(mean_HR[ i ])) < 0 )
                printf("invalid freq measure : %i\n", int(round(mean_HR[ i ])));

        }

        imagesc( "repart freq", repart, 200, 200 );

        CvPlot::clear("histo_freq");
        CvPlot::plot( "histo_freq", HR_final, size_signals );

        CvPlot::clear("1-Variance");
        CvPlot::plot( "1-Variance", variance, nb_sp );

        CvPlot::clear("SNR");
        CvPlot::plot( "SNR", SNR, nb_sp );

        for( int j=1; j<size_signals-2; j++ )
            HR_final[j] = HR_final[j+1];

        HR_final[ size_signals-2 ] = HR;
        HR_final[ 0 ] = 40;
        HR_final[ size_signals-1 ] = 140;
        CvPlot::clear("HR");
        CvPlot::plot( "HR", HR_final, size_signals );

        cv::waitKey( 1 );*/

    }

}

float* Signal_processing::get_SNR() {
    return SNR;

}


void Signal_processing::add_frame( int* parent, float* c_1, float* c_2, float* c_3, int nb_SP ) {
    nb_sp = nb_SP;

    for( int i=0; i<nb_sp; i++ ) {
        if( parent[ i ] == i ) {
            circular_data_c1[ max_SP * index_circular + i ] = c_1[ i ];
            circular_data_c2[ max_SP * index_circular + i ] = c_2[ i ];
            circular_data_c3[ max_SP * index_circular + i ] = c_3[ i ];

        }

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
    float c1, c2, c3;

//#pragma omp parallel for num_threads(4)
    for( int j=0; j<nb_sp; j++ ) {
         parent_index = j;

        for( int i=0; i<size_signals-1; i++ ) {
            temp_index = index_circular - 1 - i;
            if( temp_index < 0 )
                temp_index += size_signals;

            // construct signals
            c1 = circular_data_c1[ max_SP * temp_index + parent_index ];
            c2 = circular_data_c2[ max_SP * temp_index + parent_index ];
            c3 = circular_data_c3[ max_SP * temp_index + parent_index ];

            buff_signals_c1[ max_SP * ( size_signals - 1 - i ) + j ] = c1;
            buff_signals_c2[ max_SP * ( size_signals - 1 - i ) + j ] = c2;
            buff_signals_c3[ max_SP * ( size_signals - 1 - i ) + j ] = c3;

            parent_index = circular_parent[ max_SP * temp_index + parent_index ];

            buff_HR[ max_SP * ( size_signals - 1 - i ) + j ] = circular_fundamental[ max_SP * temp_index + parent_index ];

            // construct SNR_history
            if( circular_SNR[ max_SP * temp_index + parent_index ] != 0 ) {
                SNR[ j ] += circular_SNR[ max_SP * temp_index + parent_index ];
                mean_HR[ j ] += circular_fundamental[ max_SP * temp_index + parent_index ];

            }

        }

        // apply rPPG algo :
        double Xf[size_signals];
        double Yf[size_signals];
        double alpha;

        // C1, C2, C3
        double C1[size_signals];
        double C2[size_signals];
        double C3[size_signals];

        for(int i=0; i<size_signals; i++) {
            C1[i] = buff_signals_c1[ max_SP * i + j ];
            C2[i] = buff_signals_c2[ max_SP * i + j ];
            C3[i] = buff_signals_c3[ max_SP * i + j ];

        }

        // normalization
        double mean_C1 = gsl_stats_mean( C1, 1, size_signals );
        double sd_C1 = gsl_stats_sd( C1, 1, size_signals );

        double mean_C2 = gsl_stats_mean( C2, 1, size_signals );
        double sd_C2 = gsl_stats_sd( C2, 1, size_signals );

        double mean_C3 = gsl_stats_mean( C3, 1, size_signals );
        double sd_C3 = gsl_stats_sd( C3, 1, size_signals );

        for(int i=0; i<size_signals; i++) {
            C1[i] = ( C1[i] - mean_C1 ) / sd_C1;
            C2[i] = ( C2[i] - mean_C2 ) / sd_C2;
            C3[i] = ( C3[i] - mean_C3 ) / sd_C3;

        }

        filter( C1, size_signals, FS, output, 4, 0.4f, 4.f);
        for(int i=0; i<size_signals; i++)
            C1[i] = double(output[i]);

        filter( C2, size_signals, FS, output, 4, 0.4f, 4.f);
        for(int i=0; i<size_signals; i++)
            C2[i] = double(output[i]);

        filter( C3, size_signals, FS, output, 4, 0.4f, 4.f);
        for(int i=0; i<size_signals; i++)
            C3[i] = double(output[i]);

        //-- -- -- -- -- CHROM

        //Xs && Ys
        for(int i=0; i<size_signals; i++) {
            Xf[i] = 3*C1[i]-2*C2[i];
            Yf[i] = 1.5*C1[i]+C2[i]-1.5*C3[i];

        }

        double std_xf = gsl_stats_sd(Xf, 1, size_signals);
        double std_yf = gsl_stats_sd(Yf, 1, size_signals);

        alpha = std_xf / std_yf;

        for( int i=0; i<size_signals-1; i++ ) {
            buff_signals[ max_SP * i + j ] = float( Xf[i] - alpha*Yf[i] );

        }

        //-- -- -- -- -- iPBV
//        if( index_processed == 1 ) {
//            //Xs && Ys
//            for(int i=0; i<size_signals; i++) {
//                Xf[i] = 3*C1[i]-2*C2[i];
//                Yf[i] = 1.5*C1[i]+C2[i]-1.5*C3[i];

//            }

//            double std_xf = gsl_stats_sd(Xf, 1, size_signals);
//            double std_yf = gsl_stats_sd(Yf, 1, size_signals);

//            alpha = std_xf / std_yf;

//            cv::Mat_<double> ref( 1, size_signals );
//            cv::Mat_<double> sig( 3, size_signals ); // 3 * S

//            cv::Mat_<double> buff_3_3( 3, 3 ); // 3 * 3
//            cv::Mat_<double> buff_1_3( 1, 3 ); // 1 * 3

//            cv::Mat_<double> Pbv( 1, 3 );

////            double std_xf = gsl_stats_sd(C2, 1, size_signals);
////            double std_yf = gsl_stats_sd(C3, 1, size_signals);

////            alpha = std_xf / std_yf;

//            for( int i=0; i<size_signals-1; i++ ) {
//                ref(0, i) = -( Xf[i] - alpha*Yf[i] );
//                //ref(0, i) = C2[ i ] - alpha*C3[ i ];

//                sig(0, i) = C1[i];
//                sig(1, i) = C2[i];
//                sig(2, i) = C3[i];

//            }

//            // compute first PBV
//            Pbv = ref * sig.t();
//            buff_3_3 = ( sig * sig.t() ).inv();
//            buff_1_3 = Pbv * buff_3_3;

//            ref = buff_1_3 * sig;
//            for( int i=0; i<size_signals-1; i++ ) {
//                float a = float( ref(0, i) );
//                buff_signals[ max_SP * i + j ] = a;

//            }

//            PBV[max_SP*0 + j] = Pbv(0, 0);
//            PBV[max_SP*1 + j] = Pbv(0, 1);
//            PBV[max_SP*2 + j] = Pbv(0, 2);

//        }
//        else {
//            cv::Mat_<double> ref( 1, size_signals );
//            cv::Mat_<double> sig( 3, size_signals ); // 3 * S

//            cv::Mat_<double> buff_3_3( 3, 3 ); // 3 * 3
//            cv::Mat_<double> buff_1_3( 1, 3 ); // 1 * 3

//            cv::Mat_<double> Pbv( 1, 3 );

//            for( int i=0; i<size_signals-1; i++ ) {
//                sig(0, i) = C1[i];
//                sig(1, i) = C2[i];
//                sig(2, i) = C3[i];

//            }

//            // compute first PBV
//            Pbv(0, 0) = PBV[max_SP*0 + j];
//            Pbv(0, 1) = PBV[max_SP*1 + j];
//            Pbv(0, 2) = PBV[max_SP*2 + j];

//            buff_3_3 = ( sig * sig.t() ).inv();
//            buff_1_3 = Pbv * buff_3_3;

//            ref = buff_1_3 * sig;
//            for( int i=0; i<size_signals-1; i++ )
//                buff_signals[ max_SP * i + j ] = float( ref( 0, i ) );

//            Pbv = ref * sig.t();

//            double alpha = 0.99;
//            PBV[max_SP*0 + j] = alpha*PBV[max_SP*0 + j] + ( 1 - alpha )*Pbv(0, 0);
//            PBV[max_SP*1 + j] = alpha*PBV[max_SP*1 + j] + ( 1 - alpha )*Pbv(0, 1);
//            PBV[max_SP*2 + j] = alpha*PBV[max_SP*2 + j] + ( 1 - alpha )*Pbv(0, 2);

//            /*PBV[max_SP*0 + j] = Pbv(0, 0);
//            PBV[max_SP*1 + j] = Pbv(0, 1);
//            PBV[max_SP*2 + j] = Pbv(0, 2);*/

//        }

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

void Signal_processing::filter(double* signal, int length, float s, float* output, int n, float f1, float f2)
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
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + float(signal[j]);
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

    if( int_noise > 0 )
        SNR = float(10*log10(int_signal / int_noise));
    else
        SNR = 0.f;

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
