/* -- Serge Bobbia : serge.bobbia@u-bourgogne.fr -- Le2i 2018
 * This work is distributed for non commercial use only,
 * it implements the IBIS method as described in the ICPR 2018 paper.
 *
 * Read the ibis.h file for options and benchmark instructions
 */

#include "ibis.h"

IBIS::IBIS(int _maxSPNum, int _compacity ) {
    labels = nullptr;
    maxSPNumber = _maxSPNum;
    compacity = _compacity;
    size = 0;

    // memory allocation
    Xseeds = new float[maxSPNumber];
    Yseeds = new float[maxSPNumber];
    lseeds = new float[maxSPNumber];
    aseeds = new float[maxSPNumber];
    bseeds = new float[maxSPNumber];

    Xseeds_init = new float[maxSPNumber];
    Yseeds_init = new float[maxSPNumber];

    Xseeds_Sum = new float[maxSPNumber];
    Yseeds_Sum = new float[maxSPNumber];
    lseeds_Sum = new float[maxSPNumber];
    aseeds_Sum = new float[maxSPNumber];
    bseeds_Sum = new float[maxSPNumber];

    countPx = new float[maxSPNumber];

    inv = new float[maxSPNumber];
    adjacent_sp = new int[size_roi*maxSPNumber];
    count_adjacent = new int[maxSPNumber];
    prev_adjacent_sp = new int[size_roi*maxSPNumber];
    prev_count_adjacent = new int[maxSPNumber];

    // temporal
    count_diff = new int[maxSPNumber];
    Xseeds_prev = new float[maxSPNumber];
    Yseeds_prev = new float[maxSPNumber];
    lseeds_prev = new float[maxSPNumber];
    aseeds_prev = new float[maxSPNumber];
    bseeds_prev = new float[maxSPNumber];
    inheritance = new int[maxSPNumber];

}

IBIS::~IBIS() {
    delete[] labels;
    delete[] inv;

    // allocated in constructor
    delete[] Xseeds;
    delete[] Yseeds;
    delete[] lseeds;
    delete[] aseeds;
    delete[] bseeds;

    delete[] Xseeds_init;
    delete[] Yseeds_init;

    delete[] Xseeds_Sum;
    delete[] Yseeds_Sum;
    delete[] lseeds_Sum;
    delete[] aseeds_Sum;
    delete[] bseeds_Sum;

    delete[] countPx;

    delete[] mask_buffer;
    delete[] mask_size;
    delete[] x_vec;
    delete[] y_vec;
    delete[] vertical_index;

    delete[] adjacent_sp;
    delete[] count_adjacent;
    delete[] prev_adjacent_sp;
    delete[] prev_count_adjacent;

    delete[] initial_repartition;
    delete[] processed;

    // temporal
    delete[] avec_bis;
    delete[] lvec_bis;
    delete[] bvec_bis;

    delete[] avec_one;
    delete[] lvec_one;
    delete[] bvec_one;

    delete[] count_diff;

    delete[] Xseeds_prev;
    delete[] Yseeds_prev;
    delete[] lseeds_prev;
    delete[] aseeds_prev;
    delete[] bseeds_prev;

    delete[] inheritance;
    delete[] updated_px;
    delete[] elligible;

}

void IBIS::initSeeds() {
    int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip, yerrperstrip;
    int xoff, yoff;
    int x, y;
    int xe, ye, xe_1, ye_1;
    int start_y, final_y, start_x, final_x;
    int* tmp_adjacent_sp = new int[size_roi*maxSPNumber];
    int* tmp_count_adjacent = new int[maxSPNumber];

    xstrips = width / SPTypicalLength;
    ystrips = height / SPTypicalLength;

    xerr = width - SPTypicalLength * xstrips;
    yerr = height - SPTypicalLength * ystrips;

    xerrperstrip = (double)xerr / xstrips;
    yerrperstrip = (double)yerr / ystrips;

    xoff = SPTypicalLength / 2;
    yoff = SPTypicalLength / 2;

    n = 0;
    for (y = 0; y < ystrips; y++)
    {
        ye = (int)(y*yerrperstrip);
        ye_1 = (int)((y+1)*yerrperstrip);

        int seedy = y * SPTypicalLength + yoff + ye;

        if( y == 0 ) {
            start_y = 0;
            final_y = SPTypicalLength + ye_1;

        }
        else {
            start_y = y * SPTypicalLength + ye;
            final_y = ( (y + 1) * SPTypicalLength + ye_1 >= height ) ? height-1 : (y + 1) * SPTypicalLength + ye_1;

        }

        for (x = 0; x < xstrips; x++)
        {
            int seedx;
            xe = (int)(x*xerrperstrip);
            xe_1 = (int)((x+1)*xerrperstrip);
            seedx = x * SPTypicalLength + xoff + xe;

            if( x == 0 ) {
                start_x = 0;
                final_x = SPTypicalLength + xe_1;

            }
            else {
                start_x = x * SPTypicalLength + xe;

                final_x = ( (x + 1) * SPTypicalLength + xe_1 > width ) ? width : (x + 1) * SPTypicalLength + xe_1;

            }

            Xseeds_init[n] = (float) seedx;
            Yseeds_init[n] = (float) seedy;

            // fill line by line
            for( int index_y=start_y; index_y<=final_y; index_y++ ) {
                std::fill( initial_repartition + index_y*width + start_x, initial_repartition + index_y*width + final_x, n );

            }

            // list adjacents seeds
            tmp_count_adjacent[n] = 0;
            for( int roi_y=-(sqrt(size_roi) - 1)/2; roi_y <= (sqrt(size_roi) - 1)/2; roi_y++ ) {
                for( int roi_x=-(sqrt(size_roi) - 1)/2; roi_x <= (sqrt(size_roi) - 1)/2; roi_x++ ) {
                    if( !( y + roi_y < 0 || y + roi_y >= ystrips || x + roi_x < 0 || x + roi_x >= xstrips ) ) {
                        tmp_adjacent_sp[size_roi*n+tmp_count_adjacent[n]] = n + roi_y*xstrips + roi_x;
                        tmp_count_adjacent[n]++;

                    }

                }

            }

            n++;
        }
    }
    SPNumber = n;

    for(int i=0; i<SPNumber; i++) {
        count_adjacent[i] = 0;

        for(int j=0; j<tmp_count_adjacent[i]; j++) {
            if( tmp_adjacent_sp[size_roi*i+j] >= 0 && tmp_adjacent_sp[size_roi*i+j] < SPNumber ) {
                adjacent_sp[size_roi*i+count_adjacent[i]] = tmp_adjacent_sp[size_roi*i+j];
                count_adjacent[i]++;

            }

        }

    }

    delete[] tmp_adjacent_sp;
    delete[] tmp_count_adjacent;

}

void IBIS::mean_seeds() {

    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] > 0 ) {
            Xseeds[ i ] = Xseeds_Sum[ i ] / countPx[ i ];
            Yseeds[ i ] = Yseeds_Sum[ i ] / countPx[ i ];
            lseeds[ i ] = lseeds_Sum[ i ] / countPx[ i ];
            aseeds[ i ] = aseeds_Sum[ i ] / countPx[ i ];
            bseeds[ i ] = bseeds_Sum[ i ] / countPx[ i ];

        }

    }

}

float IBIS::get_complexity() {
    int count_px_processed = 0;
    for( int i=0; i<size; i++ ) {
        if( processed[ i ] )
            count_px_processed++;

    }

    return float(count_px_processed) / float(size);

}

double IBIS::now_ms(void) {
    double milliseconds_since_epoch = (double) (std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    return milliseconds_since_epoch;

}

// return the current actual SP number (the method may reduce the actual Sp number).
void IBIS::enforceConnectivity() {
    //local var
    int label = 0;
    int i, j, k;
    int n, c, count;
    int x, y;
    int ind;
    int oindex, adjlabel;
    int nindex;
    const int dx4[4] = { -1,  0,  1,  0 };
    const int dy4[4] = { 0, -1,  0,  1 };
    std::fill( nlabels, nlabels + size, -1 );

    oindex = 0;
    adjlabel = 0;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            if (nlabels[oindex] < 0)
            {
                label = labels[ oindex ];
                nlabels[oindex] = label;// !! labels[oindex] --> label

                x_vec[0] = k;
                y_vec[0] = j;

                for (n = 0; n < 4; n++)
                {
                    x = x_vec[0] + dx4[n];
                    y = y_vec[0] + dy4[n];

                    if ((x >= 0 && x < width) && (y >= 0 && y < height))
                    {
                        nindex = y*width + x;

                        if (nlabels[nindex] >= 0)
                            adjlabel = nlabels[nindex];
                    }
                }

                count = 1;
                for (c = 0; c < count; c++)
                {
                    for (n = 0; n < 4; n++)
                    {
                        x = x_vec[c] + dx4[n];
                        y = y_vec[c] + dy4[n];
                        if ((x >= 0 && x < width) && (y >= 0 && y < height))
                        {
                            nindex = y*width + x;

                            if (nlabels[nindex] < 0 && labels[oindex] == labels[nindex])
                            {
                                x_vec[count] = x;
                                y_vec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }

                if (count <= minSPSizeThreshold)
                {
                    for (c = 0; c < count; c++)
                    {
                        ind = y_vec[c] * width + x_vec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }

    for (i = 0; i < size; i++)
        labels[i] = nlabels[i];

}

void IBIS::init() {
    int ii, step;

    //store mean distance between 2 seeds info
    SPTypicalLength = (int)(std::sqrt((float)(size) / (float)(maxSPNumber))) + 1;

    // compacity weight
    invwt = (float)SPTypicalLength / compacity;
    invwt = 1.0f / (invwt * invwt);

    // inferior limit for superpixels size
    minSPSizeThreshold = (size / maxSPNumber) / 4;

    // set the top mask size
    index_mask = 1;
    while( SPTypicalLength > ( pow( 2.0, index_mask ) + 1 ) )
        index_mask++;
    index_mask--;

    // set mask size
    mask_size = new int[ index_mask ];
    for( int k=0; k<index_mask; k++ )
        mask_size[ k ] = pow( 2.0, k+1 ) + 1;

    // center of the first block
    start_xy = mask_size[ index_mask - 2 ] - 1;

    // enforce connectivity buffer
    x_vec = new int[ size ];
    y_vec = new int[ size ];

    // visu the processed pixels
    processed = new int[ size ];
    std::fill( processed, processed + size, 0 );

    // precalculate vertical indexes
    vertical_index = new int[ height + mask_size[ index_mask - 1 ] ];
    for( int k=0; k < height + mask_size[ index_mask - 1 ]; k++ ) {
        vertical_index[ k ] = k*width;
    }

    // output labels buffer
    labels = new int[size];

    // repartition of pixels at start
    initial_repartition = new int[size];

    // limit for masks definittion
    y_limit = height + mask_size[ index_mask-1 ];
    x_limit = width + mask_size[ index_mask-1 ];

    // create MASKs
    step = mask_size[ index_mask-1 ];
    ii=0;
    int col_count = 0;
    int row_count = 0;
    for( int y=start_xy; y<y_limit; y+=step ) {
        row_count++;
        for( int x=start_xy; x<x_limit; x+=step ) {
            if( row_count == 1 )
                col_count++;

            ii++;

        }

    }

    count_mask = ii;
    mask_buffer = new MASK[ count_mask ];
    for( int y=start_xy, yy=0; y<y_limit; y+=step, yy++ ) {
        for( int x=start_xy, xx=0; x<x_limit; x+=step, xx++ ) {
            int adjacent[9] = {-1};
            ii = yy*col_count + xx;

            int jj=0;
            for( int i=-1; i<=1; i++ ) {
                for( int j=-1; j<=1; j++ ) {
                    int pot = ii + i*col_count + j;
                    adjacent[ jj ] =  ( !( yy + i < 0 || yy + i >= row_count || xx + j < 0 || xx + j >= col_count ) ) ? pot : -1;
                    jj++;

                }

            }

            mask_buffer[ ii ].init( pow(2.0, index_mask+1), y, x, index_mask-1, this, adjacent );

        }

    }

    // TEMPORAL init
    img_bis = false;

    // image lab buffer
    avec_bis = new float[size];
    bvec_bis = new float[size];
    lvec_bis = new float[size];

    avec_one = new float[size];
    bvec_one = new float[size];
    lvec_one = new float[size];

    updated_px = new bool[size];
    elligible = new int[size];

    nlabels = new int[ size ];

    std::fill( labels, labels + size, -1 );
    std::fill( elligible, elligible + size, 0 );

}

void IBIS::reset() {
    int index_xy;

    st4 = 0;
    st2 = 0;
    st3 = 0;

    std::fill( countPx, countPx + maxSPNumber, 0 );


#if !fixed_background
    std::fill( elligible, elligible + size, 0 );
    std::fill( labels, labels + size, -1 );

#endif

    if( index_frame == 0 ) {
        for( int i=0; i < SPNumber; i++ ) {
            Xseeds[ i ] = Xseeds_init[ i ];
            Yseeds[ i ] = Yseeds_init[ i ];

            index_xy = vertical_index[ (int) Yseeds[ i ] ] + Xseeds[ i ];

            lseeds[ i ] = lvec[ index_xy ];
            aseeds[ i ] = avec[ index_xy ];
            bseeds[ i ] = bvec[ index_xy ];

        }

    }

    memset( lseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( aseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( bseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Xseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Yseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( processed, 0, sizeof( int ) * size );

    for( int i=0; i<count_mask; i++ )
        mask_buffer[ i ].reset();

}

void IBIS::getLAB( cv::Mat* img ) {
    cv::Mat lab_image;
    cv::cvtColor(*img, lab_image, CV_BGR2Lab, 0);

    int ii = 0;

    if( img_bis ) {
        for (int i = 0; i < size * 3; i += 3) {
#if MATLAB_lab
            lvec_bis[ii] = lab_image.ptr()[i] * 100 / 255;
            avec_bis[ii] = lab_image.ptr()[i + 1] - 128;
            bvec_bis[ii] = lab_image.ptr()[i + 2] - 128;
#else
            lvec_bis[ii] = lab_image.ptr()[i];
            avec_bis[ii] = lab_image.ptr()[i + 1];
            bvec_bis[ii] = lab_image.ptr()[i + 2];
#endif
            ii++;
        }

        lvec = lvec_bis;
        avec = avec_bis;
        bvec = bvec_bis;

        img_bis = false;

    }
    else {
        for (int i = 0; i < size * 3; i += 3) {
#if MATLAB_lab
            lvec_one[ii] = lab_image.ptr()[i] * 100 / 255;
            avec_one[ii] = lab_image.ptr()[i + 1] - 128;
            bvec_one[ii] = lab_image.ptr()[i + 2] - 128;
#else
            lvec_one[ii] = lab_image.ptr()[i];
            avec_one[ii] = lab_image.ptr()[i + 1];
            bvec_one[ii] = lab_image.ptr()[i + 2];
#endif
            ii++;
        }

        lvec = lvec_one;
        avec = avec_one;
        bvec = bvec_one;

        img_bis = true;

    }

}

void IBIS::global_mean_seeds() {
    memset( countPx, 0, sizeof(int) * maxSPNumber );
    memset( Xseeds_Sum, 0, sizeof(float) * maxSPNumber );
    memset( Yseeds_Sum, 0, sizeof(float) * maxSPNumber );
    memset( lseeds_Sum, 0, sizeof(float) * maxSPNumber );
    memset( aseeds_Sum, 0, sizeof(float) * maxSPNumber );
    memset( bseeds_Sum, 0, sizeof(float) * maxSPNumber );

    memset( Xseeds, 0, sizeof(float) * maxSPNumber );
    memset( Yseeds, 0, sizeof(float) * maxSPNumber );
    memset( lseeds, 0, sizeof(float) * maxSPNumber );
    memset( aseeds, 0, sizeof(float) * maxSPNumber );
    memset( bseeds, 0, sizeof(float) * maxSPNumber );

    int sp;
#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
    for( int i=0; i<height; i++ ) {
        for( int j=0; j<width; j++ ) {
            sp = labels[ vertical_index[ i ] + j ];
            countPx[ sp ]++;
            Xseeds_Sum[ sp ] += j;
            Yseeds_Sum[ sp ] += i;
            lseeds_Sum[ sp ] += lvec[ vertical_index[ i ] + j ];
            aseeds_Sum[ sp ] += avec[ vertical_index[ i ] + j ];
            bseeds_Sum[ sp ] += bvec[ vertical_index[ i ] + j ];

        }

    }

#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] > 0 ) {
            Xseeds[ i ] = Xseeds_Sum[ i ] / countPx[ i ];
            Yseeds[ i ] = Yseeds_Sum[ i ] / countPx[ i ];
            lseeds[ i ] = lseeds_Sum[ i ] / countPx[ i ];
            aseeds[ i ] = aseeds_Sum[ i ] / countPx[ i ];
            bseeds[ i ] = bseeds_Sum[ i ] / countPx[ i ];

        }

    }

}

void IBIS::process( cv::Mat* img ) {
    double lap;

    if( size == 0 ) {
        size = img->cols * img->rows;
        width = img->cols;
        height = img->rows;

        // initialise all the buffer and inner parameters
        init();

        // STEP 1 : initialize with fix grid seeds value
        initSeeds();

        for( int i=0; i < SPNumber; i++ ) {
            Xseeds[ i ] = Xseeds_init[ i ];
            Yseeds[ i ] = Yseeds_init[ i ];

            int index_xy = vertical_index[ (int) Yseeds[ i ] ] + Xseeds[ i ];

            lseeds[ i ] = lvec_one[ index_xy ];
            aseeds[ i ] = avec_one[ index_xy ];
            bseeds[ i ] = bvec_one[ index_xy ];

        }

        // convert to Lab
        getLAB( img );
        index_frame = 0;

    }
    else {
        // convert to Lab
        getLAB( img );
        index_frame++;

#if fixed_background
        diff_frame();
#endif

    }

    // STEP 2 : process IBIS
#if OUTPUT_log
    lap = now_ms();
#endif

    reset();
    mask_propagate_SP();

#if OUTPUT_log
    st3 = now_ms() - lap;
#endif

    // STEP 3 : post processing
#if OUTPUT_log
    lap = now_ms();
#endif

    enforceConnectivity();
    update_adj();
    global_mean_seeds();
    memcpy( initial_repartition, labels, sizeof(int) * size );

    if( creation_deletion() ) {
#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
        for( int i=0; i<count_mask; i++ )
            mask_buffer[ i ].reset();

        mask_propagate_SP();

        //global_mean_seeds();
        update_adj();
        memcpy( initial_repartition, labels, sizeof(int) * size );

    }

    memcpy( Xseeds_prev, Xseeds, sizeof( float ) * maxSPNumber );
    memcpy( Yseeds_prev, Yseeds, sizeof( float ) * maxSPNumber );
    memcpy( lseeds_prev, lseeds, sizeof( float ) * maxSPNumber );
    memcpy( aseeds_prev, aseeds, sizeof( float ) * maxSPNumber );
    memcpy( bseeds_prev, bseeds, sizeof( float ) * maxSPNumber );
    memcpy( prev_adjacent_sp, adjacent_sp, sizeof( float ) * size_roi * maxSPNumber );
    memcpy( prev_count_adjacent, count_adjacent, sizeof( float ) * maxSPNumber );

#if OUTPUT_log
    st4 = now_ms() - lap;
#endif

    // output log
#if OUTPUT_log
    printf("-----------------\n");
    printf("PERF_T %lf\n", st3+st4);
    printf("IBIS.process\t\t%lf\t ms\n", st3);
    printf("IBIS.post_process\t%lf\t ms\n", st4);

    #if MASK_chrono
    float chrono[4] = { 0.f };
    for( int i=0; i < count_mask; i++ )
        mask_buffer[i].get_chrono( chrono );

    float total_chrono = chrono[0] + chrono[1] + chrono[2] + st2;

    printf("-----------------------------------\n");
    printf("Pixels processed:\t%lf\t\%\n", get_complexity()*100 );
    printf("-----------------\n");
    printf("\tMASK.angular_assign()\t\t%lf \%\n", chrono[0]/total_chrono);
    printf("\tMASK.fill_mask()\t\t%lf \%\n", chrono[1]/total_chrono);
    printf("\tMASK.assign_last()\t\t%lf \%\n", chrono[2]/total_chrono);
    printf("\tIBIS.mean_seeds()\t\t%lf \%\n", st2/total_chrono);

    #if THREAD_count > 1
    printf("-----------------\n");
    printf("multi-thread accel:\t\t\t%lf times\n", total_chrono/st3);
    printf("-----------------\n");
    #endif

    #endif

#endif

}

bool IBIS::creation_deletion() {

    std::fill( elligible, elligible + size, 0 );
    std::fill( count_diff, count_diff + maxSPNumber, 0 );

    bool status = false;
    int biggest_sp = 0;
    int best_value;
    float dist, D, conti;

    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] < minSPSizeThreshold && !count_diff[ i ] ) {
            status = true;

            //get the index of the biggest sp
            best_value = 0;
            for( int j=0; j<SPNumber; j++ ) {
                if( countPx[ j ] > best_value ) {
                    best_value = countPx[ j ];
                    biggest_sp = j;

                }

            }

            // copy seeds
            Xseeds[ i ] = Xseeds[ biggest_sp ] + 1;
            Yseeds[ i ] = Yseeds[ biggest_sp ] + 1;
            lseeds[ i ] = lseeds[ biggest_sp ];
            aseeds[ i ] = aseeds[ biggest_sp ];
            bseeds[ i ] = bseeds[ biggest_sp ];

            count_diff[ i ] = 1;
            count_diff[ biggest_sp ] = 1;

            // historical inheritance
            inheritance[ i ] = biggest_sp;

            adjacent_sp[size_roi*biggest_sp+count_adjacent[biggest_sp]] = i;
            count_adjacent[biggest_sp]++;

            // reinit SP
            countPx[ i ] = 0;
            countPx[ biggest_sp ] = 0;

            Xseeds_Sum[ i ] = 0;
            Yseeds_Sum[ i ] = 0;

            lseeds_Sum[ i ] = 0;
            aseeds_Sum[ i ] = 0;
            bseeds_Sum[ i ] = 0;

            Xseeds_Sum[ biggest_sp ] = 0;
            Yseeds_Sum[ biggest_sp ] = 0;
            lseeds_Sum[ biggest_sp ] = 0;
            aseeds_Sum[ biggest_sp ] = 0;
            bseeds_Sum[ biggest_sp ] = 0;

        }
        else
            inheritance[ i ] = i;

    }

    if( index_frame > 0 ) {
        float invwt_1 = (float)SPTypicalLength / 20;
        invwt_1 = 1.0f / (invwt_1 * invwt_1);

#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
        for( int i=0; i<SPNumber; i++ ) {
            // test continuity
            if( !count_diff[ i ] ) {
                best_value = -1;
                D = 0;

                for( int j=0; j<prev_count_adjacent[ i ]; j++ ) {
                    biggest_sp = prev_adjacent_sp[ size_roi*i + j ];

                    dist = ( Xseeds[ i ] - Xseeds_prev[ biggest_sp ] ) * ( Xseeds[ i ] - Xseeds_prev[ biggest_sp ] ) +
                           ( Yseeds[ i ] - Yseeds_prev[ biggest_sp ] ) * ( Yseeds[ i ] - Yseeds_prev[ biggest_sp ] );

                    conti = ( lseeds[ i ] - lseeds_prev[ biggest_sp ] ) * ( lseeds[ i ] - lseeds_prev[ biggest_sp ] ) +
                            ( aseeds[ i ] - aseeds_prev[ biggest_sp ] ) * ( aseeds[ i ] - aseeds_prev[ biggest_sp ] ) +
                            ( bseeds[ i ] - bseeds_prev[ biggest_sp ] ) * ( bseeds[ i ] - bseeds_prev[ biggest_sp ] );

                    conti += dist * invwt_1;

                    if( best_value < 0 || conti < D ) {
                        best_value = biggest_sp;
                        D = conti;

                    }

                }

                inheritance[ i ] = best_value;

            }

        }

    }

    if( status ) {
        for( int i=0; i<size; i++ ) {
            if( count_diff[ labels[i] ] )
                elligible[ i ] = 1;

        }

    }

    return status;

}

void IBIS::mask_propagate_SP() {

    double lap;
    st2=0;

    for( int mask_index=index_mask-1; mask_index>=0; mask_index--) {

#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
        for( int i=0; i<count_mask; i++ ) {
            mask_buffer[ i ].process( mask_index );

        }

#if MASK_chrono
        lap = now_ms();
#endif
        mean_seeds();
#if MASK_chrono
        st2 += now_ms() - lap;
#endif

#if VISU
        imagesc( std::string("processed"), processed, width, height );
        imagesc( std::string("labels"), labels, width, height );
        cv::waitKey(0);
#endif

    }

}

void IBIS::update_adj() {
    std::fill( updated_px, updated_px + size, false );
    std::fill( count_adjacent, count_adjacent + maxSPNumber, 0 );
    memset( x_vec, 0, sizeof(int) * size );
    memset( y_vec, 0, sizeof(int) * size );

    int dx4[] = { -1, 0, 1, 0 };
    int dy4[] = { 0, -1, 0, 1 };
    int* adj_label = new int[ getActualSPNumber() ];

    int current_sp, prospect_sp;
    int j, k;
    int current_index, prospect_index;

    // assure that a SP is a single 4-connected cluster
    int nb_adj;
    int count_vec;

    bool inscription;

    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // current index
            current_index = vertical_index[ y ] + x;

            if( !updated_px[ current_index ] ) {
                nb_adj = 1;
                current_sp = labels[  current_index  ];
                adj_label[ 0 ] = current_sp;

                // get adjacent px for direct neighbour
                x_vec[ 0 ] = x;
                y_vec[ 0 ] = y;
                count_vec = 1;

                for( int index_vec = 0; index_vec < count_vec; index_vec++  ) {
                    for( int index_var=0; index_var<4; index_var++ ) {
                        j = y_vec[ index_vec ] + dy4[ index_var ];
                        k = x_vec[ index_vec ] + dx4[ index_var ];

                        if( j >= 0 && j < height && k >= 0 && k < width ) {
                            prospect_index = vertical_index[ j ] + k;
                            prospect_sp = labels[ prospect_index ];

                            if( !updated_px[ prospect_index ] && prospect_sp == current_sp ) {

                                count_vec++;
                                x_vec[ count_vec ] = k;
                                y_vec[ count_vec ] = j;

                                updated_px[ prospect_index ] = true;

                            }
                            else {
                                if( prospect_sp != current_sp ) {
                                    inscription = true;

                                    // check if we know already about this one
                                    for( int index_adj=0; index_adj<nb_adj; index_adj++ ) {
                                        if( adj_label[ index_adj ] == prospect_sp ) {
                                            inscription = false;
                                            break;
                                        }

                                    }

                                    if( inscription ) {
                                        adj_label[ nb_adj ] = prospect_sp;
                                        nb_adj++;

                                    }

                                }

                            }

                        }

                    }

                }

                // single pixel
                if( count_vec == 1 )
                    updated_px[ current_index ] = true;

                // update sp neighbour
                for( int index_adj=0; index_adj<nb_adj; index_adj++ ) {
                    if( count_adjacent[ current_sp ] == 0 ) {
                        adjacent_sp[ size_roi*current_sp + count_adjacent[ current_sp ] ] = adj_label[ index_adj ];
                        count_adjacent[ current_sp ]++;

                    }
                    else {

                        bool status = false;
                        for( int j=0; j<count_adjacent[ current_sp ]; j++ ) {
                            if( adjacent_sp[ size_roi*current_sp + j ] == adj_label[ index_adj ] ) {
                                status = true;

                            }

                        }

                        if( !status ) {
                            adjacent_sp[ size_roi*current_sp + count_adjacent[ current_sp ] ] = adj_label[ index_adj ];
                            count_adjacent[ current_sp ]++;

                        }

                    }

                }

            }

        }

    }

    delete[] adj_label;

}

// ------------------------------------------------------------------------------- IBIS temporal
void IBIS::diff_frame() {
    int i = 0;
    float diff_dist;
    float sensitivity_limit;

    int sensitivity = 10;
    int inertia = 5;

    std::fill( count_diff, count_diff + SPNumber, 0 );
    std::fill( elligible, elligible + size, 0 );

    for( i = 0; i < size; i++ ) {
        diff_dist = ( lvec_one[ i ] - lvec_bis[ i ] ) * ( lvec_one[ i ] - lvec_bis[ i ] ) +
                    ( avec_one[ i ] - avec_bis[ i ] ) * ( avec_one[ i ] - avec_bis[ i ] ) +
                    ( bvec_one[ i ] - bvec_bis[ i ] ) * ( bvec_one[ i ] - bvec_bis[ i ] );

        if( diff_dist > inertia*inertia ) {
            elligible[ i ] = 1;
            count_diff[ labels[ i ] ]++;

        }

    }

    int count_reset=0;
    for( i = 0; i < SPNumber; i++ ) {
        sensitivity_limit = ( ( float( sensitivity ) / 100 ) * countPx[ i ] );

        if( count_diff[ i ] > sensitivity_limit ) {
            count_diff[ i ] = 1;
            count_reset++;

        }
        else
            count_diff[ i ] = 0;

    }

}

// ------------------------------------------------------------------------------- IBIS::MASK

void IBIS::MASK::init( int size_in, int y_in, int x_in, int mask_level, IBIS* ptr_IBIS, int* adjacent) {
    IBIS_data = ptr_IBIS;

    x_var = new int[ size_in ];
    y_var = new int[ size_in ];
    xy_var = new int[ size_in ];
    limit_value_fill = new int[4];
    spotted_sp = new int[255];
    count_adjacent_mask = 0;
    count_last=0;

    size = size_in;
    mask_index = mask_level;

    x = x_in;
    y = y_in;

    top_mask = false;
    int limit_value = ( IBIS_data->mask_size[ mask_index ] - 1 )/2;
    limit_value_fill[0] = ( x - limit_value < 0 ) ? 0 : x - limit_value;
    limit_value_fill[1] = ( x + limit_value > IBIS_data->width-1 ) ? IBIS_data->width-1 : x + limit_value;
    limit_value_fill[2] = ( y - limit_value < 0 ) ? 0 : y - limit_value;
    limit_value_fill[3] = ( y + limit_value >= IBIS_data->height ) ? IBIS_data->height - 1 : y + limit_value;

    if( mask_index > 0 ) {
        //generate sub_masks
        sub_mask = new IBIS::MASK[ 4 ];

        sub_mask[0].init( pow(2.0, mask_index+1), y - pow(2.0, mask_index-1), x - pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[1].init( pow(2.0, mask_index+1), y - pow(2.0, mask_index-1), x + pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[2].init( pow(2.0, mask_index+1), y + pow(2.0, mask_index-1), x - pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[3].init( pow(2.0, mask_index+1), y + pow(2.0, mask_index-1), x + pow(2.0, mask_index-1), mask_index-1, IBIS_data );

        if( adjacent != NULL ) { // top level mask
            top_mask=true;
            adjacent_mask = new int[9];

            for( int i=0; i<9; i++ ) {
                if( adjacent[i] >= 0 ) {
                    adjacent_mask[ count_adjacent_mask ] = adjacent[i];
                    count_adjacent_mask++;

                }

            }

        }

    }
    else {
        last_parent = new int[4];
        last_px_x = new int[9];
        last_px_y = new int[9];
        last_px_xy = new int[9];

        last_parent[0] = -1;
        last_parent[1] = -1;
        last_parent[2] = -1;
        last_parent[3] = -1;

        int j, k, index_y, index_xy;
        for( int index_var_y = -1; index_var_y <= 1; index_var_y++ ) {
            j = y + index_var_y;

            if( j < IBIS_data->height ) {
                index_y = IBIS_data->vertical_index[ j ];

                for( int index_var_x = -1; index_var_x <= 1; index_var_x++ ) {
                    k = x + index_var_x;
                    index_xy = index_y + k;

                    if( k < IBIS_data->width &&
                            !(index_var_x == -1 && index_var_y == -1 ) &&
                            !(index_var_x == -1 && index_var_y == 1 ) &&
                            !(index_var_x == 1 && index_var_y == -1 ) &&
                            !(index_var_x == 1 && index_var_y == 1 ) ) {
                        last_px_x[ count_last ] = k;
                        last_px_y[ count_last ] = j;
                        last_px_xy[ count_last ] = IBIS_data->vertical_index[ j ] + k;

                        count_last++;

                    }

                }

            }

        }

    }

    // fill mask coordinates
    generate_mask();

    // reset mask
    reset();

}

IBIS::MASK::~MASK() {
    delete[] x_var;
    delete[] y_var;
    delete[] xy_var;
    delete[] limit_value_fill;

    if( mask_index > 0 ) {
        delete[] sub_mask;

        if( top_mask )
            delete[] adjacent_mask;
    }
    else {
        delete[] last_px_x;
        delete[] last_px_y;
        delete[] last_px_xy;
        delete[] last_parent;

    }

}

void IBIS::MASK::reset() {
    filled = false;
    Mt1 = 0;
    Mt2 = 0;
    Mt3 = 0;
    count_SP=0;

    if( mask_index > 0 ) {
        sub_mask[ 0 ].reset();
        sub_mask[ 1 ].reset();
        sub_mask[ 2 ].reset();
        sub_mask[ 3 ].reset();

    }
    else {
        last_parent[0] = -1;
        last_parent[1] = -1;
        last_parent[2] = -1;
        last_parent[3] = -1;

    }

}

void IBIS::MASK::get_chrono( float* chrono ) {

    chrono[0]  += Mt1;
    chrono[1]  += Mt2;
    chrono[2]  += Mt3;

    if( mask_index > 0 ) {
        sub_mask[ 0 ].get_chrono( chrono );
        sub_mask[ 1 ].get_chrono( chrono );
        sub_mask[ 2 ].get_chrono( chrono );
        sub_mask[ 3 ].get_chrono( chrono );

    }

}

void IBIS::MASK::assign_labels( int y, int x, int index_xy, int value ) {
    IBIS_data->labels[ index_xy ] = value;

    IBIS_data->countPx[ value ]++;
    IBIS_data->lseeds_Sum[ value ]+=IBIS_data->lvec[ index_xy ];
    IBIS_data->aseeds_Sum[ value ]+=IBIS_data->avec[ index_xy ];
    IBIS_data->bseeds_Sum[ value ]+=IBIS_data->bvec[ index_xy ];
    IBIS_data->Xseeds_Sum[ value ]+=x;
    IBIS_data->Yseeds_Sum[ value ]+=y;

}

int IBIS::MASK::get_adj_mask() {
    return count_adjacent_mask;

}

int IBIS::MASK::get_list_sp() {
    // construct list
    list_sp();

    // return info
    return count_spotted_sp;

}

void IBIS::MASK::list_sp() {
    count_spotted_sp = 0;
    //std::fill( spotted_sp, spotted_sp + 255, -1 );

    if( mask_index > 0 ) {
        if( filled ) {
            count_spotted_sp = 1;
            spotted_sp[0] = IBIS_data->labels[ IBIS_data->vertical_index[ y ] + x ];

        }
        else {
            int* buff;
            for(int i=0; i<4; i++) {
                int c = sub_mask[i].get_list_sp();
                buff = sub_mask[i].get_spotted();

                if( c > 0 ) {
                    for( int j=0; j<c; j++ ) {
                        bool spotted = false;

                        if( count_spotted_sp > 0 ) {
                            for( int k=0; k<count_spotted_sp; k++ ) {
                                if( spotted_sp[k] == buff[j] ) {
                                    spotted = true;
                                    break;

                                }

                            }

                            if( !spotted ) {
                                spotted_sp[ count_spotted_sp ] = buff[j];
                                count_spotted_sp++;

                            }

                        }
                        else {
                            count_spotted_sp = 1;
                            spotted_sp[0] = buff[j];

                        }

                    }

                }

            }

        }

    }
    else {
        for( int i=0; i<count_last; i++ ) {
            if( count_spotted_sp > 0 ) {
                bool spotted = false;
                for( int k=0; k<count_spotted_sp; k++ ) {
                    if( spotted_sp[ k ] == IBIS_data->labels[ last_px_xy[ i ] ] ) {
                        spotted = true;
                        break;

                    }

                }

                if( !spotted ) {
                    spotted_sp[ count_spotted_sp ] = IBIS_data->labels[ last_px_xy[ i ] ];
                    count_spotted_sp++;

                }

            }
            else {
               count_spotted_sp = 1;
               spotted_sp[ 0 ] = IBIS_data->labels[ last_px_xy[ i ] ];

            }

        }

    }

}

void IBIS::MASK::generate_mask() {
    int* tmp_x_var = new int[ size ];
    int* tmp_y_var = new int[ size ];
    int limit_val, vertical_val, value_assign;
    int k = mask_index;
    int table_index = 0;

    if( k > 0 ) {
        limit_val = IBIS_data->mask_size[ k - 1 ] - 1;
        vertical_val = -limit_val;

        for( int index_var=1; index_var <= IBIS_data->mask_size[ k ]; index_var++ ) {
            if( index_var == 1 ) { //top border

                value_assign = -limit_val;
                for( int i=table_index; i<=table_index+limit_val; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + vertical_val;

                    value_assign += 2;
                }

                table_index += limit_val + 1;
                vertical_val += 2;

            }
            else if( index_var > 1 && index_var < IBIS_data->mask_size[ k - 1 ] ) { // vertical border

                value_assign = -limit_val;
                for( int i=table_index; i<table_index+2; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + vertical_val;

                    value_assign = limit_val;
                }

                table_index += 2;
                vertical_val += 2;

            }
            else { // bot border

                value_assign = -limit_val;
                for( int i=table_index; i<=table_index+limit_val; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + limit_val;

                    value_assign += 2;
                }

            }

        }

        // eliminate impossible coordinates
        count_var = 0;
        for( int i=0; i <= table_index + limit_val; i++ ) {
            if( tmp_x_var[ i ] >= 0 && tmp_x_var[ i ] < IBIS_data->width && tmp_y_var[ i ] >= 0 && tmp_y_var[ i ] < IBIS_data->height  ) {
                x_var[ count_var ] = tmp_x_var[ i ];
                y_var[ count_var ] = tmp_y_var[ i ];
                xy_var[ count_var ] = IBIS_data->vertical_index[ tmp_y_var[ i ] ] + tmp_x_var[ i ];
                count_var++;

            }

        }

    }
    else {
        tmp_x_var[ 0 ] = x-1;
        tmp_y_var[ 0 ] = y-1;

        tmp_x_var[ 1 ] = x+1;
        tmp_y_var[ 1 ] = y-1;

        tmp_x_var[ 2 ] = x-1;
        tmp_y_var[ 2 ] = y+1;

        tmp_x_var[ 3 ] = x+1;
        tmp_y_var[ 3 ] = y+1;

        count_var = 0;
        for( int i=0; i < 4; i++ ) {
            if( tmp_x_var[ i ] >= 0 && tmp_x_var[ i ] < IBIS_data->width && tmp_y_var[ i ] >= 0 && tmp_y_var[ i ] < IBIS_data->height  ) {
                x_var[ count_var ] = tmp_x_var[ i ];
                y_var[ count_var ] = tmp_y_var[ i ];
                xy_var[ count_var ] = IBIS_data->vertical_index[ tmp_y_var[ i ] ] + tmp_x_var[ i ];
                count_var++;

            }

        }

    }

    if( count_var == 0 )
        filled = 1;

    delete[] tmp_x_var;
    delete[] tmp_y_var;

}

void IBIS::MASK::process( int mask_level) {
    bool stat;
    double lap;

    if( ! filled ) {
        if( mask_level == mask_index ) {
#if MASK_chrono
            lap = IBIS_data->now_ms();
#endif
            stat = angular_assign();
#if MASK_chrono
            Mt1 += IBIS_data->now_ms() - lap;
#endif

            if( stat ) {
#if MASK_chrono
                lap = IBIS_data->now_ms();
#endif
                fill_mask();
#if MASK_chrono
                Mt2 += IBIS_data->now_ms() - lap;
#endif

                return;
            }
            else {
                if( mask_index == 0 ) {
#if MASK_chrono
                    lap = IBIS_data->now_ms();
#endif
                    assign_last();
#if MASK_chrono
                    Mt3 += IBIS_data->now_ms() - lap;
#endif

                }

            }

        }
        else {
            if( mask_index > 0 ) {
                sub_mask[ 0 ].process( mask_level );
                sub_mask[ 1 ].process( mask_level );
                sub_mask[ 2 ].process( mask_level );
                sub_mask[ 3 ].process( mask_level );

            }

        }

    }

#if VISU_all
    imagesc( std::string("processed"), IBIS_data->processed, IBIS_data->width, IBIS_data->height );
    cv::waitKey(1);
#endif

}

int IBIS::MASK::assign_last_px( int y, int x, int index_xy ) {
    int best_sp = -1;
    float dist_xy;
    float l, a, b;
    float dist_lab;
    int index_sp;
    float D=-1.f;
    float total_dist;

    for( int i=0; i<count_last_parent; i++ ) {
        index_sp = last_parent[ i ];

        l = IBIS_data->lvec[ index_xy ];
        a = IBIS_data->avec[ index_xy ];
        b = IBIS_data->bvec[ index_xy ];

        dist_lab = ( l - IBIS_data->lseeds[ index_sp ]) * ( l - IBIS_data->lseeds[ index_sp ]) +
                   ( a - IBIS_data->aseeds[ index_sp ]) * ( a - IBIS_data->aseeds[ index_sp ]) +
                   ( b - IBIS_data->bseeds[ index_sp ]) * ( b - IBIS_data->bseeds[ index_sp ]);

        dist_xy = ( x - IBIS_data->Xseeds[ index_sp ] ) * ( x - IBIS_data->Xseeds[ index_sp ] ) +
                  ( y - IBIS_data->Yseeds[ index_sp ] ) * ( y - IBIS_data->Yseeds[ index_sp ] );


        total_dist = dist_lab + dist_xy * IBIS_data->invwt;

        if( total_dist < D || D < 0) {
            best_sp = index_sp;
            D = total_dist;

        }

    }

    return best_sp;

}

int IBIS::MASK::assign_px( int y, int x, int index_xy ) {
    int best_sp = -1;
    float dist_xy;
    float l, a, b;
    float dist_lab;
    int index_sp;
    float D=-1.f;
    float total_dist;

    for( int i=0; i<IBIS_data->count_adjacent[IBIS_data->initial_repartition[index_xy]]; i++ ) {
        index_sp = IBIS_data->adjacent_sp[ size_roi*IBIS_data->initial_repartition[index_xy] + i ];

        l = IBIS_data->lvec[ index_xy ];
        a = IBIS_data->avec[ index_xy ];
        b = IBIS_data->bvec[ index_xy ];

        dist_lab = ( l - IBIS_data->lseeds[ index_sp ]) * ( l - IBIS_data->lseeds[ index_sp ]) +
                   ( a - IBIS_data->aseeds[ index_sp ]) * ( a - IBIS_data->aseeds[ index_sp ]) +
                   ( b - IBIS_data->bseeds[ index_sp ]) * ( b - IBIS_data->bseeds[ index_sp ]);

        dist_xy = ( x - IBIS_data->Xseeds[ index_sp ] ) * ( x - IBIS_data->Xseeds[ index_sp ] ) +
                  ( y - IBIS_data->Yseeds[ index_sp ] ) * ( y - IBIS_data->Yseeds[ index_sp ] );


        total_dist = dist_lab + dist_xy * IBIS_data->invwt;

        if( total_dist < D || D < 0) {
            best_sp = index_sp;
            D = total_dist;

        }

    }

    return best_sp;
}

void IBIS::MASK::assign_last() {

    int value;

    for( int index_last=0; index_last < count_last; index_last++ ) {
        if( IBIS_data->labels[ last_px_xy[ index_last ] ] < 0  || IBIS_data->elligible[ last_px_xy[ index_last ] ] ) {
            value = assign_last_px( last_px_y[ index_last ], last_px_x[ index_last ], last_px_xy[ index_last ] );
            assign_labels( last_px_y[ index_last ], last_px_x[ index_last ], last_px_xy[ index_last ], value );

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ last_px_xy[ index_last ] ] = 1;
#endif

        }

    }

}

void IBIS::MASK::fill_mask() {

    int j;
    int index_y;

    for( int index_var_y = limit_value_fill[2]; index_var_y <= limit_value_fill[3]; index_var_y++ ) {
        for( int index_var_x = limit_value_fill[0]; index_var_x <= limit_value_fill[1]; index_var_x++ ) {
            assign_labels( index_var_y, index_var_x, IBIS_data->vertical_index[ index_var_y ] + index_var_x, angular );

        }
    }

    filled = true;

}

bool IBIS::MASK::angular_assign() {
    int j, k;
    int best_sp;
    int index_xy;
    bool first=true;
    bool output = true;
    int ref_sp;
    angular = -1;

    if( mask_index > 0 ) {
        for( int index_var=0; index_var < count_var; index_var++ ) {
            j = y_var[ index_var ];
            k = x_var[ index_var ];
            index_xy = xy_var[ index_var ];

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ index_xy ] = 1;
#endif

            if( IBIS_data->labels[ index_xy ] < 0 || IBIS_data->elligible[ index_xy ] ) {
                // assign px
                angular = assign_px( j, k, index_xy );

            }
            else
                angular = IBIS_data->labels[ index_xy ];

            if( first ) {
                first = false;
                ref_sp = angular;

            }

            if( ref_sp != angular && !first ) {
                 return false;

            }

        }

    }
    else {
        count_last_parent = 0;

        for( int index_var=0; index_var < count_var; index_var++ ) {
            j = y_var[ index_var ];
            k = x_var[ index_var ];
            index_xy = xy_var[ index_var ];

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ index_xy ] = 1;
#endif

            if( IBIS_data->labels[ index_xy ] < 0 || IBIS_data->elligible[ index_xy ] ) {
                // assign px
                best_sp = assign_px( j, k, index_xy );

                //IBIS_data->labels[ index_xy ] = best_sp;
                assign_labels( j, k, IBIS_data->vertical_index[ j ] + k, best_sp );

            }

            angular = IBIS_data->labels[ index_xy ];

            if( (last_parent[0] != angular) &&
                (last_parent[1] != angular) &&
                (last_parent[2] != angular) &&
                (last_parent[3] != angular) ) {

                last_parent[ count_last_parent ] = angular;
                count_last_parent++;

            }

            if( first ) {
                ref_sp = angular;
                first = false;

            }

            if( ref_sp != angular && !first )
                output = false;

        }

    }

    return output;

}
