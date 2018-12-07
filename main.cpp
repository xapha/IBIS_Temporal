/* -- Serge Bobbia : serge.bobbia@u-bourgogne.fr -- Le2i 2018
 * This work is distributed for non commercial use only,
 * it implements the IBIS method as described in the ICPR 2018 paper.
 * Read the ibis.h file for options and benchmark instructions
 *
 * This file show how to instanciate the IBIS class
 * You can either provide a file, or a directory, path to segment images
 */

#include <iostream>
#include "ibis.h"
#include <opencv2/opencv.hpp>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <highgui.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "utils.h"
#include "signal_processing.h"

#define SAVE_output         0
#define visu                1
#define visu_SNR	    1
#define signal_size         300
#define signal_processing   1

using namespace std;
//=================================================================================
/// DrawContoursAroundSegments
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void DrawContoursAroundSegments(
    unsigned char*&			ubuff,
    int*&					labels,
    const int&				width,
    const int&				height,
    const unsigned int&				color )
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int sz = width*height;
    vector<bool> istaken(sz, false);
    vector<int> contourx(sz);
    vector<int> contoury(sz);
    int mainindex(0);int cind(0);

    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    //if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind;//int(contourx.size());
    for( int j = 0; j < numboundpix; j++ )
    {
        int ii = contoury[j]*width + contourx[j];
        ubuff[ii] = 0xff;

        for( int n = 0; n < 8; n++ )
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            if( (x >= 0 && x < width) && (y >= 0 && y < height) )
            {
                int ind = y*width + x;
                if(!istaken[ind])
                    ubuff[ind] = 0;
            }
        }
    }
}


void write_labels(IplImage* input, const std::string& output_labels)
{
    std::ofstream file;
    file.open(output_labels.c_str());

    unsigned char* data = (unsigned char*)input->imageData;
    for (int y=0 ; y<input->height ; y++)
    {
        for (int x=0 ; x<input->width-1 ; x++)
        {
            file << (int) data[y*input->widthStep + x*input->nChannels] << " ";

        }
        file << (int) data[y*input->widthStep + (input->width -1)*input->nChannels] << std::endl;

    }

    file.close();
}

void write_traces(float* C1, float* C2, float* C3, const std::string& output_labels, IBIS* SP)
{
    std::ofstream file;
    file.open(output_labels.c_str());
    int max_sp = SP->getMaxSPNumber();

    for (int y=0 ; y<SP->getActualSPNumber() ; y++)
    {
        file << (double) C1[y] << " ";
        file << (double) C2[y] << " ";
        file << (double) C3[y] << std::endl;

    }

    file.close();
}

void execute_IBIS( int K, int compa, IBIS* Super_Pixel, Signal_processing* Signal, cv::Mat* img, std::string output_basename, int frame_index ) {

    int width = img->cols;
    int height = img->rows;
    int size = width * height;

    // process IBIS
    Super_Pixel->process( img );

    int* labels = Super_Pixel->getLabels();

    cv::Mat* output_bounds = new cv::Mat(cvSize(width, height), CV_8UC1);
    const int color = 0xFFFFFFFF;

    unsigned char* ubuff = output_bounds->ptr();
    std::fill(ubuff, ubuff + (width*height), 0);

    DrawContoursAroundSegments(ubuff, labels, width, height, color);

    cv::Mat* pImg = new cv::Mat(cvSize(width, height), CV_8UC3);
    float* sum_rgb = new float[Super_Pixel->getMaxSPNumber()*3];
    int* count_px = new int[Super_Pixel->getMaxSPNumber()];
    std::fill(sum_rgb, sum_rgb+Super_Pixel->getMaxSPNumber()*3, 0.f);
    std::fill(count_px, count_px+Super_Pixel->getMaxSPNumber(), 0);

    int ii = 0, i;
    for (i = 0; i < 3 * size; i += 3, ii++) {
        count_px[ labels[ii] ]++;
        sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 0 ] += img->ptr()[i];
        sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 1 ] += img->ptr()[i+1];
        sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 2 ] += img->ptr()[i+2];

    }

    float* R = new float[Super_Pixel->getMaxSPNumber()];
    float* G = new float[Super_Pixel->getMaxSPNumber()];
    float* B = new float[Super_Pixel->getMaxSPNumber()];

    float* R_avg = new float[Super_Pixel->getMaxSPNumber()];
    float* G_avg = new float[Super_Pixel->getMaxSPNumber()];
    float* B_avg = new float[Super_Pixel->getMaxSPNumber()];
    memset( R_avg, 0, sizeof(float) * Super_Pixel->getMaxSPNumber() );
    memset( G_avg, 0, sizeof(float) * Super_Pixel->getMaxSPNumber() );
    memset( B_avg, 0, sizeof(float) * Super_Pixel->getMaxSPNumber() );

    for (i=0; i<Super_Pixel->getMaxSPNumber(); i++) {
        sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 0 ] /= count_px[ i ];
        sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 1 ] /= count_px[ i ];
        sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 2 ] /= count_px[ i ];

        R[i] = sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 2 ];
        G[i] = sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 1 ];
        B[i] = sum_rgb[ i + Super_Pixel->getMaxSPNumber() * 0 ];

    }

    /*R = Super_Pixel->get_lseeds();
    G = Super_Pixel->get_aseeds();
    B = Super_Pixel->get_bseeds();*/

    // increase stat
    int* adj = Super_Pixel->get_adjacent_sp();
    int* nb_adj = Super_Pixel->nb_adjacent_sp();
    float dist;
    for (i=0; i<Super_Pixel->getActualSPNumber(); i++) {
        ii=0;

        for( int j=0; j<nb_adj[i]; j++ ) {
            int sp = adj[9*i+j];

            dist = ( R[ i ] - R[ sp ] ) * ( R[ i ] - R[ sp ] ) +
                   ( G[ i ] - G[ sp ] ) * ( G[ i ] - G[ sp ] ) +
                   ( B[ i ] - B[ sp ] ) * ( B[ i ] - B[ sp ] );
            dist /= 9;

            if( dist < 100.f ) {
                R_avg[ i ] += R[ sp ];
                G_avg[ i ] += G[ sp ];
                B_avg[ i ] += B[ sp ];
                ii++;

            }

        }

	R_avg[i] = R[i]/Super_Pixel->getActualSPNumber();
	G_avg[i] = G[i]/Super_Pixel->getActualSPNumber();
	B_avg[i] = B[i]/Super_Pixel->getActualSPNumber();
    }

    // signal processing
#if signal_processing
    Signal->add_frame( Super_Pixel->get_inheritance(),
                       R_avg,
                       G_avg,
                       B_avg,
                       Super_Pixel->getActualSPNumber() );

    Signal->process();
#endif
    if( frame_index % 30 == 0 ) {
        printf("-frame\t%i\n", frame_index);

    }

#if visu
        // SNR superposition
        const float* SNR;
        if( frame_index > signal_size ) {
            SNR = Signal->get_SNR();
        }

        for (i=0, ii=0; i < 3 * size; i += 3, ii++) {
            int sp = labels[ii];

            if (sp >= 0) {

                pImg->ptr()[i + 2]  = (unsigned char) img->ptr()[i+2];
                pImg->ptr()[i + 1]  = (unsigned char) img->ptr()[i+1];
                pImg->ptr()[i]      = (unsigned char) img->ptr()[i+0];

                /*pImg->ptr()[i + 2]  = (sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 2 ]);
                pImg->ptr()[i + 1]  = (sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 1 ]);
                pImg->ptr()[i]      = (sum_rgb[ labels[ii] + Super_Pixel->getMaxSPNumber() * 0 ]);*/

                if( ubuff[ ii ] == 255 ) {
                    pImg->ptr()[i + 2]  = 255;
                    pImg->ptr()[i + 1]  = 255;
                    pImg->ptr()[i]      = 255;

                }

#if signal_processing
                if( frame_index > signal_size ) {
                    int* color = Signal->get_color();

                    if( SNR[ labels[ii] ] > 0 && ubuff[ ii ] == 255 ) {
                        pImg->ptr()[i + 2]  = 0;
                        pImg->ptr()[i + 1]  = 0;
                        pImg->ptr()[i]      = 0;

                        if( SNR[ labels[ii] ] > 5 )
                            pImg->ptr()[i + color[ labels[ii] ]]  = 255;
                        else
                            pImg->ptr()[i + color[ labels[ii] ]]  = (unsigned char)(255 * SNR[ labels[ii] ] / 5 );

                    }

                }
#endif

            }

        }


#if signal_processing
        // add text
        char text[255] = "";
        sprintf( text, "HR: %i", Signal->get_HR() );
        cv::putText(*pImg, text, cv::Point(30,30),
            cv::FONT_HERSHEY_COMPLEX_SMALL, 0.8, cv::Scalar(200,200,250), 1, CV_AA);

#if visu_SNR
        if (frame_index > signal_size) {
        for( int i=0; i<Super_Pixel->getActualSPNumber(); i++ ) {
            sprintf( text, "%.1f", SNR[i] );

            cv::putText(*pImg, text, cv::Point( int(round(double(Super_Pixel->get_Xseeds()[i]))), int(round(double(Super_Pixel->get_Yseeds()[i]))) ),
                cv::FONT_HERSHEY_COMPLEX_SMALL, 0.6, cv::Scalar(0,0,250), 1, CV_AA);


        }
}
#endif
#endif

        cv::imshow("rgb mean", *pImg);
        cv::waitKey( 1 );

    //}
#endif

#if SAVE_output
    char output_labels[255] = {0};
    sprintf(output_labels, "results/%s/traces_%04i.seg", output_basename.c_str(), frame_index);
    write_traces( R, G, B, output_labels, Super_Pixel );

    sprintf(output_labels, "results/%s/parent_%04i.seg", output_basename.c_str(), frame_index);

    std::ofstream file;
    file.open(output_labels);
    int* parent = Super_Pixel->get_inheritance();

    for (int y=0; y<Super_Pixel->getActualSPNumber(); y++)
        file << parent[y] << std::endl;

    file.close();

    /*sprintf(output_labels, "results/%s/frame_%i.seg", output_basename.c_str(), frame_index);
    ofstream outfile;
    outfile.open(output_labels);
    for (int y=0 ; y<height ; y++)
    {
        for (int x=0 ; x<width-1 ; x++)
        {
            outfile << labels[y*width + x] << " ";
        }
        outfile << labels[y*width + width-1] << std::endl;
    }
    outfile.close();*/

#endif
    delete pImg;
    delete output_bounds;
    delete[] sum_rgb;
    delete[] count_px;
    delete[] R;
    delete[] G;
    delete[] B;
    delete[] R_avg;
    delete[] G_avg;
    delete[] B_avg;

}

int filter( const struct dirent *name ) {
    std::string file_name = std::string( name->d_name );
    std::size_t found = file_name.find(".avi");
    if (found!=std::string::npos) {
        return 1;

    }

    return 0;

}

int main( int argc, char* argv[] )
{
    printf(" - Temporal IBIS - \n\n");

    int K;
    int compa;

    if( argc != 4 ) {
        printf("--> usage ./IBIS_temporal SP_number Compacity File_path\n");
        printf(" |-> SP_number: user fixed number of superpixels, > 0\n");
        printf(" |-> Compacity: factor of caompacity, set to 20 for benchmark, > 0\n");
        printf(" |-> File_path: path to the file to compute\n");
        printf("  |-> if file_path is a directory, all the image within file_path/ are processed\n");
        printf("format: .png, .jpg, .ppm, .bmp, .tiff\n");
        printf("\n");
        printf("--> output file are saved in a \"./results\" directory\n");

        exit(EXIT_SUCCESS);

    }
    else {
        K = atoi( argv[ 1 ] );
        compa = atoi( argv[ 2 ] );

        if( K < 0 || compa < 0 ) {
            printf("--> usage ./IBIS_temporal SP_number Compacity File_path\n");
            printf(" |-> SP_number: user fixed number of superpixels, > 0\n");
            printf(" |-> Compacity: factor of caompacity, set to 20 for benchmark, > 0\n");
            printf(" |-> File_path: path to the file to compute\n");
            printf("  |-> if file_path is a directory, all the image within file_path/ are processed\n");
            printf("format: .png, .jpg, .ppm\n");
            printf("\n");
            printf("--> output file are saved in a \"./results\" directory\n");

            exit(EXIT_SUCCESS);

        }

    }

    // determine mode : file or path
    struct stat sb;

    if (stat(argv[3], &sb) == -1) {
        perror("stat");
        exit(EXIT_SUCCESS);
    }

    int type;
    //printf("file type : %i\n", sb.st_mode & S_IFMT);

    switch (sb.st_mode & S_IFMT) {
    case S_IFDIR:
        printf("directory processing\n");
        type=0;
        break;

    case S_IFREG:
        printf("single file processing\n");
        type=1;
        break;

    case 8192:
        printf("Device video processing\n");
        type=2;
        break;

    default:
        type=-1;
        break;

    }

    if( type == -1 )
        exit(EXIT_SUCCESS);
    else if( type >= 1 ) {
        // IBIS
        IBIS Super_Pixel( K, compa );
        Signal_processing Signal( K, signal_size );

        // get picture
        cv::VideoCapture video( argv[ 3 ] );
        if(!video.isOpened()) { // check if we succeeded
            printf("Can't open this device or video file.\n");
            exit(EXIT_SUCCESS);

        }

        cv::Mat img;
        int ii=0;
        std::string output_basename = std::string(argv[3]);

#if SAVE_output
        if( type == 1 ) {
            char command[255] = {0};
            sprintf( command, "mkdir -p results/%s\n", output_basename.c_str() );
            system( command );

        }
#endif

        while( video.read( img ) ) {
            execute_IBIS( K, compa, &Super_Pixel, &Signal, &img, output_basename, ii );
            ii++;

        }

    }
    else if( type == 0 ) {
        // get file list
        struct dirent **namelist;
        int n = scandir(argv[3], &namelist, &filter, alphasort);
        if (n == -1) {
           perror("scandir");
           exit(EXIT_FAILURE);

        }

        printf(" %i file(s) identified\n", n);
        if( n == 0 )
            exit(EXIT_SUCCESS);

        // process file list
        int width = 0;
        int height = 0;
        IBIS* Super_Pixel;
        Signal_processing Signal( K, signal_size );
        char* image_name = (char*)malloc(255);
        while (n--) {
            printf("processing %s\n", namelist[n]->d_name);

            // get picture

            sprintf(image_name, "%s/%s", argv[3], namelist[n]->d_name );
            cv::Mat img = cv::imread( image_name );

            // execute IBIS
            if( width == 0 ) {
                width = img.cols;
                height = img.rows;

                // IBIS
                Super_Pixel = new IBIS( K, compa );

            }
            else {
                if( width != img.cols ) {
                    delete Super_Pixel;
                    Super_Pixel = new IBIS( K, compa );

                    width = img.cols;
                    height = img.rows;

                }

            }

            execute_IBIS( K, compa, Super_Pixel, &Signal, &img, image_name, 0 );

            free(namelist[n]);

            printf("\n");
        }

        free( image_name );
        delete Super_Pixel;
        free( namelist );

    }

    exit(EXIT_SUCCESS);
}
