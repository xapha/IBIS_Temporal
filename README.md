# IBIS_Temporal
Temporal propagation for IBIS superpixels.
The signal processing implementation done in this code along with the visualisation tools are only provided for test purpose for users. The signal processing implementation is far from optimized and can't run the code with real time constraint in that state. Futhermore, the results obtained in the CVPM publication refers to results obtained with an offline signal processing step done in Matlab.

The output is a list of file, one per frame, giving the superpixels inheritance for the current frame regarding the previous one and the average RGB values per superpixel.

## IBIS implementation

This implementation is done in C++ and was tested on unix environment.
For benchmark comparison, please refer to the *ibis.h* file for options description.
Futhermore, the visualisation along with the signal processing flag should be deactivated in the main.cpp file as they are only provided for users test.

## Dependences :

You will need *cmake* and *openCV* ( versions 3.X ) to run this code.
For paralell execution, you will need *openMP* as well.

## Compilation

```Shell Session
git clone "https://github.com/xapha/IBIS_temporal.git"
cd IBIS_temporal
mkdir build && cd build
mkdir results
cmake -D CMAKE_BUILD_TYPE=Release ..
make
./IBIS_temporal
