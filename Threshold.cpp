#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "image.h"


void readImageHeader(char[], int &, int &, int &, bool &);

void readImage(char[], ImageType &);

void writeImage(char[], ImageType &);

int main(int argc, char *argv[]) {
    int i, j;
    int M, N, Q;
    bool type;

    RGB val,sval;
    float thresh;
    int SM, SN, SQ;
    bool Stype;
    fstream fout2;
    fout2.open("skinsample.csv", ios::out);

    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    // allocate memory for the image array
    ImageType image(N, M, Q);
    // read image
    readImage(argv[1], image);

    readImageHeader(argv[2], SN, SM, SQ, Stype);
    // allocate memory for the image array
    ImageType Simage(SN, SM, SQ);
    // read image
    readImage(argv[2], Simage);

    /*cout << "Enter threshold: ";
    cin >> thresh;*/

    // threshold image
    RGB white(255, 255, 255);
    RGB black(0, 0, 0);
   // cout << thresh << '\n';
    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++) {
            Simage.getPixelVal(i, j, sval);
            image.getPixelVal(i, j, val);
            float nr = (float) val.r / (val.r + val.g + val.b);
            float ng = (float) val.g / (val.r + val.g + val.b);
            if (nr < thresh)
                image.setPixelVal(i, j, white);
            else
                image.setPixelVal(i, j, val);

            if (sval==white)
            { image.setPixelVal(i, j, val);
                fout2<< nr<<","<<ng<<"\n";}
            else
                image.setPixelVal(i, j, white);

        }

    // write image
    writeImage(argv[3], image);

    return (1);
}
