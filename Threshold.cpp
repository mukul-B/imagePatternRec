#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

#include "image.h"
#include "MatricOP.h"
#include "DiscriminantCases.h"

#define black_ind 0
#define non_black_ind 1


void readImageHeader(char[], int &, int &, int &, bool &);

void readImage(char[], ImageType &);

void writeImage(char[], ImageType &);

void derive_distribution(char [], ImageType &);

ThresholdCase3 get_classifier();
vector<int> get_reference(char sample_file[]) {
    RGB sval;
    RGB white(255, 255, 255);
    int SM, SN, SQ;
    bool Stype;
    readImageHeader(sample_file, SN, SM, SQ, Stype);
    // allocate memory for the image array
    ImageType Simage(SN, SM, SQ);
    // read image
    readImage(sample_file, Simage);

    RGB black(0, 0, 0);
    vector<int> ref;
    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);

            if (sval == black) {
                ref.push_back(black_ind);
                // cout<< sval.r << sval.g << sval.b ;
            } else {
                ref.push_back(non_black_ind);
                //cout<< sval.r << sval.g << sval.b ;
            }
        }
    return ref;
}

vector<Matrix> get_normailizeScala(ImageType image) {
    pair<int, int> mudimen = {2, 1};
    vector<Matrix> sam;
    Matrix sample;
    RGB val;
    int M, N, Q;

    image.getImageInfo(N, M, Q);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);
            vector<vector<double> > sampler
                    = {{nr},
                       {ng}};
            //Matrix sample;
            sample.input(sampler, mudimen);
            sam.push_back(sample);
        }
            return sam;
}

void roc(vector<Matrix> sam, vector<int> ref_image, int N, int M, ThresholdCase3 classifier) {
    /*fstream fout;
    fout.open("roc.csv", ios::out);*/
    int FP = 0, FN = 0;
    for (int t = 60; t < 181; t= t+ 10) {
        FP = 0;
        FN = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (classifier.getDecision(sam[i * M + j]) > t) {
                    if (ref_image[i * M + j] == black_ind)
                        FP++;
                } else {
                    if (ref_image[i * M + j] == non_black_ind)
                        FN++;
                }
            }
        }
        cout << t << "," << FP << "," << FN <<","<<FP-FN<< "\n";
    }
}
ThresholdCase3 get_classifier(){
     vector<vector<double> > mu1 = {{0.432224},
                                    {0.295771}};
     vector<vector<double> > cov1 = {{0.00243324,-0.00111725},
                                     {-0.00111725,0.000788441}};
     //-0.00111725
     /*vector<vector<double> > mu1 = {{1},
                                    {1}};
     vector<vector<double> > cov1  = {{1, 0},
                                      {0, 1}};*/

     int noOfFfeatures = 2;
     pair<int, int> mudimen = {noOfFfeatures, 1};
     pair<int, int> covdimen = {noOfFfeatures, noOfFfeatures};
     Matrix m1;
     m1.input(mu1, mudimen);
     Matrix cv1;
     cv1.input(cov1, covdimen);
     ThresholdCase3 classifier = ThresholdCase3(m1, cv1, 0.5);
   // cout << classifier.getDecision(m1)<< endl;
    //5.88885
     return classifier;

}
int main(int argc, char *argv[]) {

    int M, N, Q;
    bool type;
    RGB black(0, 0, 0);
    RGB val;

    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    // allocate memory for the image array
    ImageType image(N, M, Q);
    // read image
    readImage(argv[1], image);
    vector<Matrix> sam=get_normailizeScala(image);
    vector<int> ref_image = get_reference(argv[2]);
    //derive_distribution(argv[2],image);
    ThresholdCase3 classifier = get_classifier();
   // roc(sam, ref_image, N, M, classifier);

   /* vector<Matrix> sami;
    Matrix sample;
    int noOfFfeatures = 2;
    pair<int, int> mudimen = {noOfFfeatures, 1};*/

   //77290,77235,55

  int FP = 0, FN = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);
            //5.99030625
            if (classifier.getDecision(sam[i * M + j]) > 127.15) {
                image.setPixelVal(i, j, val);
                if (ref_image[i * M + j] == black_ind)
                    FP++;
            } else {
                image.setPixelVal(i, j, black);
                if (ref_image[i * M + j] == non_black_ind)
                    FN++;
            }

        }
    cout << FP << "," << FN <<","<<FP-FN<< "\n";
    // write image
    writeImage(argv[3], image);


    return 0;
}

void derive_distribution(char sample_file[], ImageType &image) {
    fstream fout2;
    fout2.open("skinsample.csv", ios::out);

    RGB sval, val;
    RGB white(255, 255, 255);
    RGB black(0, 0, 0);

    int SM, SN, SQ;
    bool Stype;
    double skin_nr_mu = 0;
    double skin_ng_mu = 0;
    int skin_sample_size = 0;
    double skin_nr_sigmaSq = 0;
    double skin_ng_sigmaSq = 0;
    double skin_co_sigmaSq = 0;

    readImageHeader(sample_file, SN, SM, SQ, Stype);
    // allocate memory for the image array
    ImageType Simage(SN, SM, SQ);
    // read image
    readImage(sample_file, Simage);

    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);
            if (sval != black) {
                skin_nr_mu += nr;
                skin_ng_mu += ng;
                skin_sample_size += 1;
            }
        }
    skin_nr_mu = skin_nr_mu / skin_sample_size;
    skin_ng_mu = skin_ng_mu / skin_sample_size;

    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);
            if (sval != black) {
                skin_nr_sigmaSq += (pow(nr - skin_nr_mu, 2));
                skin_ng_sigmaSq += (pow(ng - skin_ng_mu, 2));
                skin_co_sigmaSq += ((ng - skin_ng_mu) * (nr - skin_nr_mu));
            }
        }


    skin_nr_sigmaSq = skin_nr_sigmaSq / skin_sample_size;
    skin_ng_sigmaSq = skin_ng_sigmaSq / skin_sample_size;
    skin_co_sigmaSq = skin_co_sigmaSq / skin_sample_size;

    cout << "mu nr: "<<skin_nr_mu << "," << "mu nr: " << skin_ng_mu << "\n";
    cout << "sigma square nr: "<< skin_nr_sigmaSq << "," << "sigma square nr: "<< skin_ng_sigmaSq << "\n";
    cout << "sigma square covariance: "<< skin_co_sigmaSq << "\n";
    cout << skin_sample_size<<"\n";
}