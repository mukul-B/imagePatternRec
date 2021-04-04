#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

#include "image.h"
#include "MatricOP.h"
#include "DiscriminantCases.h"


void readImageHeader(char[], int &, int &, int &, bool &);

void readImage(char[], ImageType &);

void writeImage(char[], ImageType &);

void derive_distribution(char [], ImageType &);

int main(int argc, char *argv[]) {
    int i, j;
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

    // derive_distribution(argv[2],image);


    vector<vector<double> > mu1 = {{0.430199},
                                   {0.296412}};
    vector<vector<double> > cov1 = {{0.000118533, 0},
                                    {0,           0.0000702336}};
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

    vector<vector<double> > sampler;


    Matrix sample;


    ThresholdCase3 classifier = ThresholdCase3(m1, cv1, 0.5);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {


            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);

            sampler = {{nr},
                       {ng}};
            sample.input(sampler, mudimen);
            if (classifier.getDecision(sample) > 5) {
                image.setPixelVal(i, j, val);
            } else {
                image.setPixelVal(i, j, black);
            }

        }

    // write image
    writeImage(argv[3], image);

    return 0;
}

void derive_distribution(char sample_file[], ImageType &image) {
    fstream fout2;
    fout2.open("skinsample.csv", ios::out);

    RGB sval, val;
    RGB white(255, 255, 255);

    int SM, SN, SQ;
    bool Stype;
    double skin_nr_mu = 0;
    double skin_ng_mu = 0;
    int skin_sample_size = 0;
    double skin_nr_sigmaSq = 0;
    double skin_ng_sigmaSq = 0;
    double skin_co_sigmaSq = 0;

    /*double non_skin_nr_mu=0;
    double non_skin_ng_mu=0;
    int non_skin_sample_size=0;
    double non_skin_nr_sigmaSq=0;
    double non_skin_ng_sigmaSq=0;
    double non_skin_co_sigmaSq=0;*/

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
            if (sval == white) { //image.setPixelVal(i, j, val);
                // fout2<< nr<<","<<ng<<"\n";
                skin_nr_mu += nr;
                skin_ng_mu += ng;
                skin_sample_size += 1;
            } else {
                // image.setPixelVal(i, j, white);
                /* non_skin_nr_mu+=nr;
                 non_skin_ng_mu+=ng;
                 non_skin_sample_size+=1;*/
            }

        }
    skin_nr_mu = skin_nr_mu / skin_sample_size;
    skin_ng_mu = skin_ng_mu / skin_sample_size;

    /*non_skin_nr_mu=non_skin_nr_mu/non_skin_sample_size;
    non_skin_ng_mu=non_skin_ng_mu/non_skin_sample_size;*/
    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            image.getPixelVal(i, j, val);
            double nr = (double) val.r / (val.r + val.g + val.b);
            double ng = (double) val.g / (val.r + val.g + val.b);

            if (sval == white) { //image.setPixelVal(i, j, val);
                //fout2<< nr<<","<<ng<<"\n";
                skin_nr_sigmaSq += (pow(nr - skin_nr_mu, 2));
                skin_ng_sigmaSq += (pow(ng - skin_ng_mu, 2));
                skin_co_sigmaSq += ((ng - skin_ng_mu) * (nr - skin_nr_mu));
                fout2 << nr << "," << ng << "," << skin_nr_sigmaSq << "," << skin_ng_sigmaSq << "\n";

            } else {
                //image.setPixelVal(i, j, white);
                /*non_skin_nr_sigmaSq+=((nr-skin_nr_mu)*(nr-skin_nr_mu));
                non_skin_ng_sigmaSq+=((ng-skin_ng_mu)*(ng-skin_ng_mu));
                non_skin_co_sigmaSq+=((ng-skin_ng_mu)*(nr-skin_nr_mu));*/
            }

        }

    cout << skin_nr_sigmaSq << "," << skin_ng_sigmaSq << "\n";
    skin_nr_sigmaSq = pow(skin_nr_sigmaSq, 0.5) / skin_sample_size;
    skin_ng_sigmaSq = pow(skin_ng_sigmaSq, 0.5) / skin_sample_size;
    skin_co_sigmaSq = pow(skin_co_sigmaSq, 0.5) / skin_sample_size;

    cout << skin_nr_mu << "," << skin_ng_mu << "\n";
    cout << skin_nr_sigmaSq << "," << skin_ng_sigmaSq << "\n";
    cout << skin_co_sigmaSq << "\n";
    cout << skin_sample_size;

    /*non_skin_nr_sigmaSq=pow(non_skin_nr_sigmaSq,0.5)/non_skin_sample_size;
    non_skin_ng_sigmaSq=pow(non_skin_ng_sigmaSq,0.5)/non_skin_sample_size;
    non_skin_co_sigmaSq=pow(non_skin_co_sigmaSq,0.5)/non_skin_sample_size;

    cout<<  non_skin_nr_mu <<","<<non_skin_ng_mu<<"\n";
    cout<<  non_skin_nr_sigmaSq <<","<<non_skin_ng_sigmaSq<<"\n";
    cout<<  non_skin_co_sigmaSq<<"\n";
    cout<< non_skin_sample_size ;*/
}
