#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "image.h"
#include "MatricOP.h"
#include "DiscriminantCases.h"

#define black_ind 0
#define non_black_ind 1

#define RGBColorSpace 0
#define YCbCrColorSpace 1

using namespace std;

void readImageHeader(char[], int &, int &, int &, bool &);

void readImage(char[], ImageType &);

void writeImage(char[], ImageType &);

vector<Matrix> getNormalisedPixelMatrix(ImageType, int);

vector<int> getReferenceImageMatrix(char[]);

void derive_distribution(vector<int>, ImageType, int);

ThresholdBased get_classifier();

void generateROC(vector<Matrix>, vector<int>, ThresholdBased, int, int, int);

vector<vector<double> > getImageFeature(RGB, int type);

void classsification_process(vector<int> ref_image, ThresholdBased classifier, ImageType &image);

struct sampleSpace {
    vector<vector<double> > mu;
    vector<vector<double> > cov;
    double threshold;
    int color_space;
} dist;

void samplespace1() {
    dist.mu = {{0.432224},
               {0.295771}};
    dist.cov = {{0.00243324,  -0.00111725},
                {-0.00111725, 0.000788441}};
    dist.threshold = 127.15;
    dist.color_space = RGBColorSpace;
}

void samplespace2() {
    dist.mu = {{23.6343},
               {-12.6913}};
    dist.cov = {{48.9456,  -21.7417},
                {-21.7417, 31.197}};
    dist.threshold = 0.00208;
    dist.color_space = YCbCrColorSpace;
}


int main(int argc, char *argv[]) {

    int M, N, Q;
    bool type;

    // read image header
    readImageHeader(argv[1], N, M, Q, type);
    // allocate memory for the image array
    ImageType image(N, M, Q);
    // read image
    readImage(argv[1], image);

    int prob;

    printf("Programing Assignment 2:  question 3\n");
    printf("------------------------------------------\n");
    printf("Main Menu\n");
    printf("1. 3a Training process with RGBColorSpace .\n");
    printf("2. 3a Testing process with RGBColorSpace.\n");
    printf("3. 3b Training process with YCbCrColorSpace.\n");
    printf("4. 3b Testing process with YCbCrColorSpace.\n");
    printf(" Please enter an option from the main menu: ");

    fflush(stdin);
    cin >> prob;


    int color_space = YCbCrColorSpace;
    if (prob == 1 || prob == 2)
        color_space = RGBColorSpace;

    bool testing = false;
    if (prob == 2 || prob == 4) {
        testing = true;
        if (color_space == RGBColorSpace) samplespace1();
        else samplespace2();
    }


    vector<int> ref_image = getReferenceImageMatrix(argv[2]);
    vector<Matrix> sam = getNormalisedPixelMatrix(image, color_space);

    if (!testing)
        derive_distribution(ref_image, image, color_space);

    ThresholdBased classifier = get_classifier();

    if (testing)
        classsification_process(ref_image, classifier, image);
    else
        generateROC(sam, ref_image, classifier, N, M, color_space);

    // write image
    writeImage(argv[3], image);

    return 0;
}

void classsification_process(vector<int> ref_image, ThresholdBased classifier, ImageType &image) {
    Matrix sample;
    int noOfFfeatures = 2;
    pair<int, int> mudimen = {noOfFfeatures, 1};
    RGB val;
    int N, M, Q;
    RGB black(0, 0, 0);

    image.getImageInfo(N, M, Q);
    int FP = 0, FN = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.getPixelVal(i, j, val);
            vector<vector<double> > sampler = getImageFeature(val, dist.color_space);
            sample.input(sampler, mudimen);
            if (classifier.getDecision(sample) > dist.threshold) {
                // if (classifier.getDecision(sam[i * M + j]) > 127.15) {
                image.setPixelVal(i, j, val);
                if (ref_image[i * M + j] == black_ind)
                    FP++;
            } else {
                image.setPixelVal(i, j, black);
                if (ref_image[i * M + j] == non_black_ind)
                    FN++;
            }

        }
    cout << FP << "," << FN << "," << FP - FN << "\n";
}

vector<vector<double> > getImageFeature(RGB val, int type = 0) {
    double nr;
    double ng;
    if (type == RGBColorSpace) {
        nr = (double) val.r / (val.r + val.g + val.b);
        ng = (double) val.g / (val.r + val.g + val.b);
    } else if (type == YCbCrColorSpace) {
        nr = (double) (0.500) * val.r - (0.419) * val.g - (0.081) * val.b;
        ng = (double) (-0.169) * val.r - (0.332) * val.g + (0.500) * val.b;
    }
    vector<vector<double> > sampler
            = {{nr},
               {ng}};
    return sampler;
}

vector<int> getReferenceImageMatrix(char sample_file[]) {
    RGB sval;
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
            } else {
                ref.push_back(non_black_ind);
            }
        }
    return ref;
}

vector<Matrix> getNormalisedPixelMatrix(ImageType image, int color_space) {
    pair<int, int> mudimen = {2, 1};
    vector<Matrix> sam;
    Matrix sample;
    RGB val;
    int M, N, Q;

    image.getImageInfo(N, M, Q);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.getPixelVal(i, j, val);
            vector<vector<double> > sampler = getImageFeature(val, color_space);
            sample.input(sampler, mudimen);
            sam.push_back(sample);
        }
    return sam;
}

void generateROC(vector<Matrix> sam, vector<int> ref_image, ThresholdBased classifier, int N, int M, int color_space) {
    /*fstream fout;
    fout.open("roc.csv", ios::out);*/

    cout << "Generating ROC data points";
    int FP = 0, FN = 0;

    int factor;
    double lw, up;

    if (color_space == RGBColorSpace) {
        factor = 1;
        up = 181;
        lw = 60;
    } else {
        factor = 100000;
        lw = 140;
        up = 261;
    }

    for (double t = lw; t < up; t = t + 10) {
        FP = 0;
        FN = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (classifier.getDecision(sam[i * M + j]) > (double) t / factor) {
                    if (ref_image[i * M + j] == black_ind)
                        FP++;
                } else {
                    if (ref_image[i * M + j] == non_black_ind)
                        FN++;
                }
            }
        }
        cout << (double) (t / factor) << "," << FP << "," << FN << "," << FP - FN << "\n";
    }
}

ThresholdBased get_classifier() {

    int noOfFfeatures = 2;
    pair<int, int> mudimen = {noOfFfeatures, 1};
    pair<int, int> covdimen = {noOfFfeatures, noOfFfeatures};
    Matrix m1;
    m1.input(dist.mu, mudimen);
    Matrix cv1;
    cv1.input(dist.cov, covdimen);
    ThresholdBased classifier = ThresholdBased(m1, cv1, 0.5);
    cout << classifier.getDecision(m1) << endl;

    return classifier;
}

void derive_distribution(vector<int> ref_image, ImageType image, int color_space) {
    fstream fout2;
    fout2.open("skinsample.csv", ios::out);
    int SM, SN, SQ;
    image.getImageInfo(SN, SM, SQ);
    RGB val;

    double skin_nr_mu = 0;
    double skin_ng_mu = 0;
    int skin_sample_size = 0;
    int total_sample_size = 0;
    double skin_nr_sigmaSq = 0;
    double skin_ng_sigmaSq = 0;
    double skin_co_sigmaSq = 0;

    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            image.getPixelVal(i, j, val);
            vector<vector<double> > sampler = getImageFeature(val, color_space);
            double nr = sampler[0][0];
            double ng = sampler[1][0];
            if (ref_image[i * SM + j] == non_black_ind) {
                skin_nr_mu += nr;
                skin_ng_mu += ng;
                skin_sample_size += 1;
            }
            total_sample_size++;
        }
    skin_nr_mu = skin_nr_mu / skin_sample_size;
    skin_ng_mu = skin_ng_mu / skin_sample_size;

    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            image.getPixelVal(i, j, val);
            vector<vector<double> > sampler = getImageFeature(val, color_space);
            double nr = sampler[0][0];
            double ng = sampler[1][0];
            if (ref_image[i * SM + j] == non_black_ind) {
                skin_nr_sigmaSq += (pow(nr - skin_nr_mu, 2));
                skin_ng_sigmaSq += (pow(ng - skin_ng_mu, 2));
                skin_co_sigmaSq += ((ng - skin_ng_mu) * (nr - skin_nr_mu));
            }
        }
    skin_nr_sigmaSq = skin_nr_sigmaSq / skin_sample_size;
    skin_ng_sigmaSq = skin_ng_sigmaSq / skin_sample_size;
    skin_co_sigmaSq = skin_co_sigmaSq / skin_sample_size;

    cout << "mu nr: " << skin_nr_mu << "," << "mu ng: " << skin_ng_mu << "\n";
    cout << "sigma square nr: " << skin_nr_sigmaSq << "," << "sigma square ng: " << skin_ng_sigmaSq << "\n";
    cout << "sigma square covariance: " << skin_co_sigmaSq << "\n";
    cout << "skin color prior: " << (double) skin_sample_size / total_sample_size << "\n";
    cout << skin_sample_size << "\n";

    dist.mu = {{skin_nr_mu},
               {skin_ng_mu}};
    dist.cov = {{skin_nr_sigmaSq, skin_co_sigmaSq},
                {skin_co_sigmaSq, skin_ng_sigmaSq}};
}
