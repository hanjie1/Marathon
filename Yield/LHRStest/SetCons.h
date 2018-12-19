#include <TMath.h>
struct TRI_VAR{
      Double_t value;
      Double_t err;
};
 
const int nTh=15;
const int nEp=32;
const double dBin=0.01;
const Double_t pi=TMath::Pi();
const double H1_xmin[5]={0.16,0.17,0.21,0.25,0.28};
const double H1_xmax[5]={0.23,0.25,0.29,0.33,0.38};
const double D2_xmin[11]={0.16,0.18,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.61,0.68};
const double D2_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};
const double He_xmin[11]={0.16,0.18,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.61,0.68};
const double He_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};
const double H3_xmin[11]={0.16,0.17,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.6,0.68};
const double H3_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};