#include <TMath.h>
struct TRI_VAR{
      Double_t value;
      Double_t err;
};

const Double_t dnew_gain=0.0003361;
const Double_t dnew_offset=0.0217;
const Double_t CM2toNB=1.0e33;

const int nTh=15;
const int nEp=32;
const Double_t pi=TMath::Pi();
//xbj with bin step 0.01;
//const double H1_xmin[5]={0.16,0.17,0.21,0.25,0.28};
//const double H1_xmax[5]={0.23,0.25,0.29,0.33,0.38};
//const double D2_xmin[11]={0.16,0.18,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.61,0.68};
//const double D2_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};
//const double He_xmin[11]={0.16,0.18,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.61,0.68};
//const double He_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};
//const double H3_xmin[11]={0.16,0.17,0.21,0.25,0.28,0.32,0.39,0.46,0.53,0.6,0.68};
//const double H3_xmax[11]={0.23,0.25,0.29,0.33,0.38,0.42,0.5,0.58,0.66,0.74,0.82};
//xbj with bin step 0.02
//const double H1_xmin[5]={0.16,0.18,0.2,0.24,0.28};
//const double H1_xmax[5]={0.22,0.24,0.28,0.32,0.36};
//const double D2_xmin[11]={0.16,0.18,0.2,0.24,0.28,0.32,0.38,0.46,0.54,0.6,0.68};
//const double D2_xmax[11]={0.22,0.24,0.28,0.32,0.36,0.4,0.5,0.58,0.66,0.74,0.8};
//const double He_xmin[11]={0.16,0.18,0.2,0.24,0.28,0.32,0.38,0.46,0.54,0.6,0.68};
//const double He_xmax[11]={0.22,0.24,0.28,0.32,0.36,0.4,0.5,0.58,0.66,0.74,0.82};
//const double H3_xmin[11]={0.16,0.18,0.2,0.24,0.28,0.32,0.38,0.46,0.54,0.6,0.68};
//const double H3_xmax[11]={0.22,0.24,0.28,0.32,0.36,0.4,0.5,0.58,0.66,0.74,0.82};
//xbj with specific x_min and x_max for each kinematics and each targets (above 25%)
const double H1_xmin[5]={0.159,0.175,0.210,0.248,0.283};
const double H1_xmax[5]={0.238,0.257,0.300,0.343,0.386};
const double D2_xmin[11]={0.159,0.176,0.211,0.247,0.283,0.318,0.391,0.463,0.535,0.607,0.676};
const double D2_xmax[11]={0.238,0.257,0.300,0.343,0.385,0.427,0.508,0.589,0.670,0.754,0.828};
const double He_xmin[11]={0.159,0.176,0.210,0.247,0.284,0.320,0.391,0.461,0.535,0.606,0.678};
const double He_xmax[11]={0.237,0.258,0.300,0.344,0.384,0.428,0.512,0.593,0.672,0.750,0.828};
const double H3_xmin[11]={0.159,0.175,0.210,0.247,0.283,0.318,0.390,0.462,0.534,0.606,0.678};
const double H3_xmax[11]={0.237,0.258,0.299,0.342,0.385,0.427,0.511,0.590,0.672,0.754,0.827};

const double xmin[11]={0.159,0.176,0.210,0.247,0.283,0.318,0.391,0.462,0.535,0.606,0.678};
const double xmax[11]={0.237,0.258,0.300,0.343,0.385,0.427,0.511,0.591,0.672,0.753,0.828};
//const int nBin[11]={4,4,4,5,5,5,6,6,8,7,6};



