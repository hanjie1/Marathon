#include <TMath.h>
struct TRI_VAR{
      Double_t value;
      Double_t err;
};
const Double_t Qe=TMath::Qe();
const Double_t Na=TMath::Na();

const Double_t CM2toNB=1.0e33;
const Double_t pi=TMath::Pi();

//boiling constants
const Double_t bH3_A=0.0001293;
const Double_t bH3_B=-0.007399;

const Double_t bHe3_A=0.00008686;
const Double_t bHe3_B=-0.004759;

const Double_t bD2_A=0.0001147;
const Double_t bD2_B=-0.006651;

const Double_t bH1_A=0.0001527;
const Double_t bH1_B=-0.008529;
 
//BCM constants
const Double_t dnew_gain=0.0003361;
const Double_t dnew_offset=0.0217;

//target thickness
const Double_t H3_pho=0.085099;
const Double_t H1_pho=0.0708;
const Double_t D2_pho=0.14215;
const Double_t He3_pho=0.0533752;

const int nTh=100;
const int nEp=50;
const Double_t deltaTh=0.08;
const Double_t deltaEp=0.008;

const Double_t xmin[11]={0.153,0.168,0.203,0.237,0.273,0.307,0.378,0.449,0.522,0.593,0.665}; 
const Double_t xmax[11]={0.263,0.283,0.329,0.370,0.418,0.454,0.543,0.641,0.717,0.783,0.849};
const int nBin[11]={4,3,4,4,4,4,4,4,4,4,4};
const Double_t dBin[11][4]={{0.0275,0.0275,0.0275,0.0275},
                           {0.038,0.029,0.048,0},
			   {0.029,0.029,0.029,0.039},
			   {0.033,0.03,0.03,0.040},
			   {0.037,0.030,0.028,0.05},
			   {0.040,0.028,0.029,0.05},
			   {0.05,0.0325,0.0325,0.05},
			   {0.05,0.034,0.037,0.05},
			   {0.05,0.045,0.05,0.05},
			   {0.05,0.04,0.05,0.05},
			   {0.05,0.044,0.045,0.05}};