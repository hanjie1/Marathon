#include <TMath.h>
struct TRI_VAR{
      Double_t value;
      Double_t err;
};
const Double_t Qe=TMath::Qe();
const Double_t Na=TMath::Na();

const Double_t CM2toNB=1.0e33;
const Double_t pi=TMath::Pi();
 
//BCM constants
const Double_t dnew_gain=0.0003361;
const Double_t dnew_offset=0.0217;

//target thickness
const Double_t H3_t0=0.085099;
const Double_t He3_t0=0.0000226;
const Double_t tau=4500.0;
