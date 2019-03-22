#include <TMath.h>
const Double_t pi=TMath::Pi();
const int nEp=50;
const int nTh=50;
const Double_t delep=0.008;
const Double_t delth=0.08;

Double_t SearchXS(Double_t aEp,Double_t aTheta,Double_t Theta[nTh],Double_t Eprime[nEp], Double_t xs_rad[nTh][nEp]){
           int nfound=0;
           int bbp=0,bbTh=0;
           Double_t crad=0.0;
	   
//      theta
//      ^         a3_______a4
//      |          |___|___|     the XS table should be in ep increase and theta increase
//      |----->Ep  |___|___|
//                 a1     a2  
           for(int ii=0;ii<nEp;ii++){
	       if(Eprime[ii]==0)break;
	       Double_t tmp_Eprime=(int)(Eprime[ii]/0.0001+0.5)*0.0001;
               Double_t dEp=aEp-tmp_Eprime;
               if(dEp<=delep && dEp>=0){
                  bbp=ii;
                  for(int jj=0;jj<nTh;jj++){
		     if(Theta[jj]==0)break;
		     Double_t tmp_Theta=(int)(Theta[jj]/0.0001+0.5)*0.0001;
                     Double_t dTheta=aTheta-tmp_Theta;
                     if(dTheta<=delth && dTheta>=0){
                        nfound=1;
                     }
                     if(nfound==1){bbTh=jj;break;}
                  }
                  if(nfound==0){cout<<"Theta is out of XS table range !!!!"<<"  "<<aTheta<<endl; return 0;}
               }
               if(nfound==1)break;
           }
           if(nfound==0){cout<<"Ep: "<<aEp<<" is out of XS table range !!!!"<<endl;; return 0;}

           Float_t a1,a2,a3,a4;
           Float_t tmp1=0.0,tmp2=0.0;


           a1=xs_rad[bbTh][bbp];
           a2=xs_rad[bbTh+1][bbp];
           a3=xs_rad[bbTh][bbp+1];
           a4=xs_rad[bbTh+1][bbp+1];

           tmp1=a1+(a2-a1)/(Theta[bbTh+1]-Theta[bbTh])*(aTheta-Theta[bbTh]);
           tmp2=a3+(a4-a3)/(Theta[bbTh+1]-Theta[bbTh])*(aTheta-Theta[bbTh]);

           crad=tmp1+(tmp2-tmp1)/(Eprime[bbp+1]-Eprime[bbp])*(aEp-Eprime[bbp]);
           return crad;
      }

