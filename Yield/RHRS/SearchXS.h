Double_t SearchXS(int Kin,Double_t aEp,Double_t aTheta,Double_t Theta[nTh],Double_t Eprime[nEp],Double_t xs_born[nTh][nEp],Double_t xs_rad[nTh][nEp]){
     Double_t Theta_center=0.0;
     //if(Kin==1)Theta_center=17.5717;
     //if(Kin==2)Theta_center=19.1125;
     //if(Kin==3)Theta_center=20.5751;
     //if(Kin==4)Theta_center=21.9401;
     Double_t th=aTheta*180/pi;
     int nfound=0;
     Double_t Radcor=0.0;
     Double_t aXS_born=0.0,aXS_rad=0.0;
     for(int ii=0;ii<nTh;ii++){
         Double_t dTheta=th-Theta[ii];
	 if(dTheta<0.25 && dTheta>=-0.25){
	    for(int jj=0;jj<nEp;jj++){
		Double_t dEp=aEp-Eprime[jj];
		if(abs(dEp)<0.0001){
		   Radcor=xs_born[ii][jj]/xs_rad[ii][jj];
                   nfound=1;
                   break;
		}
                if(aEp>Eprime[jj] && aEp<Eprime[jj+1]){
		   aXS_born=xs_born[ii][jj]+(xs_born[ii][jj+1]-xs_born[ii][jj])/(Eprime[jj+1]-Eprime[jj])*(aEp-Eprime[jj]);
		   aXS_rad=xs_rad[ii][jj]+(xs_rad[ii][jj+1]-xs_rad[ii][jj])/(Eprime[jj+1]-Eprime[jj])*(aEp-Eprime[jj]);
		   Radcor=aXS_born/aXS_rad;
                   nfound=1;
                   break;
		}
	    }
	 }
	 if(nfound==1)break;
         
     }
//     cout<<aEp<<"  "<<th<<"  "<<Radcor<<endl;
     return Radcor;

}
