void plot_SLAC()
{
     Double_t SLAC_x[12]={0.1722,0.1898,0.2044,0.2206,0.2407,0.2585,0.276,0.293,0.312,0.3305,0.357,0.3852};
     Double_t SLAC_D2[12]={11700,11700,11300,11000,10300,9850,9490,8970,8590,8210,7390,6900};
     Double_t SLAC_D2_err[12]={297,230,214,212,202,200,174,153,145,142,133,127};
     Double_t SLAC_H1[12]={6520,6670,6670,6420,6330,5860,5550,5340,5080,5030,4770,4510};
     Double_t SLAC_H1_err[12]={171,162,152,163,165,144,141,99,108,112,113,101};
     Double_t SLAC_N[12]={5340,5250,4890,4760,4180,4190,4150,3840,3710,3380,2810,2570}; 
     Double_t SLAC_N_err[12]={313,242,222,232,228,216,193,147,150,153,150,140};

     Double_t SLAC_r[12],SLAC_r_err[12];
     Double_t R=1.0;
     for(int ii=0;ii<12;ii++){
         R=1.0095-0.0109*SLAC_x[ii]-0.0821*SLAC_x[ii]*SLAC_x[ii];
	 //SLAC_r[ii]=SLAC_D2[ii]/SLAC_H1[ii]/R-1.0;
	 SLAC_r[ii]=SLAC_D2[ii]/(SLAC_H1[ii]+SLAC_N[ii]);
//         SLAC_r_err[ii]=SLAC_r[ii]*sqrt((SLAC_D2_err[ii]/SLAC_D2[ii])*(SLAC_D2_err[ii]/SLAC_D2[ii])+SLAC_H1_err[ii]*SLAC_H1_err[ii]/((SLAC_H1[ii]+SLAC_N[ii])*(SLAC_H1[ii]+SLAC_N[ii]))+SLAC_N_err[ii]*SLAC_N_err[ii]/((SLAC_H1[ii]+SLAC_N[ii])*(SLAC_H1[ii]+SLAC_N[ii])));
         SLAC_r_err[ii]=sqrt(pow(SLAC_D2_err[ii]/(SLAC_H1[ii]*R),2)+pow((SLAC_D2[ii]*SLAC_H1_err[ii]/(SLAC_H1[ii]*SLAC_H1[ii]*R)),2));
         cout<<SLAC_r_err[ii]<<",";
     }
     cout<<endl;


     TGraphErrors *hr=new TGraphErrors(12,SLAC_x,SLAC_r,0,SLAC_r_err);
     hr->SetMarkerStyle(8);
     hr->SetMarkerColor(1);
     hr->Draw("AP");
     hr->SetTitle("sigma_D/(sigma_p+sigma_n);xbj");
}
