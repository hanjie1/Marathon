int CalcBoilingCorr(int run_number, int kin, TString target, Double_t &Ibeam, Double_t &dbCorr, Double_t &bCorr)
{
     Ibeam=0.0, dbCorr=0.0, bCorr=1.0;
     
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);

     Double_t dnewr;
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evLeftdnew",1);
     T->SetBranchAddress("evLeftdnew_r",&dnewr);

     Int_t nentries=T->GetEntries();
     Double_t totalI=0;
     Int_t nI=0;
     for(int ii=0;ii<nentries;ii++)
      {
         T->GetEntry(ii);
         Double_t tmpI = dnew_gain*dnewr+dnew_offset;
         if(tmpI>dnew_offset && tmpI<30){
            totalI+=tmpI;
            nI++;
         }
      }

     delete T;

     Double_t avgI=0;
     if(nI!=0)avgI=totalI/nI;

     Double_t boiling_corr=1.0;
     Double_t boiling_err=0.0;
     Double_t massA=1.0;
     if(target=="H1") {
        boiling_corr=bH1_A*avgI*avgI+bH1_B*avgI+1.0;
	boiling_err=bH1_VC+pow(avgI,2)*bH1_VB+pow(avgI,4)*bH1_VA
		    +2.0*avgI*bH1_COV_BC+2.0*pow(avgI,2)*bH1_COV_AC+2.0*pow(avgI,3)*bH1_COV_AB;
     }
     if(target=="D2"){
        boiling_corr=bD2_A*avgI*avgI+bD2_B*avgI+1.0;
	boiling_err=bD2_VC+pow(avgI,2)*bD2_VB+pow(avgI,4)*bD2_VA
		    +2.0*avgI*bD2_COV_BC+2.0*pow(avgI,2)*bD2_COV_AC+2.0*pow(avgI,3)*bD2_COV_AB;
     }
     if(target=="He3"){
        boiling_corr=bHe3_A*avgI*avgI+bHe3_B*avgI+1.0;
	boiling_err=bHe3_VC+pow(avgI,2)*bHe3_VB+pow(avgI,4)*bHe3_VA
		    +2.0*avgI*bHe3_COV_BC+2.0*pow(avgI,2)*bHe3_COV_AC+2.0*pow(avgI,3)*bHe3_COV_AB;
     }
     if(target=="H3"){
        boiling_corr=bH3_A*avgI*avgI+bH3_B*avgI+1.0;
	boiling_err=bH3_VC+pow(avgI,2)*bH3_VB+pow(avgI,4)*bH3_VA
		    +2.0*avgI*bH3_COV_BC+2.0*pow(avgI,2)*bH3_COV_AC+2.0*pow(avgI,3)*bH3_COV_AB;
     }

     Ibeam=avgI;
     bCorr=boiling_corr;
     dbCorr=boiling_err; //square of the uncertainty

     return 1;
}



