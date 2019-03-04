void RunLum(int run_number,int kin,Double_t& Charge,Double_t& Ntarg,TString target)
{
     Charge=0;
     Ntarg=0;
 
     Double_t aCharge=CalcCharge(run_number,kin);
      
     Charge=aCharge/(Qe*1e6);

     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);

     Double_t dnewr;
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evRightdnew",1);

     T->SetBranchAddress("evRightdnew_r",&dnewr);
 
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
     
     Double_t avgI=0;
     if(nI!=0)avgI=totalI/nI;
    
     Double_t boiling_corr=1;
     Double_t massA=1.0;
     if(target=="H1") {
        boiling_corr=bH1_A*avgI*avgI+bH1_B*avgI+1.0;
        massA=1.0;
        Ntarg=H1_pho*boiling_corr*Na/massA;
     }
     if(target=="D2"){
        boiling_corr=bD2_A*avgI*avgI+bD2_B*avgI+1.0;
        massA=2.0;
        Ntarg=D2_pho*boiling_corr*Na/massA;
     }
     if(target=="He3"){
        boiling_corr=bHe3_A*avgI*avgI+bHe3_B*avgI+1.0;
        massA=3.0;
        Ntarg=He3_pho*boiling_corr*Na/massA;
     }
     if(target=="H3"){
        boiling_corr=bH3_A*avgI*avgI+bH3_B*avgI+1.0;
        massA=3.0;
        Ntarg=H3_pho*boiling_corr*Na/massA;
     }
     return;
}

Double_t CalcLum(int kin, TString target){
     int nrun=0;
     vector<Int_t> runList;
     nrun=GetRunList(runList,kin,target);
     if(nrun==0)exit(0);

     int run_number,success=0;
     Ssiz_t from=0;
     TString content;
     Double_t NNtarg=0.0,Ncharge=0.0;
     Double_t LUM=0.0;
     for(int ii=0;ii<runList.size();ii++){
         run_number = runList[ii];
         RunLum(run_number,kin,Ncharge,NNtarg,target);
         LUM+=Ncharge*NNtarg/CM2toNB;
    }
    
    return LUM;

}






