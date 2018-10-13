Double_t CalcCharge(int run_number,int kin)
{
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);

     Double_t gain=0.0003264;
     Double_t offset=0.1055;

     Double_t dnewr,Lclock,current_dnew;
     Double_t beamUp[5];
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evLeftdnew_r",1);
     T->SetBranchStatus("evLeftLclock",1);
     T->SetBranchStatus("LeftBCMev.BeamUp_time_v1495",1);
     T->SetBranchStatus("LeftBCMev.current_dnew",1);

     T->SetBranchAddress("evLeftdnew_r",&dnewr);
     T->SetBranchAddress("evLeftLclock",&Lclock);
     T->SetBranchAddress("LeftBCMev.BeamUp_time_v1495",beamUp);
     T->SetBranchAddress("LeftBCMev.current_dnew",&current_dnew);

     Int_t nentries=T->GetEntries();
     Double_t current,time;
     Double_t lastcount=0.0;
     Double_t charge=0.0;
     for(int ii=0;ii<nentries;ii++){
        T->GetEntry(ii);
        if(current_dnew>4){
           current = gain*dnewr+offset;
	   charge = charge+current*(Lclock-lastcount)/103700.0;
	}
	lastcount=Lclock;
     }
  
     return charge;

}
