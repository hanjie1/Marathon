Double_t CalcCharge(int run_number,int kin)
{
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);
     TChain *EndT=GetTree(run_number,kin,"EndLeft");
     int noend=0;
     if(EndT==NULL)noend=1;

     Double_t gain=0.0003361;
     Double_t offset=0.0217;

     Double_t dnewr,dnewc,Lclock,isrenewed;
     Double_t beamUp[5];
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evLeftdnew_r",1);
     T->SetBranchStatus("evLeftdnew",1);
     T->SetBranchStatus("evLeftLclock",1);
     T->SetBranchStatus("LeftBCMev.isrenewed",1);
//     T->SetBranchStatus("LeftBCMev.BeamUp_time_v1495",1);
//     T->SetBranchStatus("LeftBCMev.current_dnew",1);

     T->SetBranchAddress("evLeftdnew_r",&dnewr);
     T->SetBranchAddress("evLeftdnew",&dnewc);
     T->SetBranchAddress("evLeftLclock",&Lclock);
     T->SetBranchAddress("LeftBCMev.isrenewed",&isrenewed);
//     T->SetBranchAddress("LeftBCMev.BeamUp_time_v1495",beamUp);
//     T->SetBranchAddress("LeftBCMev.current_dnew",&current_dnew);

     Double_t dnew_max,Lclock_max;
     if(noend==0){
        dnew_max=EndT->GetMaximum("EndLeftdnew");
	Lclock_max=EndT->GetMaximum("EndLeftLclock");
     }
     Int_t nentries=T->GetEntries();
     Double_t current,time;
     Double_t last_clock=0.0;
     Double_t last_dnew=0.0;
     Double_t charge=0.0;
     for(int ii=0;ii<nentries;ii++){
        T->GetEntry(ii);
        current = gain*dnewr+offset;
        if(isrenewed){
            if(current>4){// && current<30){
	       charge = charge+current*(Lclock-last_clock)/103700.0;
             } 
	    last_clock=Lclock;
	    last_dnew=dnewc;
        }
     }
     delete T;  

     if(last_clock!=Lclock_max && noend==0){
        Double_t tmp_t=(Lclock_max-last_clock)/103700.0;
        current = (dnew_max-dnewc)/tmp_t*gain+offset;
        if(current>0 && current<30){
             charge = charge+current*tmp_t;
             cout<<"charge:  "<<charge<<"  "<<current<<"  "<<tmp_t<<endl;
        }
     }

     return charge;

}
