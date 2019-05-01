Double_t CalcCharge(int run_number,int kin)
{
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);
     if((T->GetFile())==NULL)return 0;
     TChain *EndT=GetTree(run_number,kin,"EndRight");
     int noend=0;
     if((EndT->GetFile())==NULL)noend=1;

     Double_t dnewr,dnewc,Lclock,isrenewed;
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("evRightdnew_r",1);
     T->SetBranchStatus("evRightdnew",1);
     T->SetBranchStatus("evRightLclock",1);
     T->SetBranchStatus("RightBCMev.isrenewed",1);

     T->SetBranchAddress("evRightdnew_r",&dnewr);
     T->SetBranchAddress("evRightdnew",&dnewc);
     T->SetBranchAddress("evRightLclock",&Lclock);
     T->SetBranchAddress("RightBCMev.isrenewed",&isrenewed);

     Double_t dnew_max,Lclock_max;
     if(noend==0){
        dnew_max=EndT->GetMaximum("EndRightdnew");
	Lclock_max=EndT->GetMaximum("EndRightLclock");
     }
     Int_t nentries=T->GetEntries();
     Double_t current,time;
     Double_t last_clock=0.0;
     Double_t last_dnew=0.0;
     Double_t charge=0.0;
     for(int ii=0;ii<nentries;ii++){
        T->GetEntry(ii);
        if(isrenewed){
            current = dnew_gain*dnewr+dnew_offset;
            if(current>dnew_offset && current<30){
	       charge = charge+current*(Lclock-last_clock)/103700.0;
             } 
	    last_clock=Lclock;
	    last_dnew=dnewc;
        }
     }
     delete T;  

     if(last_clock!=Lclock_max && noend==0){
        Double_t tmp_t=(Lclock_max-last_clock)/103700.0;
        current = (dnew_max-dnewc)/tmp_t*dnew_gain+dnew_offset;
        if(current>dnew_offset && current<30){
             charge = charge+current*tmp_t;
        }
     }
     cout<<"charge:  "<<charge<<endl;

     return charge;

}
