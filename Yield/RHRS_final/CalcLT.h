TRI_VAR CalcLT(int run_number,int kin,int beamcut=0)
{
     TRI_VAR LT;
     LT.value=0.0;
     LT.err=0.0;
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);
     if((T->GetFile())==NULL)return LT;
     TChain *EndT=GetTree(run_number,kin,"EndRight");
     int noend=0;
     if((EndT->GetFile())==NULL)noend=1;

     Double_t Nscaler=0;
     Double_t Nmeas=0;
     Double_t Livetime=0;
     Double_t Livetime_err=0;
     if(beamcut==0){
	if(noend)Nscaler = T->GetMaximum("evRightT5");
	else Nscaler = EndT->GetMaximum("EndRightT5");
        Nmeas = T->GetEntries(trigger2);
        if(Nscaler!=0){
           Livetime=Nmeas/Nscaler;
	   Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
	}
        cout<<Livetime<<endl;
     }

     if(beamcut==1){
        T->SetBranchStatus("*",0);
        T->SetBranchStatus("evRightT5",1);     
        T->SetBranchStatus("DR.evtypebits",1); 
        T->SetBranchStatus("evRightdnew_r",1);
        T->SetBranchStatus("RightBCMev.isrenewed",1);

	Double_t evT2,evtypebits,current_dnew,isrenewed;
        Double_t beamUp[5];
        T->SetBranchAddress("evRightT5",&evT2);
        T->SetBranchAddress("DR.evtypebits",&evtypebits);
        T->SetBranchAddress("evRightdnew_r",&current_dnew);
        T->SetBranchAddress("RightBCMev.isrenewed",&isrenewed);

        Double_t maxT2;
        if(!noend)maxT2=EndT->GetMaximum("EndRightT5");
        Double_t nentries=T->GetEntries();

	Double_t lastcount=0;
        Double_t entryT2=0;
        Nscaler=0.0;Nmeas=0.0;
        for(int ii=0;ii<nentries;ii++){
	    T->GetEntry(ii);
            Int_t trigger=(Int_t)evtypebits;                
	    if((trigger>>5)&1)entryT2=entryT2+1;

	    if(isrenewed){
               current_dnew=current_dnew*dnew_gain;
               if(current_dnew>0){
                  Nscaler=Nscaler+evT2-lastcount;
                  Nmeas=Nmeas+entryT2;
               }
               entryT2=0;
               lastcount=evT2;   
            }
        }
        
        if(maxT2!=lastcount && noend==0)Nscaler=Nscaler+maxT2-lastcount;
        Livetime=Nmeas/Nscaler;
        Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
        cout<<Livetime<<endl;

     }

     LT.value=Livetime;
     LT.err=Livetime_err;
     delete T;
     return LT;
}
