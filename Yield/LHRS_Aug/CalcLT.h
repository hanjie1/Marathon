TRI_VAR CalcLT(int run_number,int kin,int beamcut=0)
{
     TString TreeName="T";
     TChain* T=GetTree(run_number,kin,TreeName);
     TChain *EndT=GetTree(run_number,kin,"EndLeft");

     Double_t Nscaler=0;
     Double_t Nmeas=0;
     Double_t Livetime=0;
     Double_t Livetime_err=0;
     TRI_VAR LT;
     if(beamcut==0){
	Nscaler = EndT->GetMaximum("EndLeftT2");
        Nmeas = T->GetEntries(trigger2);
        if(Nscaler!=0){
           Livetime=Nmeas/Nscaler;
	   Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
	}
        cout<<Livetime<<endl;
     }

     if(beamcut==1){
        T->SetBranchStatus("*",0);
        T->SetBranchStatus("evLeftT2",1);     
        T->SetBranchStatus("LeftBCMev.BeamUp_time_v1495",1);
        T->SetBranchStatus("DL.evtypebits",1); 
        T->SetBranchStatus("LeftBCMev.current_dnew",1);
        T->SetBranchStatus("LeftBCMev.isrenewed",1);

	Double_t evT2,evtypebits,current_dnew,isrenewed;
        Double_t beamUp[5];
        T->SetBranchAddress("evLeftT2",&evT2);
        T->SetBranchAddress("DL.evtypebits",&evtypebits);
        T->SetBranchAddress("LeftBCMev.BeamUp_time_v1495",beamUp);
        T->SetBranchAddress("LeftBCMev.current_dnew",&current_dnew);
        T->SetBranchAddress("LeftBCMev.isrenewed",&isrenewed);

        Double_t maxT2=EndT->GetMaximum("EndLeftT2");
        Double_t nentries=T->GetEntries();

	Double_t lastcount=0;
        Double_t entryT2=0;
        Nscaler=0.0;Nmeas=0.0;
        for(int ii=0;ii<nentries;ii++){
	    T->GetEntry(ii);
	    if(isrenewed){
               if(current_dnew>4){
                  Nscaler=Nscaler+evT2-lastcount;
                  Nmeas=Nmeas+entryT2;
               }
               entryT2=0;
               lastcount=evT2;   
            }
	
            Int_t trigger=(Int_t)evtypebits;                
	    if((trigger>>2)&1)entryT2=entryT2+1;
        }
        
        if(maxT2!=lastcount)Nscaler=Nscaler+maxT2-lastcount;
        Livetime=Nmeas/Nscaler;
        Livetime_err=sqrt(1.0/Nmeas+1.0/Nscaler)*Livetime;
        cout<<Livetime<<endl;

     }

     LT.value=Livetime;
     LT.err=Livetime_err;
     return LT;
}
