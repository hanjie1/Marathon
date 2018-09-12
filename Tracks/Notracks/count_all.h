int count_all(Int_t run_number,TString &kin,Double_t ele[],Double_t nontrk[])
{
     TChain *T=new TChain("T");
      TString File=Form("/lustre/expphy/volatile/halla/triton/Marathon_Rootfiles/pass1_fix/tritium_%d.root",run_number);

     int index=1;
     if (!gSystem->AccessPathName(File)){
        T->Add(File);
        File=Form("/lustre/expphy/volatile/halla/triton/Marathon_Rootfiles/pass1_fix/tritium_%d_%d.root",run_number,index);
        while(!gSystem->AccessPathName(File)){
	    T->Add(File);
            index++;
            File=Form("/lustre/expphy/volatile/halla/triton/Marathon_Rootfiles/pass1_fix/tritium_%d_%d.root",run_number,index);
        }
     }
     else {
        return -1;
     }

  int LHRS=0,RHRS=0;
  if(run_number<20000)LHRS=1;
  else RHRS=1;

  TCut mTRK,CK,totalE,trigger2,noTRK,btrip;
  if(LHRS){
     totalE = "(L.prl1.e+L.prl2.e)/3100.0>0.85";
     CK = "L.cer.asum_c>1500";
     trigger2 = "(DL.evtypebits>>2)&1";
     mTRK = "L.tr.n>0";
     noTRK = "L.tr.n==0";
     btrip = "LeftBCMev.BeamUp_time_v1495[0]>3.0";
   }
  if(RHRS){
     totalE = "(R.ps.e+R.sh.e)/2900.0>0.85";
     CK = "R.cer.asum_c>2000";
     trigger2 = "(DR.evtypebits>>5)&1";
     mTRK = "R.tr.n>0";
     noTRK = "R.tr.n==0";
     btrip = "RightBCMev.BeamUp_time_v1495[0]>3.0";
   }

   Double_t goodele=0;
   Double_t notrkele=0; 

   goodele=T->GetEntries(trigger2+CK+totalE+btrip);
   notrkele=T->GetEntries(trigger2+CK+totalE+noTRK+btrip); 


  int change=0;

  if(kin=="kin1"){
     ele[0]=ele[0]+goodele;
     nontrk[0]=nontrk[0]+notrkele;      
     change=1;
   }
  
  if(kin=="kin2"){
     ele[1]=ele[1]+goodele;
     nontrk[1]=nontrk[1]+notrkele;      
     change=1;
   }

  if(kin=="kin3"){
     ele[2]=ele[2]+goodele;
     nontrk[2]=nontrk[2]+notrkele;      
     change=1;
   }

  if(kin=="kin4"){
     ele[3]=ele[3]+goodele;
     nontrk[3]=nontrk[3]+notrkele;      
     change=1;
   }

  if(kin=="kin5"){
     ele[4]=ele[4]+goodele;
     nontrk[4]=nontrk[4]+notrkele;      
     change=1;
   }

  if(kin=="kin7"){
     ele[5]=ele[5]+goodele;
     nontrk[5]=nontrk[5]+notrkele;      
     change=1;
   }

  if(kin=="kin9"){
     ele[6]=ele[6]+goodele;
     nontrk[6]=nontrk[6]+notrkele;      
     change=1;
   }

  if(kin=="kin11"){
     ele[7]=ele[7]+goodele;
     nontrk[7]=nontrk[7]+notrkele;      
     change=1;
   }

  if(kin=="kin13"){
     ele[8]=ele[8]+goodele;
     nontrk[8]=nontrk[8]+notrkele;      
     change=1;
   }

  if(kin=="kin15"){
     ele[9]=ele[9]+goodele;
     nontrk[9]=nontrk[9]+notrkele;      
     change=1;
   }

  if(kin=="kin16"){
     ele[10]=ele[10]+goodele;
     nontrk[10]=nontrk[10]+notrkele;      
     change=1;
   }


  delete T;
  if(change==0)return 0;
  else  return 1;


}
