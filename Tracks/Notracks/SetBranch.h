Double_t s2_nhit;
Double_t s2_tpads[5];
Double_t s2_time[16];
Double_t s0_time[1];
Double_t ntrk;

void SetBranch(TChain *T, int HRS)
{
 if(HRS==1){
   T->SetBranchAddress("L.tr.n",&ntrk);
   T->SetBranchAddress("L.s2.nthit",&s2_nhit);
   T->SetBranchAddress("L.s2.t_pads",s2_tpads);
   T->SetBranchAddress("L.s2.time",s2_time);
   T->SetBranchAddress("L.s0.time",s0_time);
  }
 else{
   T->SetBranchAddress("R.tr.n",&ntrk);
   T->SetBranchAddress("R.s2.nthit",&s2_nhit);
   T->SetBranchAddress("R.s2.t_pads",s2_tpads);
   T->SetBranchAddress("R.s2.time",s2_time);
   T->SetBranchAddress("R.s0.time",s0_time);
  }
}
