  Double_t th_max=0.06;
  Double_t th_min=-0.06;
  Double_t ph_max=0.03;
  Double_t ph_min=-0.03;
  Double_t dp_max=0.040;
  Double_t dp_min=-0.040;
  Double_t vz_max=0.11;
  Double_t vz_min=-0.105;

  TCut beta = "R.tr.beta>0";
  TCut Ep = "(R.ps.e+R.sh.e)/(1000*R.gold.p)>0.7";
  TCut ACC = Form("R.tr.tg_th>%f && R.tr.tg_th<%f && R.tr.tg_ph>%f && R.tr.tg_ph<%f && R.tr.tg_dp>%f && R.tr.tg_dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
  TCut VZ = Form("R.tr.vz>%f && R.tr.vz<%f",vz_min,vz_max);
  TCut CK = "R.cer.asum_c>2000";
  TCut trigger2 = "(DR.evtypebits>>5)&1";
  TCut trigger1 = "(DR.evtypebits>>4)&1";
  TCut TRK = "R.tr.n==1";
  TCut W2="EKRxe.W2>3";
