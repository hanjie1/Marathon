  Double_t th_max=0.06;
  Double_t th_min=-0.06;
  Double_t ph_max=0.03;
  Double_t ph_min=-0.03;
  Double_t dp_max=0.040;
  Double_t dp_min=-0.040;
  Double_t vz_max[11]={0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.105,0.105};
  Double_t vz_min[11]={-0.08,-0.08,-0.08,-0.08,-0.08,-0.09,-0.09,-0.095,-0.095,-0.10,-0.10};

  TCut beta = "L.tr.beta>0";
  TCut beta_R = "R.tr.beta>0";
  TCut Ep = "(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>0.7";
  TCut ACC = Form("L.tr.tg_th>%f && L.tr.tg_th<%f && L.tr.tg_ph>%f && L.tr.tg_ph<%f && L.tr.tg_dp>%f && L.tr.tg_dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);
//  TCut CK = "L.cer.asum_c>1500";
//  TCut CK_R = "R.cer.asum_c>2000";
  TCut trigger2 = "(DL.evtypebits>>2)&1";
  TCut trigger1 = "(DL.evtypebits>>1)&1";
  TCut trigger3 = "(DL.evtypebits>>3)&1";
  TCut trigger4 = "(DR.evtypebits>>4)&1";
  TCut trigger5 = "(DR.evtypebits>>5)&1";
  TCut trigger6 = "(DR.evtypebits>>6)&1";
  TCut TRK = "L.tr.n==1";
  TCut W2 = "EKLxe.W2>3.0";

  TCut ACC_R = Form("R.tr.tg_th>%f && R.tr.tg_th<%f && R.tr.tg_ph>%f && R.tr.tg_ph<%f && R.tr.tg_dp>%f && R.tr.tg_dp<%f",th_min,th_max,ph_min,ph_max,dp_min,dp_max);

