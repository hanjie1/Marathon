void Plot()
{
     TCut TRK="L.tr.n==1";
     TCut ACC="abs(L.tr.tg_th)<0.06 && abs(L.tr.tg_ph)<0.03 && abs(L.tr.tg_dp)<0.04";
     TCut VZ="L.tr.vz<0.10 && L.tr.vz>-0.08";
     TCut T1="(DL.evtypebits>>1)&1";
     TCut T2="(DL.evtypebits>>2)&1";
     TCut beta="L.tr.beta>0";

     TCut Ep="(L.prl1.e+L.prl2.e)/(L.gold.p*1000)>0.8 && (L.prl1.e+L.prl2.e)/(L.gold.p*1000)<1.2";
     TCut Ep_pi="(L.prl1.e+L.prl2.e)/(L.gold.p*1000)<0.2";

     TCut CK="L.cer.asum_c>2200 && L.cer.asum_c<8000";
     TCut CK_pi="L.cer.asum_c<200";

     TString rootpath="/lustre19/expphy/cache/halla/triton/prod/marathon/pass2/kin1";
     TChain *T=new TChain("T");
     T->Add(Form("%s/tritium_1213*",rootpath.Data()));

     TH1F *h1=new TH1F("h1","Cherenkov sum",500,0,10000);
     TH1F *h2=new TH1F("h2","Cherenkov sum",500,0,10000);
     TH1F *h3=new TH1F("h3","Cherenkov sum",500,0,10000);

     TH1F *h4=new TH1F("h4","E/p",300,0,1.5);
     TH1F *h5=new TH1F("h5","E/p",300,0,1.5);
     TH1F *h6=new TH1F("h6","E/p",300,0,1.5);

     TCanvas *c1=new TCanvas("c1");
     T->Draw("L.cer.asum_c>>h1",TRK+beta+ACC+VZ+T1);
     T->Draw("L.cer.asum_c>>h2",TRK+beta+ACC+VZ+T2+Ep,"same");
     T->Draw("L.cer.asum_c>>h3",TRK+beta+ACC+VZ+T1+Ep_pi,"same");

     TCanvas *c2=new TCanvas("c2");
     T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>h4",TRK+beta+ACC+VZ+T1);
     T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>h5",TRK+beta+ACC+VZ+T2+CK,"same");
     T->Draw("(L.prl1.e+L.prl2.e)/(1000*L.gold.p)>>h6",TRK+beta+ACC+VZ+T1+CK_pi,"same");

     h2->SetLineColor(2);
     h3->SetLineColor(8);

     h5->SetLineColor(2);
     h6->SetLineColor(8);




}


