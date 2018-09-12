void plot_s2Lped()
{
     ifstream file1;
     file1.open("../OUT/Ls2L_ped.dat");

     Double_t nrun[1400]={0.0};
     Double_t s2L1[1400]={0.0};
     Double_t s2L2[1400]={0.0};
     Double_t s2L3[1400]={0.0};
     Double_t s2L4[1400]={0.0};
     Double_t s2L5[1400]={0.0};
     Double_t s2L6[1400]={0.0};
     Double_t s2L7[1400]={0.0};
     Double_t s2L8[1400]={0.0};
     Double_t s2L9[1400]={0.0};
     Double_t s2L10[1400]={0.0};
     Double_t s2L11[1400]={0.0};
     Double_t s2L12[1400]={0.0};
     Double_t s2L13[1400]={0.0};
     Double_t s2L14[1400]={0.0};
     Double_t s2L15[1400]={0.0};
     Double_t s2L16[1400]={0.0};

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0,total=0;
     while(tmp.ReadLine(file1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from,"  ");
          nrun[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L3[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L4[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L5[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L6[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L7[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L8[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L9[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L10[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L11[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L12[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L13[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L14[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L15[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          s2L16[nn-1]=atof(content.Data());
          nn++;
          total++;
          from=0;
     }
     file1.close();

     TGraph *hs2L0=new TGraph();
     TGraph *hs2L1=new TGraph();
     TGraph *hs2L2=new TGraph();
     TGraph *hs2L3=new TGraph();
     TGraph *hs2L4=new TGraph();
     TGraph *hs2L5=new TGraph();
     TGraph *hs2L6=new TGraph();
     TGraph *hs2L7=new TGraph();
     TGraph *hs2L8=new TGraph();
     TGraph *hs2L9=new TGraph();
     TGraph *hs2L10=new TGraph();
     TGraph *hs2L11=new TGraph();
     TGraph *hs2L12=new TGraph();
     TGraph *hs2L13=new TGraph();
     TGraph *hs2L14=new TGraph();
     TGraph *hs2L15=new TGraph();

     int mm=0;
     for(int ii=0;ii<total;ii++){
         hs2L0->SetPoint(mm,nrun[ii],s2L1[ii]);
         hs2L1->SetPoint(mm,nrun[ii],s2L2[ii]);
         hs2L2->SetPoint(mm,nrun[ii],s2L3[ii]);
         hs2L3->SetPoint(mm,nrun[ii],s2L4[ii]);
         hs2L4->SetPoint(mm,nrun[ii],s2L5[ii]);
         hs2L5->SetPoint(mm,nrun[ii],s2L6[ii]);
         hs2L6->SetPoint(mm,nrun[ii],s2L7[ii]);
         hs2L7->SetPoint(mm,nrun[ii],s2L8[ii]);
         hs2L8->SetPoint(mm,nrun[ii],s2L9[ii]);
         hs2L9->SetPoint(mm,nrun[ii],s2L10[ii]);
         hs2L10->SetPoint(mm,nrun[ii],s2L11[ii]);
         hs2L11->SetPoint(mm,nrun[ii],s2L12[ii]);
         hs2L12->SetPoint(mm,nrun[ii],s2L13[ii]);
         hs2L13->SetPoint(mm,nrun[ii],s2L14[ii]);
         hs2L14->SetPoint(mm,nrun[ii],s2L15[ii]);
         hs2L15->SetPoint(mm,nrun[ii],s2L16[ii]);
         mm++;
     }
     TCanvas *c1=new TCanvas("c1");
     hs2L0->SetMarkerStyle(8);
     hs2L0->SetMarkerColor(4);
     hs2L0->Draw("AP");
     hs2L0->SetTitle("LHRS s2L1;run number;pedestal");

     TCanvas *c2=new TCanvas("c2");
     hs2L1->SetMarkerStyle(8);
     hs2L1->SetMarkerColor(4);
     hs2L1->Draw("AP");
     hs2L1->SetTitle("LHRS s2L2;run number;pedestal");

     TCanvas *c3=new TCanvas("c3");
     hs2L2->SetMarkerStyle(8);
     hs2L2->SetMarkerColor(4);
     hs2L2->Draw("AP");
     hs2L2->SetTitle("LHRS s2L3;run number;pedestal");

     TCanvas *c4=new TCanvas("c4");
     hs2L3->SetMarkerStyle(8);
     hs2L3->SetMarkerColor(4);
     hs2L3->Draw("AP");
     hs2L3->SetTitle("LHRS s2L4;run number;pedestal");

     TCanvas *c5=new TCanvas("c5");
     hs2L4->SetMarkerStyle(8);
     hs2L4->SetMarkerColor(4);
     hs2L4->Draw("AP");
     hs2L4->SetTitle("LHRS s2L5;run number;pedestal");

     TCanvas *c6=new TCanvas("c6");
     hs2L5->SetMarkerStyle(8);
     hs2L5->SetMarkerColor(4);
     hs2L5->Draw("AP");
     hs2L5->SetTitle("LHRS s2L6;run number;pedestal");

     TCanvas *c7=new TCanvas("c7");
     hs2L6->SetMarkerStyle(8);
     hs2L6->SetMarkerColor(4);
     hs2L6->Draw("AP");
     hs2L6->SetTitle("LHRS s2L7;run number;pedestal");

     TCanvas *c8=new TCanvas("c8");
     hs2L7->SetMarkerStyle(8);
     hs2L7->SetMarkerColor(4);
     hs2L7->Draw("AP");
     hs2L7->SetTitle("LHRS s2L8;run number;pedestal");

     TCanvas *c9=new TCanvas("c9");
     hs2L8->SetMarkerStyle(8);
     hs2L8->SetMarkerColor(4);
     hs2L8->Draw("AP");
     hs2L8->SetTitle("LHRS s2L9;run number;pedestal");

     TCanvas *c10=new TCanvas("c10");
     hs2L9->SetMarkerStyle(8);
     hs2L9->SetMarkerColor(4);
     hs2L9->Draw("AP");
     hs2L9->SetTitle("LHRS s2L10;run number;pedestal");

     TCanvas *c11=new TCanvas("c11");
     hs2L10->SetMarkerStyle(8);
     hs2L10->SetMarkerColor(4);
     hs2L10->Draw("AP");
     hs2L10->SetTitle("LHRS s2L11;run number;pedestal");

     TCanvas *c12=new TCanvas("c12");
     hs2L11->SetMarkerStyle(8);
     hs2L11->SetMarkerColor(4);
     hs2L11->Draw("AP");
     hs2L11->SetTitle("LHRS s2L12;run number;pedestal");

     TCanvas *c13=new TCanvas("c13");
     hs2L12->SetMarkerStyle(8);
     hs2L12->SetMarkerColor(4);
     hs2L12->Draw("AP");
     hs2L12->SetTitle("LHRS s2L13;run number;pedestal");

     TCanvas *c14=new TCanvas("c14");
     hs2L13->SetMarkerStyle(8);
     hs2L13->SetMarkerColor(4);
     hs2L13->Draw("AP");
     hs2L13->SetTitle("LHRS s2L14;run number;pedestal");

     TCanvas *c15=new TCanvas("c15");
     hs2L14->SetMarkerStyle(8);
     hs2L14->SetMarkerColor(4);
     hs2L14->Draw("AP");
     hs2L14->SetTitle("LHRS s2L15;run number;pedestal");

     TCanvas *c16=new TCanvas("c16");
     hs2L15->SetMarkerStyle(8);
     hs2L15->SetMarkerColor(4);
     hs2L15->Draw("AP");
     hs2L15->SetTitle("LHRS s2L16;run number;pedestal");

     c1->Print("Ls2L_ped.pdf[");
     c1->Print("Ls2L_ped.pdf");
     c2->Print("Ls2L_ped.pdf");
     c3->Print("Ls2L_ped.pdf");
     c4->Print("Ls2L_ped.pdf");
     c5->Print("Ls2L_ped.pdf");
     c6->Print("Ls2L_ped.pdf");
     c7->Print("Ls2L_ped.pdf");
     c8->Print("Ls2L_ped.pdf");
     c9->Print("Ls2L_ped.pdf");
     c10->Print("Ls2L_ped.pdf");
     c11->Print("Ls2L_ped.pdf");
     c12->Print("Ls2L_ped.pdf");
     c13->Print("Ls2L_ped.pdf");
     c14->Print("Ls2L_ped.pdf");
     c15->Print("Ls2L_ped.pdf");
     c16->Print("Ls2L_ped.pdf");
     c16->Print("Ls2L_ped.pdf]");






}


