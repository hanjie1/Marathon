void plot_RC()
{
     ifstream H1_file,D2_file;
     ifstream H1_fileR,D2_fileR;

     Double_t H1_yield[17]={0.0},D2_yield[17]={0.0};
     Double_t H1_Q2[17]={0.0},D2_Q2[17]={0.0};
     Double_t H1_nQ2[17]={0.0},D2_nQ2[17]={0.0};
     Double_t H1_RC[17]={0.0},D2_RC[17]={0.0};
     Double_t H1_yield_err[17]={0.0},D2_yield_err[17]={0.0};
     Double_t H1_RC_err[17]={0.0},D2_RC_err[17]={0.0};
     Double_t xbj[17]={0.0};
     for(int ii=0;ii<17;ii++)
         xbj[ii]=0.17+ii*0.02;
     
     for(int ii=0;ii<4;ii++)
      {
         H1_file.open(Form("../RawYield/H1_kin%d.txt",ii+1));
         Ssiz_t from=0;
         TString content,tmp;
         int nn=0,nnx=0; 
         while(tmp.ReadLine(H1_file)){
	     if(nn==0){nn++;continue;}	
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             H1_yield[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_yield_err[nn-1]+=atof(content.Data());
             nn++;
             from=0;
         }
	H1_file.close();

         H1_fileR.open(Form("../RC_Yield/H1_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(H1_fileR)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             H1_RC[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_RC_err[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             H1_Q2[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)H1_nQ2[nn-1]+=1.0;
             nn++;
             from=0;
         }
        H1_fileR.close();


         D2_file.open(Form("../RawYield/D2_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(D2_file)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             D2_yield[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_yield_err[nn-1]+=atof(content.Data());
             nn++;
             from=0;
         }
        D2_file.close();

         D2_fileR.open(Form("../RC_Yield/D2_kin%d.txt",ii+1));
         from=0;
         nn=0,nnx=0;
         while(tmp.ReadLine(D2_fileR)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from," ");
             tmp.Tokenize(content,from," ");
             D2_RC[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_RC_err[nn-1]+=atof(content.Data());
             tmp.Tokenize(content,from," ");
             D2_Q2[nn-1]+=atof(content.Data());
             if(atof(content.Data())!=0)D2_nQ2[nn-1]+=1.0;
             nn++;
             from=0;
         }
        D2_fileR.close();


      }

     Double_t ratio[17]={0.0},RC_ratio[17]={0.0};
     Double_t ratio_err[17]={0.0},RC_ratio_err[17]={0.0};
     for(int ii=0;ii<17;ii++){
	 if(H1_yield[ii]!=0){
            ratio[ii]=D2_yield[ii]/H1_yield[ii];
	    ratio_err[ii]=ratio[ii]*sqrt((D2_yield_err[ii]/D2_yield[ii])*(D2_yield_err[ii]/D2_yield[ii])+(H1_yield_err[ii]/H1_yield[ii])*(H1_yield_err[ii]/H1_yield[ii]));
         }
         if(H1_RC[ii]!=0){
            RC_ratio[ii]=D2_RC[ii]/H1_RC[ii];
            RC_ratio_err[ii]=ratio[ii]*sqrt((D2_RC_err[ii]/D2_RC[ii])*(D2_RC_err[ii]/D2_RC[ii])+(H1_RC_err[ii]/H1_RC[ii])*(H1_RC_err[ii]/H1_RC[ii]));
         }
         if(H1_nQ2[ii]!=0)H1_Q2[ii]=H1_Q2[ii]/H1_nQ2[ii];
         if(D2_nQ2[ii]!=0)D2_Q2[ii]=D2_Q2[ii]/D2_nQ2[ii];
         //cout<<fixed<<setprecision(3)<<endl;
         //cout<<xbj[ii]<<" &  "<<H1_Q2[ii]<<" &  "<<D2_Q2[ii]<<" &  "<<ratio[ii]<<" & "<<RC_ratio[ii]<<endl;
         cout<<D2_Q2[ii]<<",";
      }
     cout<<endl;

     Double_t xbj_err[17]={0.0};
     TCanvas *c1=new TCanvas();
     TGraphErrors *hH1=new TGraphErrors(17,xbj,H1_yield,xbj_err,H1_yield_err);
     TGraphErrors *hH1_RC=new TGraphErrors(17,xbj,H1_RC,xbj_err,H1_RC_err);
     TMultiGraph *mg_H1=new TMultiGraph();
     hH1->SetMarkerStyle(8);
     hH1->SetMarkerColor(1);
     hH1_RC->SetMarkerStyle(8);
     hH1_RC->SetMarkerColor(2);
     mg_H1->Add(hH1);
     mg_H1->Add(hH1_RC);
     mg_H1->Draw("AP");
     mg_H1->SetTitle("Hydrogen Yield;xbj;nb");

     auto leg1 = new TLegend(0.7,0.7,0.85,0.85);
     leg1->AddEntry(hH1,"before RC","P");
     leg1->AddEntry(hH1_RC,"after RC","P");
     leg1->Draw(); 

     TCanvas *c3=new TCanvas();
     TGraphErrors *hD2=new TGraphErrors(17,xbj,D2_yield,xbj_err,D2_yield_err);
     TGraphErrors *hD2_RC=new TGraphErrors(17,xbj,D2_RC,xbj_err,D2_RC_err);
     TMultiGraph *mg_D2=new TMultiGraph();
     hD2->SetMarkerStyle(8);
     hD2->SetMarkerColor(1);
     hD2_RC->SetMarkerStyle(8);
     hD2_RC->SetMarkerColor(2);
     mg_D2->Add(hD2);
     mg_D2->Add(hD2_RC);
     mg_D2->Draw("AP");
     mg_D2->SetTitle("Deuterium Yield;xbj;nb");

     auto leg2 = new TLegend(0.7,0.7,0.85,0.85);
     leg2->AddEntry(hD2,"before RC","P");
     leg2->AddEntry(hD2_RC,"after RC","P");
     leg2->Draw(); 

     Double_t SLAC_x[13]={0.1722,0.1898,0.2044,0.2206,0.2407,0.2585,0.276,0.293,0.312,0.3305,0.357,0.3852,0.4151};
     Double_t SLAC_Q2[13]={2.5482,2.7419,2.8948,3.0579,3.2515,3.4146,3.5675,3.7102,3.8631,4.0058,4.1995,4.3931,4.5868};
     Double_t SLAC_D2[13]={11700,11700,11300,11000,10300,9850,9490,8970,8590,8210,7390,6900,6480};
     Double_t SLAC_H1[13]={6520,6670,6670,6420,6330,5860,5550,5340,5080,5030,4770,4510,4020};

     Double_t SLAC_ratio[13];
     for(int ii=0;ii<13;ii++){
	 SLAC_ratio[ii]=SLAC_D2[ii]/SLAC_H1[ii];
     //    cout<<fixed<<setprecision(3);
      //   cout<<SLAC_x[ii]<<" & "<<SLAC_Q2[ii]<<" & "<<SLAC_ratio[ii]<<endl;
     }


     TCanvas *c2=new TCanvas();
     TGraphErrors *hratio=new TGraphErrors(17,xbj,ratio,xbj_err,ratio_err);
     TGraphErrors *hratio_RC=new TGraphErrors(17,xbj,RC_ratio,xbj_err,RC_ratio_err);
     TGraph *hSLAC=new TGraph(13,SLAC_x,SLAC_ratio);
     TMultiGraph *mg_ratio=new TMultiGraph();
     hratio->SetMarkerStyle(8);
     hratio->SetMarkerColor(1);
     hratio_RC->SetMarkerStyle(8);
     hratio_RC->SetMarkerColor(2);
     hSLAC->SetMarkerStyle(8);
     hSLAC->SetMarkerColor(4);
     mg_ratio->Add(hratio);
     mg_ratio->Add(hratio_RC);
     mg_ratio->Add(hSLAC);
     mg_ratio->Draw("AP");
     mg_ratio->SetTitle("Deuterium/proton Data Yield;xbj");

     auto leg3 = new TLegend(0.7,0.7,0.85,0.85);
     leg3->AddEntry(hratio,"before RC","P");
     leg3->AddEntry(hratio_RC,"after RC","P");
     leg3->AddEntry(hSLAC,"SLAC","P");
     leg3->Draw(); 

}
