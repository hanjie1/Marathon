void plot_Yield()
{
     ifstream H1_file,D2_file;

     Double_t H1_yield[17]={0.0},D2_yield[17]={0.0};
     Double_t H1_yield_err[17]={0.0},D2_yield_err[17]={0.0};
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

      }

     Double_t ratio[17]={0.0};
     Double_t ratio_err[17]={0.0};
     for(int ii=0;ii<17;ii++){
	 if(H1_yield[ii]!=0){
            ratio[ii]=D2_yield[ii]/H1_yield[ii];
	    ratio_err[ii]=ratio[ii]*sqrt((D2_yield_err[ii]/D2_yield[ii])*(D2_yield_err[ii]/D2_yield[ii])+(H1_yield_err[ii]/H1_yield[ii])*(H1_yield_err[ii]/H1_yield[ii]));
         }
      }

     Double_t xbj_err[17]={0.0};
     TCanvas *c1=new TCanvas();
     TGraphErrors *hH1=new TGraphErrors(17,xbj,H1_yield,xbj_err,H1_yield_err);
     hH1->SetMarkerStyle(8);
     hH1->SetMarkerColor(1);
     hH1->Draw("AP");
     hH1->SetTitle("H1 Data Yield;xbj;nb");

     TCanvas *c2=new TCanvas();
     TGraphErrors *hD2=new TGraphErrors(17,xbj,D2_yield,xbj_err,D2_yield_err);
     hD2->SetMarkerStyle(8);
     hD2->SetMarkerColor(1);
     hD2->Draw("AP");
     hD2->SetTitle("D2 Data Yield;xbj;nb");

     TCanvas *c3=new TCanvas();
     TGraphErrors *hratio=new TGraphErrors(17,xbj,ratio,xbj_err,ratio_err);
     hratio->SetMarkerStyle(8);
     hratio->SetMarkerColor(1);
     hratio->Draw("AP");
     hratio->SetTitle("Deuterium/proton Data Yield;xbj");


}
