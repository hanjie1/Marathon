#include "count_all.h"
void eff_all()
{
     TString He_file,H3_file;

     Double_t He_gele[11],He_notrack[11];
     Double_t H3_gele[11],H3_notrack[11];
     Double_t xbj[11]={1,2,3,4,5,7,9,11,13,15,16};

     for(int jj=0;jj<11;jj++)
      {
         He_gele[jj]=0; He_notrack[jj]=0;
         H3_gele[jj]=0; H3_notrack[jj]=0;
      }

     for(int ii=0;ii<11;ii++){
         He_file=Form("./Runlist/He3_kin%d.dat",(int)xbj[ii]);     
         H3_file=Form("./Runlist/H3_kin%d.dat",(int)xbj[ii]);  

         ifstream infile1,infile2;
	 infile1.open(He_file);
         infile2.open(H3_file);

         TString line;
         int run_number;
         int Check=0;
         int nrun=0;  
 
         TString kin=Form("kin%d",(int)xbj[ii]);
         while(line.ReadToken(infile1) ){
              run_number=line.Atoi();
         //     cout<<run_number<<endl;
              Check=count_all(run_number,kin,He_gele,He_notrack); 
              if(Check==-1){
                 cout<<run_number<<" rootfile couldn't find"<<endl;
                 continue;
              }
             if(Check==0){
                cout<<run_number<<" can't find Kinematics setting"<<endl;
                continue;
              }
             nrun++;
	 }
         infile1.close();
         cout<<"He3 "<<kin<<" finished"<<endl;

         nrun=0;
         while(line.ReadToken(infile2) ){
              run_number=line.Atoi();
              Check=count_all(run_number,kin,H3_gele,H3_notrack); 
              if(Check==-1){
                 cout<<run_number<<" rootfile couldn't find"<<endl;
                 continue;
              }
             if(Check==0){
                cout<<run_number<<" can't find Kinematics setting"<<endl;
                continue;
              }
             nrun++;
	 }
         infile2.close();
         cout<<"H3 "<<kin<<" finished"<<endl;

     }


    Double_t He_eff[11],H3_eff[11];
    Double_t He_eff_err[11],H3_eff_err[11];

    for(int jj=0;jj<11;jj++){
	He_eff[jj]=0;        H3_eff[jj]=0;
        He_eff_err[jj]=0;    H3_eff_err[jj]=0;
    }

    for(int jj=0;jj<11;jj++){
        if(He_gele[jj]==0)He_eff[jj]=0;
        else{
	     He_eff[jj]=He_notrack[jj]/He_gele[jj];
             He_eff_err[jj]=He_eff[jj]*sqrt(1.0/He_notrack[jj]-1.0/He_gele[jj]);
        }

        if(H3_gele[jj]==0)H3_eff[jj]=0;
        else{
             H3_eff[jj]=H3_notrack[jj]/H3_gele[jj];
             H3_eff_err[jj]=H3_eff[jj]*sqrt(1.0/H3_notrack[jj]-1.0/H3_gele[jj]);
        }
    }

   Double_t xbj_err[11]={0};
   TGraphErrors *hHe=new TGraphErrors(11,xbj,He_eff,xbj_err,He_eff_err);
   TGraphErrors *hH3=new TGraphErrors(11,xbj,H3_eff,xbj_err,H3_eff_err);

   TCanvas *c1=new TCanvas();
   TMultiGraph *mg=new TMultiGraph();
   mg->Add(hHe);
   hHe->SetMarkerColor(1);
   hHe->SetMarkerStyle(8);
   mg->Add(hH3);
   hH3->SetMarkerColor(2);
   hH3->SetMarkerStyle(8);
   mg->Draw("AP");

   auto leg = new TLegend(0.1,0.75,0.25,0.90);
   leg->AddEntry(hHe,"He3","P");
   leg->AddEntry(hH3,"Tritium","P");
   leg->Draw();

   c1->Print("H3_He_withcosmic.png");

   ofstream myfile;
   myfile.open("H3_He_withcosmic.txt");
   for(int jj=0;jj<11;jj++){
      myfile<<xbj[jj]<<"  "<<H3_eff[jj]<<"  "<<H3_eff_err[jj]<<endl;
    }
   for(int jj=0;jj<11;jj++){
      myfile<<xbj[jj]<<"  "<<He_eff[jj]<<"  "<<He_eff_err[jj]<<endl;
    }
   myfile.close();

}
