#include "GetTrees.h"
#include "SetCut.h"

void CutFlow()
{
     TString TreeName="T";
     TString filename;
     cout<<"Input filename: ";
     cin>>filename;
     TChain* T=GetFiles(filename,TreeName);

     Double_t Nele=T->GetEntries();    
     Double_t Ngood_1=T->GetEntries(trigger2);
     Double_t Ngood_2=T->GetEntries(trigger2+CK);
     Double_t Ngood_3=T->GetEntries(trigger2+CK+Ep);
     Double_t Ngood_4=T->GetEntries(trigger2+CK+Ep+beta);
     Double_t Ngood_5=T->GetEntries(trigger2+CK+Ep+beta+ACC);
     Double_t Ngood_6=T->GetEntries(trigger2+CK+Ep+beta+ACC+VZ);
     Double_t Ngood_7=T->GetEntries(trigger2+CK+Ep+beta+ACC+VZ+TRK);


     Double_t pcut1=Ngood_1/Nele;
     Double_t pcut2=Ngood_2/Ngood_1;
     Double_t pcut3=Ngood_3/Ngood_2;
     Double_t pcut4=Ngood_4/Ngood_3;
     Double_t pcut5=Ngood_5/Ngood_4;
     Double_t pcut6=Ngood_6/Ngood_5;
     Double_t pcut7=Ngood_7/Ngood_6;

     ofstream ofile1;
     TString outfile="./Output_cutflow/"+filename+".txt";
     ofile1.open(outfile);
     ofile1<<filename<<endl;
     ofile1<<"trigger2    "<<pcut1<<endl;
     ofile1<<"CK          "<<pcut2<<endl;
     ofile1<<"E/p         "<<pcut3<<endl;
     ofile1<<"beta        "<<pcut4<<endl;
     ofile1<<"ACC         "<<pcut5<<endl;
     ofile1<<"VZ          "<<pcut6<<endl;
     ofile1<<"TRK         "<<pcut7<<endl;

     ofile1.close();


}
