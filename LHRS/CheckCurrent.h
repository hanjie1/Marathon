int CheckCurrent(const int run_number,int kin)
{
     TString TreeName="T";
     TChain *T=GetTree(run_number,kin,TreeName);
     cout<<run_number<<endl;
     Double_t Current[10];
     for(int ii=0;ii<10;ii++){
	 Current[ii]=0;
     }
     TCut timeup="LeftBCMev.BeamUp_time_v1495[0]>60";
     T->Draw(">>beamup",timeup);
     TEventList *beamup=(TEventList*)gDirectory->Get("beamup");
     T->SetEventList(beamup);
  
     T->SetBranchStatus("*",0);
     T->SetBranchStatus("LeftBCMev.current_dnew",1);

     Double_t c_dnew;
     T->SetBranchAddress("LeftBCMev.current_dnew",&c_dnew);

     Double_t lastcurrent=0;
     int index=0;
     Int_t nentries=beamup->GetN();
     for(int ii=0;ii<nentries;ii++){
        T->GetEntry(beamup->GetEntry(ii));
        Double_t tmp=c_dnew-lastcurrent;
        if(abs(tmp)>2.0){
           Double_t curr=c_dnew;
           int mm=0;
           for(int jj=0;jj<10;jj++){
		if(abs(Current[jj]-curr)<=2){
                   mm=1;
                   break;}
           }            
           lastcurrent=c_dnew;
           if(mm==1)continue;
           if(index<10)Current[index++]=curr;
           else{cout<<"Run "<<run_number<<" more than 10 currents!"<<endl;
                break;}
        }
     }
     delete T;

     ofstream myfile;
     myfile.open("current.txt",fstream::app); 
     myfile<<run_number<<"  "<<index<<"  ";
     for(int ii=0;ii<10;ii++){
	 if(Current[ii]==0)continue;
         myfile<<fixed<<setprecision(1)<<Current[ii]<<"  ";
     }
     myfile<<endl;
     myfile.close();
     return index;
}
