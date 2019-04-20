int SearchACC(Double_t dTheta, Double_t dEp, Double_t ACC_th[nTh], Double_t ACC_table[nTh][nEp], int ACC_matrix[nTh][nEp]){
    int aACC=0;
    int found=0;

    for(int ii=0;ii<nTh;ii++){
	Double_t tmpTh=dTheta-ACC_th[ii];
        if(tmpTh>=deltaTh || tmpTh<0)continue;
	for(int jj=0;jj<nEp;jj++){
           Double_t tmpEp=dEp-ACC_table[ii][jj];
	   if(tmpEp>=deltaEp || tmpEp<0)continue;
	   found=1;
	   aACC=ACC_matrix[ii][jj];
	   break;
	}
        if(found==1)break;
    }

    if(found==0)cout<<"No ACC matrix found !!  "<<dTheta<<"  "<<dEp<<endl;
    return aACC;

}
