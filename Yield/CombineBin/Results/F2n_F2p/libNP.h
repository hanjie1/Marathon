Double_t SLAC_NP(Double_t x){
	Double_t np=1.0-0.8*x;
	return np;
} 

Double_t NMC_NP(Double_t x, Double_t Q2){
	Double_t Ax=0.979-1.692*x+2.797*pow(x,2)-4.313*pow(x,3)+3.075*pow(x,4);
	Double_t Bx=-0.171*x+0.244*pow(x,2);
	Double_t np=Ax*pow(Q2/20.0,Bx)*(1.0+x*x/Q2);
	return np;
} 
