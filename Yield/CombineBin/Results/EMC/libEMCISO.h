Double_t H3_ISO(Double_t x){
        Double_t R=1.0;
        R=1.09527-1.73012*x+14.9157*pow(x,2)-72.0011*pow(x,3)+209.873*pow(x,4)-384.816*pow(x,5)+434.79*pow(x,6)
          -276.485*pow(x,7)+75.879*pow(x,8);
        return R;
}

Double_t He_ISO(Double_t x){
        Double_t R=1.0;
        R=1.02927-0.336485*x+2.5182*pow(x,2)-11.9978*pow(x,3)+35.3705*pow(x,4)-69.2632*pow(x,5)+87.672*pow(x,6)
          -64.3681*pow(x,7)+20.8184*pow(x,8);
        return R;
}
 
Double_t SLAC_EMC( Double_t x){
	 Double_t A=3.0;
         Double_t Cx=0.017+0.018*log(x)+0.005*log(x)*log(x);
	 Double_t alphax=-0.07+2.189*x-24.667*pow(x,2)+145.291*pow(x,3)-497.237*pow(x,4)+1013.129*pow(x,5)
			-1208.393*pow(x,6)+775.767*pow(x,7)-205.872*pow(x,8);
	 Double_t R=exp(Cx)*pow(A,alphax);
	 return R;
}

Double_t SLAC_EMC_Den(Double_t x, Double_t A, Double_t Z){
	Double_t Betax=17.9419*pow(x,0)-435.651*pow(x,1)+4657.98*pow(x,2)-28480.2*pow(x,3)+109172*pow(x,4)-270877*pow(x,5)+432968*pow(x,6)-428761*pow(x,7)+238534*pow(x,8)-56827.9*pow(x,9);
	Double_t Dx=1.08253-3.46447*x+45.0686*pow(x,2)-281.017*pow(x,3)+969.875*pow(x,4)-1948.29*pow(x,5)
                       +2265.94*pow(x,6)-1412.72*pow(x,7)+364.558*pow(x,8);

	Double_t r_H3=1.68;
	Double_t r_He=1.88;
	Double_t Re=1.0;
	if(Z==1.0)Re=sqrt(5.0/3.0)*r_H3;	
	if(Z==2.0)Re=sqrt(5.0/3.0)*r_He;	
	Double_t rho=3.0*A/(4.0*TMath::Pi()*pow(Re,3));
	//cout<<Z<<"  "<<rho<<endl;
	Double_t EMC=Dx*(1.0+Betax*rho);
	return EMC;

}
