void tria2()
{
    //bool isNormalHierarchy = true; // ly luan
    double sinsq12=0.307; //sin squared 12
    
    double sinsq23N=0.386; //sin squared 23 Normal 
    double sinsq23I=0.392; //sin squared 23 invert
    //double sinsq23 = isNormalHierarchy?sinsq23N:sinsq23I; // 
    
    double sinsq13N=0.0241; 
    double sinsq13I=0.0244;
    //double sinsq13 = isNormalHierarchy?sinsq13N:sinsq13I;
    
    double deltaN=1.08*TMath::Pi(); //delta_{CP}
    double deltaI=1.09*TMath::Pi();
    //double delta = isNormalHierarchy?deltaN:deltaI;   //  
    
cout<<"The global analysis input" <<endl; //xuat du lieu
cout<<"sinsq12: "<<sinsq12<<endl;
cout<<"sinsq23: Normal "<<sinsq23N<<" Invert "<<sinsq23I<<endl;
cout<<"sinsq13: Normal "<<sinsq13N<<" Invert "<<sinsq13I<<endl;
cout<<"deltacp: Normal "<<deltaN<<" Invert "<<deltaI<<endl;
    
    
    //Draw a triangle;
    TCanvas* c1 = new TCanvas("c1","",3000,5000,5000,3000);
  //  c1 = new TCanvas("c1");
    c1-> Divide(2,1);
    gStyle->SetOptStat(0); //tat bang thong ke 
    c1->Range(-1.0,-1.0,1.0,1.0);
    
    for (int i = 1; i<3; i++)
  {  
    if(i==1) 
    { 
		cout<<"Normal hierarchy" <<endl; //if is normal, output value of s12,s23,s13,... and vice versa
		// else cout<<"Invert hierarchy" <<endl;
		c1-> cd(1);
    
		//write in sin and cos
		//calculate s & c
		double s12=TMath::Sqrt(sinsq12); 
		double c12=TMath::Sqrt(1-sinsq12);
    
		double s23=TMath::Sqrt(sinsq23N);
		double c23=TMath::Sqrt(1-sinsq23N);
    
		double s13=TMath::Sqrt(sinsq13N);
		double c13=TMath::Sqrt(1-sinsq13N);
    
		//using TComplex for CP phase 
		double scp=TMath::Sin(deltaN);
		double ccp=TMath::Cos(deltaN); 
		cpPhasePos= TComplex(ccp,scp); //xd complex number e^{i*delta}
		cpPhaseNeg= TComplex(ccp,-scp);  //negative e^{-i*delta}			
	}
    
   
    
	if(i==2) 
	{
		cout<<"Invert hierarchy" <<endl;
		c1-> cd(2);
		//write in sin and cos
		//calculate s & c
		double s12=TMath::Sqrt(sinsq12); 
		double c12=TMath::Sqrt(1-sinsq12);
    
		double s23=TMath::Sqrt(sinsq23I);
		double c23=TMath::Sqrt(1-sinsq23I);
		
		double s13=TMath::Sqrt(sinsq13I);
		double c13=TMath::Sqrt(1-sinsq13I);
    
		//using TComplex for CP phase 
		double scp=TMath::Sin(deltaI);
		double ccp=TMath::Cos(deltaI);
		cpPhasePos= TComplex(ccp,scp); //xd complex number e^{i*delta}
		cpPhaseNeg= TComplex(ccp,-scp);  //negative e^{-i*delta}	
	}  
  
    cout<<"s12: "<<s12<<" c12 "<<c12<<endl;
    cout<<"s23: "<<s23<<" c23 "<<c23<<endl;
    cout<<"s13: "<<s13<<" c13 "<<c13<<endl;
	cout<<"scp: "<<cpPhasePos.Re()<<" ccp "<<cpPhasePos.Im()<<endl; //output Real and Im in e^delta_{CP} = cos + i*sin 			 
   
    
    //PMNS element, the value
    //electron
    double Ue1=c12*c13; //Ue1 = c12*(1-s13)
    cout<<"Ue1 "<<Ue1<<endl; 
    double Ue2=s12*c13; //Ue2 = s12*(1-s13)
    cout<<"Ue2 "<<Ue2<<endl;
    TComplex Ue3 = cpPhaseNeg*s13; // s13 * e^{-i*delta}  
    
    //muon 
    cout<<"Ue3 Re "<<Ue3.Re() <<" Im "<<Ue3.Im()<<endl; //take Re & Im of Ue3
    TComplex Umu1 = -s12 * c23 - c12 * s23 * s13 * cpPhasePos;   
    cout<<"Umu1 Re "<<Umu1.Re() <<" Im "<<Umu1.Im()<<endl;
    TComplex Umu2 = c12*c23-s12*s23*s13*cpPhasePos;
    cout<<"Umu2 Re "<<Umu2.Re()<<" Im "<<Umu2.Im()<<endl;   
    double Umu3 = s23*c13;
    cout<<"Umu3  "<<Umu3<<endl;
    
    //tau 
    TComplex Utau1=s12*s23-c12*c23*s13*cpPhasePos;
    cout<<"Utau1 Re "<<Utau1.Re()<<" Im "<<Utau1.Im()<<endl;
    TComplex Utau2=-c12*s23-s12*c23*s13*cpPhasePos;
    cout<<"Utau2 Re "<<Utau2.Re()<<" Im "<<Utau2.Im()<<endl;
    double Utau3=c23*c13;
    cout<<"Utau3  "<<Utau3<<endl;
    
    //test unitary condition (=1)
    double norm_e = Ue1**2 + Ue2**2 + Ue3.Rho2();
    double norm_mu = Umu1.Rho2() + Umu2.Rho2() + Umu3**2;
    double norm_tau = Utau1.Rho2() + Utau2.Rho2() + Utau3**2;  //Rho2 x=a+bi ->|x|^2 = a^2+b^2
    cout<<"norm e "<<norm_e<<endl;
    cout<<"norm mu "<<norm_mu<<endl;
    cout<<"norm tau "<<norm_tau<<endl;
    
    //calculate the triangle
    TComplex Ue1Umu1=Ue1*Umu1;//Ue1 is real 
    cout<<"Ue1Umu1 Re "<<Ue1Umu1.Re()<<" Im "<<Ue1Umu1.Im()<<" Amp "<<Ue1Umu1.Rho()<<endl; //Rho x=a+bi -> |x| = sqrt(a^2+b^2)
    TComplex Ue2Umu2=Ue2*Umu2;//Ue2 is real
    cout<<"Ue2Umu2 Re "<<Ue2Umu2.Re()<<" Im "<<Ue2Umu2.Im()<<" Amp "<<Ue2Umu2.Rho()<<endl;
    TComplex Ue3Umu3 = (TComplex::Conjugate(Ue3)) *Umu3;  //Umu3 is real (Conjugate ??)
    cout<<"Ue3Umu3 Re "<<Ue3Umu3.Re()<<" Im "<<Ue3Umu3.Im()<<" Amp "<<Ue3Umu3.Rho()<<endl;
    TComplex sumCheck = Ue1Umu1 + Ue2Umu2 + Ue3Umu3; //equal 0
    cout<<"sumCheck Re "<<sumCheck.Re()<<" Im "<<sumCheck.Im()<<endl;
    
    // sides
    TComplex sideA=Ue1Umu1/Ue2Umu2;
    TComplex sideB=Ue3Umu3/Ue2Umu2;
    cout<<"sideA Re "<<sideA.Re()<<" Im "<<sideA.Im()<<" Amp "<<sideA.Rho()<<" angle "<<sideA.Theta()<<endl; //theta??
    cout<<"sideB Re "<<sideB.Re()<<" Im "<<sideB.Im()<<" Amp "<<sideB.Rho()<<" angle "<<sideB.Theta()<<endl;
    
    
	TH2D *h2 = new TH2D("h2","",200,-0.5,1.2,200, -0.1,0.1); 
    h2->GetXaxis()->SetTitle("#rho (Real)"); 
    h2->GetYaxis()->SetTitle("#eta (Imagine)");
    h2->GetYaxis()->SetTitleOffset(1.2); //
    h2->GetXaxis()->CenterTitle(); //dat tieu de o vi tri o giua
    h2->GetYaxis()->CenterTitle();
    gPad->SetGrid();  //ve luoi cho truc toa do trong Pad 
   
	//draw lines
    TLine *linec= new TLine(0,0,1,0);
    TLine *linea= new TLine(1,0,1+sideA.Re(),sideA.Im()); //
    TLine *lineb= new TLine(0,0,-sideB.Re(),-sideB.Im()); //
    linec->SetLineWidth(1);
    linea->SetLineWidth(1);
    lineb->SetLineWidth(1);
    h2->Draw("AXIS");
    linec->Draw("same");
    linea->Draw("same");
    lineb->Draw("same"); 
    
    //ky hieu latex , fig2
    TLatex vertexA(-0.05,-0.015,"A(0,0)");  //dinh
    TLatex vertexB(0.9,-0.015,"B(0,1)");
    TLatex texsideA(0.3,0.2,"#left| #frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}} #right|");
    TLatex texsideB(0.8,0.035,"#left| #frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}} #right|");
    vertexA.SetTextColor(kSpring);
    vertexB.SetTextColor(kSpring);
     texsideA.SetTextSize(0.035);
    texsideB.SetTextSize(0.035);
    vertexA.Draw("same");
    vertexB.Draw("same");
    texsideA.Draw("same");
    texsideB.Draw("same"); 
    
    vertexA.DrawLatex(-0.05,-0.015,"A(0,0)");  //dinh  fig1
    vertexB.DrawLatex(0.9,-0.015,"B(0,0)");
    texsideA.DrawLatex(-0.35,0.03,"#left| #frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}} #right|");
    texsideB.DrawLatex(0.8,0.035,"#left| #frac{U_{e3}^{*}U_{#mu3}}{U_{e2}^{*}U_{#mu2}} #right|");
    vertexA.SetTextColor(kMagenta);
    vertexB.SetTextColor(kMagenta); 
    
    
    if(i==1)
    {
		TLatex vertexC(0.22,0.07,"C(0.241,0.057)");
		TLatex angleA(0.1,0.005,"#gamma =13^{o}");
		TLatex angleB(0.6,0.005,"#beta =4^{o}");
		vertexC.Draw("same");
		angleA.Draw("same");
		angleB.Draw("same");
    }
    
    if(i==2)
    {    
		vertexC.DrawLatex(0.22,0.07,"C(0.243,0.065)");
		angleA.DrawLatex(0.1,0.005,"#gamma =15^{o}");
		angleB.DrawLatex(0.6,0.005,"#beta =5^{o}");
	}
	
	angleA.SetTextColor(kAzure);
    angleB.SetTextColor(kAzure);
	vertexC.SetTextColor(kRed);
    c1->Print("triangle_normal_and_invert.eps");
}

}
