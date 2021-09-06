using namespace std;
using namespace TMath;

Double_t Tb_func(Double_t theta, Double_t ma, Double_t mb, Double_t my, Double_t Ta, Double_t Q)
{
    return Power((Cos(theta)*Sqrt(ma*mb*Ta)+Sqrt(ma*mb*Ta*Power(Cos(theta),2)+(my+mb)*(my*(Q)+(my-ma)*Ta)))/(my+mb),2);;
}


void kinematics(){
    
    //Levels in 15N
//    Double_t Ex = 0.; // MeV
//    Double_t Ex = 5.269;//5/2+
    Double_t Ex = 5.298;//1/2+
//    Double_t Ex = (5.269+5.298)/2;
//    Double_t Ex = 6.32235;
//    Double_t Ex = 7.15505; //  for the 1.88477 line
//    Double_t Ex = 7.29892;
    
    //Levels in 16O
//    Double_t Ex = 0.; // MeV
//    Double_t Ex = 6.12863;//3-
//    Double_t Ex = 6.9155;//2+
    
    //Levels in 12C
//    Double_t Ex = 0.; // MeV
//    Double_t Ex = 4.43894;//2+
    
    //Levels in 19F
//    Double_t Ex = 0.; // MeV
//    Double_t Ex = 0.110;
//    Double_t Ex = 0.197;
//    Double_t Ex = 1.554;// 1.3 MeV line
//    Double_t Ex = 2.779849;// 2.5MeV line
    
    Double_t Siri_angle[16] = {140.0,138.0,136.0,134.0,132.0,130.0,128.0,126.0,54.0,52.0,50.0,48.0,46.0,44.0,42.0,40.0};
    Double_t Siri_angle_rad[16];
    Double_t Siri_angle_cos[16];
    
    for (Int_t i=0; i<16; i++) {
        Siri_angle_rad[i]=Pi()*Siri_angle[i]/180;
        Siri_angle_cos[i]=Cos(Siri_angle_rad[i]);
    }
    
    Double_t theta = Siri_angle_rad[0];

//    Double_t theta = Pi()*175.0/180;
//    Double_t theta = Pi()*10.0/180;
    
    
    Double_t Q;
    
    Double_t Ta = 26.;//Beam energy
    Double_t ma = 3728.400916;//alpha
    Double_t pa2 = Ta*Ta + 2*ma*Ta;
    
    Double_t mb = 938.7829706;//proton
//    Double_t mb = 3728.400916;//alpha
    
    Double_t mx = 11177.928;//12C
//    Double_t mx = 14899.167;//16O
    
    Double_t my = 13972.51144+Ex;//15N
//    Double_t my = 14899.167+Ex;//16O
//    Double_t my = 11177.928+Ex;//12C
//    Double_t my = 17696.89856+Ex;//19F
	
    Double_t Tb = 0., Ty = 0.;
    Double_t k0, k1, a, b, c;
    
    //---------------------------------------------------------------------
    
    Q=mx+ma-my-mb;
    cout << "Q0: " << setprecision(7) << mx+ma-my-mb+Ex  << " MeV" << endl;
    
    Tb = Tb_func(theta,ma,mb,my,Ta,Q);
    Ty = Q - Tb +Ta;
    
//    cout << "SiRi angle : " << theta*180/Pi() << endl;
//    cout << "outgoing particle energy: " << Tb  << " MeV" << endl;
//    cout << "recoil energy: " << Ty  << " MeV" << endl;
//    cout << "Excitation energy: " << Ex  << " MeV" << endl;
//    cout << "Q: " << Q  << " MeV" << endl;
//    cout << "recoil velocity: " << 100*Sqrt(2*Ty/my)  << " %" << endl;
//    cout << "Beam energy (test): " << Ty + Tb - Q  << endl;
//    cout << "Doppler shift for the 5.2 MeV 15N line: " << 5283.489*(1+Sqrt(2*Ty/my)*Cos(Pi()*37.0/180))-5283.489  << " keV" << endl;
//    cout << "Recoil angle : " << ASin(Sqrt(mb*Tb)*Sin(theta)/Sqrt(my*Ty))*180/Pi() << endl;
//    cout << "Doppler shifted energy for 6128.63 keV 16O line: " << 6128.63*(1+Sqrt(2*Ty/my)*Cos(Pi()*79.0/180))  << " MeV" << endl;
//    cout << "Doppler shifted energy for 2582.7 keV 19F line: " << 2582.7*(1+Sqrt(2*Ty/my)*Cos(Pi()*79.0/180))  << " MeV" << endl;
    
    
    
//    Double_t temp;
//    temp = Sqrt(2*Ty/my);
//    Ex = 7.29892;
//    my = 13972.51144+Ex;
//    Q=mx+ma-my-mb;
//    Tb = Power((Cos(theta)*Sqrt(ma*mb*Ta)+Sqrt(ma*mb*Ta*Power(Cos(theta),2)+(my+mb)*(my*(Q)+(my-ma)*Ta)))/(my+mb),2);
//    Ty = Q - Tb +Ta;
//    cout << "Excitation energy: " << Ex  << " MeV" << endl<< endl;
//    cout << "recoil velocity: " << 100*Sqrt(2*Ty/my)  << " %" << endl<< endl;
//    cout << setprecision(10) << "(1+b1cos)/(1+b2cos): " << (1+temp*Cos(Pi()*37.0/180))/(1+Sqrt(2*Ty/my)*Cos(Pi()*37.0/180)) << endl<< endl;
    
    Double_t factor[16], velocity[16];
    
    for (Int_t i=0; i<16; i++) {
        factor[i]=1/(1+Power(Sin(Siri_angle_rad[i])/(Sqrt(Ta/Tb_func(Siri_angle_rad[i],ma,mb,my,Ta,Q))-Cos(Siri_angle_rad[i])),2));
        velocity[i]=100*Sqrt(2*(Q-Tb_func(Siri_angle_rad[i],ma,mb,my,Ta,Q)+Ta)/my);
    }
    
    TCanvas *c22 = new TCanvas("c22", "factor vs angle", 1200, 600);
    c22->Divide(2,1);
    
    c22->cd(1);
    TGraph *gr_factor = new TGraphErrors(16,Siri_angle_cos,factor);
    gr_factor->SetMarkerStyle(20.);
    gr_factor->SetMarkerSize(0.7);
    gr_factor->Draw("ap");
    gr_factor->GetXaxis()->SetRangeUser(-1.,1.);
    gr_factor->GetYaxis()->SetRangeUser(0.,1.);
    gr_factor->SetTitle("Factor f");
    TF1 *fit_factor= new TF1("fit_factor","pol0", -1.,0);
    gr_factor->Fit("fit_factor","RQ");
    
    Double_t factor_error = (factor[0]-factor[7])/2;
    cout << "f: " << fit_factor->GetParameter(0) << " ± " << factor_error << " (" << 100*factor_error/fit_factor->GetParameter(0) << "%)" << endl;
    
    c22->cd(2);
    TGraph *gr_velocity = new TGraphErrors(16,Siri_angle_cos,velocity);
    gr_velocity->SetMarkerStyle(20.);
    gr_velocity->SetMarkerSize(0.7);
    gr_velocity->Draw("ap");
    gr_velocity->GetXaxis()->SetRangeUser(-1.,1.);
    gr_velocity->GetYaxis()->SetRangeUser(0.,4.);
    gr_velocity->SetTitle("Velocity #beta [%]");
    TF1 *fit_velocity= new TF1("fit_velocity","pol0", -1.,0);
    gr_velocity->Fit("fit_velocity","RQ");
    
    Double_t velocity_error = (velocity[0]-velocity[7])/2;
    cout << "beta: " << fit_velocity->GetParameter(0)/100 << " ± " << velocity_error/100 << " (" << 100*velocity_error/fit_velocity->GetParameter(0) << "%)" << endl;
    
}
