#include "particle.h"
#include "constants.h"

//Unità Naturali (c=1)


Particle::Particle() {
	X=0.0; // cm
	Y=0.0; // cm
	Z=0.0; // cm
	ThetaX=0.0; // rad
	ThetaY=0.0; // rad
	Momentum=0.0; // GeV/c
	Charge=0;
}

Particle::Particle(double x,double y, double z, double tx, double ty, double m, int c){
	X=x; // cm
	Y=y; // cm
	Z=z; // cm
	ThetaX=tx; // rad
	ThetaY=ty; // rad
	Momentum=m; // GeV/c
	Charge=c;
}


//Getter

double Particle::getX(){
    return X;
}

double Particle::getY(){
    return Y;
}

double Particle::getZ(){
    return Z;
}

double Particle::getTX(){
    return ThetaX;
}

double Particle::getTY(){
    return ThetaY;
}

double Particle::getMomentum(){
    return Momentum;
}

int Particle::getCharge(){
	return Charge;
}


//Metodi di movimento della particella

void Particle::move(double r) {
	X+=(r/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Cos(ThetaY);
	Y+=(r/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Sin(ThetaY);
	Z+=r;
}

void Particle::move_BGO(double r, double &energy_inside_calorimeter, TH1D *histogram_energy) {
	double step=length_of_calorimeter/13.0;
	double save=r;
	//In questa parte del movimento viene valutato se la particella deve compiere uno spostamento più lungo o meno lungo
	//di una lunghezza di radiazione
	if ((Momentum-r*minimum_ionization) > critical_energy) {
		while (true)  {
			if(save > step) {
				X+=(step/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Cos(ThetaY);
				Y+=(step/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Sin(ThetaY);
				Z+=step;
				energy_inside_calorimeter+=minimum_ionization*step;
				Momentum-=minimum_ionization*step;
				histogram_energy->Fill((Z-start_of_calorimeter-step)/radiation_length_BGO_cm,minimum_ionization*step);
				save-=step;
			} else if (save <= step) {
				X+=(save/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Cos(ThetaY);	
				Y+=(save/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Sin(ThetaY);
				Z+=save;
				energy_inside_calorimeter+=minimum_ionization*save;
				Momentum-=minimum_ionization*save;
				histogram_energy->Fill((Z-start_of_calorimeter)/radiation_length_BGO_cm,minimum_ionization*save);
				break;				
			}
		}
	//Quando la particella non ha abbastanza momento dopo lo spostamento viene suddiviso lo spostamento in dieci tratti più piccoli
	//in modo da non far arrivare la particella ad energie negative
	} else if ((Momentum-r*minimum_ionization) <= critical_energy && (Momentum-r*minimum_ionization) > 0.0) {
		while (true) {
			if (Momentum > critical_energy) {
				X+=(0.1*r/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Cos(ThetaY);
				Y+=(0.1*r/TMath::Cos(ThetaX))*TMath::Sin(ThetaX)*TMath::Sin(ThetaY);
				Z+=0.1*r;
				energy_inside_calorimeter+=minimum_ionization*0.1*r;
				Momentum-=minimum_ionization*0.1*r;
				histogram_energy->Fill((Z-start_of_calorimeter)/radiation_length_BGO_cm,minimum_ionization*0.1*r);			
			} else if (Momentum < critical_energy) {
				energy_inside_calorimeter+=Momentum;
				histogram_energy->Fill((Z-start_of_calorimeter)/radiation_length_BGO_cm,Momentum);
				Momentum=0.0;
				break;
			}
		}

	} else if (Momentum < 0.0) {
		Momentum=0.0;
	}
}


//Metodi di interazione

void Particle::tracker(TF1 *tf, TRandom3 *rand, std::vector<double> &values_from_landau, TH1D *histogram_landau) {
	double land=0.0;
	double Theta_0=0.0;
	land=tf->GetRandom();
	Theta_0=(13.6*TMath::Power(10,-3)/Momentum) * TMath::Abs(Charge) * TMath::Sqrt(layer_thickness/radiation_length_SI) * (1+0.038 * TMath::Log(layer_thickness/radiation_length_SI)); 
	ThetaX=rand->Gaus(0,Theta_0);
	ThetaY=rand->Gaus(0,Theta_0);
	Momentum-=land;
	//Inefficienza
	if (rand->Binomial(1,0.95) == 1) {
		values_from_landau.push_back(land);
	}
}

void Particle::pair_production(TF1 *pdf_pp, TRandom3 *rand, std::vector<Particle *> &par, double &energy_inside_calorimeter) {
	if (Momentum > 0.0) {
		double x=0.0;
		double energy_pp=0.0;
		double sigma=0.0;
		x=pdf_pp->GetRandom();
		energy_pp=x*Momentum;
		sigma=electron_mass/Momentum;
    	Momentum=Momentum-energy_pp;
    	ThetaX=(rand->Gaus(0,sigma))/TMath::Sqrt(2);
    	ThetaY=(rand->Gaus(0,sigma))/TMath::Sqrt(2);
    	Charge=-1;
		par.push_back(new Particle(X,Y,Z,rand->Gaus(0,sigma)/TMath::Sqrt(2),rand->Gaus(0,sigma)/TMath::Sqrt(2),energy_pp, 1));
	} else {
		Momentum=0.0;
	}
}

void Particle::bremsstrahlung (TF1 *pdf_bremm, TRandom3 *rand, std::vector<Particle *> &par, double &energy_inside_calorimeter, TH1D *histogram_energy) {
	if (Momentum > critical_energy) {
		double y=0.0;
		while(y<1.0/3.0 && Momentum>critical_energy) {
			double sigma=0.0;
			double energy_bremm=0.0;
			y=pdf_bremm->GetRandom();
			energy_bremm=y*Momentum;
			sigma=electron_mass/Momentum;
			if (y<1.0/3.0) {
				if (energy_bremm>critical_energy) {
					Momentum-=energy_bremm;
					par.push_back(new Particle(X,Y,Z,rand->Gaus(0,sigma)/TMath::Sqrt(2),rand->Gaus(0,sigma)/TMath::Sqrt(2),energy_bremm, 0));
					
				} else {
					Momentum-=energy_bremm;
					energy_inside_calorimeter+=energy_bremm;
					histogram_energy->Fill((Z-start_of_calorimeter)/radiation_length_BGO_cm,energy_bremm);						
				}
			} else {
				ThetaX=(rand->Gaus(0,sigma))/TMath::Sqrt(2);
    			ThetaY=(rand->Gaus(0,sigma))/TMath::Sqrt(2);
				Momentum-=energy_bremm;
				par.push_back(new Particle(X,Y,Z,rand->Gaus(0,sigma)/TMath::Sqrt(2),rand->Gaus(0,sigma)/TMath::Sqrt(2),energy_bremm, 0));
				break;
			}
		}
	} else if (Momentum < critical_energy && Momentum>0.0) {
		energy_inside_calorimeter+=Momentum;
		histogram_energy->Fill((Z-start_of_calorimeter)/radiation_length_BGO_cm,Momentum);
		Momentum=0.0;
	} else if (Momentum < 0.0) {
		Momentum=0.0;
	}
}