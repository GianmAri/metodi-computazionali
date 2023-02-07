#include <iostream>
#include <vector>
#include <ctime>
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1D.h"

class Particle {
private:
	double X, Y, Z, ThetaX, ThetaY, Momentum;
	int Charge;
public:
    Particle();
	Particle(double x,double y, double z, double tx, double ty, double m, int c);
   	double getX();
    double getY();
	double getZ();
	double getTX();
	double getTY();
	double getMomentum();
	int getCharge();
	void move(double r);
	void move_BGO(double r, double &energy_inside_calorimeter, TH1D *histogram_energy);
	void tracker(TF1 *tf, TRandom3 *rand, std::vector <double> &values_from_landau, TH1D *histogram_landau);
	void pair_production(TF1 *pdf_pp, TRandom3 *rand, std::vector<Particle *> &par, double &energy_inside_calorimeter);
	void bremsstrahlung(TF1 *pdf_bremm, TRandom3 *rand, std::vector<Particle *> &par, double &energy_inside_calorimeter, TH1D *histogram_energy);
};