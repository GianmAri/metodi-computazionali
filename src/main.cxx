#include "particle.h"
#include "constants.h"
#include "functions.cxx"
#include <chrono>


int main(int argc,char *argv[]) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        TApplication *myApp = new TApplication("myApp",&argc,argv);

        //Generatore di numeri casuali Mersenne-Twister
        int seed=time(0);

        TRandom3 *rand=new TRandom3(seed);  

        //Canvas per i grafici
        TCanvas *canvas_1 = new TCanvas("canvas_1", "canvas_1");
        TCanvas *canvas_2= new TCanvas("canvas_2", "canvas_2");

        //Istogramma carica ricostruita
        TH1D *histogram_landau= new TH1D("histogram_landau","Rebuilt charge",TMath::Sqrt(initial_electrons),(most_probable_value_landau-5*2*xi_landau),(most_probable_value_landau+15*4*xi_landau));
        histogram_landau->GetXaxis()->SetTitle("Energy (GeV)");
        histogram_landau->GetYaxis()->SetTitle("Entries");

        //Istogramma posizione particella in relazione all'energia
        TH1D *histogram_energy= new TH1D("histogram_energy","Stored energy",int(length_of_calorimeter/radiation_length_BGO_cm),0,length_of_calorimeter); 
        histogram_energy->GetXaxis()->SetTitle("Depth in radiation length");
        histogram_energy->GetYaxis()->SetTitle("Energy (GeV)");
        

        //Funzione Landau per tracker
	TF1 *landau = new TF1("tf","TMath::Landau(x, [0],[1])",(most_probable_value_landau-5*2*xi_landau),(most_probable_value_landau+15*2*xi_landau));
	landau->SetParameter(0,most_probable_value_landau);
	landau->SetParameter(1,2.0*xi_landau);

        //Funzione di estrazione di probabilità di interazione Pair Production
        TF1 *P_int_pp = new TF1("P_int_pp","TMath::Exp([0]*x)",0,2*length_of_calorimeter);
	P_int_pp->SetParameter(0,-mu_pp);

        //Funzione di estrazione di probabilità di interazione bremsstrahlung
        TF1 *P_int_bremm = new TF1("P_int_bremm","TMath::Exp([0]*x)",0,2*length_of_calorimeter);
	P_int_bremm->SetParameter(0,-mu_bremm);

        //Funzione di densità di probabilità Pair Production
        TF1 *pdf_pp=new TF1("pdf_pp", "(1-(4/3)*x+(4/3)*TMath::Power(x,2))",0,1);

        //Funzione di densità di probabilità bremsstrahlung
        TF1 *pdf_bremm=new TF1("pdf_pp", "((4/3)*-(4/3)*x+(4/3)*TMath::Power(x,2))",0,1); 

        //Energia al di fuori del calorimetro
        double energy_outside_calorimeter=0.0;

        //Energia depositata nel calorimetro
        double energy_inside_calorimeter=0.0;

        double x_move=0.0;

        //Vector contenente le particelle
        std::vector <Particle *> par;
        
        std::vector <double> values_from_landau;
        std::vector <double> trimmed_value;
        std::vector <double> energy_before_calorimeter;
        std::cout << std::endl;
        std::cout << "PARAMETRI SIMULAZIONE" << std::endl;
        std::cout << "Numero di elettroni iniziali N0: " << initial_electrons << std::endl;
        std::cout << "Valore seed generatore di numeri casuali: " << seed << std::endl;
        std::cout << std::endl;

        for (int i=0; i<initial_electrons; i++) {
                std::cout <<"Percentuale di completamento: " << (i+1)/(initial_electrons*TMath::Power(10,-2)) <<  "% \r";
                std::cout.flush();
                //CREAZIONE PARTICELLA
                //Particle P(rand->Gaus(0,1),rand->Gaus(0,1),0.0,rand->Gaus(0,10)*TMath::Power(10,-6),rand->Gaus(0,10)*TMath::Power(10,-6),initial_energy,-1); //cm, cm, cm, rad, rad, GeV/c, numero
                //par.push_back(P);
                par.push_back(new Particle(0,rand->Gaus(0,1),rand->Gaus(0,1),rand->Gaus(0,10)*TMath::Power(10,-6),rand->Gaus(0,10)*TMath::Power(10,-6),initial_energy,-1));
                par[0]->move(1000);

                //TRACKER
                for(int q=1;q<=number_of_layers;q++) {        
                        par[0]->tracker(landau, rand, values_from_landau, histogram_landau);
                        par[0]->move(layer_thickness);
                        if (q < number_of_layers) {
                                par[0]->move(distance_of_layers);
                        }
                }
                trimmed_value.push_back(trimmed_mean(values_from_landau));
                values_from_landau.clear();
                par[0]->move(10.0);
                energy_before_calorimeter.push_back(par[0]->getMomentum());

                //CALORIMETRO
                for (unsigned int j=0; j<par.size(); j++) {
                        while (par[j]->getZ() < end_of_calorimeter && par[j]->getMomentum() > 0) {
                                if (par[j]->getCharge() == 0) {
                                        x_move=P_int_pp->GetRandom();
                                        if ((par[j]->getZ()+x_move) < end_of_calorimeter) {
                                                par[j]->move(x_move); //cm
                                                par[j]->pair_production(pdf_pp,rand,par,energy_inside_calorimeter);
                                        } else if ((par[j]->getZ()+x_move) > end_of_calorimeter) {
                                                par[j]->move(end_of_calorimeter-par[j]->getZ());
                                        }       
                                } else {
                                        x_move=P_int_bremm->GetRandom();
                                        if ((par[j]->getZ()+x_move) < end_of_calorimeter) {
                                                par[j]->move_BGO(x_move,energy_inside_calorimeter, histogram_energy);
                                                par[j]->bremsstrahlung(pdf_bremm,rand,par,energy_inside_calorimeter, histogram_energy);
                                        } else if ((par[j]->getZ()+x_move) > end_of_calorimeter) {
                                                par[j]->move_BGO((end_of_calorimeter-par[j]->getZ()), energy_inside_calorimeter, histogram_energy);
                                        }
                                }
                        }
                }
                //Somma momenti esterni al calorimetro
                for (unsigned int j=0; j<par.size(); j++) {
                        if (par[j]->getZ() >= end_of_calorimeter) {
                                energy_outside_calorimeter+=par[j]->getMomentum();
                        }
                }  
                for (auto p : par) {
                        delete p;
                }
                par.clear();    
        }
        double mean_energy_before_calorimeter=mean(energy_before_calorimeter);
        energy_before_calorimeter.clear();
        double b=0.5;
        double a=b*(TMath::Log(mean_energy_before_calorimeter/critical_energy)-0.5)+1.0;
        TF1 *montecarlo_function=new TF1("montecarlo_function","([0]*[1]*((TMath::Power([1]*x,[2])*TMath::Exp(-[1]*x))/TMath::Gamma([3])))", 0, 15);
        montecarlo_function->SetParameter(0, mean_energy_before_calorimeter);
        montecarlo_function->SetParameter(1,b);
        montecarlo_function->SetParameter(2,a-1);
        montecarlo_function->SetParameter(3,a);
        double mean_landau=mean(trimmed_value);
        double standard_deviation_landau=standard_deviation(trimmed_value);
        for (unsigned int i=0;i<trimmed_value.size();i++) {
                histogram_landau->Fill(trimmed_value[i]);
        }
        double value_simpson=simpson(0,length_of_calorimeter/radiation_length_BGO_cm,mean_energy_before_calorimeter,montecarlo_function);
        trimmed_value.clear();
        std::cout << std::endl;
        std::cout << std::endl;
        
        std::cout << "MISURA CARICA" << std::endl;
        std::cout << "Media: " << mean_landau << " GeV " << std::endl;
        std::cout << "Deviazione standard: " << standard_deviation_landau << " GeV" << std::endl;
        std::cout << "Stimatore di carica: " << standard_deviation_landau/mean_landau << std::endl;
        std::cout << std::endl;
        
        std::cout << "RISOLUZIONE EQUAZIONE DIFFERENZIALE CON METODO DI SIMPSON" << std::endl;
        std::cout << "Calcolo energia contenuta nel calorimetro: " << value_simpson << " GeV  Calcolo energia contenuta moltiplicato per N0: " << value_simpson*initial_electrons << " GeV" << std::endl;
        std::cout << "Calcolo energia fuoriuscita dal calorimetro: " << initial_energy - value_simpson << " GeV  Calcolo energia fuoriuscita dal calorimetro moltiplicato per N0: " <<(initial_energy - value_simpson)*initial_electrons << " GeV" << std::endl;
        std::cout << "Percentuale energia contenuta nel calorimetro secondo il valore ottenuto dalla risoluzione dell'equazione differenziale: " << value_simpson/initial_energy *100 << "%" <<std::endl;
        std::cout << std::endl;

        std::cout << "CALORIMETRO" << std::endl;
        std::cout << "Energia media rispetto a N0 contenuta nel calorimetro: " << energy_inside_calorimeter/initial_electrons << " GeV Energia totale contenuta nel calorimetro: " << energy_inside_calorimeter << " GeV" << std::endl;
        std::cout << "Energia media rispetto a N0 fuoriuscita dal calorimetro: " << energy_outside_calorimeter/initial_electrons << " GeV Energia totale fuoriuscita dal calorimetro: " <<energy_outside_calorimeter << " GeV" << std::endl;
        std::cout << "Percentuale energia contenuta nel calorimetro: " << energy_inside_calorimeter/(initial_energy*initial_electrons) *100 << "%" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        canvas_1->cd(); histogram_landau->Draw();
        double normalize=1./double(initial_electrons);
        histogram_energy->Scale(normalize);
        canvas_2->cd(); montecarlo_function->Draw();  histogram_energy->Draw("HISTO,SAME");
        canvas_1->SaveAs("landau_histogram.ps");
        canvas_2->SaveAs("calorimeter_histogram.ps");
        end = std::chrono::system_clock::now(); 
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Tempo impiegato dalla simulazione: " << elapsed_seconds.count() << " s" << std::endl;
        myApp->Run();
        delete rand;
        delete landau;
        delete P_int_pp; delete P_int_bremm;
        delete pdf_bremm; delete pdf_pp;
        delete histogram_landau; delete histogram_energy;
        delete canvas_1; delete canvas_2;
        delete myApp;
        return 0;
}