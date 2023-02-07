//Costanti iniziali
const int initial_electrons=1*TMath::Power(10,3); //number
const double initial_energy=180.0; //GeV/c


//Costanti Landau e tracker
const int number_of_layers=10;
const double distance_of_layers=3.0; //cm
const double layer_thickness=300.0*TMath::Power(10,-4); //cm
const double couple_energy=3.6*TMath::Power(10,-9); //GeV
const double couple_generation=80.0*TMath::Power(10,4); //number/cm
const double most_probable_value_landau=(couple_generation*couple_energy*layer_thickness); //GeV
const double K=0.307075*TMath::Power(10,-3); //GeV mol^-1 cm^2
const double atomic_number_SI=14.0;
const double mass_number_SI=28.0855; //g/mol
const double SI_density=2.329; //g/cm^3
const double mass_thickness=layer_thickness*SI_density; //g/cm^2
const double xi_landau=(K/2.0)*(atomic_number_SI/mass_number_SI)*mass_thickness; //GeV
const double radiation_length_SI=9.37; //cm

//Costanti Pair Production e Bremmstrahlung
const double BGO_density=7.130; //g/cm^3
const double radiation_length_BGO_mass=7.97; //g/cm^2
const double radiation_length_BGO_cm=1.118; //cm
const double mean_free_path_pp=(9.0/7.0)*(radiation_length_BGO_mass/BGO_density); //cm
const double mu_pp=TMath::Power(mean_free_path_pp,-1); //1/cm
const double mean_free_path_bremm=(radiation_length_BGO_mass/BGO_density); //cm
const double mu_bremm=TMath::Power(mean_free_path_bremm,-1); //1/cm
const int atomic_number_O=8; //number
const int atomic_number_Ge=32; //number
const int atomic_number_Bi=83; //number
const double fraction_atomic_number_mass_number_BGO=0.42065; //number
const int atomic_number_BGO=(atomic_number_Bi*2+atomic_number_O*3)*2+(atomic_number_Ge+atomic_number_O*2)*3; //number
const double mass_number_BGO=(double(atomic_number_BGO)/fraction_atomic_number_mass_number_BGO); //number
const double minimum_ionization=8.918*TMath::Power(10,-3); //GeV/cm
const double critical_energy=10.32*TMath::Power(10,-3); //GeV, average value for e- and e+
const double electron_mass=0.51099895000*TMath::Power(10,-3); //GeV/c^2
const double length_of_calorimeter=15.0; //cm
const double start_of_calorimeter=1037.3; //cm
const double end_of_calorimeter=start_of_calorimeter+length_of_calorimeter; //cm
