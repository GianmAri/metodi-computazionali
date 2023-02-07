double mean(std::vector <double> &ver) {
    double somma=0.0;
    for (unsigned int i=0;i<ver.size();i++) {
        somma+=ver[i];
    }
    return somma/ver.size();
}

double trimmed_mean(std::vector<double> &ver) {
    double somma=0.0;
    double high=ver[0];
    for (unsigned int i=0;i<ver.size();i++) {
        if(ver[i]>high) {
            high=ver[i];
        }
    }
    double j=0.0;
    for (unsigned int i=0;i<ver.size();i++) {
        if(ver[i]!=high) {
            somma+=ver[i];
            j++;
        }
    }
    return somma/double(j);
}

double standard_deviation (std::vector <double> &ver) {
    double sum=0.0;
    double media=mean(ver);
    for (unsigned int i=0;i<ver.size();i++) {
        sum=sum+TMath::Power(ver[i]-media,2);
    }
    return TMath::Sqrt(sum/ver.size());
}


double simpson(double a, double b, double  E, TF1 *montecarlo_function) {
    double N=10000.0;
    double h=(b-a)/N;
    double value=0.0;
    double x=0.0;
    
    for (int i=1;i<N-1;i++) {
        x=(a+i*h)/radiation_length_BGO_cm;
        if (i%2 == 0) {
            value+=montecarlo_function->Eval(x)*2;
        }
        if (i%2 != 0) {
            value+=montecarlo_function->Eval(x)*4;
        }
    }
    value=(value+montecarlo_function->Eval(a)+montecarlo_function->Eval(b))*h/3.;
    return value;
}