//
//  LinearInterpolator.cpp
//  LinearInterpolator
//
//  Created by Hugh Dickinson on 9/18/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#include "LinearInterpolator.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>

#include <TTextToData.h>
#include <TF2.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>

LinearInterpolator<3> * gli = NULL;

double processEnergy(double energy){
    return energy > 0.0 ? std::log10(energy) : 0.0;
}

double processFlux(double flux, double energy){
    return flux > 0.0 ? std::log10(flux/(energy*energy)) : 0.0;
}

double calcBField(int step, double minBField = 1e-17, double maxBField = 1e-15, int nSteps = 6){
    double bFieldStep = (maxBField - minBField)/static_cast<double>(nSteps-1);
    return minBField + step*bFieldStep;
}

double calcRedshift(int step, double minRedShift = 0.3, double maxRedshift = 1, int nSteps = 6){
    double redshiftStep = (maxRedshift - minRedShift)/static_cast<double>(nSteps-1);
    return minRedShift + step*redshiftStep;
}

void LoadFile(LinearInterpolator<3> & li, std::string const & prefix){
	TFieldSelector cols("0, 1");
    TFieldSelector rows;
    
    for(int i = 1; i <= 6; i++){ // b-field from 1e-17 to 1e-15 Gauss
        for(int j = 1; j <= 6; j++){ // redshift from 0.3 to 1
            std::stringstream filename;
            filename << prefix << i << "_" << j;
            TTextToData<double> tttd(filename.str(), cols, rows);
            std::vector<double> *x, *y;
            x = tttd.GetCol(0);
            y = tttd.GetCol(1);
            std::transform(x->begin(), x->end(), x->begin(), processEnergy);
            std::transform(y->begin(), y->end(), x->begin(), y->begin(), processFlux);

			std::vector<double>::iterator energyIt, fluxIt;
			
			for(energyIt = x->begin(), fluxIt = y->begin(); energyIt != x->end(); ++energyIt, ++fluxIt){
				std::cout << "Inserting (" << calcRedshift(j-1) << ", " << calcBField(i-1) << ", " << *energyIt << ") = " << *fluxIt << std::endl;
				li[calcRedshift(j-1)][calcBField(i-1)][*energyIt] = *fluxIt;
			}
				
            delete x;
            delete y;
        }
    }

}

double InterpolatedWrapper(double * x, double * pars){
	double energy = x[0];
	double bField = pars[0];
	double redshift = x[1];
	
	return (*gli)(energy)(bField)(redshift);
}

int main(int argc, char * argv[]){
    
	TApplication app("app", &argc, argv);
	
	gStyle->SetCanvasPreferGL(true);
	
    LinearInterpolator<3> li;
	gli = &li;
	
	LoadFile(li, app.Argv(1));
	
	TF2 * func = new TF2("func", InterpolatedWrapper, 6.5, 14.0, 0.3, 1.0, 1);
	
	TCanvas * canv = new TCanvas("gl_canv", "canv", 1);
	gPad->SetPhi(330);

	for(int i = 0; i < 20; i++){
		std::stringstream title;
		double bVal = 1e-17+i*((1e-15 - 1e-17)/20.0);
		title << "B = " << bVal << ";log_{10}(E_{#gamma})[eV];Redshift;log_{10}(E_{#gamma}F_{E_{#gamma}})[eV s^{-1}]";
		
		func->SetParameter(0, bVal);
		func->SetTitle(title.str().c_str());
		func->Draw("surf1");
		
		func->GetXaxis()->SetTitleOffset(1.8);
		func->GetYaxis()->SetTitleOffset(1.8);
		func->GetZaxis()->SetTitleOffset(1.2);
		
		canv->Modified();
		canv->Update();
		
		canv->Print("distro.gif+50");
	}
	
	app.Run();
	
    return 0;
}
