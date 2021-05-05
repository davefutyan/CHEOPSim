/*
 * Star.cxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

#include <stdexcept>
#include <iostream>
#include <fstream>

#include "Star.hxx"
#include "source/include/TransitFluxModulator.hxx"

const string Star::spectralTypeString[64] = {"Undefined","O7","O8","O9",
														 "B0","B1","B2","B3","B4","B5","B6","B7","B8","B9",
														 "A0","A1","A2","A3","A4","A5","A6","A7","A8","A9",
														 "F0","F1","F2","F3","F4","F5","F6","F7","F8","F9",
														 "G0","G1","G2","G3","G4","G5","G6","G7","G8","G9",
														 "K0","K1","K2","K3","K4","K5","K6","K7","K8","K9",
														 "M0","M1","M2","M3","M4","M5","M6","M7","M8","M9"};

double Star::getEffectiveTemperature(Star::spectralType specType) {

	switch (specType) {
	case Star::O7:
		return 36500.; break;
	case Star::O8:
		return 34500.; break;
	case Star::O9:
		return 32500.; break;
	case Star::B0:
		return 31500.; break;
	case Star::B1:
		return 26000.; break;
	case Star::B2:
		return 20600.; break;
	case Star::B3:
		return 17000.; break;
	case Star::B4:
		return 16700.; break;
	case Star::B5:
		return 15700.; break;
	case Star::B6:
		return 14500.; break;
	case Star::B7:
		return 14000.; break;
	case Star::B8:
		return 12500.; break;
	case Star::B9:
		return 10700.; break;
	case Star::A0:
		return 9700.; break;
	case Star::A1:
		return 9200.; break;
	case Star::A2:
		return 8840.; break;
	case Star::A3:
		return 8550.; break;
	case Star::A4:
		return 8270.; break;
	case Star::A5:
		return 8080.; break;
	case Star::A6:
		return 8000.; break;
	case Star::A7:
		return 7800.; break;
	case Star::A8:
		return 7500.; break;
	case Star::A9:
		return 7440.; break;
	case Star::F0:
		return 7220.; break;
	case Star::F1:
		return 7030.; break;
	case Star::F2:
		return 6810.; break;
	case Star::F3:
		return 6720.; break;
	case Star::F4:
		return 6640.; break;
	case Star::F5:
		return 6510.; break;
	case Star::F6:
		return 6340.; break;
	case Star::F7:
		return 6240.; break;
	case Star::F8:
		return 6170.; break;
	case Star::F9:
		return 6040.; break;
	case Star::G0:
		return 5920.; break;
	case Star::G1:
		return 5880.; break;
	case Star::G2:
		return 5770.; break;
	case Star::G3:
		return 5720.; break;
	case Star::G4:
		return 5680.; break;
	case Star::G5:
		return 5660.; break;
	case Star::G6:
		return 5590.; break;
	case Star::G7:
		return 5530.; break;
	case Star::G8:
		return 5490.; break;
	case Star::G9:
		return 5340.; break;
	case Star::K0:
		return 5280.; break;
	case Star::K1:
		return 5170.; break;
	case Star::K2:
		return 5040.; break;
	case Star::K3:
		return 4840.; break;
	case Star::K4:
		return 4620.; break;
	case Star::K5:
		return 4450.; break;
	case Star::K6:
		return 4200.; break;
	case Star::K7:
		return 4050.; break;
	case Star::K8:
		return 3970.; break;
	case Star::K9:
		return 3880.; break;
	case Star::M0:
		return 3850.; break;
	case Star::M1:
		return 3700.; break;
	case Star::M2:
		return 3550.; break;
	case Star::M3:
		return 3400.; break;
	case Star::M4:
		return 3200.; break;
	case Star::M5:
		return 3050.; break;
	case Star::M6:
		return 2800.; break;
	case Star::M7:
		return 2650.; break;
	case Star::M8:
		return 2500.; break;
	case Star::M9:
		return 2450.; break;
	default:
		throw runtime_error("ERROR in Star constructor: Invalid spectral type");
	}

}

Star::Star(SkyPosition pos, double VbandMag, double cheopsMag, double gaiaMag, spectralType specType, double effectiveTemperature, string id) :
	m_pos(pos),m_Vmag(VbandMag),m_cheopsMag(cheopsMag),m_gaiaMag(gaiaMag),m_VmagErr(0.),m_cheopsMagErr(0.),m_gaiaMagErr(0.),m_spectralType(specType),m_id(id),m_radius(0.),m_mass(0.),
	m_effectiveTemperature(effectiveTemperature),m_effectiveTemperatureErr(0.),m_transitTime(BJD()),m_transitPeriod(0.),m_haloFlux(0.) {

	if (m_spectralType==Undefined && m_effectiveTemperature==0.){
		throw runtime_error("ERROR in Star constructor: Either spectral type or effective temperature must be defined");
	}

	if (m_spectralType==Undefined && m_effectiveTemperature>0.) {
		//Define the spectral type based on the effective temperature

		if (m_effectiveTemperature > (36500.+34500.)/2.) m_spectralType = O7;
		else if (m_effectiveTemperature > (34500.+32500.)/2.) m_spectralType = O8;
		else if (m_effectiveTemperature > (32500.+31500.)/2.) m_spectralType = O9;
		else if (m_effectiveTemperature > (31500.+26000.)/2.) m_spectralType = B0;
		else if (m_effectiveTemperature > (26000.+20600.)/2.) m_spectralType = B1;
		else if (m_effectiveTemperature > (20600.+17000.)/2.) m_spectralType = B2;
		else if (m_effectiveTemperature > (17000.+16700.)/2.) m_spectralType = B3;
		else if (m_effectiveTemperature > (16700.+15700.)/2.) m_spectralType = B4;
		else if (m_effectiveTemperature > (15700.+14500.)/2.) m_spectralType = B5;
		else if (m_effectiveTemperature > (14500.+14000.)/2.) m_spectralType = B6;
		else if (m_effectiveTemperature > (14000.+12500.)/2.) m_spectralType = B7;
		else if (m_effectiveTemperature > (12500.+10700.)/2.) m_spectralType = B8;
		else if (m_effectiveTemperature > (10700.+9700.)/2.) m_spectralType = B9;
		else if (m_effectiveTemperature > (9700.+9200.)/2.) m_spectralType = A0;
		else if (m_effectiveTemperature > (9200.+8840.)/2.) m_spectralType = A1;
		else if (m_effectiveTemperature > (8840.+8550.)/2.) m_spectralType = A2;
		else if (m_effectiveTemperature > (8550.+8270.)/2.) m_spectralType = A3;
		else if (m_effectiveTemperature > (8270.+8080.)/2.) m_spectralType = A4;
		else if (m_effectiveTemperature > (8080.+8000.)/2.) m_spectralType = A5;
		else if (m_effectiveTemperature > (8000.+7800.)/2.) m_spectralType = A6;
		else if (m_effectiveTemperature > (7800.+7500.)/2.) m_spectralType = A7;
		else if (m_effectiveTemperature > (7500.+7440.)/2.) m_spectralType = A8;
		else if (m_effectiveTemperature > (7440.+7220.)/2.) m_spectralType = A9;
		else if (m_effectiveTemperature > (7220.+7030.)/2.) m_spectralType = F0;
		else if (m_effectiveTemperature > (7030.+6810.)/2.) m_spectralType = F1;
		else if (m_effectiveTemperature > (6810.+6720.)/2.) m_spectralType = F2;
		else if (m_effectiveTemperature > (6720.+6640.)/2.) m_spectralType = F3;
		else if (m_effectiveTemperature > (6640.+6510.)/2.) m_spectralType = F4;
		else if (m_effectiveTemperature > (6510.+6340.)/2.) m_spectralType = F5;
		else if (m_effectiveTemperature > (6340.+6240.)/2.) m_spectralType = F6;
		else if (m_effectiveTemperature > (6240.+6170.)/2.) m_spectralType = F7;
		else if (m_effectiveTemperature > (6170.+6040.)/2.) m_spectralType = F8;
		else if (m_effectiveTemperature > (6040.+5920.)/2.) m_spectralType = F9;
		else if (m_effectiveTemperature > (5920.+5880.)/2.) m_spectralType = G0;
		else if (m_effectiveTemperature > (5880.+5770.)/2.) m_spectralType = G1;
		else if (m_effectiveTemperature > (5770.+5720.)/2.) m_spectralType = G2;
		else if (m_effectiveTemperature > (5720.+5680.)/2.) m_spectralType = G3;
		else if (m_effectiveTemperature > (5680.+5660.)/2.) m_spectralType = G4;
		else if (m_effectiveTemperature > (5660.+5590.)/2.) m_spectralType = G5;
		else if (m_effectiveTemperature > (5590.+5530.)/2.) m_spectralType = G6;
		else if (m_effectiveTemperature > (5530.+5490.)/2.) m_spectralType = G7;
		else if (m_effectiveTemperature > (5490.+5340.)/2.) m_spectralType = G8;
		else if (m_effectiveTemperature > (5340.+5280.)/2.) m_spectralType = G9;
		else if (m_effectiveTemperature > (5280.+5170.)/2.) m_spectralType = K0;
		else if (m_effectiveTemperature > (5170.+5040.)/2.) m_spectralType = K1;
		else if (m_effectiveTemperature > (5040.+4840.)/2.) m_spectralType = K2;
		else if (m_effectiveTemperature > (4840.+4620.)/2.) m_spectralType = K3;
		else if (m_effectiveTemperature > (4620.+4450.)/2.) m_spectralType = K4;
		else if (m_effectiveTemperature > (4450.+4200.)/2.) m_spectralType = K5;
		else if (m_effectiveTemperature > (4200.+4050.)/2.) m_spectralType = K6;
		else if (m_effectiveTemperature > (4050.+3970.)/2.) m_spectralType = K7;
		else if (m_effectiveTemperature > (3970.+3880.)/2.) m_spectralType = K8;
		else if (m_effectiveTemperature > (3880.+3850.)/2.) m_spectralType = K9;
		else if (m_effectiveTemperature > (3850.+3700.)/2.) m_spectralType = M0;
		else if (m_effectiveTemperature > (3700.+3550.)/2.) m_spectralType = M1;
		else if (m_effectiveTemperature > (3550.+3400.)/2.) m_spectralType = M2;
		else if (m_effectiveTemperature > (3400.+3200.)/2.) m_spectralType = M3;
		else if (m_effectiveTemperature > (3200.+3050.)/2.) m_spectralType = M4;
		else if (m_effectiveTemperature > (3050.+2800.)/2.) m_spectralType = M5;
		else if (m_effectiveTemperature > (2800.+2650.)/2.) m_spectralType = M6;
		else if (m_effectiveTemperature > (2650.+2500.)/2.) m_spectralType = M7;
		else if (m_effectiveTemperature > (2500.+2450.)/2.) m_spectralType = M8;
		else m_spectralType = M9;

	} else if (m_spectralType!=Undefined && m_effectiveTemperature==0.) {
		//Define the effective temperature based on the spectral type
		m_effectiveTemperature = getEffectiveTemperature(m_spectralType);
	}

	//Define mass and radius based on the spectral type
	double logL = 0.;
	switch (m_spectralType) {
	case Star::O7:
		m_mass = 28.0;
		logL = 5.17;
		break;
	case Star::O8:
		m_mass = 22.9;
		logL = 4.99;
		break;
	case Star::O9:
		m_mass = 19.7;
		logL = 4.82;
		break;
	case Star::B0:
		m_mass = 17.5;
		logL = 4.70;
		break;
	case Star::B1:
		m_mass = 11;
		logL = 4.13;
		break;
	case Star::B2:
		m_mass = 7.3;
		logL = 3.38;
		break;
	case Star::B3:
		m_mass = 5.4;
		logL = 2.96;
		break;
	case Star::B4:
		m_mass = 5.0;
		logL = 2.91;
		break;
	case Star::B5:
		m_mass = 4.6;
		logL = 2.80;
		break;
	case Star::B6:
		m_mass = 4.0;
		logL = 2.54;
		break;
	case Star::B7:
		m_mass = 3.9;
		logL = 2.49;
		break;
	case Star::B8:
		m_mass = 3.4;
		logL = 2.27;
		break;
	case Star::B9:
		m_mass = 2.8;
		logL = 1.79;
		break;
	case Star::A0:
		m_mass = 2.3;
		logL = 1.54;
		break;
	case Star::A1:
		m_mass = 2.15;
		logL = 1.41;
		break;
	case Star::A2:
		m_mass = 2.05;
		logL = 1.33;
		break;
	case Star::A3:
		m_mass = 2.00;
		logL = 1.29;
		break;
	case Star::A4:
		m_mass = 1.90;
		logL = 1.20;
		break;
	case Star::A5:
		m_mass = 1.85;
		logL = 1.16;
		break;
	case Star::A6:
		m_mass = 1.83;
		logL = 1.14;
		break;
	case Star::A7:
		m_mass = 1.76;
		logL = 1.06;
		break;
	case Star::A8:
		m_mass = 1.67;
		logL = 0.97;
		break;
	case Star::A9:
		m_mass = 1.67;
		logL = 0.97;
		break;
	case Star::F0:
		m_mass = 1.59;
		logL = 0.89;
		break;
	case Star::F1:
		m_mass = 1.50;
		logL = 0.77;
		break;
	case Star::F2:
		m_mass = 1.44;
		logL = 0.70;
		break;
	case Star::F3:
		m_mass = 1.43;
		logL = 0.67;
		break;
	case Star::F4:
		m_mass = 1.39;
		logL = 0.61;
		break;
	case Star::F5:
		m_mass = 1.33;
		logL = 0.54;
		break;
	case Star::F6:
		m_mass = 1.25;
		logL = 0.43;
		break;
	case Star::F7:
		m_mass = 1.21;
		logL = 0.36;
		break;
	case Star::F8:
		m_mass = 1.18;
		logL = 0.31;
		break;
	case Star::F9:
		m_mass = 1.14;
		logL = 0.26;
		break;
	case Star::G0:
		m_mass = 1.08;
		logL = 0.14;
		break;
	case Star::G1:
		m_mass = 1.07;
		logL = 0.13;
		break;
	case Star::G2:
		m_mass = 1.02;
		logL = 0.01;
		break;
	case Star::G3:
		m_mass = 1.00;
		logL = -0.01;
		break;
	case Star::G4:
		m_mass = 0.99;
		logL = -0.04;
		break;
	case Star::G5:
		m_mass = 0.98;
		logL = -0.05;
		break;
	case Star::G6:
		m_mass = 0.97;
		logL = -0.11;
		break;
	case Star::G7:
		m_mass = 0.96;
		logL = -0.12;
		break;
	case Star::G8:
		m_mass = 0.94;
		logL = -0.17;
		break;
	case Star::G9:
		m_mass = 0.90;
		logL = -0.25;
		break;
	case Star::K0:
		m_mass = 0.87;
		logL = -0.33;
		break;
	case Star::K1:
		m_mass = 0.85;
		logL = -0.37;
		break;
	case Star::K2:
		m_mass = 0.82;
		logL = -0.47;
		break;
	case Star::K3:
		m_mass = 0.78;
		logL = -0.58;
		break;
	case Star::K4:
		m_mass = 0.73;
		logL = -0.71;
		break;
	case Star::K5:
		m_mass = 0.72;
		logL = -0.75;
		break;
	case Star::K6:
		m_mass = 0.70;
		logL = -0.92;
		break;
	case Star::K7:
		m_mass = 0.64;
		logL = -1.04;
		break;
	case Star::K8:
		m_mass = 0.61;
		logL = -1.08;
		break;
	case Star::K9:
		m_mass = 0.58;
		logL = -1.22;
		break;
	case Star::M0:
		m_mass = 0.58;
		logL = -1.26;
		break;
	case Star::M1:
		m_mass = 0.52;
		logL = -1.42;
		break;
	case Star::M2:
		m_mass = 0.48;
		logL = -1.57;
		break;
	case Star::M3:
		m_mass = 0.43;
		logL = -1.78;
		break;
	case Star::M4:
		m_mass = 0.24;
		logL = -2.20;
		break;
	case Star::M5:
		m_mass = 0.15;
		logL = -2.52;
		break;
	case Star::M6:
		m_mass = 0.10;
		logL = -3.02;
		break;
	case Star::M7:
		m_mass = 0.098;
		logL = -3.21;
		break;
	case Star::M8:
		m_mass = 0.082;
		logL = -3.34;
		break;
	case Star::M9:
		m_mass = 0.065;
		logL = -3.53;
		break;
	default:
		throw runtime_error("ERROR in Star constructor: Invalid spectral type");
	}

	double luminosity = kLSun*pow(10,logL);
	m_radius = sqrt(luminosity/(4.*M_PI*kSigma*pow(m_effectiveTemperature,4)));
	m_radius /= 1000.*TransitFluxModulator::kSunRadius;

//	int prec = m_mass>=1. ? 2:1;
//	cout << "cp " << " M\\=" << setprecision(prec)
//		 << m_mass << "0_R\\=" << m_radius << "0_Teff\\=" << int(round(m_effectiveTemperature/500.)*500) << ".00_Tstep\\=15.00.dat "
//		 << " ../granulation_" << setprecision(prec) << spectralTypeString[m_spectralType]
//		 << "_M" << m_mass << "_R" << m_radius << "_Teff" << int(round(m_effectiveTemperature/500.)*500) << ".00_Tstep15.dat"<< endl;
//	ofstream stellar_params("stellar_params.txt", ios::app);
//	stellar_params << spectralTypeString[m_spectralType] << " \& " << setprecision(5) << m_effectiveTemperature << " \& " << m_mass << " \& " << setprecision(3) << m_radius << " \\\\" << endl;
//	stellar_params.close();

	//Define flux in ergs s^-1 cm^-2 nm^-1, using flux for 0 mag star in V band = 3.61E-11 W m^-2 nm^-1
	//Taken form http://www.vikdhillon.staff.shef.ac.uk/teaching/phy217/instruments/phy217_inst_photsys.html#table1
	double fluxPerNm = pow(10,-0.4*m_Vmag) * 3.61E-8;

	//Set flux integrated over wavelength for each wavelength bin (assuming uniform flux vs wavelength),
	//for bins of width 10nm from 330nm to 1100nm
	//NOTE: this is just to initialize to reasonable values.
	//Flux vs wavelength should be set by calling StarProducer::evaluateFlux(Star * star);
	for (int i=0; i<kNWavelength; i++) {
		m_wavelengthBinLowEdge[i] = kWavelengthLowBound + double(i)*10.;
		m_meanFlux[i] = fluxPerNm*10.;
	}

	//Read in limb darkening coefficients calculated using https://github.com/nespinoza/limb-darkening
	string limbDarkening_filename = "limb_darkening_apr2018.dat";
	ifstream limbDarkCoeffs(string(getenv("CHEOPS_SW"))+"/resources/"+limbDarkening_filename,ios::in);
    if (!limbDarkCoeffs.good()) throw runtime_error("Error in Star constructor: Error opening file "+limbDarkening_filename);
	string specType_LD,dummy;
	double coeff1,coeff2;
	bool found = false;
	while(limbDarkCoeffs.good())
	{
		limbDarkCoeffs >> specType_LD;
		if (specType_LD == getSpectralTypeString()) {
			limbDarkCoeffs >> coeff1 >> coeff2;
			found=true;
			break;
		} else {
			limbDarkCoeffs >> dummy >> dummy;
		}
	}
	limbDarkCoeffs.close();
	if (!found) {
		throw runtime_error("Error in Star constructor: No match found for spectral type "+getSpectralTypeString()+" in "+limbDarkening_filename);
	}
	m_limbDarkeningCoeffs = make_pair(coeff1,coeff2);

	//Write input file for calculation of limb darkening coefficients using https://github.com/nespinoza/limb-darkening
	//calculating log(g) = 4.43812 +log(M/R^2), using g_sun=4.43812, from https://sites.google.com/site/mamajeksstarnotes/basic-astronomical-data-for-the-sun
	//Apr 2018: this calculation is now replaced my calculation in spreadsheet /Users/futyand/CHEOPS/limb-darkening-master/mamajek.xlsx
//	ofstream limbDarkInputs("limb_darkening_inputs.dat",ios::app);
//	limbDarkInputs << getSpectralTypeString() << "\t" << m_effectiveTemperature << "\t" << (4.43812+log10(m_mass/(m_radius*m_radius))) << "\t0.0\t2.0\tthroughputQE_BoL_april2017.dat\tA100\t-1\t-1" << endl;
//	limbDarkInputs.close();
	//The file throughputQE_BoL_april2017.dat is generated using commented code in GlobalThroughputGenerator::doBegin()

}

Star::~Star() {
    for (vector<StarData*>::const_iterator it = m_timeData.begin(); it!=m_timeData.end(); ++it) delete (*it);
}

vector<StarData*> Star::getTimeSeries() const {
	if (m_timeData.size()>0) {
		return m_timeData;
	} else {
		throw runtime_error("Error in Star::getTimeSeries : time series is empty");
	}
}

double Star::getMeanFlux(double wavelength) const {

	int wavelengthBin=-1;
	for (int i=0; i<kNWavelength-1; i++) {
		if (wavelength>=m_wavelengthBinLowEdge[i] && wavelength<m_wavelengthBinLowEdge[i+1]) {
			wavelengthBin=i;
			break;
		}
	}
	if (wavelength>=m_wavelengthBinLowEdge[kNWavelength-1]) wavelengthBin=kNWavelength-1;
	if (wavelengthBin==-1) {
		return 0.;
	} else {
		return m_meanFlux[wavelengthBin];
	}

}

double Star::getMeanPhotonFlux(double wavelength) const {

	int wavelengthBin=-1;
	for (int i=0; i<kNWavelength-1; i++) {
		if (wavelength>=m_wavelengthBinLowEdge[i] && wavelength<m_wavelengthBinLowEdge[i+1]) {
			wavelengthBin=i;
			break;
		}
	}
	if (wavelength>=m_wavelengthBinLowEdge[kNWavelength-1]) wavelengthBin=kNWavelength-1;
	if (wavelengthBin==-1) {
		return 0.;
	} else {
		return m_meanPhotonFlux[wavelengthBin];
	}

}

double Star::getMeanFlux() const {

	double flux=0.;
	for (int i=0; i<kNWavelength; i++) flux += m_meanFlux[i];
	return flux;

}

double Star::getMeanPhotonFlux() const {

	double flux=0.;
	for (int i=0; i<kNWavelength; i++) flux += m_meanPhotonFlux[i];
	return flux;

}

void Star::initializeTimeSeries(int nTimeSteps) {
	for (int i=0; i<nTimeSteps; i++) {
		StarData * starData = new StarData();
		m_timeData.push_back(starData);
	}
}
