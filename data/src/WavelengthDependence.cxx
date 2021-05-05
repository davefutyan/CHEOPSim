/*
 * WavelengthDependence.cxx
 *
 *  Created on: 4 Feb 2016
 *      Author: futyand
 */

#include "REF_APP_Throughput.hxx"
#include "REF_APP_QE.hxx"
#include "REF_APP_SEDTeff.hxx"

#include "source/include/StarProducer.hxx"

#include "WavelengthDependence.hxx"

WavelengthDependence::WavelengthDependence(bool applyThroughput, string throughputFilename,
										   bool applyQE, string QEFilename, string SEDFilename, string targetSEDFilename) :
				 	 	 		 m_applyThroughput(applyThroughput), m_throughputFilename(throughputFilename),
								 m_applyQE(applyQE), m_QEFilename(QEFilename), m_SEDFilename(SEDFilename), m_targetSEDFilename(targetSEDFilename), m_VegaSpectrum(false) {

	for (int i=0; i<kNWavelength; i++) {
		m_wavelength[i] = Star::kWavelengthLowBound + 0.5*i;
		m_throughput[i] = 0.;
		m_transmission_Vband[i] = 0.;
		m_transmission_GaiaBand[i] = 0.;
		m_flux_Vega[i] = 0.;
		m_qe[i] = 0.;
		m_qeVsTempSlope[i] = 0.;
		m_weight[i] = 0.;
	}

    m_blackBody = (m_SEDFilename.find("REF_APP_SEDTeff") == string::npos); //Use black body if name of SED file does not include "REF_APP_SEDTeff" (e.g. empty string or "null")
    if (!m_blackBody) readSEDs();
	if (m_targetSEDFilename != "") readUserSED();

	readThroughput();
	readQuantumEfficiency();
    generateVbandTransmission();
    generateGaiaBandTransmission();
	generateVegaSpectrum();

}

double WavelengthDependence::getIntegratedRefBandTransmission(REFERENCE_BAND refBand, double effectiveTemperature, bool target) const {

	double photonSpectrum[kNWavelength]; //Array storing the stellar spectrum for each of the wavelength values in m_wavelength
    generatePhotonSpectrum(photonSpectrum,effectiveTemperature,target);

	double integratedRefBandTransmission = 0;
	for (unsigned i=0; i<kNWavelength; i++) {
		double weightedRefBandTransmission = photonSpectrum[i];
		//cout << m_wavelength[i] << "	" << photonSpectrum[i] << "	" << m_transmission_GaiaBand[i] << endl;
		weightedRefBandTransmission *= (refBand==GaiaBand ? m_transmission_GaiaBand[i] : m_transmission_Vband[i]);
		integratedRefBandTransmission += weightedRefBandTransmission;
	}
  	//cout << "Integrated reference band transmission: " << integratedRefBandTransmission << " for TEff = " << effectiveTemperature << endl;

	return integratedRefBandTransmission;

}

double WavelengthDependence::getIntegratedThroughputQE(double effectiveTemperature, double ccdTemperature, bool target) const {

	double photonSpectrum[kNWavelength]; //Array storing the stellar spectrum for each of the wavelength values in m_wavelength
    generatePhotonSpectrum(photonSpectrum,effectiveTemperature,target);

	double integratedThroughputQE = 0;
	for (unsigned i=0; i<kNWavelength; i++) {
		double weightedThroughputQE = photonSpectrum[i];
		if (m_applyThroughput) {
			weightedThroughputQE *= m_throughput[i];
		}
		if (m_applyQE) {
			//m_qeVsTempSlope value is ppm per milliKelvin
			double qe = m_qe[i] * (1 + (m_qeVsTempSlope[i]/1000.)*(ccdTemperature-SatelliteData::kDefaultCcdTemperature));
			//cout << ccdTemperature << " " << i << " " <<  m_qe[i] << " " << m_qeVsTempSlope[i] << " " << qe << endl;
			weightedThroughputQE *= qe;
		}
		integratedThroughputQE += weightedThroughputQE;
	}
  	//cout << "Integrated throughput * QE: " << integratedThroughputQE << " for TEff = " << effectiveTemperature << " and CCD temperature = " << ccdTemperature << endl;

	return integratedThroughputQE;

}

double WavelengthDependence::getIntegratedThroughput(double effectiveTemperature) const {

	double photonSpectrum[kNWavelength]; //Array storing the stellar spectrum for each of the wavelength values in m_wavelength
    generatePhotonSpectrum(photonSpectrum,effectiveTemperature);

	double integratedThroughput = 0;
	for (unsigned i=0; i<kNWavelength; i++) {
		double weightedThroughput = photonSpectrum[i];
		if (m_applyThroughput) weightedThroughput *= m_throughput[i];
		integratedThroughput += weightedThroughput;
	}
  	//cout << "Integrated throughput: " << integratedThroughput << " for TEff = " << effectiveTemperature << endl;

	return integratedThroughput;

}

void WavelengthDependence::generateWeights(double effectiveTemperature) {

	double photonSpectrum[kNWavelength]; //Array storing the stellar spectrum for each of the wavelength values in m_wavelength
    generatePhotonSpectrum(photonSpectrum,effectiveTemperature);

    double integral = 0.;
	for (int i=0; i<kNWavelength; i++) {
		m_weight[i] = photonSpectrum[i];
		if (m_applyThroughput) m_weight[i] *= m_throughput[i];
		if (m_applyQE) m_weight[i] *= m_qe[i];
		integral += m_weight[i];
	}

	//Normalization
	for (int i=0; i<kNWavelength; i++) {
		m_weight[i]/=integral;
		//cout << m_wavelength[i] << " " << m_weight[i]*20 << " " << integral << " " << effectiveTemperature << endl;
	}

}

double WavelengthDependence::getWeight(double wavelengthLow, double wavelengthHigh) const {

	if (wavelengthLow < m_wavelength[0] || wavelengthHigh > m_wavelength[kNWavelength-1]) {
		cout << m_wavelength[0] << " " << m_wavelength[kNWavelength-1] << endl;
		throw runtime_error("Error in WavelengthDependence::getWeight: wavelengths must be in range 330 to 1100, values are "+
				to_string(wavelengthLow)+","+to_string(wavelengthHigh));
	}

	double weight  = 0.;
	for (int i=0; i<kNWavelength; i++) {
		if (m_wavelength[i] > wavelengthLow &&  m_wavelength[i] <= wavelengthHigh) weight += m_weight[i];
	}
	return weight;

}

int WavelengthDependence::getBin(double wavelength) const {

	int bin = round((wavelength - Star::kWavelengthLowBound)*2.);
	if (bin<0 || bin>kNWavelength) {
		throw runtime_error("Error in WavelengthDependence::getBin: wavelength must be in range 330 to 1100, value is "+to_string(wavelength));
	}
	return bin;

}

void WavelengthDependence::readThroughput() {

	double throughput_fromFile[2000], wavelength_fromFile[2000];
	unsigned size_fromFile = 0;

    if (m_throughputFilename.substr(m_throughputFilename.find_last_of(".")+1) == "txt") {

		//user uploaded throughput file, to be read from execution directory
		ifstream throughput_in(m_throughputFilename, ios::in);
	    if (!throughput_in.good()) throw runtime_error("Error in WavelengthDependence::readThroughput: Error opening file "+m_throughputFilename);
	    string line;
	    while (getline(throughput_in, line)) {
	    	if(line[0] != '#') {
				if (size_fromFile>=2000) throw runtime_error("ERROR: Throughput file size greater than maximum 2000 rows");
	    		stringstream(line) >> wavelength_fromFile[size_fromFile] >> throughput_fromFile[size_fromFile];
	    		if (throughput_fromFile[size_fromFile] > 1. || throughput_fromFile[size_fromFile] < 0.) {
	    			throw runtime_error("Error in WavelengthDependence::readThroughput: Read invalid value from "+m_throughputFilename
	    								+": Values for throughput must be between 0 an 1");
	    		}
				size_fromFile++;
	    	}
		}
		throughput_in.close();

	} else {

		//Read in the telescope throughput
		string filename = string(getenv("CHEOPS_SW"))+"/resources/"+m_throughputFilename;
		RefAppThroughput * throughput_file = new RefAppThroughput(filename,"READONLY");
		while (throughput_file->ReadRow()) {
			if (size_fromFile>=2000) throw runtime_error("ERROR: throughput file size greater than maximum 2000 rows");
			wavelength_fromFile[size_fromFile] = throughput_file->getCellWavelength();
			throughput_fromFile[size_fromFile] = throughput_file->getCellThroughput();
			if (throughput_fromFile[size_fromFile]>1.1) throw runtime_error("ERROR: invalid value read from throughput file");
			size_fromFile++;
		}
		delete throughput_file;

		if (m_applyThroughput) {
			ofstream referenceFilesList("reference_files.txt", ios::app);
			if (m_throughputFilename == "throughput_BoL_april2017.fits") {
				referenceFilesList << "CH_TU1970-01-01T00-00-00_REF_APP_Throughput-BeginOfLife_V0001.fits" << endl;
			} else if (m_throughputFilename == "throughput_EoL_april2017.fits") {
				referenceFilesList << "CH_TU1970-01-01T00-00-00_REF_APP_Throughput-EndOfLife_V0002.fits" << endl;
			} else {
				referenceFilesList << m_throughputFilename << endl;
			}
			referenceFilesList.close();
		}

	}

	linearInterpolation(wavelength_fromFile,throughput_fromFile,size_fromFile,m_wavelength,m_throughput,kNWavelength);
	//for (unsigned i=0; i<size_fromFile; i++)  cout << wavelength_fromFile[i] << " " << throughput_fromFile[i] << endl;
	//for (unsigned i=0; i<kNWavelength; i++)  cout << m_wavelength[i] << " " << m_throughput[i] << endl;

}

void WavelengthDependence::generateVbandTransmission() {

	//Define the V band transmission curve, as defined in Buser & Kurucz (1978)
	double wavelength_Vband_ref[54], transmission_Vband_ref[54];
	for (int i=0; i<54; i++) {
		wavelength_Vband_ref[i] = 475. + 5.*i;  //Units: nm
		transmission_Vband_ref[i] = VbandTransmission(i);
	}

	//	ifstream vband_in("/Users/futyand/CHEOPS/V_CHEOPS/Bessel_V.txt", ios::in);
	//    if (!vband_in.good()) throw runtime_error("Error in WavelengthDependence::generateVbandTransmission Error opening Bessel_V.txt");
	//	double wavelength_Vband_ref[1801], transmission_Vband_ref[1801];
	//	unsigned i=0;
	//	double dummy;
	//	while (vband_in.good()) {
	//		vband_in >> dummy;
	//		if(vband_in.eof()) break;
	//		wavelength_Vband_ref[i] = dummy;
	//		vband_in >> transmission_Vband_ref[i];
	//		i++;
	//	}
	//	for (unsigned k=0;k<1801;k++) cout << wavelength_Vband_ref[k] << " " << transmission_Vband_ref[k] << endl;

	//Interpolate V band transmission spectrum to a 0.5 nm fixed step array
	StarProducer::interpolateFlux(wavelength_Vband_ref,transmission_Vband_ref,54,m_wavelength,m_transmission_Vband,kNWavelength);
	for (int i=0;i<kNWavelength;i++) {
		if (m_wavelength[i]<475. || m_wavelength[i]>=740.) m_transmission_Vband[i]=0.; //Exclude bad extrapolations outside the V band limits
	}

}

void WavelengthDependence::generateGaiaBandTransmission() {

	//Read in the Gaia passband
	double wavelength_GaiaBand_ref[3195], transmission_GaiaBand_ref[3195];
	ifstream gaia_passband_in(string(getenv("CHEOPS_SW"))+"/resources/GaiaDR2_RevisedPassbands.dat", ios::in);
	if (!gaia_passband_in.good()) throw runtime_error("Error in WavelengthDependence::generateGaiaBandTransmission: Error opening Gaia passband file");
	string line;
	unsigned size_fromFile = 0;
	while (getline(gaia_passband_in, line)) {
		if(line[0] != '#') {
			stringstream(line) >> wavelength_GaiaBand_ref[size_fromFile] >> transmission_GaiaBand_ref[size_fromFile];
			if (transmission_GaiaBand_ref[size_fromFile] > 1. || transmission_GaiaBand_ref[size_fromFile] < 0.) {
				throw runtime_error("Error in StarProducer::initialize: Read invalid value "+to_string(size_fromFile)+" "+to_string(transmission_GaiaBand_ref[size_fromFile])+" from Gaia passband file: Values for must be between 0 and 1");
			}
			size_fromFile++;
			if (size_fromFile>3194) break;
		}
	}
	gaia_passband_in.close();

	//Interpolate Gaia band transmission spectrum to a 0.5 nm fixed step array
	StarProducer::interpolateFlux(wavelength_GaiaBand_ref,transmission_GaiaBand_ref,3195,m_wavelength,m_transmission_GaiaBand,kNWavelength);

}

void WavelengthDependence::generateVegaSpectrum() {

	//Read in Vega spectrum from http://dls.physics.ucdavis.edu/calib/Vega.sed   Units: ergs/s/cm2/A
    double wavelength_Vega_fromFile[1141], flux_Vega_fromFile[1141];
    StarProducer::readVegaSpectrum(wavelength_Vega_fromFile,flux_Vega_fromFile);

    //convert from Angstroms to nm
	for (int i=0; i<1141; i++) wavelength_Vega_fromFile[i] /= 10.;

	//Interpolate Vega flux to a 0.5nm fixed step array
	StarProducer::interpolateFlux(wavelength_Vega_fromFile,flux_Vega_fromFile,1141,m_wavelength,m_flux_Vega,kNWavelength);

}

void WavelengthDependence::readQuantumEfficiency() {

    double wavelength_fromFile[2000], qe_fromFile[2000], qeVsTempSlope_fromFile[2000];
    unsigned size_fromFile = 0;

    if (m_QEFilename.substr(m_QEFilename.find_last_of(".")+1) == "txt") {

		//user uploaded qe file, to be read from execution directory
		ifstream qe_in(m_QEFilename, ios::in);
	    if (!qe_in.good()) throw runtime_error("Error in WavelengthDependence::readQuantumEfficiency: Error opening file "+m_QEFilename);
	    string line;
	    while (getline(qe_in, line)) {
	    	if(line[0] != '#') {
				if (size_fromFile>=2000) throw runtime_error("ERROR: QE file size greater than maximum 2000 rows");
	    		stringstream(line) >> wavelength_fromFile[size_fromFile] >> qe_fromFile[size_fromFile];
	    		if (qe_fromFile[size_fromFile] > 1. || qe_fromFile[size_fromFile] < 0.) {
	    			throw runtime_error("Error in WavelengthDependence::readQuantumEfficiency: Read invalid value from "+m_QEFilename
	    								+": Values for qe must be between 0 an 1");
	    		}
				qeVsTempSlope_fromFile[size_fromFile] = 0.;
				size_fromFile++;
	    	}
	    }
		qe_in.close();

	} else {

		//Read in the QE
		string filename = string(getenv("CHEOPS_SW"))+"/resources/"+m_QEFilename;
		RefAppQe * qe_file = new RefAppQe(filename,"READONLY");
		while (qe_file->ReadRow()) {
			if (size_fromFile>=2000) throw runtime_error("ERROR: QE file size greater than maximum 2000 rows");
			wavelength_fromFile[size_fromFile] = qe_file->getCellWavelength();
			qe_fromFile[size_fromFile] = qe_file->getCellQe();
			if (!qe_file->isNullQeVsTempSlope()) {
				qeVsTempSlope_fromFile[size_fromFile] = qe_file->getCellQeVsTempSlope();
			} else {
				qeVsTempSlope_fromFile[size_fromFile] = 0.;
			}
			if (qe_fromFile[size_fromFile]>1.1) throw runtime_error("ERROR: invalid value read from qe file");
			size_fromFile++;
		}
		delete qe_file;

		if (m_applyQE) {
			ofstream referenceFilesList("reference_files.txt", ios::app);
			if (m_QEFilename == "qe_may2016.fits") {
				referenceFilesList << "CH_TU1970-01-01T00-00-00_REF_APP_QE_V0000.fits" << endl;
			} else {
				referenceFilesList << m_QEFilename << endl;
			}
			referenceFilesList.close();
		}

	}

	linearInterpolation(wavelength_fromFile,qe_fromFile,size_fromFile,m_wavelength,m_qe,kNWavelength);
	linearInterpolation(wavelength_fromFile,qeVsTempSlope_fromFile,size_fromFile,m_wavelength,m_qeVsTempSlope,kNWavelength);
	//for (unsigned i=0; i<size_fromFile; i++)  cout << wavelength_fromFile[i] << " " << qe_fromFile[i] << " " << qeVsTempSlope_fromFile[i] << endl;
	//for (unsigned i=0; i<kNWavelength; i++)  cout << m_wavelength[i] << " " << m_qe[i] << " " << m_qeVsTempSlope[i] << endl;

}

void WavelengthDependence::readSEDs() {

	double wavelength_fromFile[2000];

	//Read in the spectral energy distributions
	string filename = string(getenv("CHEOPS_SW"))+"/resources/"+m_SEDFilename;
	RefAppSedteff * sed_file = new RefAppSedteff(filename,"READONLY");
	unsigned row = 0;
	while (sed_file->ReadRow()) {
		if (row==0) {
			std::vector<float> wavelengths = sed_file->getCellVectorWavelength();
			if (wavelengths.size() >= 2000) throw runtime_error("ERROR: REF_APP_SEDTeff file number of wavelength bins greater than maximum 2000");
			for (unsigned i=0; i<wavelengths.size(); i++) wavelength_fromFile[i] = wavelengths[i];
		}
		m_SEDTeff[row] = sed_file->getCellTemperatur();
		double flux_fromFile[2000], flux[kNWavelength];
		std::vector<float> fluxes = sed_file->getCellVectorFlux();
		for (unsigned i=0; i<fluxes.size(); i++) flux_fromFile[i] = fluxes[i];
		linearInterpolation(wavelength_fromFile,flux_fromFile,fluxes.size(),m_wavelength,flux,kNWavelength);
		for (unsigned i=0; i<kNWavelength; i++) m_SEDFlux[i][row] = flux[i];
		row++;
		if (row>kNTEff) throw runtime_error("ERROR: Number of rows in REF_APP_SEDTeff file exceeds the expected "+to_string(kNTEff));
	}
	delete sed_file;
	if (row!=kNTEff) throw runtime_error("ERROR: Number of rows in REF_APP_SEDTeff file is fewer than the expected "+to_string(kNTEff));

//	for (unsigned teff=0; teff<kNTEff; teff++) {
//		fstream sed_file("sed_"+to_string(m_SEDTeff[teff])+".dat", ios::app);
//		cout << "Effective temperature: " << m_SEDTeff[teff] << endl;
//		for (unsigned i=0; i<kNWavelength; i++) sed_file << m_SEDFlux[i][teff] << endl;
//		sed_file.close();
//	}

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << m_SEDFilename << endl;
	referenceFilesList.close();

}

void WavelengthDependence::readUserSED() {

	double flux[kNWavelength];
	double wavelength_fromFile[2000],flux_fromFile[2000];
	double wavelength_in,flux_in;
	vector<double> wavelengths,fluxes;

	//Read in the spectral energy distribution
	FitsDalTable * sed_file = new FitsDalTable(m_targetSEDFilename,"READONLY");
	int tableLength = sed_file->GetAttr<int>("NAXIS2");
	sed_file->Assign("WAVELENGTH", &wavelength_in, tableLength);
	sed_file->Assign("FLUX", &flux_in, tableLength);

	while (sed_file->ReadRow()) {
		if (wavelength_in>=3200. && wavelength_in<11100.) {
			wavelengths.push_back(wavelength_in/10.);
			fluxes.push_back(flux_in);
		}
	}

	if (wavelengths.size() >= 2000) throw runtime_error("ERROR: User SED file number of wavelength bins between 320 and 1110nm greater than maximum 2000: "+to_string(wavelengths.size()));
	if (wavelengths.front()>330. || wavelengths.back()<1100.) throw runtime_error("ERROR: User SED file only covers "+to_string(wavelengths.front())+" to "+to_string(wavelengths.back())+"nm");
	//cout << wavelengths.size() << " " << to_string(wavelengths.front()) << " " << to_string(wavelengths.back()) << endl;
	for (unsigned i=0; i<wavelengths.size(); i++) {
		wavelength_fromFile[i] = wavelengths[i];
		flux_fromFile[i] = fluxes[i];
	}

	linearInterpolation(wavelength_fromFile,flux_fromFile,2000,m_wavelength,flux,kNWavelength);
	for (unsigned i=0; i<kNWavelength; i++) m_targetSEDFlux[i] = flux[i];
	delete sed_file;

//	fstream sed_out("userSED.dat", ios::app);
//	for (unsigned i=0; i<kNWavelength; i++) sed_out << m_targetSEDFlux[i] << endl;
//	sed_out.close();

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << m_targetSEDFilename << endl;
	referenceFilesList.close();

}

void WavelengthDependence::generatePhotonSpectrum(double (&photonSpectrum)[kNWavelength], double effectiveTemperature, bool target) const {

	//Get the photon spectrum for the star in 0.5nm steps
	double spectrumIntegral = 0.;
	for (int i=0; i<kNWavelength; i++) {
		//Energy per photon is inversely proportional to wavelength, so number or photons is proportional to the spectrum energy mulitplied by the wavelength
		if (m_VegaSpectrum || effectiveTemperature==0.) { //If setVegaFluxSpectrum has been called by StarProducer, or if TEff is 0., use Vega spectrum
			photonSpectrum[i] = m_flux_Vega[i]*m_wavelength[i];
		} else if (effectiveTemperature<0.) { //If TEff is negative, use a flat weight for all wavelengths
			photonSpectrum[i] = 1.;
		} else if (target && m_targetSEDFilename != "") { //User provided target SED spectrum
			photonSpectrum[i] = m_targetSEDFlux[i]*m_wavelength[i];
		} else if (m_blackBody || effectiveTemperature > m_SEDTeff[kNTEff-1]) { //Black body spectrum (SEDs only available up to 7200K so use black body for higher Teffs)
			double blackBody = planckRadiance(m_wavelength[i]*10.,effectiveTemperature);
			photonSpectrum[i] = blackBody*m_wavelength[i];
		} else if (m_SEDFilename.find("REF_APP_SEDTeff") != string::npos) { //Reference SED spectrum
//			unsigned nearest_iTEff = 0;
//			double min_temp_diff = 1.E9;
//			for (unsigned j=0; j<kNTEff; j++) {
//				double temp_diff =  fabs(effectiveTemperature-m_SEDTeff[j]);
//				if (temp_diff < min_temp_diff) {
//					min_temp_diff = temp_diff;
//					nearest_iTEff = j;
//				}
//			}
//			photonSpectrum[i] = m_SEDFlux[i][nearest_iTEff]*m_wavelength[i];
			unsigned iTeff = 0;
			double mu = 0.;
			for (unsigned j=0; j<kNTEff; j++) {
				if (m_SEDTeff[j] >= effectiveTemperature) {
					if (j==0) { //For Teff equal to or lower than the lowest Teff for which there is a SED available, use the lowest available SED
						iTeff = j+1;
						mu = 0;
					} else {
						iTeff = j;
						mu = (effectiveTemperature - m_SEDTeff[j-1]) / (m_SEDTeff[j] - m_SEDTeff[j-1]);
					}
					//if (i==0) cout << effectiveTemperature << " " << iTeff << " "  << m_SEDTeff[iTeff-1] << " "  << m_SEDTeff[iTeff] << " " << mu << endl;
					break;
				}
			}
			if (iTeff==0) throw runtime_error("Error in WavelengthDependence::generatePhotonSpectrum: cannot find SED corresponding to effective temperature"+to_string(effectiveTemperature));
			photonSpectrum[i] = (m_SEDFlux[i][iTeff-1]*(1.-mu)+m_SEDFlux[i][iTeff]*mu)*m_wavelength[i];
		}
		spectrumIntegral += photonSpectrum[i];
	}

	//Normalization
	for (int i=0; i<kNWavelength; i++) {
		photonSpectrum[i]/=spectrumIntegral;
		//cout << m_wavelength[i] << "	" << photonSpectrum[i] << endl;
	}

}

void WavelengthDependence::getStellarEnergySpectrum_1nm(double (&energySpectrum)[Star::kNWavelength*10], double effectiveTemperature, bool target) const {

	//Get the energy spectrum for the star in 1nm steps as expected by StarProducer. Since the arrays are in 0.5nm steps, take every second element [i*2]
	double spectrumIntegral = 0.;
	for (int i=0; i<Star::kNWavelength*10; i++) {
		if (m_VegaSpectrum || effectiveTemperature==0.) { //If setVegaFluxSpectrum has been called by StarProducer, or if TEff is 0., use Vega spectrum
			if (i==0) cout << "   Using Vega spectrum" << endl;
			energySpectrum[i] = m_flux_Vega[i*2];
		} else if (effectiveTemperature<0.) { //If TEff is negative, use a flat weight for all wavelengths
			if (i==0) cout << "   Using flat spectrum" << endl;
			energySpectrum[i] = 1.;
		} else if (target && m_targetSEDFilename != "") { //User provided target SED spectrum
			if (i==0) cout << "   Using user SED spectrum" << endl;
			energySpectrum[i] = m_targetSEDFlux[i*2];
		} else if (m_blackBody || effectiveTemperature > m_SEDTeff[kNTEff-1]) { //Black body spectrum (SEDs only available up to 7200K so use black body for higher Teffs)
			if (i==0) cout << "   Using Planck spectrum" << endl;
			energySpectrum[i] = planckRadiance(m_wavelength[i*2]*10.,effectiveTemperature);
		} else if (m_SEDFilename.find("REF_APP_SEDTeff") != string::npos) { //Reference SED spectrum
			if (i==0) cout << "   Using reference SED spectrum" << endl;
			unsigned iTeff = 0;
			double mu = 0.;
			for (unsigned j=0; j<kNTEff; j++) {
				if (m_SEDTeff[j] >= effectiveTemperature) {
					if (j==0) { //For Teff equal to or lower than the lowest Teff for which there is a SED available, use the lowest available SED
						iTeff = j+1;
						mu = 0;
					} else {
						iTeff = j;
						mu = (effectiveTemperature - m_SEDTeff[j-1]) / (m_SEDTeff[j] - m_SEDTeff[j-1]);
					}
					break;
				}
			}
			if (iTeff==0) throw runtime_error("Error in WavelengthDependence::getStellarEnergySpectrum_1nm: cannot find SED corresponding to effective temperature"+to_string(effectiveTemperature));
			energySpectrum[i] = m_SEDFlux[i*2][iTeff-1]*(1.-mu)+m_SEDFlux[i*2][iTeff]*mu;
		}
		spectrumIntegral += energySpectrum[i];
	}

	//Normalization
	for (int i=0; i<Star::kNWavelength*10; i++) energySpectrum[i]/=spectrumIntegral;

}

void WavelengthDependence::linearInterpolation(const double (&x_ref)[2000], const double (&y_ref)[2000], unsigned size_ref, double* x, double* y, unsigned size) const {

	for (unsigned i=0; i<size; i++) {
		for (unsigned j=1; j<size_ref; j++) {
			if (x_ref[j]>=x[i]) {
				double mu = (x[i]-x_ref[j-1])/(x_ref[j]-x_ref[j-1]);
				y[i] = y_ref[j-1]*(1.-mu)+y_ref[j]*mu;
				break;
			}
		}
	}

}

double WavelengthDependence::planckRadiance(double wavelength, double effectiveTemperature) {
	auto h = kPlanck;
	auto c = kSpeedOfLight;
	auto k = kBoltzmann;
    return (2.*h*c*c/pow(wavelength*1.e-10,5)) / (exp(h*c/(wavelength*1.e-10*k*effectiveTemperature))-1.);
}

void WavelengthDependence::setVegaFluxSpectrum(double flux_Vega[kNWavelength]) {

	for (int i=0; i<kNWavelength; i++) {
		m_flux_Vega[i] = flux_Vega[i];
		m_VegaSpectrum=true;
	}

}
