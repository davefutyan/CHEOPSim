/*
 * StarProducer.cxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

#include <boost/algorithm/string.hpp>
#include "boost/date_time/gregorian/gregorian_types.hpp"

#include "CreateFitsFile.hxx"
#include "EXT_PRE_StarCatalogue.hxx"
#include "EXT_DRFT_StarCatalogue.hxx"
#include "REF_APP_GainCorrection.hxx"

#include "StarProducer.hxx"
#include "telescope/include/PSFGenerator.hxx"

using namespace std;

void StarProducer::initialize(const ModuleParams & params) {

	m_VegaSpectrum = false; //Set to true to use the spectrum of Vega rather than a black body

	m_FOVradius = params.GetAsDouble("FOVradius");
	m_minMagnitude = params.GetAsDouble("minMagnitude");
	m_starFilename = params.GetAsString("starFilename");
	string targetMagnitude = params.GetAsString("targetMagnitude");
	if (targetMagnitude=="") {
		m_overrideTargetMagnitude = false;
		m_targetMagnitude = -999.;
	} else {
		m_overrideTargetMagnitude = true;
		m_targetMagnitude = atof(targetMagnitude.c_str());
	}
	m_gaiaBand = params.GetAsBool("gaiaBand");
	m_useTargetStarForPointing = params.GetAsBool("useTargetStarForPointing");
	m_doBackgroundStars = params.GetAsBool("doBackgroundStars");
	string fieldCrowding = params.GetAsString("fieldCrowding");

	if (fieldCrowding == "low") {
		m_fieldCrowding = FIELD_CROWDING::low;
	} else if (fieldCrowding == "medium") {
		m_fieldCrowding = FIELD_CROWDING::medium;
	} else if (fieldCrowding == "extreme") {
		m_fieldCrowding = FIELD_CROWDING::extreme;
	} else {
		throw runtime_error("Error in StarProducer::initialize: invalid value for fieldCrowding parameter: "+fieldCrowding);
	}

	cout << "Field of view radius (arcseconds)                              :" << m_FOVradius << endl;
	cout << "Star list filename                                             :" << m_starFilename << endl;
	if (m_doBackgroundStars) cout << "Background star field crowding                                 :" << fieldCrowding << endl;

	//Define wavelength bins in Angstroms
	for (int i=0; i<Star::kNWavelength*10; i++) m_wavelength[i] = Star::kWavelengthLowBound*10. + 10.*i; //10 Angstrom steps
	for (int i=0; i<WavelengthDependence::kNWavelength; i++) m_wavelength_fine[i] = Star::kWavelengthLowBound*10 + 5*i; //0.5nm steps

	//Read in Vega spectrum from http://dls.physics.ucdavis.edu/calib/Vega.sed   Units: ergs/s/cm2/A
    double wavelength_Vega_fromFile[1141], flux_Vega_fromFile[1141];
    readVegaSpectrum(wavelength_Vega_fromFile,flux_Vega_fromFile);

	//Interpolate Vega flux to a 10 Angstrom fixed step array
	interpolateFlux(wavelength_Vega_fromFile,flux_Vega_fromFile,1141,m_wavelength,m_flux_Vega,Star::kNWavelength*10);
	//Interpolate Vega flux to a 0.5nm fixed step array
	interpolateFlux(wavelength_Vega_fromFile,flux_Vega_fromFile,1141,m_wavelength_fine,m_flux_Vega_fine,WavelengthDependence::kNWavelength);

	//Initialize the FluxConverter class to convert between Gaia magnitude and CHEOPS magnitude.
	//The parameters are determined using CHEOPSim/resources/FluxConversion.ipynb
	m_fluxConv = new FluxConverter();

}

void StarProducer::doBegin(Data * data, bool fullFrame) {

	// Read in the V-band or Gaia band transmission curve
	if (!m_gaiaBand) {

		//Define the V band transmission curve, as defined in Buser & Kurucz (1978)
		double wavelength_Vband[54], transmission_Vband[54];
		for (int i=0; i<54; i++) {
			wavelength_Vband[i] = 4750. + 50.*i;  //Units: Angstroms
			transmission_Vband[i] = WavelengthDependence::VbandTransmission(i);
		}

		//Interpolate V band transmission spectrum to a 10 Angstrom fixed step array
		interpolateFlux(wavelength_Vband,transmission_Vband,54,m_wavelength,m_transmission_refBand,Star::kNWavelength*10);

		for (int i=0;i<Star::kNWavelength*10;i++) {
			if (m_wavelength[i]<4750 || m_wavelength[i]>=7400) m_transmission_refBand[i]=0.; //Exclude bad extrapolations outside the V band limits
			//cout << "m_wavelength[" << i << "]=" << m_wavelength[i] << ", m_flux_Vega[" << i << "]=" << m_flux_Vega[i]
			//     << ", m_transmission_refBand[" << i << "]=" << m_transmission_refBand[i] << endl;
		}

	} else {

		double wavelength_Gband[3195], transmission_Gband[3195];
		ifstream gaia_passband_in(string(getenv("CHEOPS_SW"))+"/resources/GaiaDR2_RevisedPassbands.dat", ios::in);
		if (!gaia_passband_in.good()) throw runtime_error("Error in StarProducer::doBegin: Error opening Gaia passband file");
		string line;
		unsigned size_fromFile = 0;
		while (getline(gaia_passband_in, line)) {
			if(line[0] != '#') {
				stringstream(line) >> wavelength_Gband[size_fromFile] >> transmission_Gband[size_fromFile];
				wavelength_Gband[size_fromFile] *= 10.;
				if (transmission_Gband[size_fromFile] > 1. || transmission_Gband[size_fromFile] < 0.) {
					throw runtime_error("Error in StarProducer::initialize: Read invalid value "+to_string(size_fromFile)+" "+to_string(transmission_Gband[size_fromFile])+" from Gaia passband file: Values for must be between 0 and 1");
				}
				size_fromFile++;
				if (size_fromFile>3194) break;
			}
		}
		gaia_passband_in.close();
		interpolateFlux(wavelength_Gband,transmission_Gband,3195,m_wavelength,m_transmission_refBand,Star::kNWavelength*10);

	}

//	ofstream file_out;
//	file_out.open("fluxVsWavelength.csv");
//	for (int i=0;i<Star::kNWavelength*10;i++) {
//		file_out << m_wavelength[i] << "," << m_flux_Vega[i] << "," << m_transmission_refBand[i] << endl;
//	}
//	file_out.close();

	//Override the configuration file FOV radius for the full frame case
	//to radius of circle enclosing 1024x1024 plus a 25 pixel margin
	double fovRadius = m_FOVradius;
	if (fullFrame) fovRadius = 759.;

	//Set the field of view radius and magnitude threshold
	SkyFieldOfView * fov = data->getFieldOfView();
	fov->setFOVradius(fovRadius);
	fov->setMinMagnitude(m_minMagnitude);

	if (m_VegaSpectrum) {
		//Set the Vega flux spectrum in the WavelengthDependence class
		data->getWavelengthDependence()->setVegaFluxSpectrum(m_flux_Vega_fine);
	}

	unsigned nTimeSteps = data->getTimeConfiguration().getNumberOfTimeSteps();

	//Read in the foreground stars
	bool writeDrftStarCatalogue = true;
	fov->resetStars();
	const WavelengthDependence * wavelengthDependence = data->getWavelengthDependence();
	MpsPreVisits * mpsVisits = data->getMpsPreVisits();
	if (m_starFilename.substr(m_starFilename.find_last_of(".")+1) == "fits") {
		if (m_starFilename.find("EXT_PRE_StarCatalogue") != string::npos) {
			readStarCatalogue<ExtPreStarcatalogue>(fov,wavelengthDependence,nTimeSteps,mpsVisits);
		} else if (m_starFilename.find("EXT_DRFT_StarCatalogue") != string::npos) {
			readStarCatalogue<ExtDrftStarcatalogue>(fov,wavelengthDependence,nTimeSteps,mpsVisits);
			writeDrftStarCatalogue = false; //No need to write an EXT_DRFT_StarCatalogue output if it has been provided as input
		} else {
			throw runtime_error("Error in StarProducer::doBegin: Name of star catalogue FITS file must contain EXT_PRE_StarCatalogue or EXT_DRFT_StarCatalogue.");
		}
	} else {
		readStars(fov,wavelengthDependence,nTimeSteps,mpsVisits);
	}

	if (m_useTargetStarForPointing && fov->getStars().size() == 0) {
		throw runtime_error("Error in StarProducer::doBegin: the parameter useTargetStarForPointing has been set to true, but no target star has been defined.");
	}

	//Read in the background stars
	if (m_doBackgroundStars) readBackgroundStars(fov,wavelengthDependence,nTimeSteps);

	if (fov->getStars().size() == 0) {
		throw runtime_error("Error in StarProducer::doBegin: There are no stars in the field of view. Please check that the RA and Dec of the target star is consistent with the RA and Dec of the pointing direction.");
	}

	//Set the Gaia magnitude error for the target star, if a non-zero value was read from MPS_PRE_Visits in Simulator::initialize()
	if (data->getVisit().m_gaiaMagError>0.) {
		fov->getStars()[0]->setGaiaMagnitudeError(data->getVisit().m_gaiaMagError);
	}

	//Write the list of stars to an EXT_PRE_StarCatalog file
	writeStars<ExtPreStarcatalogue>(data,data->getVisit().m_visitCtr);
	//If there was no EXT_DRFT_StarCatalogue input file , also write the list of stars to an EXT_DRFT_StarCatalog file
	if (writeDrftStarCatalogue) writeStars<ExtDrftStarcatalogue>(data,0);
}

template <class T> void StarProducer::readStarCatalogue(SkyFieldOfView* fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps, MpsPreVisits * mpsVisits) const {

	T * catalogue = new T(m_starFilename,"READONLY");

	unsigned iStar = 0;
	while (catalogue->ReadRow()) {
		double ra = catalogue->getCellRa()*3600.;
		double raErr = catalogue->isNullRaErr() ? 0. : catalogue->getCellRaErr()*3600.;
		double dec = catalogue->getCellDec()*3600.;
		double decErr = catalogue->isNullDecErr() ? 0. : catalogue->getCellDecErr()*3600.;
		double effectiveTemperature = catalogue->getCellTEff();
		double effectiveTemperatureErr = catalogue->isNullTEffErr() ? 0. : catalogue->getCellTEffErr();
		double gaiaMag = catalogue->getCellMagGaia();
		double gaiaMagErr = catalogue->isNullMagGaiaErr() ? 0. : catalogue->getCellMagGaiaErr();
		string id = catalogue->getCellId();
		if (std::isnan(effectiveTemperature)) effectiveTemperature = 5660.; //Default to effective temperature for G5 star
		if (std::isnan(ra) || std::isnan(dec) || std::isnan(gaiaMag)) {
			throw runtime_error("ERROR reading star catalogue file: value is nan for RA, DEC or MAG");
		}

		double VbandMag,cheopsMag;
		if (iStar == 0 && m_overrideTargetMagnitude) {
			if (m_gaiaBand) {
				gaiaMag = m_targetMagnitude;
				VbandMag = gaiaMag - GaiaMinusVbandMagnitude(effectiveTemperature);
				cheopsMag = gaiaMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,effectiveTemperature,wavelengthDependence,iStar==0);
			} else{
				VbandMag = m_targetMagnitude;
				gaiaMag = VbandMag + GaiaMinusVbandMagnitude(effectiveTemperature);
				cheopsMag = VbandMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::Vband,effectiveTemperature,wavelengthDependence,iStar==0);
			}
		} else {
			VbandMag = gaiaMag - GaiaMinusVbandMagnitude(effectiveTemperature);
			cheopsMag = gaiaMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,effectiveTemperature,wavelengthDependence,iStar==0);
		}

		Star * star = new Star(SkyPosition(ra,dec,raErr,decErr),VbandMag,cheopsMag,gaiaMag,Star::Undefined,effectiveTemperature,id);

		//Set the pointing direction to match the target star if requested
		if (iStar == 0 && m_useTargetStarForPointing) {
			fov->setPointingDirection(star->getSkyPosition());
			mpsVisits->setCellRaTarg(ra/3600.);
			mpsVisits->setCellDecTarg(dec/3600.);
		}

		if (fov->isInside(star->getSkyPosition()) && gaiaMag<m_minMagnitude) {
			cout << "star " << iStar << (iStar==0 ? " (target star): ":": ") << endl;
			cout << "   V-band magnitude " << VbandMag << ", CHEOPS magnitude " << cheopsMag << ", Gaia magnitude " << gaiaMag << ", spectral type " << star->getSpectralTypeString() << ", RA=" << ra/3600. << ", dec=" << dec/3600. << endl;
			evaluateFlux(star,wavelengthDependence,iStar==0);
			star->initializeTimeSeries(nTimeSteps);
			star->setGaiaMagnitudeError(gaiaMagErr);
			if (iStar == 0) star->setCheopsMagnitudeError(catalogue->getKeyMagCerr());
			star->setEffectiveTemperatureError(effectiveTemperatureErr);
			fov->addStar(star);
			iStar++;
		}

	}

	delete catalogue;

}

void StarProducer::readStars(SkyFieldOfView* fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps, MpsPreVisits * mpsVisits) const {

	string star_ra,star_dec;
	double gaiaMag;
	string star_specType_str;

	ifstream stars_in;
	if (m_starFilename=="stars.txt" || m_starFilename=="onestar.txt" || m_starFilename=="twostars.txt" ||
	    m_starFilename=="fivestars.txt" || m_starFilename=="nostars.txt") {
		//default star file to be read from resources directory
		stars_in.open(string(getenv("CHEOPS_SW"))+"/resources/"+m_starFilename, ios::in);
	} else {
		//user defined star file, to be read from execution directory
		stars_in.open(m_starFilename, ios::in);
	}
    if (!stars_in.good()) throw runtime_error("Error in StarProducer: Error opening file "+m_starFilename);

    int iStar=-1;
	while (stars_in.good()) {
		iStar++;

		stars_in >> star_ra;
		if(stars_in.eof()) break;
		stars_in >> star_dec;
		stars_in >> gaiaMag;
		stars_in >> star_specType_str;

		Star::spectralType star_specType;
		if (star_specType_str=="O7") star_specType=Star::O7;
		else if (star_specType_str=="O8") star_specType=Star::O8;
		else if (star_specType_str=="O9") star_specType=Star::O9;
		else if (star_specType_str=="B0") star_specType=Star::B0;
		else if (star_specType_str=="B1") star_specType=Star::B1;
		else if (star_specType_str=="B2") star_specType=Star::B2;
		else if (star_specType_str=="B3") star_specType=Star::B3;
		else if (star_specType_str=="B4") star_specType=Star::B4;
		else if (star_specType_str=="B5") star_specType=Star::B5;
		else if (star_specType_str=="B6") star_specType=Star::B6;
		else if (star_specType_str=="B7") star_specType=Star::B7;
		else if (star_specType_str=="B8") star_specType=Star::B8;
		else if (star_specType_str=="B9") star_specType=Star::B9;
		else if (star_specType_str=="A0") star_specType=Star::A0;
		else if (star_specType_str=="A1") star_specType=Star::A1;
		else if (star_specType_str=="A2") star_specType=Star::A2;
		else if (star_specType_str=="A3") star_specType=Star::A3;
		else if (star_specType_str=="A4") star_specType=Star::A4;
		else if (star_specType_str=="A5") star_specType=Star::A5;
		else if (star_specType_str=="A6") star_specType=Star::A6;
		else if (star_specType_str=="A7") star_specType=Star::A7;
		else if (star_specType_str=="A8") star_specType=Star::A8;
		else if (star_specType_str=="A9") star_specType=Star::A9;
		else if (star_specType_str=="F0") star_specType=Star::F0;
		else if (star_specType_str=="F1") star_specType=Star::F1;
		else if (star_specType_str=="F2") star_specType=Star::F2;
		else if (star_specType_str=="F3") star_specType=Star::F3;
		else if (star_specType_str=="F4") star_specType=Star::F4;
		else if (star_specType_str=="F5") star_specType=Star::F5;
		else if (star_specType_str=="F6") star_specType=Star::F6;
		else if (star_specType_str=="F7") star_specType=Star::F7;
		else if (star_specType_str=="F8") star_specType=Star::F8;
		else if (star_specType_str=="F9") star_specType=Star::F9;
		else if (star_specType_str=="G0") star_specType=Star::G0;
		else if (star_specType_str=="G1") star_specType=Star::G1;
		else if (star_specType_str=="G2") star_specType=Star::G2;
		else if (star_specType_str=="G3") star_specType=Star::G3;
		else if (star_specType_str=="G4") star_specType=Star::G4;
		else if (star_specType_str=="G5") star_specType=Star::G5;
		else if (star_specType_str=="G6") star_specType=Star::G6;
		else if (star_specType_str=="G7") star_specType=Star::G7;
		else if (star_specType_str=="G8") star_specType=Star::G8;
		else if (star_specType_str=="G9") star_specType=Star::G9;
		else if (star_specType_str=="K0") star_specType=Star::K0;
		else if (star_specType_str=="K1") star_specType=Star::K1;
		else if (star_specType_str=="K2") star_specType=Star::K2;
		else if (star_specType_str=="K3") star_specType=Star::K3;
		else if (star_specType_str=="K4") star_specType=Star::K4;
		else if (star_specType_str=="K5") star_specType=Star::K5;
		else if (star_specType_str=="K6") star_specType=Star::K6;
		else if (star_specType_str=="K7") star_specType=Star::K7;
		else if (star_specType_str=="K8") star_specType=Star::K8;
		else if (star_specType_str=="K9") star_specType=Star::K9;
		else if (star_specType_str=="M0") star_specType=Star::M0;
		else if (star_specType_str=="M1") star_specType=Star::M1;
		else if (star_specType_str=="M2") star_specType=Star::M2;
		else if (star_specType_str=="M3") star_specType=Star::M3;
		else if (star_specType_str=="M4") star_specType=Star::M4;
		else if (star_specType_str=="M5") star_specType=Star::M5;
		else if (star_specType_str=="M6") star_specType=Star::M6;
		else if (star_specType_str=="M7") star_specType=Star::M7;
		else if (star_specType_str=="M8") star_specType=Star::M8;
		else if (star_specType_str=="M9") star_specType=Star::M9;
		else if (star_specType_str.substr(0,1)=="O") {
        	cout << "   Omitting star in with unsupported spectral type " << star_specType_str << endl;
        	cerr << "StarProducer: [W] WARNING: Star in " << m_starFilename << " with unsupported spectral type "
        		 << star_specType_str << " is being omitted. Valid spectral types are O7, O8, O9 and [B,A,F,G,K,M][0-9], e.g. K5" << endl;
        	continue;
		} else {
			throw runtime_error("Invalid spectral type in "+m_starFilename+": "+star_specType_str+" - must be O7, O8, O9 or [B,A,F,G,K,M][0-9], e.g. K5");
		}

		double effectiveTemperature = Star::getEffectiveTemperature(star_specType);

		double VbandMag,cheopsMag;
		if (iStar == 0 && m_overrideTargetMagnitude) {
			if (m_gaiaBand) {
				gaiaMag = m_targetMagnitude;
				VbandMag = gaiaMag - GaiaMinusVbandMagnitude(effectiveTemperature);
				cheopsMag = gaiaMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,effectiveTemperature,wavelengthDependence,iStar==0);
			} else{
				VbandMag = m_targetMagnitude;
				gaiaMag = VbandMag + GaiaMinusVbandMagnitude(effectiveTemperature);
				cheopsMag = VbandMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::Vband,effectiveTemperature,wavelengthDependence,iStar==0);
			}
		} else {
			VbandMag = gaiaMag - GaiaMinusVbandMagnitude(effectiveTemperature);
			cheopsMag = gaiaMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,effectiveTemperature,wavelengthDependence,iStar==0);
		}

		Star * star = new Star(SkyPosition(convertToArcseconds(star_ra,true),convertToArcseconds(star_dec,false)),VbandMag,cheopsMag,gaiaMag,star_specType);

		//Set the pointing direction to match the target star if requested
		if (iStar == 0 && m_useTargetStarForPointing) {
			fov->setPointingDirection(star->getSkyPosition());
			mpsVisits->setCellRaTarg(star->getSkyPosition().getRightAscension()/3600.);
			mpsVisits->setCellDecTarg(star->getSkyPosition().getDeclination()/3600.);
		}

		if (fov->isInside(star->getSkyPosition()) && gaiaMag<m_minMagnitude) {
			cout << "star " << iStar << (iStar==0 ? " (target star): ":": ") << endl;
			cout << "   Gaia magnitude " << gaiaMag << ", CHEOPS magnitude " << cheopsMag << ", V-band magnitude " << VbandMag << ", spectral type " << star_specType_str << ", RA=" << star_ra << ", dec=" << star_dec << endl;
			evaluateFlux(star,wavelengthDependence,iStar==0);
			star->initializeTimeSeries(nTimeSteps);
			fov->addStar(star);
//		} else {
//			if (gaiaMag>m_minMagnitude) {
//				cout << "Omitting star due to Gaia magnitude less than minMagnitude: " << gaiaMag << endl;
//			} else {
//				double dec_rad = fov->getPointingDirection().getDeclination()/3600.*M_PI/180.;
//				double delta_RA = (star->getSkyPosition().getRightAscension()-fov->getPointingDirection().getRightAscension())*cos(dec_rad);
//				double delta_dec = star->getSkyPosition().getDeclination()-fov->getPointingDirection().getDeclination();
//				cout << "Omitting star (outside field of view), delta_RA*cos(dec) = " << delta_RA << ", delta_dec=" << delta_dec << endl;
//			}
		}

	}
	stars_in.close();

}

void StarProducer::readBackgroundStars(SkyFieldOfView *fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps) const {

	string star_ra,star_dec;
	double VbandMag;
	string star_specType_str;

	//Open the input file
	string field_filename;
	double ra_field,dec_field;
	if (m_fieldCrowding == FIELD_CROWDING::low) {
		field_filename = "Low_field.txt";
		ra_field = convertToArcseconds("14:12:35.507",true); //target star position in Low_field.txt
		dec_field = convertToArcseconds("04:18:19.147",false); //target star position in Low_field.txt
	} else if (m_fieldCrowding == FIELD_CROWDING::medium) {
		field_filename = "Medium_field.txt";
		ra_field = convertToArcseconds("14:10:01.662",true); //target star position in Medium_field.txt
		dec_field = convertToArcseconds("07:35:40.666",false); //target star position in Medium_field.txt
	} else if (m_fieldCrowding == FIELD_CROWDING::extreme) {
		field_filename = "Extreme_field.txt";
		ra_field = convertToArcseconds("17:04:28.376",true); //target star position in Extreme_field.txt
		dec_field = convertToArcseconds("-28:50:00.303",false); //target star position in Extreme_field.txt
	}
	ifstream stars_in(string(getenv("CHEOPS_SW"))+"/resources/"+field_filename, ios::in);
    if (!stars_in.good()) throw runtime_error("Error in StarProducer: Error opening file "+field_filename);

    //Open the output file to contain the list of background stars with positions relative to the pointing direction
	//ofstream stars_out("background_stars.txt", ios::out);

    unsigned nStars=0;
    bool targetStar=true; //first star is a target star, to be skipped
	while (stars_in.good()) {

		stars_in >> star_ra;
		if(stars_in.eof()) break;
		stars_in >> star_dec;
		stars_in >> VbandMag; //Low_field.txt, Medium_field.txt and Extreme_field.txt were generated using V-band magnitudes
		stars_in >> star_specType_str;

		if (targetStar) { //skip the first star in the list (target star)
			targetStar=false;
			continue;
		}

		SkyPosition star_position(convertToArcseconds(star_ra,true)   - ra_field  + fov->getPointingDirection().getRightAscension(),
								  convertToArcseconds(star_dec,false) - dec_field + fov->getPointingDirection().getDeclination());

		//cout << star_ra << " " << convertToArcseconds(star_ra,true) << " " << ra_field << " " << star_position.getRightAscension() << endl;
		//cout << star_dec << " " << convertToArcseconds(star_dec,false) << " " << dec_field << " " << star_position.getDeclination() << endl << endl;

		Star::spectralType star_specType = Star::G0;
		if (star_specType_str=="M") star_specType = Star::M0;
		else if (star_specType_str=="K") star_specType = Star::K0;
		double effectiveTemperature = Star::getEffectiveTemperature(star_specType);

		double cheopsMag,gaiaMag;
		if (!m_gaiaBand) {
			cheopsMag = VbandMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::Vband,effectiveTemperature,wavelengthDependence);
			gaiaMag = VbandMag + GaiaMinusVbandMagnitude(effectiveTemperature);
		} else {
			gaiaMag = VbandMag + GaiaMinusVbandMagnitude(effectiveTemperature);
			cheopsMag = gaiaMag + m_fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,effectiveTemperature,wavelengthDependence);
		}

		if (fov->isInside(star_position) && gaiaMag<m_minMagnitude) {

			cout << "Background star " << ++nStars << ":" << endl;
			cout << "   Gaia magnitude " << gaiaMag << ", V-band magnitude " << VbandMag << ", spectral type " << star_specType_str << endl;
			cout << "   (RA,dec) (degrees): (" << star_position.getRightAscension()/3600. << "," << star_position.getDeclination()/3600. << ")" << endl;

			//stars_out << setw(7) << star_position.getRightAscension()/3600. << setw(9) << star_position.getDeclination()/3600. << setw(7) << VbandMag << " " << star_specType_str << endl;

			Star * star = new Star(SkyPosition(star_position),VbandMag,cheopsMag,gaiaMag,star_specType);
			evaluateFlux(star,wavelengthDependence);
			star->initializeTimeSeries(nTimeSteps);
			fov->addStar(star);

		}

	}

	stars_in.close();
	//stars_out.close();

	cout << endl << "Number of background stars generated: " << nStars << endl;

}

double StarProducer::GaiaMinusVbandMagnitude(double effectiveTemperature) {

	double GMinusV = 0.;
	//Read in the conversion between V-band magnitude and Gaia magnitude
	ifstream gminusv_in(string(getenv("CHEOPS_SW"))+"/resources/G-V.txt", ios::in);
	if (!gminusv_in.good()) throw runtime_error("Error in StarProducer::initialize: Error opening G-V file");
	string line;
	double GMinusV_prev,Teff_prev;
	while (getline(gminusv_in, line)) {
		if(line[0] != '#') {
			double GMinusV_line,Teff_GMinusV;
			stringstream(line) >> Teff_GMinusV >> GMinusV_line;
			if (Teff_GMinusV <= effectiveTemperature) {
				double mu = (effectiveTemperature - Teff_prev) / (Teff_GMinusV - Teff_prev);
				GMinusV = (GMinusV_prev*(1.-mu)+GMinusV_line*mu);
				break;
			} else {
				GMinusV_prev = GMinusV_line;
				Teff_prev = Teff_GMinusV;
			}
		}
	}
	gminusv_in.close();
	//cout << star->getSpectralTypeString() << " " << star->getEffectiveTemperature() << " " << GMinusV << endl;
	return GMinusV;

}

void StarProducer::evaluateFlux(Star * star, const WavelengthDependence * wavelengthDependence, bool target) const {

	//Get the energy spectrum for the star in 10 Angstrom steps
	double stellarEnergySpectrum[Star::kNWavelength*10];
	wavelengthDependence->getStellarEnergySpectrum_1nm(stellarEnergySpectrum,star->getEffectiveTemperature(),target);

	//Calculate fluxes for Vega and for the star in the reference frequency band (V-band or Gaia band)
	//Factor 10 in the calculations below is due to fact that wavelength bins are 10 Angstroms wide
	double flux_Vega_refBand = 0.;
	double photonFlux_Vega_refBand = 0.;
	double flux_star_refBand = 0.;
	double photonFlux_star_refBand = 0.;
	double photonEnergy[Star::kNWavelength*10];
	for (int i=0; i<Star::kNWavelength*10; i++) {
		//Calculate photon energy
		photonEnergy[i] = WavelengthDependence::kPlanck*WavelengthDependence::kSpeedOfLight/(m_wavelength[i]*1.e-10);  //Units: ergs
		//Get the flux for Vega in the reference band
		flux_Vega_refBand += 10.*m_flux_Vega[i]*m_transmission_refBand[i];
		//Get the photon flux for Vega in the reference band
		photonFlux_Vega_refBand += 10.*m_flux_Vega[i]*m_transmission_refBand[i]/photonEnergy[i];
		//Get the flux for the star in the reference band
		flux_star_refBand += 10.*stellarEnergySpectrum[i]*m_transmission_refBand[i];
		//Get the photon flux for the star in the reference band
		photonFlux_star_refBand += 10.*stellarEnergySpectrum[i]*m_transmission_refBand[i]/photonEnergy[i];
	}

	//cout << "flux_Vega_refBand: " << flux_Vega_refBand << ", photonFlux_Vega_refBand: " << photonFlux_Vega_refBand
	//     << ", flux_star_refBand: " << flux_star_refBand << ", photonFlux_star_refBand: " << photonFlux_star_refBand << endl;

	//Calculate star flux for the full spectrum normalized to that of a zero magnitude star
	double flux_star[Star::kNWavelength*10];
	double photonFlux_star[Star::kNWavelength*10];
	double integratedPhotonFlux_Star = 0.;
	for (int i=0; i<Star::kNWavelength*10; i++) {
		if (m_VegaSpectrum) {
			flux_star[i] = 10.*m_flux_Vega[i] ;
			photonFlux_star[i] = 10.*m_flux_Vega[i]/photonEnergy[i];
		} else {
			flux_star[i] = 10.*stellarEnergySpectrum[i] * flux_Vega_refBand/flux_star_refBand;
			photonFlux_star[i] = 10.*stellarEnergySpectrum[i]/photonEnergy[i] * photonFlux_Vega_refBand/photonFlux_star_refBand;
			integratedPhotonFlux_Star += photonFlux_star[i];
		}
	}

	//cout << "photonFlux_Vega_refBand: " << photonFlux_Vega_refBand*PSFGenerator::kTelescopeArea << ", photonFlux_star_refBand: " << photonFlux_star_refBand*PSFGenerator::kTelescopeArea
	//	 << ", photonFlux_star: " << integratedPhotonFlux_Star*PSFGenerator::kTelescopeArea << endl;

	//Use the magnitude of the star to convert the zero magnitude flux to the true flux
	double magnitude_Vega_refBand = m_gaiaBand ? 0.0 : 0.035;
	double deltaMagnitude = 0.;
	if (m_gaiaBand) {
		deltaMagnitude = star->getGaiaMagnitude() - magnitude_Vega_refBand;
	} else {
		deltaMagnitude = star->getVbandMagnitude() - magnitude_Vega_refBand;
	}
	double fluxRatio = pow(10,-0.4*deltaMagnitude);
	for (int i=0; i<Star::kNWavelength*10; i++) {
		flux_star[i] *= fluxRatio;
		photonFlux_star[i] *= fluxRatio;
	    //cout << "flux_star["<< i <<"]: " << flux_star[i] << ", photonFlux_star["<< i <<"]: " << photonFlux_star[i] << endl;
	}

//	ofstream spectrum_file("spectrum.dat", ios::app);
//	for (int i=0; i<Star::kNWavelength*10; i++) {
//		//stellarEnergySpectrum[i] = WavelengthDependence::planckRadiance(m_wavelength[i],star->getEffectiveTemperature());
//		//spectrum_file << stellarEnergySpectrum[i] << endl;
//		spectrum_file << photonFlux_star[i] << endl;
//	}

	//Rebin the flux arrays by factor 10
	double integratedFlux = 0.;
	double integratedPhotonFlux = 0.;
	double flux[Star::kNWavelength];
	double photonFlux[Star::kNWavelength];
	for (int i=0; i<Star::kNWavelength; i++) {
		flux[i]=0.;
		photonFlux[i]=0.;
		for (int j=0; j<10; j++) {
			flux[i]+=flux_star[10*i+j];
			photonFlux[i]+=photonFlux_star[10*i+j];
		}
		integratedFlux += flux[i];
		integratedPhotonFlux += photonFlux[i];
	    //cout << "flux["<< i <<"]: " << flux[i] << ", photonFlux["<< i <<"]: " << photonFlux[i] << endl;
	}

	//cout << "   integrated flux (" << Star::kWavelengthLowBound << " to " << Star::kWavelengthLowBound+10*Star::kNWavelength << " nm): "
	//		<< integratedFlux*PSFGenerator::kTelescopeArea << " ergs per second" << endl;
	cout << "   integrated flux on telescope aperture (" << Star::kWavelengthLowBound << " to " << Star::kWavelengthLowBound+10*Star::kNWavelength << " nm, Teff = " << star->getEffectiveTemperature() << "K): "
			<< setprecision(10) << integratedPhotonFlux*PSFGenerator::kTelescopeArea << " photons per second" << endl;

	//Set the Star flux arrays
	star->setMeanFlux(flux);
	star->setMeanPhotonFlux(photonFlux);

}

template <class T> void StarProducer::writeStars(Data * data, unsigned visitCounter) const {

	//Get the validity interval
	TimeConfiguration timeConf = data->getTimeConfiguration();
	UTC startTime_utc = timeConf.getVisitStartTimeUTC();
	UTC endTime_utc = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());

	//Open the output file
	Data::Visit visit = data->getVisit();
	string dirname = data->getOutputDirectory()+"/data";
	boost::filesystem::create_directory(dirname);
	char hstr[100];
	sprintf(hstr, "OBSID%010u",visit.m_obsId);
	string dataname = hstr;
	string filename = buildFileName(dirname, startTime_utc.getUtc(),
				VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visitCounter), PassId(),
				dataname, T::getExtName());
	if (boost::filesystem::exists(filename)) return;
	T * stars_out = createFitsFile<T>(dirname,startTime_utc.getUtc(),
			VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visitCounter), PassId(), dataname);

	vector<Star*> stars = data->getFieldOfView()->getStars();
	stars_out->setKeyDataname(dataname);
	stars_out->setKeyArchRev(0);
	stars_out->setKeyProcNum(visit.m_versionNum);
	stars_out->setKeyVStrtU(startTime_utc);
	stars_out->setKeyVStopU(endTime_utc);
	stars_out->setKeyPiName(visit.m_piName);
	stars_out->setKeyPiUid(visit.m_piUid);
	stars_out->setKeyObsCat(visit.m_obsCat);
	stars_out->setKeyProgtype(visit.m_progType);
	stars_out->setKeyProgId(visit.m_progId);
	stars_out->setKeyReqId(visit.m_reqId);
	stars_out->setKeyVisitctr(visitCounter);
	stars_out->setKeyObsid(visit.m_obsId);
	stars_out->setKeyPrpVst1(visit.m_prpFirst);
	stars_out->setKeyPrpVstn(visit.m_prpLast);

	boost::posix_time::ptime startTime = timeConf.getStartTime();
	boost::gregorian::greg_year year = startTime.date().year();
	boost::posix_time::ptime startOfYear = boost::posix_time::ptime(boost::gregorian::date(year,1,1));
	double secondsSinceYearStart = (startTime - startOfYear).total_seconds();
	double epoch = year + secondsSinceYearStart/(3600*24*365);

	stars_out->setKeyObsepoch(epoch);
	stars_out->setKeyCentRa((*stars.begin())->getSkyPosition().getRightAscension()/3600.);
	stars_out->setKeyCentDec((*stars.begin())->getSkyPosition().getDeclination()/3600.);
	stars_out->setKeyFov(data->getFieldOfView()->getFOVradius());
	stars_out->setKeyPitl(true);
	stars_out->setKeyTEff((*stars.begin())->getEffectiveTemperature());

	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+data->getGainFilename(),"READONLY");
	stars_out->setKeyGain(gainCorrection_file->getKeyGainNom());
	delete gainCorrection_file;

	for (vector<Star*>::const_iterator star = stars.begin(); star!=stars.end(); ++star) {
		if (star==stars.begin()) {
			stars_out->setKeyTargname(visit.m_targName);
			stars_out->setKeySpectype((*star)->getSpectralTypeString());
			stars_out->setKeyMagG((*star)->getGaiaMagnitude());
			stars_out->setKeyMagGerr((*star)->getGaiaMagnitudeError());
			stars_out->setKeyMagChps((*star)->getCheopsMagnitude());
			stars_out->setKeyMagCerr((*star)->getCheopsMagnitudeError());
			stars_out->setKeyRaTarg((*star)->getSkyPosition().getRightAscension()/3600.);
			stars_out->setKeyDecTarg((*star)->getSkyPosition().getDeclination()/3600.);
		} else { //Leave these as NULL for the target star for consistency with MPS provided EXT_DRFT_StarCatalogue files
			stars_out->setCellRaErr((*star)->getSkyPosition().getRightAscensionError()/3600.);
			stars_out->setCellDecErr((*star)->getSkyPosition().getDeclinationError()/3600.);
		}
		stars_out->setCellRa((*star)->getSkyPosition().getRightAscension()/3600.);
		stars_out->setCellDec((*star)->getSkyPosition().getDeclination()/3600.);
		stars_out->setCellDistance((*star)->getSkyPosition().separation((*(stars.begin()))->getSkyPosition()));
		stars_out->setCellTEff((*star)->getEffectiveTemperature());
		if ((*star)->getEffectiveTemperatureError() == 0.) {
			stars_out->setNullTEffErr();
		} else {
			stars_out->setCellTEffErr((*star)->getEffectiveTemperatureError());
		}
		stars_out->setCellMagGaia((*star)->getGaiaMagnitude());
		if ((*star)->getGaiaMagnitudeError()>0.) {
			stars_out->setCellMagGaiaErr((*star)->getGaiaMagnitudeError());
		} else {
			stars_out->setNullMagGaiaErr();
		}
		stars_out->setCellId((*star)->getCatalogueId());
		stars_out->setCellTarget(star==stars.begin());
		stars_out->WriteRow();
	}

	delete stars_out;

	cout << endl << "Number of stars within (radius " << data->getFieldOfView()->getFOVradius()
		 << " arcseconds) and with Gaia magnitude<" << m_minMagnitude << ": " << data->getFieldOfView()->getStars().size() << endl << endl;

}

double StarProducer::convertToArcseconds(string hoursMinsSecs, bool isRA) {
	if (hoursMinsSecs == "As target star") return 0.;
	if (hoursMinsSecs.find(":") == string::npos) {
		//already in decimal format: convert from degrees to arcseconds
		return boost::lexical_cast<double>(hoursMinsSecs)*3600.;
	} else {
		vector<string> substrings;
		boost::split(substrings,hoursMinsSecs,boost::is_any_of(":"));
		int hours = boost::lexical_cast<int>(substrings[0]);
		int minutes = boost::lexical_cast<int>(substrings[1]);
		double seconds = boost::lexical_cast<double>(substrings[2]);
		double return_value = double(3600*abs(hours)+60*minutes)+seconds;
		if (isRA) return_value*=15.;
		if (hoursMinsSecs.substr(0,1)=="-") return_value *= -1.;
		//cout << hoursMinsSecs << " " << hours << ":" << minutes << ":" << seconds << " " <<  return_value << endl;
		return return_value;
	}
}
