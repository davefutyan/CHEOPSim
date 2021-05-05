/*
 * ReferenceFileConverter.cxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#include <fstream>

#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian_types.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"

#include "REF_APP_DarkFrame.hxx"
#include "REF_APP_WhitePSF.hxx"
#include "REF_APP_ColouredPSF.hxx"
#include "REF_APP_OversampledWhitePSF.hxx"
#include "REF_APP_OversampledColouredPSF.hxx"
#include "REF_APP_WhitePSFMetadata.hxx"
#include "REF_APP_ColouredPSFMetadata.hxx"
#include "REF_APP_WhiteCCDLocationPSF.hxx"
#include "REF_APP_WhiteCCDLocationPSFMetadata.hxx"
#include "REF_APP_Jitter.hxx"
#include "AUX_RES_Orbit.hxx"
#include "REF_APP_Throughput.hxx"
#include "REF_APP_QE.hxx"
#include "REF_APP_Temperature.hxx"
#include "REF_APP_StrayLight.hxx"
#include "SCI_RAW_SubArray.hxx"
#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "CreateFitsFile.hxx"

#include "detector/include/FlatFieldGenerator.hxx"
#include "satellite/include/JitterProducer.hxx"
#include "telescope/include/PSFGenerator.hxx"

#include "ReferenceFileConverter.hxx"

ReferenceFileConverter::ReferenceFileConverter(const ParamsPtr params): m_params(params) {

	m_validityStart = params->GetAsString("validityStart");
	m_validityStop = params->GetAsString("validityStop");
	m_inputFilename = params->GetAsString("inputFilename");
	m_provider = params->GetAsString("provider");
	m_description = params->GetAsString("description");

}

void ReferenceFileConverter::process() const {

	bool singleTemp = true;
	double psfTemp = -999.;
	string type = m_params->GetAsString("referenceFileType");
	if (type.find("PSF") != std::string::npos) {
		singleTemp = m_params->GetAsBool("singleTemperaturePSF");
		psfTemp = m_params->GetAsDouble("psfTemperature");
	}

	if (type == "REF_APP_QE") {
		convertQE(m_params->GetAsDouble("ccdTemperature"));
	} else if (type == "REF_APP_Throughput") {
		convertThroughput(m_params->GetAsBool("endOfLife"));
	} else if (type == "REF_APP_WhitePSF") {
		convertWhitePSF(singleTemp,psfTemp);
	} else if (type == "REF_APP_WhiteCCDLocationPSF") {
		convertCcdLocationPSF();
	} else if (type == "REF_APP_DarkFrame") {
		convertDarkFrame();
	} else if (type == "REF_APP_FlatField") {
		convertFlatField();
	} else if (type == "AUX_RES_Orbit") {
		convertOrbit();
	} else if (type == "REF_APP_Jitter") {
		convertJitter(m_params->GetAsDouble("jitterOffset"));
	} else if (type == "REF_APP_StrayLight") {
		convertStrayLight(m_params->GetAsString("ltan"),
				 m_params->GetAsDouble("altitude"),
				 m_params->GetAsDouble("pointingRA"),
				 m_params->GetAsDouble("pointingDec"));
	} else if (type == "REF_APP_Temperature") {
		convertTemperature();
	} else if (type == "REF_APP_OversampledWhitePSF") {
		convertWhitePSF(singleTemp,psfTemp,true);
	} else if (type == "REF_APP_ColouredPSF") {
		convertColouredPSF(singleTemp,psfTemp);
	} else if (type == "REF_APP_OversampledColouredPSF") {
		convertColouredPSF(singleTemp,psfTemp, true);
	}

}

void ReferenceFileConverter::convertQE(double ccdTemp) const {

	bool slopeFromInput = true;

	//Open output fits file
	boost::filesystem::path file_path(m_inputFilename);
	RefAppQe * qe_out = new RefAppQe(buildFileName(".", m_validityStart,
			VisitId(), PassId(), std::string(), RefAppQe::getExtName()),"CREATE");
	qe_out->setKeyVStrtU(UTC(m_validityStart));
	qe_out->setKeyVStopU(UTC(m_validityStop));
	qe_out->setKeyProvider(m_provider);
	qe_out->setKeyDescrip(m_description);
	//qe_out->setKeyTempCcd(ccdTemp);

	ifstream qe_in(m_inputFilename, ios::in);
	ifstream qe_vs_temp_in(string(getenv("CHEOPS_SW"))+"/resources/qe_may2016.dat", ios::in); //provided in case slope is not provided in input file
    if (!qe_in.good()) throw runtime_error("Error opening quantum efficiency file");
    string line,dummy;
    int i=0;
    while (getline(qe_in, line)) {
    	if(line[0] != '#') {
    	    double qe, wavelength, err, slopeVsTemp;
    	    if (slopeFromInput) {
    	    	stringstream(line) >> wavelength >> qe >> err >> slopeVsTemp;
    	    } else {
    	    	err = 0.;
    	    	stringstream(line) >> wavelength >> qe;
        		getline(qe_vs_temp_in, line);
        		stringstream(line) >> dummy >> dummy >> slopeVsTemp;
        		getline(qe_vs_temp_in, line);
        		stringstream(line) >> dummy >> dummy >> dummy; //comment this line if the input file has the same cadence as qe_may2016.dat
    	    }
    		if (wavelength>2000.) wavelength/=10.; //convert angstroms to nm
    		//qe/=100.; //convert from per cent
    		cout << "wavelength[" << i << "]=" << wavelength << ", qe[" << i << "]=" << qe << ", err[" << i << "]=" << err << ", slopeVsTemp[" << i << "]=" << slopeVsTemp << endl;
    		i++;
    		qe_out->setCellWavelength(wavelength);

    		qe_out->setCellQe(qe);
    		if (err>0.) {
    			qe_out->setCellQeError(err);
    		} else {
    			qe_out->setNullQeError();
    		}
    		if (slopeFromInput && slopeVsTemp>0.) {
    			qe_out->setCellQeVsTempSlope(slopeVsTemp);
    		} else {
    			qe_out->setNullQeVsTempSlope();
    		}
    		qe_out->WriteRow();
    	}
	}
	qe_in.close();
	delete qe_out;

}

void ReferenceFileConverter::convertThroughput(bool endOfLife) const {

	//Open output fits file
	RefAppThroughput * throughput_out = new RefAppThroughput(buildFileName(".", m_validityStart,
			VisitId(), PassId(), std::string(), RefAppThroughput::getExtName()+(endOfLife ? "-EndOfLife" : "-BeginOfLife")),"CREATE");
	throughput_out->setKeyVStrtU(UTC(m_validityStart));
	throughput_out->setKeyVStopU(UTC(m_validityStop));
	throughput_out->setKeyProvider(m_provider);
	throughput_out->setKeyDescrip(m_description);

	//Read in the telescope throughput and fill the output fits file
	ifstream throughput_in(m_inputFilename, ios::in);
    if (!throughput_in.good()) throw runtime_error("Error opening throughput file");
    string line;
    int i = 0;
    while (getline(throughput_in, line)) {
    	if(line[0] != '#') {
    	    double throughput_bol, throughput_eol, wavelength;
    		stringstream(line) >> wavelength >> throughput_bol >> throughput_eol;
    		double throughput = endOfLife ? throughput_eol : throughput_bol;
    		if (wavelength>2000.) wavelength/=10.; //convert angstroms to nm
    		if (wavelength<2.) wavelength*=1000.; //convert microns to nm
			if (throughput>2.) throughput/=100.; //convert from per cent
    		cout << "throughput[" << i << "]=" << throughput << ", wavelength[" << i << "]=" << wavelength << endl;
    		i++;
    		throughput_out->setCellWavelength(wavelength);
    		throughput_out->setCellThroughput(throughput);
    		throughput_out->WriteRow();
    	}
	}

    throughput_in.close();
	delete throughput_out;

}

void ReferenceFileConverter::convertTemperature() const {

	//Open output fits file
	RefAppTemperature * temperature_out = new RefAppTemperature("temperature_telescope.fits", "CREATE");
	temperature_out->setKeyProvider(m_provider);
	temperature_out->setKeyDescrip(m_description);

	//Read in the telescope temperature and fill the output fits file
	string filename = m_inputFilename.substr(0,m_inputFilename.find_last_of("."));
	ifstream temperature_in(string(getenv("CHEOPS_SW"))+"/resources/"+filename+".txt", ios::in);
    if (!temperature_in.good()) throw runtime_error("Error opening temperature file");
    string line;
    int i = 0;
    while (getline(temperature_in, line)) {
    	if(line[0] != '#') {
    	    double temperature, time;
    		stringstream(line) >> time >> temperature;
    		//cout << "temperature[" << i << "]=" << temperature[i] << ", time[" << i << "]=" << time[i] << endl;
    		i++;
    		temperature_out->setCellTime(time);
    		temperature_out->setCellTemperature(temperature);
    		temperature_out->WriteRow();
    	}
	}
	temperature_in.close();
	delete temperature_out;

}

//The content of this method is commented because the data structure definition for the output reference file has since been updated,
//breaking compilation, and the method has not subsequently been needed and therefore has not been updated accordingly
void ReferenceFileConverter::convertDarkFrame() const {
//
//	// Output should be in e-/s, whereas input file is in ADU/s, so the values read from the file are all divided by the gain in the code below.
//	// Future measured dark frames will be provided in e-/s.
//	// Gain value in ADU/e- corresponding to the FM2 measurement. Confirmed by email from Adrien Deline 14 Nov 2017.
//	double gain = 0.4402;
//
//	//Open the input dark frame fits file
//	FitsDalImage<double> * dark_frame_in = new FitsDalImage<double>(string(getenv("CHEOPS_SW"))+"/resources/"+m_inputFilename, "READONLY",{1020,1020});
//
//	UTC currentTime = UTC(to_iso_extended_string(boost::posix_time::second_clock::universal_time()));
//	UTC endValidity = UTC(to_iso_extended_string(boost::posix_time::ptime(boost::gregorian::date(2040,1,1))));
//
//	//Open the output dark frame fits files including margins
//	RefAppDarkframe * dark_frame_out = new RefAppDarkframe(buildFileName(".", currentTime, VisitId(),
//			PassId(), std::string(), RefAppDarkframe::getExtName()), "CREATE",
//			{Image::kXDim,Image::kYDim});
//	RefAppDarkleftimage * dark_left_out =  Append_RefAppDarkleftimage(dark_frame_out,{Image::kNDarkCols,Image::kYDim});
//	RefAppDarkrightimage * dark_right_out =  Append_RefAppDarkrightimage(dark_frame_out,{Image::kNDarkCols,Image::kYDim});
//	RefAppDarktopimage * dark_top_out =  Append_RefAppDarktopimage(dark_frame_out,{Image::kXDim,Image::kNDarkRows});
//	RefAppDarkbottomimage * dark_bottoout =  Append_RefAppDarkbottomimage(dark_frame_out,{Image::kXDim,Image::kNDarkRows});
//	RefAppBlankleftimage * blank_left_out =  Append_RefAppBlankleftimage(dark_frame_out,{Image::kNBlankCols,Image::kYDim});
//	RefAppBlankrightimage * blank_right_out =  Append_RefAppBlankrightimage(dark_frame_out,{Image::kNBlankCols,Image::kYDim});
//
//	//Set header keywords
//	dark_frame_out->setKeyProvider(m_provider);
//	dark_frame_out->setKeyDescrip(m_description);
//	dark_frame_out->setKeyTemp(233);
//	dark_frame_out->setKeyVStrtU(currentTime);
//	dark_frame_out->setKeyVStopU(endValidity);
//	dark_left_out->setKeyVStrtU(currentTime);
//	dark_left_out->setKeyVStopU(endValidity);
//	dark_right_out->setKeyVStrtU(currentTime);
//	dark_right_out->setKeyVStopU(endValidity);
//	dark_top_out->setKeyVStrtU(currentTime);
//	dark_top_out->setKeyVStopU(endValidity);
//	dark_bottoout->setKeyVStrtU(currentTime);
//	dark_bottoout->setKeyVStopU(endValidity);
//	blank_left_out->setKeyVStrtU(currentTime);
//	blank_left_out->setKeyVStopU(endValidity);
//	blank_right_out->setKeyVStrtU(currentTime);
//	blank_right_out->setKeyVStopU(endValidity);
//
//	//Initialize random number generator to draw randomly from the 1020x1020 pixels in the input dark frame
//	int seed = 1;
//	boost::mt19937 randomNumberEngine(seed);
//    boost::uniforint<int> uniforrandom(0,1020);
//    RANDOUNIFORM * uniformRandomGenerator = new RANDOUNIFORM(randomNumberEngine,uniforrandom);
//
//    //Read the input dark frame and copy the values to the output FITS structure
//	for (int ix=0; ix<Image::kXDim; ix++) {
//		for (int iy=0; iy<Image::kYDim; iy++) {
//			if (ix>=2 && ix<Image::kXDim-2 && iy>=2 && iy<Image::kYDim-2) {
//				//FM1 dark frame covers only 1020x1020 so excludes dark reference columns and first and last two rows and columns
//				(*dark_frame_out)[iy][ix] = double((*dark_frame_in)[iy-2][ix-2])/gain;
//			} else {
//				//Assign the dark frame values for physical pixels for which measurements are not available by setting
//				//their values equal to that of a randomly selected pixel for which a measured value exists
//				(*dark_frame_out)[iy][ix] = double((*dark_frame_in)[(*uniformRandomGenerator)()][(*uniformRandomGenerator)()])/gain;
//			}
//		}
//	}
//
//	for (int ix=0; ix<Image::kNDarkCols; ix++) {
//		for (int iy=0; iy<Image::kYDim; iy++) {
//			(*dark_left_out)[iy][ix] = double((*dark_frame_in)[(*uniformRandomGenerator)()][(*uniformRandomGenerator)()])/gain;
//			(*dark_right_out)[iy][ix] = double((*dark_frame_in)[(*uniformRandomGenerator)()][(*uniformRandomGenerator)()])/gain;
//		}
//	}
//	for (int ix=0; ix<Image::kXDim; ix++) {
//		for (int iy=0; iy<Image::kNDarkRows; iy++) {
//			(*dark_top_out)[iy][ix] = double((*dark_frame_in)[(*uniformRandomGenerator)()][(*uniformRandomGenerator)()])/gain;
//			(*dark_bottoout)[iy][ix] = double((*dark_frame_in)[(*uniformRandomGenerator)()][(*uniformRandomGenerator)()])/gain;
//		}
//	}
//	for (int ix=0; ix<Image::kNBlankCols; ix++) {
//		for (int iy=0; iy<Image::kYDim; iy++) {
//			(*blank_left_out)[iy][ix] = 0.;
//			(*blank_right_out)[iy][ix] = 0.;
//		}
//	}
//
//	delete dark_frame_out;
//	delete dark_left_out;
//	delete dark_right_out;
//	delete dark_top_out;
//	delete dark_bottoout;
//	delete blank_left_out;
//	delete blank_right_out;
//	delete dark_frame_in;

}

//The content of this method is commented because the data structure definition for the output reference file has since been updated,
//breaking compilation, and the method has not subsequently been needed and therefore has not been updated accordingly
void ReferenceFileConverter::convertFlatField() const {

//	double flatFieldWavelength[FlatFieldGenerator::kNFlatFieldWavelength] =
//		{400.0,450.0,500.0,550.0,600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.0,1050.0,1100.0};
//	string flatFieldWavelengthString[FlatFieldGenerator::kNFlatFieldWavelength] =
//		{"400","450","500","550","600","650","700","750","800","850","900","950","1000","1050","1100"};
//	double exposureTime[FlatFieldGenerator::kNFlatFieldWavelength] =
//	{50.,10.,5.,2.5,2.5,2.5,2.5,5.,5.,5.,5.,10.,20.,120.,500.};
//
//	UTC currentTime = UTC(to_iso_extended_string(boost::posix_time::second_clock::universal_time()));
//	UTC endValidity = UTC(to_iso_extended_string(boost::posix_time::ptime(boost::gregorian::date(2040,1,1))));
//
//	//Open the output flat field fits files including margins
//	string filename = buildFileName(".", currentTime, VisitId(), PassId(), std::string(), RefAppFlatfield::getExtName());
//	RefAppFlatfield * flat_field_out = new RefAppFlatfield(filename, "CREATE",{Image::kXDim,Image::kYDim,15});
//	RefAppFlatfieldmetadata * flat_field_out_metadata = new RefAppFlatfieldmetadata(filename, "APPEND");
//
//	flat_field_out->setKeyProvider(m_provider);
//	flat_field_out->setKeyDescrip(m_description);
//	flat_field_out->setKeyBandwid(50.);
//	flat_field_out->setKeyTemp(233.);
//	flat_field_out->setKeyVStrtU(currentTime);
//	flat_field_out->setKeyVStopU(endValidity);
//	flat_field_out_metadata->setKeyVStrtU(currentTime);
//	flat_field_out_metadata->setKeyVStopU(endValidity);
//
//	string topDirectory = string(getenv("CHEOPS_SW"))+"/resources/"+m_inputFilename;
//
//	for (unsigned wavelengthBin=0; wavelengthBin<FlatFieldGenerator::kNFlatFieldWavelength; wavelengthBin++) {
//
//		for (int j = 0; j<Image::kXDim; j++) {
//			for (int i = 0; i<Image::kYDim; i++) {
//				(*flat_field_out)[wavelengthBin][j][i] = 0.;
//			}
//		}
//
//		double bias_offset = 0.;
//		int nBiasOffset = 0;
//		int nImages = 0;
//
//		//Fill the data
//		string inputDirectory = topDirectory+"/PRNU_"+flatFieldWavelengthString[wavelengthBin];
//		boost::filesystem::directory_iterator it(inputDirectory), eod;
//		BOOST_FOREACH(boost::filesystem::path const &fs_path, std::make_pair(it, eod)) {
//			if(is_regular_file(fs_path)) {
//
//				nImages++;
//
//				//See email Adrien Deline 3/7/2017: there are 3 dark rows top and bottom in the measurement
//				SciRawSubarray flat_field_in(fs_path.string(),"READONLY");
//				for (int j = 0; j<Image::kXDim; j++) {
//					for (int i = 0; i<Image::kYDim; i++) {
//						(*flat_field_out)[wavelengthBin][j][i] += static_cast<double>(flat_field_in[j+3][i+24]);
//					}
//				}
//				//Use the blank columns to determine the bias offset
//				for (int i=0; i<Image::kNBlankCols; i++) {
//					for (int j=0; j<Image::kYDim; j++) {
//						bias_offset += static_cast<double>(flat_field_in[j+3][i]);
//						nBiasOffset++;
//						bias_offset += static_cast<double>(flat_field_in[j+3][i+1064]);
//						nBiasOffset++;
//					}
//				}
//
//			}
//
//		}
//
//		bias_offset /= nBiasOffset;
//		cout << "Bias offset =  " << bias_offset << endl;
//
//		//Subtract the bias offset
//		for (int j = 0; j<Image::kXDim; j++) {
//			for (int i = 0; i<Image::kYDim; i++) {
//				(*flat_field_out)[wavelengthBin][j][i] -= bias_offset*nImages;
//			}
//		}
//
//		//Calculate the mean within the well illuminated region (central 1020x1020)
//		double mean = 0.;
//		int N_mean = 0;
//		for (int j = 2; j<Image::kXDim-2; j++) {
//			for (int i = 2; i<Image::kYDim-2; i++) {
//				mean += (*flat_field_out)[wavelengthBin][j][i];
//				N_mean++;
//			}
//		}
//		mean /= N_mean;
//		cout << "Mean in central 1020x1020 =  " << mean << endl;
//
//		//Normalize the pixel values to have mean 1 in the well illuminated region
//		for (int j = 0; j<Image::kXDim; j++) {
//			for (int i = 0; i<Image::kYDim; i++) {
//				(*flat_field_out)[wavelengthBin][j][i] /= mean;
//			}
//		}
//
//		flat_field_out_metadata->setCellStatus(1);
//		flat_field_out_metadata->setCellWavelength(flatFieldWavelength[wavelengthBin]);
//		flat_field_out_metadata->setCellExptime(exposureTime[wavelengthBin]);
//		flat_field_out_metadata->WriteRow();
//
//	}
//
//	delete flat_field_out;
//	delete flat_field_out_metadata;

}

void ReferenceFileConverter::convertJitter(double jitterOffset) const {

	//Open output fits file
	string filename = m_inputFilename.substr(0,m_inputFilename.find_last_of("."));
	RefAppJitter * jitter_out = new RefAppJitter(filename+".fits", "CREATE");
	jitter_out->setKeyProvider(m_provider);
	jitter_out->setKeyDescrip(m_description);

	//Read in jitter data and fill the output fits file
	ifstream jitter_in(string(getenv("CHEOPS_SW"))+"/resources/"+filename+".txt", ios::in);
	if (!jitter_in.good()) throw runtime_error("Error opening jitter file");
	string line;
	double jitterX_mean = 0.;
	double jitterY_mean = 0.;
	double jitterZ_mean = 0.;
	double time_mean = 0.;
	bool pitlFlag1_1s = true;
	bool pitlFlag2_1s = true;
	int N = 0;
	unsigned NRows = 0;
	while (getline(jitter_in, line)) {
	  if(line[0] != '#') {
	    double time,jitterX,jitterY,jitterZ;
	    if (m_inputFilename.find("1s") != string::npos || m_inputFilename.find("1Hz") != string::npos) {
	    	//averaging per 1s already done
	    	int pitlFlag1,pitlFlag2;
	    	stringstream(line) >> time >> pitlFlag1 >> pitlFlag2 >> jitterX >> jitterY >> jitterZ;
	    	if (time<jitterOffset) continue;
	    	time-=jitterOffset;
	    	if (time>48.*3600.) continue;
	    	jitter_out->setCellTime(time);
	    	jitter_out->setCellValidAocs(pitlFlag1);
	    	jitter_out->setCellValidScience(pitlFlag2!=255);
	    	jitter_out->setCellRoll(jitterX);
	    	jitter_out->setCellPitch(jitterY);
	    	jitter_out->setCellYaw(jitterZ);
	    	jitter_out->WriteRow();
	    	if (time < 600.) cout << time << " " << pitlFlag1_1s << " " << pitlFlag2_1s << " " << jitterX << " " << jitterY << " " << jitterZ << endl;
	    } else {
	    	//do averaging per 1s
	    	int itime_old = -1;
		    int pitlFlag1,pitlFlag2;
	    	stringstream(line) >> time >> pitlFlag1 >> pitlFlag2 >> jitterX >> jitterY >> jitterZ;
	    	if (time<jitterOffset) continue;
	    	time-=jitterOffset;
	    	if (time>48.*3600.) continue;
			int itime = static_cast<int>(time+1E-9);
	    	if (itime!=0 && itime > itime_old) {
	    		itime_old = itime;
	    		jitter_out->setCellTime(static_cast<int>(time_mean/N+1E-9));
		    	jitter_out->setCellValidAocs(pitlFlag1_1s);
		    	jitter_out->setCellValidScience(pitlFlag2_1s);
	    		jitter_out->setCellRoll(jitterX_mean/N);
	    		jitter_out->setCellPitch(jitterY_mean/N);
	    		jitter_out->setCellYaw(jitterZ_mean/N);
	    		jitter_out->WriteRow();
	    		if (itime < 100) cout << "*** " << time << " " << time_mean/N << " " << pitlFlag1_1s << " " << pitlFlag2_1s << " " << jitterX_mean/N << " " << jitterY_mean/N << " " << jitterZ_mean/N << endl;
	    		pitlFlag1_1s = true;
	    		pitlFlag2_1s = true;
	    		jitterX_mean = 0.;
	    		jitterY_mean = 0.;
	    		jitterZ_mean = 0.;
	    		time_mean = 0.;
	    		N = 0;
	    		NRows++;
	    		if (NRows==JitterProducer::kNJitter) break;
	    	}
	    	if (!pitlFlag1) pitlFlag1_1s = false;
	    	if (pitlFlag2==255) pitlFlag2_1s = false;
	    	jitterX_mean += jitterX;
	    	jitterY_mean += jitterY;
	    	jitterZ_mean += jitterZ;
	    	time_mean += time;
	    	N++;
	    	if (time < 100.) cout << N << " " << time << " " << pitlFlag1 << " " << pitlFlag2 << " " << jitterX << " " << jitterY << " " << jitterZ << endl;
	    }
	  }
	}
	jitter_in.close();
	delete jitter_out;

}

void ReferenceFileConverter::convertOrbit() const {

	//Set the leap second file
	list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0101.fits"};
	LeapSeconds::setLeapSecondsFileNames(leapSecondFile);

	//Open output fits file
	string filename = m_inputFilename.substr(0,m_inputFilename.find_last_of("."));
	AuxResOrbit * orbit_out = new AuxResOrbit(filename+".fits", "CREATE");

	//Read in the orbit data and fill the output fits file
	ifstream orbit_in(string(getenv("CHEOPS_SW"))+"/resources/"+filename+".txt", ios::in);
    if (!orbit_in.good()) throw runtime_error("Error opening orbit file");
    string line;
    UTC utc(2018,1,1,0,0,0);
    while (getline(orbit_in, line)) {
    	if(line[0] != '#') {
    	    unsigned time;
    	    double x,y,z,latitude,longitude;
    		stringstream(line) >> time >> x >> y >> z >> longitude >> latitude;
    		//stringstream(line) >> time >> x >> y >> z;
    		orbit_out->setCellEpoch(utc+DeltaTime(time*60));
    		orbit_out->setCellX(x);
    		orbit_out->setCellY(y);
    		orbit_out->setCellZ(z);
    		orbit_out->setCellLatitude(latitude);
    		orbit_out->setCellLongitude(longitude>180.?longitude-360.:longitude);
    		orbit_out->WriteRow();
    	}
	}
	orbit_in.close();
	delete orbit_out;

}

void ReferenceFileConverter::convertStrayLight(string ltan, double altitude, double pointingRA, double pointingDec) const {

	//Open output fits file
	string filename = m_inputFilename.substr(0,m_inputFilename.find_last_of("."));
	RefAppStraylight * sl_out = new RefAppStraylight(filename+".fits", "CREATE");
	sl_out->setKeyProvider(m_provider);
	sl_out->setKeyDescrip(m_description);
	sl_out->setKeyLtan(ltan);
	sl_out->setKeyAltitude(altitude);
	sl_out->setKeyPointra(pointingRA);
	sl_out->setKeyPointdec(pointingDec);

	ifstream sl_in(string(getenv("CHEOPS_SW"))+"/resources/"+filename+".dat", ios::in);
    if (!sl_in.good()) throw runtime_error("Error opening stray light file");
    string line;
    int i=0;
    while (getline(sl_in, line)) {
    	if(line[0] != '#') {
    	    double time, flux;
    		stringstream(line) >> time >> flux;
    		cout << "time[" << i << "]=" << time << ", flux[" << i << "]=" << flux << endl;
    		i++;
    		sl_out->setCellTime(time);
    		sl_out->setCellFlux(flux);
    		sl_out->WriteRow();
    	}
	}
	sl_in.close();
	delete sl_out;

}

void ReferenceFileConverter::convertColouredPSF(bool singleTemperaturePSF, double psfTemperature, bool oversampledPsf) const {

	int startWavelength = 375;
	vector<string> substrings;
	boost::split(substrings,m_inputFilename,boost::is_any_of(","));
	string field = substrings[0];
	string shape = substrings[1];

	// Define the temperatures corresponding to the available PSFs
	string thermalMap[PSFGenerator::kNpsfTemp] = {"Nominal_-10C","Nominal_CC","Nominal_HC1","Nominal_HC2"};
	if (shape != "Nominal") {
		thermalMap[0] = "Worst_-10C";
		thermalMap[1] = "Worst_CC";
		thermalMap[2] = "Worst_HC1";
		thermalMap[3] = "Worst_HC2";
	}
	string thermalMap_out[PSFGenerator::kNpsfTemp] = {"fixed","cold","hot1","hot2"};

	array4D * psf_50nm = new array4D(boost::extents[2000][2000][PSFGenerator::kNpsfWavelength][PSFGenerator::kNpsfTemp]);

	//Read in the PSF(s) from the ascii file(s), looping over wavelengths and temperatures
	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {

		int wavelength = startWavelength;
		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {

			string filename;
			if (shape == "Nominal") {
				filename = "/Volumes/My Passport/PSF_12.0px/Matlab/"+thermalMap[iTemp]+"/"+field+"/PSF_"+field+"_def0.000m"+to_string(wavelength)+"nm.txt";
			} else {
				filename = "/Volumes/My Passport/PSF_12.0px/Matlab/"+thermalMap[iTemp]+"/"+field+"/"+shape+"/PSF_"+field+"_"+shape+"_"+to_string(wavelength)+"nm.txt";
			}
			cout << "opening PSF file: " << filename << endl;
			ifstream psf_stream(filename, ios::in);
			if (!psf_stream.good()) throw runtime_error("Error opening PSF file");

			string line;
			int j=0;
			while (getline(psf_stream, line)) {
				stringstream ss;
				ss << line;
				int jj = j-2048+1000;
				if (jj>=0) {
					for (int i=0; i<4096; i++) {
						double val;
						ss >> val;
						int ii = i-2048+1000;
						if (ii>=0) (*psf_50nm)[ii][jj][iWavelength][iTemp] = val;
						if(ii==1999) break;
					}
				}
				if (jj==1999) break;
				j++;
			}

			psf_stream.close();
			wavelength += 50;
		}

	}

	//		cout << "writing 5nm PSF Fits cube..." << endl;
	//		RefAppOversampledcolouredpsf * psf_out = new RefAppOversampledcolouredpsf(m_inputFilename+"_5nm.fits", "APPEND",{2000,2000,141});
	//		RefAppColouredpsfmetadata * psf_out_metadata = new RefAppColouredpsfmetadata(m_inputFilename+"_5nm.fits", "APPEND");
	//		psf_out->setKeyProvider(m_provider);
	//		psf_out->setKeyDescrip(m_description);
	//		wavelength = 400;
	//		for (int iWavelength=0; iWavelength<141; iWavelength++) {
	//			cout << wavelength << endl;
	//			for (int j = 0; j<2000; j++) {
	//				for (int i = 0; i<2000; i++) {
	//					(*psf_out)[iWavelength][j][i] = (*psf_5nm)[i][j][iWavelength];
	//				}
	//			}
	//			psf_out_metadata->setCellStatus(1);
	//			psf_out_metadata->setCellIndexWavelength(iWavelength);
	//			psf_out_metadata->setCellWavelength(wavelength);
	//			psf_out_metadata->WriteRow();
	//			wavelength +=5;
	//		}
	//		delete psf_out;
	//		delete psf_out_metadata;

//	array4D * psf_50nm = new array4D(boost::extents[2000][2000][PSFGenerator::kNpsfWavelength][PSFGenerator::kNpsfTemp]);
//
//	cout << "initializing 50nm array..." << endl;
//	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {
//		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {
//			cout << iTemp << " " << iWavelength << endl;
//			for (int j = 0; j<2000; j++) {
//				for (int i = 0; i<2000; i++) {
//					(*psf_50nm)[i][j][iWavelength][iTemp] = 0.;
//				}
//			}
//		}
//	}
//
//	cout << "Filling 50nm array..." << endl;
//	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {
//		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {
//			cout << iTemp << " " << iWavelength << endl;
//			for (int j = 0; j<2000; j++) {
//				for (int i = 0; i<2000; i++) {
//					for (int ii = 0; ii<10; ii++) {
//						(*psf_50nm)[i][j][iWavelength][iTemp] += (*psf_5nm)[i][j][iWavelength*10+ii][iTemp]/10.;
//					}
//				}
//			}
//		}
//	}
//
//	delete psf_5nm;

	//Ensure PSFs are normalized to 1
	cout << "Checking PSF normalization..." << endl;
	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {
		int wavelength = startWavelength;
		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {
			double psf_sum = 0.;
			for (int j = 0; j<2000; j++) {
				for (int i = 0; i<2000; i++) {
					psf_sum += (*psf_50nm)[i][j][iWavelength][iTemp];
				}
			}
			if (fabs(psf_sum-1.)>1.0E-6) {
				cout << " Normalizing PSF for "+thermalMap[iTemp]+", wavelength " << to_string(wavelength) << " (old normalization " << psf_sum << ")" << endl;
				for (int j = 0; j<2000; j++) {
					for (int i = 0; i<2000; i++) {
						(*psf_50nm)[i][j][iWavelength][iTemp] /= psf_sum;
					}
				}
			}
			wavelength += 50;
		}

	}

	cout << "writing 50nm PSF Fits cube..." << endl;
	string filename;
	if (shape == "Nominal") {
		filename = "PSF_"+field+"_def0.000mW_350_1100nDwave5nradius12px_nominal_monochromatic";
	} else {
		filename = "PSF_"+field+"_"+shape+"_W_350_1100nDwave5nradius12px_worst_monochromatic";
	}
	RefAppOversampledcolouredpsf * psf_out = new RefAppOversampledcolouredpsf(filename+".fits", "APPEND",{2000,2000,PSFGenerator::kNpsfTemp,PSFGenerator::kNpsfWavelength});
	RefAppColouredpsfmetadata * psf_out_metadata = new RefAppColouredpsfmetadata(filename+".fits", "APPEND");
	psf_out->setKeyProvider(m_provider);
	psf_out->setKeyDescrip(m_description);
	for (int iTemp=0; iTemp<PSFGenerator::kNpsfTemp; iTemp++) {
		double wavelength = startWavelength;
		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {
			cout << wavelength << endl;
			for (int j = 0; j<2000; j++) {
				for (int i = 0; i<2000; i++) {
					(*psf_out)[iTemp][iWavelength][j][i] = (*psf_50nm)[i][j][iWavelength][singleTemperaturePSF ? 0 : iTemp];
				}
			}
			psf_out_metadata->setCellStatus(1);
			psf_out_metadata->setCellWavelength(wavelength);
			psf_out_metadata->setCellIndexWavelength(iWavelength);
			psf_out_metadata->setCellTempTel(psfTemperature);
			psf_out_metadata->setCellThermalMap(thermalMap_out[iTemp]);
			psf_out_metadata->setCellIndexTemp(iTemp);
			psf_out_metadata->WriteRow();
			wavelength += 50.;
		}
	}
	delete psf_out;
	delete psf_out_metadata;

	cout << "writing rebinned PSF Fits cube..." << endl;
	RefAppColouredpsf * psf_out_rebin = new RefAppColouredpsf(filename+"_rebin.fits", "APPEND",{200,200,PSFGenerator::kNpsfTemp,PSFGenerator::kNpsfWavelength});
	RefAppColouredpsfmetadata * psf_out_rebin_metadata = new RefAppColouredpsfmetadata(filename+"_rebin.fits", "APPEND");
	psf_out_rebin->setKeyProvider(m_provider);
	psf_out_rebin->setKeyDescrip(m_description);
	for (int iTemp=0; iTemp<PSFGenerator::kNpsfTemp; iTemp++) {
		double wavelength = startWavelength;
		for (int iWavelength=0; iWavelength<PSFGenerator::kNpsfWavelength; iWavelength++) {
			cout << wavelength << " " << iTemp << endl;
			for (int j = 0; j<200; j++) {
				for (int i = 0; i<200; i++) {
					for (int ii=0; ii<10; ii++) {
						for (int jj=0; jj<10; jj++) {
							(*psf_out_rebin)[iTemp][iWavelength][j][i] += (*psf_50nm)[i*10+ii][j*10+jj][iWavelength][singleTemperaturePSF ? 0 : iTemp];
						}
					}
				}
			}
			psf_out_rebin_metadata->setCellStatus(1);
			psf_out_rebin_metadata->setCellWavelength(wavelength);
			psf_out_rebin_metadata->setCellIndexWavelength(iWavelength);
			psf_out_rebin_metadata->setCellTempTel(psfTemperature);
			psf_out_rebin_metadata->setCellThermalMap(thermalMap_out[iTemp]);
			psf_out_rebin_metadata->setCellIndexTemp(iTemp);
			psf_out_rebin_metadata->WriteRow();
			wavelength += 50.;
		}
	}
	delete psf_50nm;
	delete psf_out_rebin;
	delete psf_out_rebin_metadata;

}

void ReferenceFileConverter::convertWhitePSF(bool singleTemperaturePSF, double psfTemperature, bool oversampledPsf) const {

	const int kDim = oversampledPsf ? 2000 : 200;

	// Extract the PSF filenames
	string prefix = m_inputFilename.substr(0,m_inputFilename.find_last_of("."));
	string suffix = m_inputFilename.substr(m_inputFilename.find_last_of("."));

	// Define the temperatures corresponding to the available PSFs
	string thermalMap[PSFGenerator::kNpsfTemp] = {"T-10C","ColdMap","HotMap1","HotMap2"};
	string thermalMap_out[PSFGenerator::kNpsfTemp] = {"fixed","cold","hot1","hot2"};

	array3D * psf = new array3D(boost::extents[kDim][kDim][PSFGenerator::kNpsfTemp]);

	//Read in the PSF(s) from the ascii file(s), looping over temperatures

	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {

		string filename;
		if (singleTemperaturePSF) {
			filename = prefix+suffix;
		} else {
			filename = prefix+"_"+thermalMap[iTemp]+suffix;
		}
		cout << "opening PSF file: " << filename << endl;
		if (suffix == ".fits") {
	    	SimRawDoubleprecisionsubarray psf_in(filename,"READONLY");
			std::vector<long> dimension = psf_in.GetSize();
			if (dimension.size() !=2) throw runtime_error("Error in ReferenceFileConverter::convertWhitePSF: PSF FITS file must contain a single image");
			if (dimension[0] != 200 || dimension[1] != 200) {
				throw runtime_error("Error in ReferenceFileConverter::convertWhitePSF: dimensions of PSF image in uploaded flat field FITS file must be 200x200. Image dimensions in uploaded file are "+
									to_string(dimension[0])+"x"+to_string(dimension[1]));
			}
			for (int j = 0; j<200; j++) {
				for (int i = 0; i<200; i++) {
					(*psf)[i][j][iTemp] = static_cast<double>(psf_in[j][i]);
				}
			}
		} else {
			ifstream psf_stream(string(getenv("CHEOPS_SW"))+"/resources/"+filename, ios::in);
			if (!psf_stream.good()) throw runtime_error("Error opening PSF file");

			string line;
			int j=0;
			while (getline(psf_stream, line)) {
				stringstream ss;
				ss << line;
				for (int i=0; i<kDim; i++) {
					ss >> (*psf)[i][j][iTemp];
				}
				j++;
			}

			psf_stream.close();
		}

	}

//	//Truncate the PSF at radius 70 pixels
//	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {
//		for (int j = 0; j<kDim; j++) {
//			for (int i = 0; i<kDim; i++) {
//				int x = i-kDim/2;
//				int y = j-kDim/2;
//				if (sqrt(double(x*x+y*y)) > oversampeldPsf ? 700. : 70.) {
//					(*psf)[i][j][iTemp] = 0.;
//				}
//			}
//		}
//	}

	//Ensure PSFs are normalized to 1
	cout << "Checking PSF normalization..." << endl;
	for (int iTemp=0; iTemp<(singleTemperaturePSF ? 1 : PSFGenerator::kNpsfTemp); iTemp++) {
		double psf_sum = 0.;
		for (int j = 0; j<kDim; j++) {
			for (int i = 0; i<kDim; i++) {
				psf_sum += (*psf)[i][j][iTemp];
			}
		}
		if (fabs(psf_sum-1.)>1.0E-6) {
			cout << "Normalizing PSF for "+thermalMap[iTemp] << " (old normalization " << psf_sum << ")" << endl;
			for (int j = 0; j<kDim; j++) {
				for (int i = 0; i<kDim; i++) {
					(*psf)[i][j][iTemp] /= psf_sum;
				}
			}
		}
	}

	if (oversampledPsf) {
		cout << "writing oversampled PSF Fits cube..." << endl;
		RefAppOversampledwhitepsf * psf_out = new RefAppOversampledwhitepsf(prefix+".fits", "APPEND",{2000,2000,PSFGenerator::kNpsfTemp});
		RefAppWhitepsfmetadata * psf_out_metadata = new RefAppWhitepsfmetadata(prefix+".fits", "APPEND");
		psf_out->setKeyProvider(m_provider);
		psf_out->setKeyDescrip(m_description);
		psf_out->setKeyVStrtU(UTC(m_validityStart));
		psf_out->setKeyVStopU(UTC(m_validityStop));
		psf_out_metadata->setKeyVStrtU(UTC(m_validityStart));
		psf_out_metadata->setKeyVStopU(UTC(m_validityStop));
		for (int iTemp=0; iTemp<PSFGenerator::kNpsfTemp; iTemp++) {
			for (int j = 0; j<2000; j++) {
				for (int i = 0; i<2000; i++) {
					(*psf_out)[iTemp][j][i] = (*psf)[i][j][singleTemperaturePSF ? 0 : iTemp];
				}
			}
			psf_out_metadata->setCellStatus(1);
			psf_out_metadata->setCellTempTel(psfTemperature);
			psf_out_metadata->setCellThermalMap(thermalMap_out[iTemp]);
			psf_out_metadata->setCellIndexTemp(iTemp);
			psf_out_metadata->WriteRow();
		}
		delete psf_out;
		delete psf_out_metadata;
	}

	cout << "writing non-oversampled PSF Fits cube..." << endl;
	string outputFilename = buildFileName(".", m_validityStart,VisitId(), PassId(), std::string(), RefAppWhitepsf::getExtName());
	if (oversampledPsf) outputFilename = prefix+"_rebin.fits";
	RefAppWhitepsf * psf_out_rebin = new RefAppWhitepsf(outputFilename, "APPEND",{200,200,PSFGenerator::kNpsfTemp});
	RefAppWhitepsfmetadata * psf_out_rebin_metadata = new RefAppWhitepsfmetadata(outputFilename, "APPEND");
	psf_out_rebin->setKeyProvider(m_provider);
	psf_out_rebin->setKeyDescrip(m_description);
	psf_out_rebin->setKeyVStrtU(UTC(m_validityStart));
	psf_out_rebin->setKeyVStopU(UTC(m_validityStop));
	psf_out_rebin_metadata->setKeyVStrtU(UTC(m_validityStart));
	psf_out_rebin_metadata->setKeyVStopU(UTC(m_validityStop));
	for (int iTemp=0; iTemp<PSFGenerator::kNpsfTemp; iTemp++) {
		for (int j = 0; j<200; j++) {
			for (int i = 0; i<200; i++) {
				if (oversampledPsf) {
					for (int ii=0; ii<10; ii++) {
						for (int jj=0; jj<10; jj++) {
							(*psf_out_rebin)[iTemp][j][i] += (*psf)[i*10+ii][j*10+jj][singleTemperaturePSF ? 0 : iTemp];
						}
					}
				} else {
					(*psf_out_rebin)[iTemp][j][i] = (*psf)[i][j][singleTemperaturePSF ? 0 : iTemp];
				}
			}
		}
		psf_out_rebin_metadata->setCellStatus(1);
		psf_out_rebin_metadata->setCellTempTel(psfTemperature);
		psf_out_rebin_metadata->setCellThermalMap(thermalMap_out[iTemp]);
		psf_out_rebin_metadata->setCellIndexTemp(iTemp);
		psf_out_rebin_metadata->WriteRow();
	}
	delete psf_out_rebin;
	delete psf_out_rebin_metadata;

	delete psf;

}


void ReferenceFileConverter::convertCcdLocationPSF() const {

	const int kNLayers = 10;
	const int kDim = 200;
	int xcentres[kNLayers] = {537,195,195,195,537,537,877,877,877,263};
	int ycentres[kNLayers] = {537,195,537,877,195,877,195,537,877,842};

	if (!boost::filesystem::is_directory(m_inputFilename)) {
		throw runtime_error("Error in ReferenceFileConverter::convertCcdLocationPSF: directory "+m_inputFilename+" does not exist");
	}

	//Obtain a sorted list of fits image files in the directory
	vector<string> psfFiles;
	boost::filesystem::directory_iterator it(m_inputFilename), eod;
	BOOST_FOREACH(boost::filesystem::path const &fs_path, std::make_pair(it, eod)) {
		if(is_regular_file(fs_path)) {
			string filename = fs_path.string();
			if(filename.substr(filename.find_last_of(".")+1) == "fits") {
				psfFiles.push_back(filename);
			}
		}
	}
	std::sort(psfFiles.begin(), psfFiles.end());

	array3D * psf = new array3D(boost::extents[kDim][kDim][kNLayers]);

	for (int iPos=0; iPos<kNLayers; iPos++) {

		cout << "opening PSF file: " << psfFiles[iPos] << endl;
		SimRawDoubleprecisionsubarray psf_in(psfFiles[iPos],"READONLY");
		std::vector<long> dimension = psf_in.GetSize();
		if (dimension.size() !=2) throw runtime_error("Error in ReferenceFileConverter::convertCcdLocationPSF: PSF FITS file must contain a single image");
		if (dimension[0] != 200 || dimension[1] != 200) {
			throw runtime_error("Error in ReferenceFileConverter::convertCcdLocationPSF: dimensions of PSF image in uploaded flat field FITS file must be 200x200. Image dimensions in uploaded file are "+
					to_string(dimension[0])+"x"+to_string(dimension[1]));
		}
		for (int j = 0; j<200; j++) {
			for (int i = 0; i<200; i++) {
				(*psf)[i][j][iPos] = static_cast<double>(psf_in[j][i]);
			}
		}

	}

	//Ensure PSFs are normalized to 1
	cout << "Checking PSF normalization..." << endl;
	for (int iPos=0; iPos<kNLayers; iPos++) {
		double psf_sum = 0.;
		for (int j = 0; j<kDim; j++) {
			for (int i = 0; i<kDim; i++) {
				psf_sum += (*psf)[i][j][iPos];
			}
		}
		if (fabs(psf_sum-1.)>1.0E-6) {
			cout << "Normalizing PSF for position " << iPos << " (old normalization " << psf_sum << ")" << endl;
			for (int j = 0; j<kDim; j++) {
				for (int i = 0; i<kDim; i++) {
					(*psf)[i][j][iPos] /= psf_sum;
				}
			}
		}
	}

	cout << "writing PSF Fits cube..." << endl;
	string outputFilename = buildFileName(".", m_validityStart,
			VisitId(), PassId(), std::string(), RefAppWhiteccdlocationpsf::getExtName());
	RefAppWhiteccdlocationpsf * psf_out = new RefAppWhiteccdlocationpsf(outputFilename, "APPEND",{200,200,kNLayers});
	RefAppWhiteccdlocationpsfmetadata * psf_out_metadata = new RefAppWhiteccdlocationpsfmetadata(outputFilename, "APPEND");
	psf_out->setKeyProvider(m_provider);
	psf_out->setKeyDescrip(m_description);
	psf_out->setKeyTemp(-10.);
	psf_out->setKeyVStrtU(UTC(m_validityStart));
	psf_out->setKeyVStopU(UTC(m_validityStop));
	psf_out_metadata->setKeyVStrtU(UTC(m_validityStart));
	psf_out_metadata->setKeyVStopU(UTC(m_validityStop));
	for (int iPos=0; iPos<kNLayers; iPos++) {
		for (int j = 0; j<200; j++) {
			for (int i = 0; i<200; i++) {
				(*psf_out)[iPos][j][i] = (*psf)[i][j][iPos];
			}
		}
		psf_out_metadata->setCellXOffFullArray(xcentres[iPos]-100);
		psf_out_metadata->setCellYOffFullArray(ycentres[iPos]-100);
		psf_out_metadata->WriteRow();
	}
	delete psf_out;
	delete psf_out_metadata;

	delete psf;

}
