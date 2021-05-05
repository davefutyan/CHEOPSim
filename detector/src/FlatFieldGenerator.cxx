/*
 * FlatFieldGenerator.cxx
 *
 *  Created on: 4 Mar 2014
 *      Author: futyand
 */

#include <fstream>

#include "boost/filesystem.hpp"

#include "REF_APP_FlatFieldTeff.hxx"
#include "REF_APP_FlatFieldTeffMetadata.hxx"
#include "REF_APP_WhiteFlatField.hxx"

#include "source/include/StarProducer.hxx"
#include "detector/include/FlatFieldGenerator.hxx"

FlatFieldGenerator::FlatFieldGenerator() : Module("FlatFieldGenerator",timeLoop) {

	for (int ix=0; ix<Image::kXDim; ix++) {
		for (int iy=0; iy<Image::kYDim; iy++) {
			m_deadPixelQE[ix][iy] = 1.;
			m_flatField[ix][iy] = 0.;
		}
	}

	//Default to effective temperature for G5 star if there are no stars in the FOV.
	m_flatFieldTeff = 5660.;

}

FlatFieldGenerator::~FlatFieldGenerator() {
	delete m_gaussianNoiseGenerator;
	delete m_gaussianNoiseGenerator2;
}

void FlatFieldGenerator::initialize(const ModuleParams& params) {

	m_gaussianFlatField = params.GetAsBool("gaussianFlatField");
	m_flatFieldSigma = params.GetAsDouble("flatFieldSigma");
	unsigned flatFieldSeed = params.GetAsInt("flatFieldSeed");

	m_applyFlatField = params.GetAsBool("applyFlatField");
	m_flatFieldFilename = params.GetAsString("flatFieldFilename");
	m_flatFieldScaleFactor = params.GetAsDouble("flatFieldScaleFactor");

	m_flatFieldTeffOffset = params.GetAsDouble("flatFieldTeffOffset");
	m_flatFieldSmearSigma = params.GetAsDouble("flatFieldSmearSigma");

	m_flatFieldToSubtractFilename = params.GetAsString("flatFieldToSubtractFilename");

	m_writeTruthFlatField = params.GetAsBool("writeTruthFlatField");

	m_doDeadPixels = params.GetAsBool("doDeadPixels");
	m_deadPixelPositionSeed = params.GetAsInt("deadPixelPositionSeed");
	double fracDeadPixels = params.GetAsDouble("fracDeadPixels");
	m_deadPixelRelativeQE = params.GetAsDouble("deadPixelRelativeQE");

	//Convert the fraction to the number of dead pixels corresponding to full frame
	int nPixels = Image::kXDim * Image::kYDim;
	m_nDeadPixels = static_cast<int>(lround(fracDeadPixels * nPixels));

	boost::mt19937 randomNumberEngine(flatFieldSeed);
	boost::mt19937 randomNumberEngine2(flatFieldSeed+1);
    boost::normal_distribution<double> normalDistribution(1.,1.);
    boost::normal_distribution<double> normalDistribution2(1.,m_flatFieldSmearSigma);
    m_gaussianNoiseGenerator = new RANDOM_GAUSS(randomNumberEngine,normalDistribution);
    m_gaussianNoiseGenerator2 = new RANDOM_GAUSS(randomNumberEngine2,normalDistribution2);

}

void FlatFieldGenerator::doBegin(Data* data, bool fullFrame) {

	data->initializeBadPixelMap(m_flatFieldFilename);

	if (m_doDeadPixels) setDeadPixels(data);

	if (data->getFieldOfView()->getStars().size()>0) {
		m_flatFieldTeff = (data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
	}

	if (m_applyFlatField && data->getFlatField() == nullptr) {//May have already been generated if doFullFrame is on
		if (m_gaussianFlatField) {
			setGaussianFlatField(data);
		} else {
			readFlatField(data);
		}
		setEmpiricalDeadPixels(data);
		if (m_writeTruthFlatField) writeTruthFlatField(data);
		//writeWhiteFlatField(data); //Generate single image in runCHEOPSim.xml, running only FocalPlaneGenerator,FlatFieldGenerator; set m_flatFieldTeff to desired value in constructor (negative for flat incident spectrum)
	}

}

void FlatFieldGenerator::setGaussianFlatField(Data* data) {

	//define the width of the Gaussian distribution to draw from
	m_gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(1.,m_flatFieldSigma));

	boost::math::tools::stats<double> flatFieldStats;

    //Fill the flat field array
    for (int i=0; i<Image::kXDim; i++) {
    	for (int j=0; j<Image::kYDim; j++) {
    		m_flatField[i][j] = (*m_gaussianNoiseGenerator)();
    		flatFieldStats.add(m_flatField[i][j]);
    	}
    }

    cout << "Gaussian flat field: mean=" << sqrt(flatFieldStats.mean())
         << ", standard deviation=" << sqrt(flatFieldStats.variance1()) << endl;

	data->setFlatField(generateFlatField(flatFieldStats));
	data->setFlatFieldToSubtract(generateFlatField(flatFieldStats));

}

void FlatFieldGenerator::readFlatField(Data* data) {

	RefAppFlatfieldteff * flatField;
	RefAppFlatfieldteffmetadata * flatFieldMetadata = nullptr;

	boost::math::tools::stats<double> dummyStats;
	boost::math::tools::stats<double> flatFieldStats;

    //If input flat field file is user uploaded, convert it to RefAppFlatfield format
    if (m_flatFieldFilename == "flatField.fits") {

    	flatField = new RefAppFlatfieldteff("REF_APP_FlatFieldTeff.fits", "CREATE",{Image::kXDim,Image::kYDim,1});
    	flatFieldMetadata = new RefAppFlatfieldteffmetadata("REF_APP_FlatFieldTeff.fits", "APPEND");

    	UTC endValidity = UTC(to_iso_extended_string(boost::posix_time::ptime(boost::gregorian::date(2050,1,1))));
    	string email = "";
    	ifstream email_in("email.txt", ios::in);
    	if (email_in.good()) {
    	    string line;
            getline(email_in, line);
            stringstream(line) >> email;
    	}
    	email_in.close();
    	flatField->setKeyProvider(email);
    	flatField->setKeyDescrip("User uploaded flat field");
    	flatField->setKeyVStopU(endValidity);
    	flatField->setKeyVStopU(endValidity);

    	//Copy the data from the user provided flat field file to the RefAppFlatfield object
    	//The user provided file is interpreted as SimRawDoubleprecisionsubarray, but requires only that the user file
    	//contains an image in the primary extension with NAXIS=2 and NAXIS1=NAXIS2=1024
    	SimRawDoubleprecisionsubarray flatField_in(m_flatFieldFilename,"READONLY");
		std::vector<long> dimension = flatField_in.GetSize();
		if (dimension.size() !=2) throw runtime_error("Error in FlatFieldGenerator::readFlatField: uploaded flat field FITS file must contain a single image");
		if (dimension[0] != Image::kXDim || dimension[1] != Image::kYDim) {
			throw runtime_error("Error in FlatFieldGenerator::readFlatField: dimensions of image in uploaded flat field FITS file must be 1024x1024. Image dimensions in uploaded file are "+
								to_string(dimension[0])+"x"+to_string(dimension[1]));
		}
		for (int j = 0; j<1024; j++) {
			for (int i = 0; i<1024; i++) {
				(*flatField)[0][j][i] = static_cast<double>(flatField_in[j][i]);
			}
		}
		flatFieldMetadata->setCellStatus(1);
		flatFieldMetadata->setCellTEff(m_flatFieldTeff);
		flatFieldMetadata->WriteRow();

	    //Fill the flat field array
	    for (int i=0; i<Image::kXDim; i++) {
	    	for (int j=0; j<Image::kYDim; j++) {
	    		m_flatField[i][j] = (*flatField)[j][i];
	    	}
	    }

	    //Define the flat field to optionally be used for flat field correction in DataReduction
	    //(corresponds to the input file with no Gaussian smearing applied)
		data->setFlatFieldToSubtract(generateFlatField(dummyStats));

	    //Define the flat field with optional smearing
		data->setFlatField(generateFlatField(flatFieldStats,m_flatFieldSmearSigma>1.E-10));

    } else {

    	//Open the flat field fits file
    	flatField = new RefAppFlatfieldteff(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_flatFieldFilename,"READONLY");
    	flatFieldMetadata = Open_RefAppFlatfieldteffmetadata(flatField);

    	ofstream referenceFilesList("reference_files.txt", ios::app);
    	referenceFilesList << m_flatFieldFilename << endl;
    	referenceFilesList.close();

    	//Determine the row with Teff closest to the target star Teff
	    vector<double> TeffsFromFile;
	    while (flatFieldMetadata->ReadRow()) {
	    	TeffsFromFile.push_back(flatFieldMetadata->getCellTEff());
	    }
	    unsigned iTeff = 0;
	    double minDiff = 1.E8;
	    for (unsigned i=0; i<TeffsFromFile.size(); i++) {
	    	double diff = fabs(m_flatFieldTeff - TeffsFromFile[i]);
	    	if (diff < minDiff) {
	    		minDiff = diff;
	    		iTeff = i;
	    	}
	    }

	    //Determine the adjacent row in the direction of m_flatFieldTeffOffset (positive or negative) relative to the target Teff
	    unsigned iTeff_offset = iTeff;
	    if (m_flatFieldTeffOffset < 0. && iTeff>0) {
	    	iTeff_offset = iTeff-1;
	    } else if (m_flatFieldTeffOffset > 0. && iTeff < TeffsFromFile.size()-1) {
	    	iTeff_offset = iTeff+1;
	    }
		double mu = 0.;
		if (iTeff_offset != iTeff) {
			mu = fabs(m_flatFieldTeffOffset) / fabs(TeffsFromFile[iTeff_offset] - TeffsFromFile[iTeff]);
			if (mu>1.) mu=1.;
			cout << "Using offset TEff weighting for flat field. Interpolating between TEff=" << TeffsFromFile[iTeff] << " (nominal) and TEff=" << TeffsFromFile[iTeff_offset] << " using mu=" << mu << endl;
		}

	    m_flatFieldTeff = TeffsFromFile[iTeff];

	    //Fill the flat field array with optional Teff offset and smearing and add it to the data
	    for (int i=0; i<Image::kXDim; i++) {
	    	for (int j=0; j<Image::kYDim; j++) {
	    		m_flatField[i][j] = (*flatField)[iTeff][j][i]*(1.-mu)+(*flatField)[iTeff_offset][j][i]*mu;
	    	}
	    }
		//cout << "applied flat field: ";
		data->setFlatField(generateFlatField(flatFieldStats,m_flatFieldSmearSigma>1.E-10));

	    if (m_flatFieldToSubtractFilename == "" || m_flatFieldToSubtractFilename == m_flatFieldFilename) {

	    	//Fill the flat field array to optionally be used for flat field correction in DataReduction
	    	//(corresponds to the reference file with no TEff offset or Gaussian smearing applied)
	    	for (int i=0; i<Image::kXDim; i++) {
	    		for (int j=0; j<Image::kYDim; j++) {
	    			m_flatField[i][j] = (*flatField)[iTeff][j][i];
	    		}
	    	}
			//cout << "flat field to subtract: ";
	    	data->setFlatFieldToSubtract(generateFlatField(dummyStats));

	    } else {

	    	//Open the flat field fits file to optionally be used for flat field correction in DataReduction
	    	delete flatField;
	    	if (flatFieldMetadata != nullptr) delete flatFieldMetadata;
	    	flatField = new RefAppFlatfieldteff(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_flatFieldToSubtractFilename,"READONLY");
	    	flatFieldMetadata = Open_RefAppFlatfieldteffmetadata(flatField);

	    	//Determine the row with Teff closest to the target star Teff
		    TeffsFromFile.clear();
		    while (flatFieldMetadata->ReadRow()) {
		    	TeffsFromFile.push_back(flatFieldMetadata->getCellTEff());
		    }
		    iTeff = 0;
		    minDiff = 1.E8;
		    for (unsigned i=0; i<TeffsFromFile.size(); i++) {
		    	double diff = fabs(m_flatFieldTeff - TeffsFromFile[i]);
		    	if (diff < minDiff) {
		    		minDiff = diff;
		    		iTeff = i;
		    	}
		    }

	    	//Fill the flat field array to optionally be used for flat field correction in DataReduction
	    	for (int i=0; i<Image::kXDim; i++) {
	    		for (int j=0; j<Image::kYDim; j++) {
	    			m_flatField[i][j] = (*flatField)[iTeff][j][i];
	    		}
	    	}
	    	data->setFlatFieldToSubtract(generateFlatField(dummyStats));

	    }

    }

    cout << "Flat field mean=" << flatFieldStats.mean() << ", standard deviation=" << sqrt(flatFieldStats.variance1()) << endl;

	delete flatField;
	if (flatFieldMetadata != nullptr) delete flatFieldMetadata;

}

Image * FlatFieldGenerator::generateFlatField(boost::math::tools::stats<double> & stats, bool withRandomNoise) {

	if (!m_gaussianFlatField) {
		//Calculate the mean within the well illuminated region (omitting 100 columns and rows around the edges of the CCD)
		double flatFieldMean = 0.;
		int nPixels = 0;
		for (int i=100; i<Image::kXDim-100; i++) {
			for (int j=100; j<Image::kYDim-100; j++) {
				flatFieldMean += m_flatField[i][j];
				nPixels +=1;
			}
		}
		flatFieldMean /= nPixels;

		//Normalize the flat field to have mean 1.0 and scale the deviations w.r.t. the mean according to the user provided scale factor
		for (int i=0; i<Image::kXDim; i++) {
			for (int j=0; j<Image::kYDim; j++) {
				double flatField = m_flatField[i][j];
				flatField /= flatFieldMean;
				flatField -= 1.;
				flatField *= m_flatFieldScaleFactor;
				flatField += 1.;
				if (flatField < 0.) flatField=0.;
				if (withRandomNoise) flatField *= (*m_gaussianNoiseGenerator2)();
				m_flatField[i][j] = flatField;
				//if (i==100) cout << i << " " << j << " " << m_flatField[i][j] << endl;
				if (i>=100 && i<Image::kXDim-100 && j>=100 && j<Image::kYDim-100) stats.add(m_flatField[i][j]);
			}
		}
	}

	//Full frame image
	Image * flatFieldImage = new Image(Image::kXDim,Image::kYDim,0,0);

	for (int ix=flatFieldImage->getXOffset(); ix<flatFieldImage->getXOffset()+flatFieldImage->getXDim(); ix++) {
		for (int iy=flatFieldImage->getYOffset(); iy<flatFieldImage->getYOffset()+flatFieldImage->getYDim(); iy++) {

			double flatField = 1.;
			if (m_applyFlatField) {
				flatField *= m_flatField[ix][iy];
			}
			if (m_doDeadPixels) flatField *= m_deadPixelQE[ix][iy];
			flatFieldImage->incrementPixelValue(ix,iy,flatField);
			//if (ix==512 && iy==512) cout << flatFieldImage->getPixelValue(ix,iy) << endl;

		}
	}

	return flatFieldImage;

}

void FlatFieldGenerator::setDeadPixels(Data* data) {

	boost::mt19937 randomNumberEngine(m_deadPixelPositionSeed);
	boost::mt19937 randomNumberEngine2(m_deadPixelPositionSeed+1);
	boost::uniform_int<int> uniform_randomX(0,Image::kXDim-1);
	boost::uniform_int<int> uniform_randomY(0,Image::kYDim-1);
	RANDOM_UNIFORM uniformRandomGeneratorX(randomNumberEngine,uniform_randomX);
	RANDOM_UNIFORM uniformRandomGeneratorY(randomNumberEngine2,uniform_randomY);
	for (int ibp=0; ibp<m_nDeadPixels; ibp++) {
		int ix_bp = uniformRandomGeneratorX();
		int iy_bp = uniformRandomGeneratorY();
		m_deadPixelQE[ix_bp][iy_bp] = m_deadPixelRelativeQE;
		//cout << ix_bp << " " << iy_bp << " " << m_deadPixelQE[ix_bp][iy_bp] << endl;

		//Add the dead pixels to the bad pixel reference file
		data->setBadPixel(Image::kLeftMargin+ix_bp, iy_bp, m_deadPixelRelativeQE<=0 ? -2 : -1);
	}

}

void FlatFieldGenerator::setEmpiricalDeadPixels(Data* data) {

	//Add the dead pixels to the bad pixel reference file
	for (int ix=1; ix<1023; ix++) { //Exclude first and last columns which are not fully illuminated
		for (int iy=0; iy<1024; iy++) {
			//Exclude corner pixels which were not fully illuminated
			double radius = sqrt(pow(ix-512,2)+pow(iy-512,2));
			if (radius<685.) {
				double sensitivity = data->getFlatField()->getPixelValue(ix,iy);
				if (m_deadPixelQE[ix][iy]==1. && sensitivity < 0.8) {
					data->setBadPixel(Image::kLeftMargin+ix, iy, -1);
					cout << "empirical dead pixel at (" << ix-Image::kLeftMargin << "," << iy << "): sensitivity=" << sensitivity << endl;
				}
			}
		}
	}

}

void FlatFieldGenerator::writeTruthFlatField(Data * data) const {

	//Get the validity interval
	TimeConfiguration timeConf = data->getTimeConfiguration();
	UTC startTime_utc = timeConf.getVisitStartTimeUTC();
	UTC endTime_utc = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());

	//Open the output file
	Data::Visit visit = data->getVisit();
	string dirname = data->getOutputDirectory()+"/truth";
	if (!boost::filesystem::exists(dirname)) boost::filesystem::create_directory(dirname);
	SimTruFlatfield * truthFlatFieldFile = createFitsFile<SimTruFlatfield>(dirname,startTime_utc.getUtc(),
			VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visit.m_visitCtr), PassId());

	const WavelengthDependence * wavelengthDependence = data->getWavelengthDependence();
	string throughputFilename = wavelengthDependence->applyThroughput() ? wavelengthDependence->getThroughputFilename() : "OFF";
	string qeFilename = wavelengthDependence->applyQE() ? wavelengthDependence->getQEFilename() : "OFF";
	string flatFieldFilename = m_gaussianFlatField ? "Gaussian" : m_flatFieldFilename;
	double flatFieldScaleFactor = m_gaussianFlatField ? m_flatFieldSigma : m_flatFieldScaleFactor;

	//Set the keywords
	truthFlatFieldFile->setKeyTeff(m_flatFieldTeff);
	truthFlatFieldFile->setKeyFfref(flatFieldFilename);
	truthFlatFieldFile->setKeyFfscale(flatFieldScaleFactor);
	truthFlatFieldFile->setKeyThrptref(throughputFilename);
	truthFlatFieldFile->setKeyQeref(qeFilename);
	truthFlatFieldFile->setKeyPiName(visit.m_piName);
	truthFlatFieldFile->setKeyPiUid(visit.m_piUid);
	truthFlatFieldFile->setKeyProcNum(visit.m_versionNum);
	truthFlatFieldFile->setKeyArchRev(0);
	truthFlatFieldFile->setKeyVStrtU(startTime_utc);
	truthFlatFieldFile->setKeyVStrtM(startTime_utc.getMjd());
	truthFlatFieldFile->setKeyVStopU(endTime_utc);
	truthFlatFieldFile->setKeyVStopM(endTime_utc.getMjd());
	truthFlatFieldFile->setKeyObsCat("flat field");
	truthFlatFieldFile->setKeyProgtype(visit.m_progType);
	truthFlatFieldFile->setKeyProgId(visit.m_progId);
	truthFlatFieldFile->setKeyReqId(visit.m_reqId);
	truthFlatFieldFile->setKeyVisitctr(visit.m_visitCtr);
	truthFlatFieldFile->setKeyObsid(visit.m_obsId);
	truthFlatFieldFile->setKeyPrpVst1(visit.m_prpFirst);
	truthFlatFieldFile->setKeyPrpVstn(visit.m_prpLast);

	//Write the data
	for (int j = 0; j<Image::kYDim; j++) {
		for (int i = 0; i<Image::kXDim; i++) {
			(*truthFlatFieldFile)[j][i] = data->getFlatField()->getPixelValue(i,j);
		}
	}

	delete truthFlatFieldFile;

}

void FlatFieldGenerator::writeWhiteFlatField(Data* data) const {

	//Open the output file
	UTC currentTime = UTC(to_iso_extended_string(boost::posix_time::second_clock::universal_time()));
	UTC endValidity = UTC(to_iso_extended_string(boost::posix_time::ptime(boost::gregorian::date(2040,1,1))));
	RefAppWhiteflatfield * whiteFlatFieldFile = new RefAppWhiteflatfield(buildFileName(".", currentTime,
			VisitId(), PassId(), std::string(), RefAppWhiteflatfield::getExtName()), "CREATE",{Image::kXDim,Image::kYDim});

	//Set the keywords
	whiteFlatFieldFile->setKeyProvider("CHEOPSim");
	whiteFlatFieldFile->setKeyDescrip("Wavelength weighting according to global throughput only. Assumes flat spectrum for incident light.");
	whiteFlatFieldFile->setKeyProcNum(data->getVisit().m_versionNum);
	whiteFlatFieldFile->setKeyArchRev(0);
	whiteFlatFieldFile->setKeyVStrtU(currentTime);
	whiteFlatFieldFile->setKeyVStopU(endValidity);
	whiteFlatFieldFile->setKeyTemp(SatelliteData::kDefaultCcdTemperature);
	whiteFlatFieldFile->setKeyExptime(data->getTimeConfiguration().getExposureTimeAsDouble());

	//Write the data
	for (int j = 0; j<Image::kYDim; j++) {
		for (int i = 0; i<Image::kXDim; i++) {
			(*whiteFlatFieldFile)[j][i] = data->getFlatField()->getPixelValue(i,j);
		}
	}

	delete whiteFlatFieldFile;

}

void FlatFieldGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	Image * flatFieldImage = data->getFlatField();
	Image * image = data->getImages().back();

	if (m_applyFlatField) {
		if (flatFieldImage == nullptr) throw runtime_error("Error in FlatFieldGenerator::process: flat field has not been initialized");
		applyFlatField(image,flatFieldImage);
		if (data->doFrameTransferSmearing() && data->getImagesToSmearUp().size() == data->getImages().size()) {
			applyFlatField(data->getImagesToSmearUp().back(),flatFieldImage);
			applyFlatField(data->getImagesToSmearDown().back(),flatFieldImage);
		}
	}

	if (m_doDeadPixels) storeDeadPixelTruthData(image,flatFieldImage);

}

void FlatFieldGenerator::applyFlatField(Image* image, const Image* flatFieldImage) const {

	//Apply flat field for each pixel
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=image->getYOffset(); iy<image->getYOffset()+image->getYDim(); iy++) {

			double nElectrons = image->getPixelValue(ix,iy) * flatFieldImage->getPixelValue(ix,iy);
			//if (ix==512 && iy==512) cout << image->getPixelValue(ix,iy) << " " << nElectrons << " " << flatFieldImage->getPixelValue(ix,iy) << endl;

			image->setPixelValue(ix,iy,nElectrons);

		}
	}

}

void FlatFieldGenerator::storeDeadPixelTruthData(Image* image, Image * flatField) const {

	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=image->getYOffset(); iy<image->getYOffset()+image->getYDim(); iy++) {
			double deadPixelQE = m_deadPixelQE[ix][iy];
			if (deadPixelQE<1. || (flatField != nullptr && flatField->getPixelValue(ix,iy)<0.8)) image->getTruthData()->addDeadPixel(ix,iy,deadPixelQE);
		}
	}

}
