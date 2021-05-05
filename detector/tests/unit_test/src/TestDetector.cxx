#define BOOST_TEST_MODULE TestCHEOPSim
#include "boost/test/unit_test.hpp"
#include "boost/test/unit_test_suite.hpp"
#include "boost/test/unit_test_log.hpp"
#include "boost/test/test_tools.hpp"
#include "boost/test/detail/global_typedef.hpp"
#include "boost/filesystem.hpp"
#include "boost/math/tools/stats.hpp"

#include <cmath>
#include <pqxx/pqxx>
#include <stdlib.h>

#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "source/include/StarProducer.hxx"
#include "telescope/include/FocalPlaneGenerator.hxx"
#include "telescope/include/PSFGenerator.hxx"
#include "detector/include/FlatFieldGenerator.hxx"
#include "detector/include/DarkCurrentGenerator.hxx"
#include "detector/include/FrameTransferSmearer.hxx"
#include "detector/include/PhotonNoiseGenerator.hxx"
#include "detector/include/CosmicRayGenerator.hxx"
#include "detector/include/ChargeTransferSimulator.hxx"
#include "detector/include/BiasGenerator.hxx"

#define private public
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct DetectorFixture {
	DetectorFixture() {
		m_params = CheopsInit(framework::master_test_suite().argc,
						      framework::master_test_suite().argv);

		// Get the time configuration
		boost::posix_time::ptime startTime = boost::posix_time::time_from_string(m_params->GetAsString("startTime"));
		unsigned numberOfStackedImages = m_params->GetAsInt("numberOfStackedImages");
		unsigned exposuresPerStack = m_params->GetAsInt("exposuresPerStack");
		double exposureTime = m_params->GetAsDouble("exposureTime");
		unsigned maxImagesPerCube = m_params->GetAsInt("maxImagesPerCube");

		// Initialize the time configuration
		TimeConfiguration timeConf(startTime,numberOfStackedImages,exposuresPerStack,exposureTime,maxImagesPerCube,m_params->GetAsBool("doFullFrame"));
		// Initialize the data
		m_data = new Data(timeConf,m_params->GetAsString("modulesToRun"));

	    //Create the output directory
		if (boost::filesystem::exists(m_data->getOutputDirectory())) {
			unlink((string(m_data->getOutputDirectory()+"/data/CH_TU2018-09-01T12-00-30_REF_APP_BadPixelMap_V0000.fits")).c_str());
			unlink((string(m_data->getOutputDirectory()+"/truth/CH_PR"+boost::lexical_cast<string>(Data::kProgramType)+"0001_TG000001_TU2018-09-01T12-00-00_SIM_TRU_FlatField_V0000.fits")).c_str());
		} else {
			boost::filesystem::create_directory(m_data->getOutputDirectory());
		}

		//Initialize the field of view by defining the pointing direction
		string pointingRA_str = m_params->GetAsString("pointingRA");
		string pointingDec_str = m_params->GetAsString("pointingDec");
		double pointingRA = StarProducer::convertToArcseconds(pointingRA_str,true);
		double pointingDec = StarProducer::convertToArcseconds(pointingDec_str,false);
		SkyPosition pointingDirection(pointingRA,pointingDec);
		SkyFieldOfView *fov = new SkyFieldOfView(pointingDirection);
		m_data->setFieldOfView(fov);

		//Set the nominal values for each of the bias voltages, extracting them from the REF_APP_GainCorrection file
		m_data->setNominalVoltages(m_params->GetAsString("gainCorrectionFilename"));

		//Initialize the wavelength dependence
		bool applyThroughput = m_params->GetAsBool("applyThroughput");
		string throughputFilename = m_params->GetAsString("throughputFilename");
		bool applyQE = m_params->GetAsBool("applyQE");
		string QEFilename = m_params->GetAsString("QEFilename");
		string SEDFilename = m_params->GetAsString("SEDFilename");
		string targetSEDFilename = m_params->GetAsString("targetSEDFilename");
		m_data->setWavelengthDependence(new WavelengthDependence(applyThroughput,throughputFilename,applyQE,QEFilename,SEDFilename,targetSEDFilename));

		//Set the leap second file and the OBT2UTC correlation file needed to define onboard time
		list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0002.fits"};
		LeapSeconds::setLeapSecondsFileNames(leapSecondFile);
		list<string> correlationFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU2015-01-01T00-00-00_AUX_RES_ObtUtcCorrelation_V000.fits"};
		OBTUTCCorrelation::SetCorrelationFileNames(correlationFile);

	}
	~DetectorFixture() {}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestDetector, DetectorFixture )

BOOST_AUTO_TEST_CASE( testFlatField )
{
	Module * starProducer = new StarProducer();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();
	Module * flatFieldGenerator = new FlatFieldGenerator();

	starProducer->initialize(m_params->Module("StarProducer"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));
	flatFieldGenerator->initialize(m_params->Module("FlatFieldGenerator"));

	starProducer->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);
	flatFieldGenerator->doBegin(m_data);

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);
	flatFieldGenerator->process(m_data,0);

	Image * image = m_data->getImages().back();
	boost::math::tools::stats<double> flatFieldStats;
	for (int i=image->getXOffset()+1; i<image->getXOffset()+image->getXDim()-1; i++) {
		for (int j=image->getYOffset()+1; j<image->getYOffset()+image->getYDim()-1; j++) {
			flatFieldStats.add(image->getPixelValue(i,j));
		}
	}

	BOOST_CHECK_CLOSE( sqrt(flatFieldStats.variance1())/flatFieldStats.mean(), 0.0131758871, 1E-3);

	delete focalPlaneGenerator;
	delete psfGenerator;
	delete flatFieldGenerator;
}

BOOST_AUTO_TEST_CASE( testDarkCurrent )
{
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * darkCurrentGenerator = new DarkCurrentGenerator();

	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	darkCurrentGenerator->initialize(m_params->Module("DarkCurrentGenerator"));

	focalPlaneGenerator->doBegin(m_data);
	darkCurrentGenerator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	darkCurrentGenerator->process(m_data,0);
	focalPlaneGenerator->process(m_data,1);
	darkCurrentGenerator->process(m_data,1);

	delete focalPlaneGenerator;
	delete darkCurrentGenerator;
}

BOOST_AUTO_TEST_CASE( testPhotonNoise )
{
	Module * starProducer = new StarProducer();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();
	Module * shotNoiseGenerator = new PhotonNoiseGenerator();

	starProducer->initialize(m_params->Module("StarProducer"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));
	shotNoiseGenerator->initialize(m_params->Module("PhotonNoiseGenerator"));

	starProducer->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);
	shotNoiseGenerator->doBegin(m_data);

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);
	shotNoiseGenerator->process(m_data,0);
	focalPlaneGenerator->process(m_data,1);
	psfGenerator->process(m_data,1);
	shotNoiseGenerator->process(m_data,1);

	Image * image = m_data->getImages().back();
	boost::math::tools::stats<double> shotNoiseStats;
	for (int i=image->getXOffset()+1; i<image->getXOffset()+image->getXDim()-1; i++) {
		for (int j=image->getYOffset()+1; j<image->getYOffset()+image->getYDim()-1; j++) {
			shotNoiseStats.add(image->getPixelValue(i,j));
		}
	}

	BOOST_CHECK_CLOSE( shotNoiseStats.variance1(), shotNoiseStats.mean(), 1.);

	delete focalPlaneGenerator;
	delete psfGenerator;
	delete shotNoiseGenerator;
}

BOOST_AUTO_TEST_CASE( testCosmics )
{
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * cosmicRayGenerator = new CosmicRayGenerator();

	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	cosmicRayGenerator->initialize(m_params->Module("CosmicRayGenerator"));

	focalPlaneGenerator->doBegin(m_data);
	cosmicRayGenerator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	cosmicRayGenerator->process(m_data,0);

	delete focalPlaneGenerator;
	delete cosmicRayGenerator;
}

BOOST_AUTO_TEST_CASE( testChargeTransfer )
{
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * chargeTransferSimulator = new ChargeTransferSimulator();

	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	chargeTransferSimulator->initialize(m_params->Module("ChargeTransferSimulator"));

	focalPlaneGenerator->doBegin(m_data);
	chargeTransferSimulator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	chargeTransferSimulator->process(m_data,0);

	delete focalPlaneGenerator;
	delete chargeTransferSimulator;
}

BOOST_AUTO_TEST_CASE( testBias )
{
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * biasGenerator = new BiasGenerator();

	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	biasGenerator->initialize(m_params->Module("BiasGenerator"));

	focalPlaneGenerator->doBegin(m_data);
	biasGenerator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	biasGenerator->process(m_data,0);
	focalPlaneGenerator->process(m_data,1);
	biasGenerator->process(m_data,1);

	Image * image = m_data->getImages().back();
	boost::math::tools::stats<double> stats;
	for (int i=image->getXOffset()+1; i<image->getXOffset()+image->getXDim()-1; i++) {
		for (int j=image->getYOffset()+1; j<image->getYOffset()+image->getYDim()-1; j++) {
			stats.add(image->getPixelValue(i,j));
		}
	}

	double biasMean = m_params->Module("BiasGenerator").GetAsDouble("biasMean");
	double biasWidth = m_params->Module("BiasGenerator").GetAsDouble("biasWidth");
	BOOST_CHECK_CLOSE( stats.mean(), biasMean , 0.2);
	BOOST_CHECK_CLOSE( sqrt(stats.variance1()), biasWidth , 0.5);

	delete focalPlaneGenerator;
	delete biasGenerator;
}

BOOST_AUTO_TEST_CASE( testCcdNonLinearity )
{
	BiasGenerator * biasGenerator = new BiasGenerator();

	biasGenerator->initialize(m_params->Module("BiasGenerator"));
	biasGenerator->doBegin(m_data);

	for (double nElectrons = 1.; nElectrons<160000.; nElectrons+=500.) {
		double nElectrons_before = nElectrons;
		biasGenerator->applyCcdNonLinearity(nElectrons);
		biasGenerator->correctCcdNonLinearity(nElectrons);
		BOOST_CHECK_CLOSE( nElectrons_before, nElectrons, 0.001);
	}

	delete biasGenerator;
}

BOOST_AUTO_TEST_SUITE_END()
