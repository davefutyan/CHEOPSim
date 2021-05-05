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
#include "detector/include/FullWellSimulator.hxx"

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

		//Initialize the field of view by defining the pointing direction
		string pointingRA_str = m_params->GetAsString("pointingRA");
		string pointingDec_str = m_params->GetAsString("pointingDec");
		double pointingRA = StarProducer::convertToArcseconds(pointingRA_str,true);
		double pointingDec = StarProducer::convertToArcseconds(pointingDec_str,false);
		SkyPosition pointingDirection(pointingRA,pointingDec);
		SkyFieldOfView *fov = new SkyFieldOfView(pointingDirection);
		m_data->setFieldOfView(fov);

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

BOOST_AUTO_TEST_CASE( testFrameTransfer )
{
	Module * starProducer = new StarProducer();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();
	Module * frameTransferSmearer = new FrameTransferSmearer();

	starProducer->initialize(m_params->Module("StarProducer"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));
	frameTransferSmearer->initialize(m_params->Module("FrameTransferSmearer"));

	starProducer->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);
	frameTransferSmearer->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);
	frameTransferSmearer->process(m_data,0);

	Image * image = m_data->getImages().back();
	double trailBrightnessRatio = image->getPixelValue(512,600)/image->getPixelValue(512,512);
	BOOST_CHECK_CLOSE( trailBrightnessRatio, 1.780653e-05, 1E-3);

	delete focalPlaneGenerator;
	delete psfGenerator;
	delete frameTransferSmearer;
}

BOOST_AUTO_TEST_CASE( testFullWell )
{
	Module * starProducer = new StarProducer();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();
	Module * fullWellSimulator = new FullWellSimulator();

	starProducer->initialize(m_params->Module("StarProducer"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));
	fullWellSimulator->initialize(m_params->Module("FullWellSimulator"));

	starProducer->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);
	fullWellSimulator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);
	fullWellSimulator->process(m_data,0);

	Image * image = m_data->getImages().back();
	double fullWellCapacity = m_params->Module("FullWellSimulator").GetAsInt("fullWellCapacity");
	BOOST_CHECK_EQUAL( image->getPixelValue(512,512), fullWellCapacity);

	delete focalPlaneGenerator;
	delete psfGenerator;
	delete fullWellSimulator;
}

BOOST_AUTO_TEST_SUITE_END()
