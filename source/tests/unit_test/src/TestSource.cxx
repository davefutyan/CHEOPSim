#define BOOST_TEST_MODULE TestCHEOPSim
#include "boost/test/unit_test.hpp"
#include "boost/test/unit_test_suite.hpp"
#include "boost/test/unit_test_log.hpp"
#include "boost/test/test_tools.hpp"
#include "boost/test/detail/global_typedef.hpp"
#include "boost/filesystem.hpp"

#include <cmath>
#include <pqxx/pqxx>
#include <stdlib.h>

#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "BarycentricOffset.hxx"
#include "REF_APP_FluxConversion.hxx"

#include "source/include/StarProducer.hxx"
#include "source/include/TransitFluxModulator.hxx"
#include "source/include/StellarNoiseFluxModulator.hxx"
#include "source/include/StellarVariationFluxModulator.hxx"
#include "source/include/ZodiacalLightGenerator.hxx"
#include "source/include/StrayLightGenerator.hxx"
#include "telescope/include/FocalPlaneGenerator.hxx"
#include "telescope/include/PSFGenerator.hxx"

#define private public
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct SourceFixture {
	SourceFixture() {
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
		boost::filesystem::create_directory(m_data->getOutputDirectory());

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

		//Set the reference files needed to calculate the transit time in BJD
		list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0002.fits"};
		LeapSeconds::setLeapSecondsFileNames(leapSecondFile);
		list<string> correlationFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU2015-01-01T00-00-00_AUX_RES_ObtUtcCorrelation_V000.fits"};
		OBTUTCCorrelation::SetCorrelationFileNames(correlationFile);
		BarycentricOffset::setEphemerisFileName(string(getenv("CHEOPS_TESTDATA"))+"/resources/CH_TU1949-12-14T00-00-00_EXT_APP_DE1_V0100.fits");

	}
	~SourceFixture() {
		delete m_data;
	}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestCHEOPSim, SourceFixture )

BOOST_AUTO_TEST_CASE( testStarProducer )
{
	Module * starProducer = new StarProducer();
	starProducer->initialize(m_params->Module("StarProducer"));
	starProducer->doBegin(m_data);
	vector<Star*> stars = m_data->getFieldOfView()->getStars();

	BOOST_CHECK_EQUAL( stars.size(), 137 );

	for (unsigned iStar=0; iStar<stars.size(); iStar++) {
		double photonFlux = stars[iStar]->getMeanPhotonFlux()*PSFGenerator::kTelescopeArea;
		if (iStar==0) {
			BOOST_CHECK_CLOSE( photonFlux, 970812., 0.0001);
		} else {
			BOOST_CHECK( photonFlux > 100. && photonFlux < 1.E6);
		}
	}

	delete starProducer;
}

BOOST_AUTO_TEST_CASE( testMagnitudeConversion )
{

	string filename = "resources/CH_TU2020-01-29T00-00-00_REF_APP_FluxConversion_V0104.fits";
	RefAppFluxconversion * conversion_file = new RefAppFluxconversion(filename,"READONLY");

	FluxConverter * fluxConv = new FluxConverter(-0.004795,0.94604,-2.0074e-06);

	for (int TEff=2300; TEff<=40000; TEff+=100) {
		double deltaMag = fluxConv->cheopsMinusRefBandMagnitude(WavelengthDependence::GaiaBand,TEff,m_data->getWavelengthDependence());

		double deltaMag_ref = -999.;
		while (conversion_file->ReadRow()) {
			if (conversion_file->getCellTEff() == TEff) {
				deltaMag_ref = conversion_file->getCellCheopsmagMinusGmag();
				break;
			}
		}

		BOOST_CHECK_CLOSE( deltaMag, deltaMag_ref, 1E-6 );

	}

	delete fluxConv;
	delete conversion_file;
}

BOOST_AUTO_TEST_CASE( testTransit )
{
	Module * starProducer = new StarProducer();
	starProducer->initialize(m_params->Module("StarProducer"));
	starProducer->doBegin(m_data);

	Module * transitFluxModulator = new TransitFluxModulator(0);
	transitFluxModulator->initialize(m_params->Module("TransitFluxModulator_star0"));
	transitFluxModulator->doBegin(m_data);

	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
		transitFluxModulator->process(m_data,iTimeStep);
	}

	Star * star = m_data->getFieldOfView()->getStars()[0];
	BOOST_CHECK_CLOSE( star->getTimeSeries()[4]->getTransitFluxFactor(), 0.998783765, 1E-6);

	delete starProducer;
	delete transitFluxModulator;
}

BOOST_AUTO_TEST_CASE( testStellarNoise )
{
	Module * starProducer = new StarProducer();
	starProducer->initialize(m_params->Module("StarProducer"));
	starProducer->doBegin(m_data);

	Module * stellarNoiseFluxModulator = new StellarNoiseFluxModulator(0);
	stellarNoiseFluxModulator->initialize(m_params->Module("StellarNoiseFluxModulator_star0"));
	stellarNoiseFluxModulator->doBegin(m_data);

	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
		stellarNoiseFluxModulator->process(m_data,iTimeStep);
	}

	Star * star = m_data->getFieldOfView()->getStars()[0];
	BOOST_CHECK_CLOSE( star->getTimeSeries()[4]->getNoiseFluxFactor(), 1.0000331557123645, 1E-6);

	delete starProducer;
	delete stellarNoiseFluxModulator;
}

//BOOST_AUTO_TEST_CASE( testStellarVariation )
//{
//	Module * starProducer = new StarProducer();
//	starProducer->initialize(m_params->Module("StarProducer"));
//	starProducer->doBegin(m_data);
//
//	Module * stellarVariationFluxModulator = new StellarVariationFluxModulator(0);
//	stellarVariationFluxModulator->initialize(m_params->Module("StellarVariationFluxModulator_star0"));
//	stellarVariationFluxModulator->doBegin(m_data);
//
//	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
//	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
//		stellarVariationFluxModulator->process(m_data,iTimeStep);
//	}
//
//	Star * star = m_data->getFieldOfView()->getStars()[0];
//	BOOST_CHECK_CLOSE( star->getTimeSeries()[4]->getVariationFluxFactor(), 0.998786714, 1E-6);
//
//	delete starProducer;
//	delete stellarVariationFluxModulator;
//}

BOOST_AUTO_TEST_CASE( testZodiacalLight )
{
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	focalPlaneGenerator->doBegin(m_data);

	Module * zodiacalLightGenerator = new ZodiacalLightGenerator();
	zodiacalLightGenerator->initialize(m_params->Module("ZodiacalLightGenerator"));
	zodiacalLightGenerator->doBegin(m_data);

	focalPlaneGenerator->process(m_data,0);
	zodiacalLightGenerator->process(m_data,0);

	BOOST_CHECK_CLOSE( (*m_data->getImages().begin())->getPixelValue(512,512), 114.357, 1E-3);

	delete focalPlaneGenerator;
	delete zodiacalLightGenerator;
}

BOOST_AUTO_TEST_CASE( testCoordinateConversion )
{

	string ra = "220:00:00.0";
	string dec = "47:00:00.0";
	double ra_arcsec = StarProducer::convertToArcseconds(ra,true);
	double dec_arcsec = StarProducer::convertToArcseconds(dec,false);
	BOOST_CHECK_CLOSE( ra_arcsec, 1.188e+07, 0.1);
	BOOST_CHECK_CLOSE( dec_arcsec, 169200, 0.1);
	SkyPosition skypos(ra_arcsec,dec_arcsec);

	double eclipticLatitude = skypos.getEclipticLatitude();
	double eclipticLongitude = skypos.getEclipticLongitude();
	BOOST_CHECK_CLOSE( eclipticLatitude, 25.8531, .0001);
	BOOST_CHECK_CLOSE( eclipticLongitude, 67.733, .0001);

}

BOOST_AUTO_TEST_SUITE_END()
