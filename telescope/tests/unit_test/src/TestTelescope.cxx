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

#include "CreateFitsFile.hxx"
#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "source/include/StarProducer.hxx"
#include "satellite/include/OrbitSimulator.hxx"
#include "telescope/include/FocalPlaneGenerator.hxx"
#include "telescope/include/PSFGenerator.hxx"
#include "telescope/include/HaloGenerator.hxx"
#include "telescope/include/GlobalThroughputGenerator.hxx"

#define private public
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct TelescopeFixture {
	TelescopeFixture() {
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
			unlink((string(m_data->getOutputDirectory()+"/aux/CH_TU2018-09-01T11-55-00_AUX_RES_Orbit_V0000.fits")).c_str());
			unlink((string(m_data->getOutputDirectory()+"/data/CH_PR"+boost::lexical_cast<string>(Data::kProgramType)+"0001_TG000001_TU2018-09-01T12-00-00_SCI_RAW_Attitude_V0000.fits")).c_str());
			unlink((string(m_data->getOutputDirectory()+"/data/CH_PR"+boost::lexical_cast<string>(Data::kProgramType)+"0001_TG000001_TU2018-09-01T12-00-00_MPS_PRE_Visits_V0000.fits")).c_str());
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

		//Initialize the wavelength dependence
		bool applyThroughput = m_params->GetAsBool("applyThroughput");
		string throughputFilename = m_params->GetAsString("throughputFilename");
		bool applyQE = m_params->GetAsBool("applyQE");
		string QEFilename = m_params->GetAsString("QEFilename");
		string SEDFilename = m_params->GetAsString("SEDFilename");
		string targetSEDFilename = m_params->GetAsString("targetSEDFilename");
		m_data->setWavelengthDependence(new WavelengthDependence(applyThroughput,throughputFilename,applyQE,QEFilename,SEDFilename,targetSEDFilename));

		//Set the nominal values for each of the bias voltages, extracting them from the REF_APP_GainCorrection file
		m_data->setNominalVoltages(m_params->GetAsString("gainCorrectionFilename"));

		//Set the leap second file and the OBT2UTC correlation file needed to define onboard time
		list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0002.fits"};
		LeapSeconds::setLeapSecondsFileNames(leapSecondFile);
		list<string> correlationFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU2015-01-01T00-00-00_AUX_RES_ObtUtcCorrelation_V000.fits"};
		OBTUTCCorrelation::SetCorrelationFileNames(correlationFile);

		//Open the output MPS_PRE_Visits file
		string dirname = m_data->getOutputDirectory()+"/data";
		boost::filesystem::create_directory(dirname);
		Data::Visit visit = m_data->getVisit();
		MpsPreVisits * mpsVisits_out = new MpsPreVisits(buildFileName(dirname, timeConf.getVisitStartTimeUTC(),
				VisitId(Data::kProgramType,visit.m_progId,visit.m_reqId,visit.m_visitCtr),
				PassId(), std::string(), MpsPreVisits::getExtName()),"CREATE");
		m_data->setMpsPreVisits(mpsVisits_out);

	}
	~TelescopeFixture() {}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestTelescope, TelescopeFixture )

BOOST_AUTO_TEST_CASE( testPSFGenerator )
{
	Module * starProducer = new StarProducer();
	Module * orbitSimulator = new OrbitSimulator();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();

	starProducer->initialize(m_params->Module("StarProducer"));
	orbitSimulator->initialize(m_params->Module("OrbitSimulator"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));

	starProducer->doBegin(m_data);
	orbitSimulator->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);

	m_data->setFrameTransferSmearing();

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);

	BOOST_CHECK_CLOSE( (m_data->getImages()[0])->getPixelValue(512,512), 22713.913, 1E-3);

	delete starProducer;
	delete orbitSimulator;
	delete focalPlaneGenerator;
	delete psfGenerator;
}

BOOST_AUTO_TEST_CASE( testGlobalThroughput )
{
	Module * starProducer = new StarProducer();
	Module * orbitSimulator = new OrbitSimulator();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * psfGenerator = new PSFGenerator();
	Module * globalThroughputGenerator = new GlobalThroughputGenerator();

	starProducer->initialize(m_params->Module("StarProducer"));
	orbitSimulator->initialize(m_params->Module("OrbitSimulator"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	psfGenerator->initialize(m_params->Module("PSFGenerator"));
	globalThroughputGenerator->initialize(m_params->Module("GlobalThroughputGenerator"));

	starProducer->doBegin(m_data);
	orbitSimulator->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	psfGenerator->doBegin(m_data);
	globalThroughputGenerator->doBegin(m_data);

	focalPlaneGenerator->process(m_data,0);
	psfGenerator->process(m_data,0);
	globalThroughputGenerator->process(m_data,0);

	BOOST_CHECK_CLOSE( (m_data->getImages()[0])->getPixelValue(512,512), 10586.088, 1E-3);

	double effectiveTemperature = (m_data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
	double integratedThroughputQE = m_data->getWavelengthDependence()->getIntegratedThroughputQE(effectiveTemperature,SatelliteData::kDefaultCcdTemperature,true);
	double integratedThroughputQE_CHEOPS_Vega = m_data->getWavelengthDependence()->getIntegratedThroughputQE(0.,SatelliteData::kDefaultCcdTemperature);
	double F0_X = StarProducer::kNPhoVega * integratedThroughputQE_CHEOPS_Vega;
	double nElectrons_CHEOPSim = (m_data->getFieldOfView()->getStars()[0])->getMeanPhotonFlux()*integratedThroughputQE*PSFGenerator::kTelescopeArea;
	double nElectrons_ETC = F0_X*pow(10,-0.4*(m_data->getFieldOfView()->getStars()[0])->getCheopsMagnitude());
	BOOST_CHECK_CLOSE( nElectrons_CHEOPSim, nElectrons_ETC, 1.);

	delete starProducer;
	delete orbitSimulator;
	delete focalPlaneGenerator;
	delete psfGenerator;
	delete globalThroughputGenerator;
}

BOOST_AUTO_TEST_CASE( testHaloGenerator )
{
	Module * starProducer = new StarProducer();
	Module * orbitSimulator = new OrbitSimulator();
	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	Module * haloGenerator = new HaloGenerator();

	starProducer->initialize(m_params->Module("StarProducer"));
	orbitSimulator->initialize(m_params->Module("OrbitSimulator"));
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	haloGenerator->initialize(m_params->Module("HaloGenerator"));

	starProducer->doBegin(m_data);
	orbitSimulator->doBegin(m_data);
	focalPlaneGenerator->doBegin(m_data);
	haloGenerator->doBegin(m_data);

	focalPlaneGenerator->process(m_data,0);
	haloGenerator->process(m_data,0);

	BOOST_CHECK_CLOSE( (m_data->getImages()[0])->getPixelValue(512,512), 18.725986, 1E-3);

	delete starProducer;
	delete focalPlaneGenerator;
	delete haloGenerator;
}

BOOST_AUTO_TEST_SUITE_END()
