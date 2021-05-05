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

#include "AUX_RES_Orbit.hxx"

#include "CreateFitsFile.hxx"
#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "source/include/StarProducer.hxx"
#include "satellite/include/OrbitSimulator.hxx"
#include "satellite/include/JitterProducer.hxx"

#define private public
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct SatelliteFixture {
	SatelliteFixture() {
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
	~SatelliteFixture() {}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestSatellite, SatelliteFixture )

BOOST_AUTO_TEST_CASE( testOrbitModel )
{
	Module * orbitSimulator = new OrbitSimulator();
	orbitSimulator->initialize(m_params->Module("OrbitSimulator"));
	orbitSimulator->doBegin(m_data);

	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
		orbitSimulator->process(m_data,iTimeStep);
	}
	delete orbitSimulator;

	SatelliteData * satData = m_data->getSatelliteData(100-m_data->getTimeConfiguration().getNumberOfPreSubArrayTimeSteps());
	BOOST_CHECK_CLOSE( satData->getTelescopeTemperature(), 262.542, 01E-3);
	BOOST_CHECK_CLOSE( satData->getCcdTemperature(), 233.008606, 01E-6);
	BOOST_CHECK_CLOSE( satData->getEarthLimbAngle(), 6.145456, 01E-3);
	BOOST_CHECK_CLOSE( satData->getRollAngles()[0], 85.399807, 01E-3);
	BOOST_CHECK_EQUAL( satData->getSAAFlag(), false);

	//Open the fits file containing the orbit data
	string orbitfilename = buildFileName(m_data->getOutputDirectory()+"/aux", m_data->getTimeConfiguration().getVisitStartTimeUTC()-DeltaTime(300),
										 VisitId(), PassId(), std::string(), AuxResOrbit::getExtName());
	AuxResOrbit * orbitData = new AuxResOrbit(string(getenv("PWD"))+"/"+orbitfilename,"READONLY");

	while (orbitData->ReadRow()) {

		double x = orbitData->getCellX();
		double y = orbitData->getCellY();
		double z = orbitData->getCellZ();
		double vx = orbitData->getCellXDot();
		double vy = orbitData->getCellYDot();
		double vz = orbitData->getCellZDot();

		double orbitRadius = sqrt(x*x + y*y + z*z);
		double orbitSpeed = sqrt(vx*vx + vy*vy + vz*vz);
		double nominalRadius = 6371.+650.;
		double nominalSpeed = 2*M_PI*nominalRadius/5850.;

		BOOST_CHECK_CLOSE( orbitRadius, nominalRadius, 0.5);
		BOOST_CHECK_CLOSE( orbitSpeed, nominalSpeed, 0.5);

	}


}

BOOST_AUTO_TEST_SUITE_END()
