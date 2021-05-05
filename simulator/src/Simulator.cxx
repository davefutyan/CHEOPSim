/*
 * Simulator.cxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

#include <iostream>
#include <fstream>
#include <vector>

#include "boost/filesystem.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian_types.hpp"
#include <pqxx/pqxx>

#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "BarycentricOffset.hxx"
#include "OrbitInterpolation.hxx"

#include "CreateFitsFile.hxx"
#include "MPS_PRE_Visits.hxx"
#include "MPS_PRE_VisitConstraints.hxx"
#include "REF_APP_FluxConversion.hxx"
#include "REF_APP_SEDTeff.hxx"

#include "source/include/StarProducer.hxx"
#include "EXT_APP_ObservationRequests_schema.hxx"
#include "ModuleRegistry.hxx"
#include "Simulator.hxx"

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"

typedef boost::variate_generator<boost::mt19937,boost::normal_distribution<double> > RANDOM_GAUSS;

using namespace std;

Simulator::Simulator(const ParamsPtr params) :
		m_params(params),m_data(nullptr),m_database(nullptr),m_modulesInitialized(false),
		m_doFullFrame(false),m_jobid("") {initialize();}


void Simulator::initialize() {

	cout << endl << "=========================================================================================================" << endl;
	cout << "|                                          Welcome to CHEOPSim!                                         |" << endl;
	cout << "=========================================================================================================" << endl;

	//Open the connection to the job submission database
	if (m_params->GetAsBool("updateDatabase")) connectToDatabase(m_params->GetAsString("databaseServer"));

	// Get the visit information
	boost::posix_time::ptime startTime = boost::posix_time::time_from_string(m_params->GetAsString("startTime"));
	unsigned numberOfStackedImages = m_params->GetAsInt("numberOfStackedImages");
	unsigned exposuresPerStack = m_params->GetAsInt("exposuresPerStack");
	double exposureTimeAsDouble = m_params->GetAsDouble("exposureTime");
	m_doFullFrame = m_params->GetAsBool("doFullFrame");
	string pointingRA_str = m_params->GetAsString("pointingRA");
	string pointingDec_str = m_params->GetAsString("pointingDec");
	double pointingRA = StarProducer::convertToArcseconds(pointingRA_str,true);
	double pointingDec = StarProducer::convertToArcseconds(pointingDec_str,false);
	unsigned maxImagesPerCube = m_params->GetAsInt("maxImagesPerCube");

	//Initialize some visit related keywords with default values
	int progType=Data::kProgramType;
	try {
		progType = m_params->GetAsInt("programType");
	} catch (const std::exception& e) {
		cerr << "Warning: programType parameter not found in configuration file, using default value " << progType << endl;
	}
	int progId = 1;
	int reqId = 0;
	int visitCtr = 1;
	string obsCat = "time critical";
	string piName = "CHEOPSim";
	try {
		piName = m_params->GetAsString("piName");
	} catch (const std::exception& e) {
		cerr << "Warning: piName parameter not found in configuration file, using default value " << piName << endl;
	}
	int piUid = 0;
	try {
		piUid = m_params->GetAsInt("piUid");
	} catch (const std::exception& e) {
		cerr << "Warning: piUid parameter not found in configuration file, using default value " << piUid << endl;
	}
	int prpFirst = 365;
	int prpLast = 548;
	string targName = "simulation";
	string readMode = "ultrabright";
	if (exposureTimeAsDouble>1.05) readMode = "bright";
	if (exposureTimeAsDouble>2.326) readMode = "faint fast";
	if (exposureTimeAsDouble>11.296) readMode = "faint";
	string marginMode = "undefined";
	double gaiaMagError = 0.;
	double strayLightThreshold = m_params->GetAsDouble("strayLightThreshold");

	// If an MPS_PRE_Visits or EXT_APP_ObservationRequests file is provided, read information from the file
	string visitFilename = m_params->GetAsString("MPSFilename");
	unsigned obsId = m_params->GetAsInt("obsId");

	MpsPreVisitconstraints * visitConstraints = nullptr;
	if (visitFilename.find("MPS_PRE_Visits") != string::npos) {

		MpsPreVisits * mpsVisits = new MpsPreVisits(visitFilename,"READONLY");
		visitConstraints = Open_MpsPreVisitconstraints(mpsVisits);
		if (visitConstraints == nullptr) {
			throw runtime_error("Error in Simulator::initialize: Error opening MPS_PRE_VisitConstraints extension of the input MPS_PRE_Visits fits file");
		}

		bool found = false;
		while (mpsVisits->ReadRow()) {
			if (mpsVisits->getCellObsid() == obsId) {
				found = true;
				progType = mpsVisits->getCellProgrammeType();
				progId = mpsVisits->getCellProgrammeId();
				reqId = mpsVisits->getCellRequestId();
				obsCat = mpsVisits->getCellObsCategory();
				piName = mpsVisits->getCellPiName();
				piUid = mpsVisits->getCellPiUid();
				prpFirst = mpsVisits->getCellPropFirstVisit();
				prpLast = mpsVisits->getCellPropLastVisit();
				targName = mpsVisits->getCellTargetName();
				if (mpsVisits->getCellStrayLightThreshold() > 0.) {
					strayLightThreshold = mpsVisits->getCellStrayLightThreshold();
				}
				readMode = mpsVisits->getCellReadoutMode();
				if (readMode != "faint" && readMode != "faint fast" && readMode != "bright" && readMode != "ultrabright") {
					throw runtime_error("ERROR in Simulator::initialize: READOUT_MODE read from MPS_PRE_Visits file must be faint, faint fast, bright or ultrabright. Value found: "+readMode);
				}
				marginMode = mpsVisits->getCellMarginMode();
				cout << "reading margin mode as " << marginMode << endl;
				if (marginMode != "image" && marginMode != "reduced" && marginMode != "total collapsed") {
					throw runtime_error("ERROR in Simulator::initialize: MARGIN_MODE read from MPS_PRE_Visits file must be image, reduced or total collapsed. Value found: "+marginMode);
				}
				break;
			}
		}
		if (!found) throw runtime_error("ERROR in Simulator::initialize: No entry found in MPS_PRE_Visits file for requested ObsId "+to_string(obsId));
		delete mpsVisits;

	} else if (visitFilename.find("EXT_APP_ObservationRequests") != string::npos) {

		// Get the schema file name of the observation requests file
		string schemaFileName = string(getenv("CHEOPS_SW"))+"/resources/EXT_APP_ObservationRequests_schema.xsd";
		if (!boost::filesystem::exists(schemaFileName)) {
			throw runtime_error("Failed to find schema of observation request file: " + schemaFileName);
		}

		// Get the xml
		xml_schema::properties props;
	    props.no_namespace_schema_location(schemaFileName);
	    props.schema_location ("http://www.w3.org/XML/1998/namespace", "xml.xsd");
	    auto_ptr<Earth_Explorer_File_type> observationRequest_xml;
	    try {
	    	observationRequest_xml = auto_ptr<Earth_Explorer_File_type>(Earth_Explorer_File(visitFilename, 0, props));
	    }
	    // Error parsing the param file. Need to pass the error object to a stream
	    // in order to get the complete error message that includes information
	    // about the line on which the parse error occurred.
	    catch (const xml_schema::exception& e) {
	    	std::stringstream message;
	    	message << e;
	    	throw runtime_error(message.str());
	    }

	    stringstream strValue;
	    strValue << observationRequest_xml->Earth_Explorer_Header().Variable_Header().Programme_Type();
	    strValue >> progType;


	    if (observationRequest_xml->Data_Block().List_of_Time_Critical_Requests().present()) {

	    	// Get the time critical request parameters
	    	List_of_Time_Critical_Requests_type time_critical_requests = observationRequest_xml->Data_Block().List_of_Time_Critical_Requests().get();
	    	List_of_Time_Critical_Requests_type::Time_Critical_Request_sequence & requests = time_critical_requests.Time_Critical_Request();
	    	List_of_Time_Critical_Requests_type::Time_Critical_Request_iterator observationRequest = requests.begin();
	    	//while (observationRequest != requests.end()) {observationRequest++;}

	    	// Read in the parameters
	    	progId = observationRequest->Programme_ID();
	    	reqId = observationRequest->Observation_Request_ID();
	    	obsCat = observationRequest->Observation_Category();
	    	prpFirst = observationRequest->Proprietary_Period_First_Visit();
	    	prpLast = observationRequest->Proprietary_Period_Last_Visit();
	    	targName = observationRequest->Target_Name();
	    	gaiaMagError = observationRequest->Target_Magnitude_Error();
	    	readMode = observationRequest->Readout_Mode();

	    } else if (observationRequest_xml->Data_Block().List_of_Non_Time_Critical_Requests().present()) {

	    	// Get the time critical request parameters
	    	List_of_Non_Time_Critical_Requests_type non_time_critical_requests = observationRequest_xml->Data_Block().List_of_Non_Time_Critical_Requests().get();
	    	List_of_Non_Time_Critical_Requests_type::Non_Time_Critical_Request_sequence & requests = non_time_critical_requests.Non_Time_Critical_Request();
	    	List_of_Non_Time_Critical_Requests_type::Non_Time_Critical_Request_iterator observationRequest = requests.begin();
	    	//while (observationRequest != requests.end()) {observationRequest++;}

	    	// Read in the parameters
	    	progId = observationRequest->Programme_ID();
	    	reqId = observationRequest->Observation_Request_ID();
	    	obsCat = observationRequest->Observation_Category();
	    	prpFirst = observationRequest->Proprietary_Period_First_Visit();
	    	prpLast = observationRequest->Proprietary_Period_Last_Visit();
	    	targName = observationRequest->Target_Name();
	    	gaiaMagError = observationRequest->Target_Magnitude_Error();
	    	readMode = observationRequest->Readout_Mode();

	    } else {
			throw runtime_error("ERROR in Simulator::initialize: Observation request xml file contains neither Time_Critical_Requests nor Non_Time_Critical_Requests");
	    }

	    if (readMode != "faint" && readMode != "faint fast" && readMode != "bright" && readMode != "ultrabright") {
			throw runtime_error("ERROR in Simulator::initialize: READOUT_MODE read from EXT_APP_ObservationRequests file must be faint, faint fast, bright or ultrabright. Value found: "+readMode);
		}

	} else if (visitFilename!="null" && visitFilename!="") {

		throw runtime_error("Visit filename must contain the string MPS_PRE_Visits or EXP_APP_ObservationRequests");

	} else if (m_database != nullptr) {

		// If progId and reqId are not provided in an MPS_PRE_Visits or EXT_APP_ObservationRequests file, define values based on the job submission database
		pqxx::work transaction(*m_database);
		pqxx::result res = transaction.exec("SELECT submitter_id,obs_id FROM cheopsimjobs WHERE id="+m_jobid+";");
		if (res.size() == 1) {
			progId = res[0]["submitter_id"].as<int>();
			reqId = res[0]["obs_id"].as<int>();
		} else {
			throw runtime_error("Error reading program ID and request ID from the database");
		}
		transaction.commit();

	}

	//Check that exposure times are within allowed ranges
	if (exposureTimeAsDouble<0. || exposureTimeAsDouble>600.) {
		throw runtime_error("ERROR in Simulator::initialize: exposure duration must be in range 0 to 600 seconds. Input value is"+to_string(exposureTimeAsDouble));
	}
	if (exposuresPerStack<1 || exposuresPerStack>60) {
		throw runtime_error("ERROR in Simulator::initialize: number of exposures per stack must be in range 0 to 60. Input value is"+to_string(exposuresPerStack));
	}
	if (exposureTimeAsDouble*exposuresPerStack>600) {
		throw runtime_error("ERROR in Simulator::initialize: product of exposure duration and number of exposures per stack cannot exceed 600 s. Input value is"+
							to_string(exposuresPerStack*exposureTimeAsDouble));
	}

	//Update the execution start time in the database
	if (m_database != nullptr) {
		string execution_start_time = boost::posix_time::to_iso_extended_string(boost::posix_time::second_clock::local_time());
		cout << "Job " << m_jobid << " started at " << execution_start_time << endl;
		pqxx::work transaction(*m_database);
		transaction.exec("UPDATE cheopsimjobs SET start_time='"+execution_start_time+"', status='RUN' WHERE id="+m_jobid+";");
		transaction.commit();
	}

	//Set the leap second file and the OBT2UTC correlation file needed to define onboard time
	list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0101.fits"};
	LeapSeconds::setLeapSecondsFileNames(leapSecondFile);
	list<string> correlationFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU2019-06-01T14-33-26_AUX_RES_ObtUtcCorrelation_V0004.fits"};
	OBTUTCCorrelation::SetCorrelationFileNames(correlationFile);
	string ephemerisFileName = string(getenv("CHEOPS_TESTDATA"))+"/resources/CH_TU1949-12-14T00-00-00_EXT_APP_DE1_V0101.fits";
	BarycentricOffset::setEphemerisFileName(ephemerisFileName);

	// Write the list of input reference files to reference_files.txt
	ofstream referenceFilesList("reference_files.txt", ios::out);
	referenceFilesList << "The following is a list of all reference files that are used as input to CHEOPSim and that are used elsewhere in the SOC software (i.e. not exclusive to CHEOPSim)." << endl << endl;
	referenceFilesList << "The reference files are provided in a tar file linked from https://redmine.isdc.unige.ch/projects/cheops/wiki/Wiki#Reference-Data" << endl;
	referenceFilesList << "There are also direct links here: https://docs.google.com/spreadsheets/d/1-ZUPEUMi7gE430UUTbTmx-xHdv1uG3BnUY_8CWzdhoA" << endl << endl;
	referenceFilesList << leapSecondFile.begin()->substr(leapSecondFile.begin()->find_last_of("/") + 1) << endl;
	referenceFilesList << correlationFile.begin()->substr(correlationFile.begin()->find_last_of("/") + 1) << endl;
	referenceFilesList << ephemerisFileName.substr(ephemerisFileName.find_last_of("/") + 1) << endl;
	referenceFilesList.close();

	// Initialize the time configuration
	TimeConfiguration timeConf(startTime,numberOfStackedImages,exposuresPerStack,exposureTimeAsDouble,maxImagesPerCube,m_doFullFrame,m_params->GetAsDouble("repetitionPeriod"));
	// Initialize the data
	m_data = new Data(timeConf,m_params->GetAsString("modulesToRun"));

	// Define the visit
	Data::Visit visit(progType,progId,reqId,obsId,visitCtr,obsCat);
	visit.m_progType = progType;
	visit.m_progId = progId;
	visit.m_reqId = reqId;
	visit.m_obsCat = obsCat;
	visit.m_piName = piName;
	visit.m_piUid = piUid;
	visit.m_prpFirst = prpFirst;
	visit.m_prpLast = prpLast;
	visit.m_targName = targName;
	visit.m_readMode = readMode;
	visit.m_gaiaMagError = gaiaMagError;
	visit.m_strayLightThreshold = strayLightThreshold;
	visit.m_marginMode = marginMode;
	m_data->setVisit(visit);
	if (visitConstraints != nullptr) m_data->setMpsVisitConstraints(visitConstraints);

	cout << endl << "=================== Time configuration ===========================" << endl;
    cout << "Simulation start time                    :" << timeConf.getStartTime() << endl;
    cout << "Simulation duration (minutes)            :" << timeConf.getDuration()/60. << endl;
    cout << "Exposure time (s)                        :" << timeConf.getExposureTimeAsDouble() << endl;
    cout << "Number of time steps (i.e. exposures)    :" << timeConf.getNumberOfTimeSteps() << endl;
    cout << "Number of exposures per stack            :" << timeConf.getExposuresPerStack() << endl;
    cout << "Number of stacked images to be output    :" << timeConf.getNumberOfStackedImages() << endl << endl;

    //Create the output directory
	boost::filesystem::create_directory(m_data->getOutputDirectory());

	//Initialize the field of view by defining the pointing direction
	SkyPosition pointingDirection(pointingRA,pointingDec);
	SkyFieldOfView *fov = new SkyFieldOfView(pointingDirection);
	m_data->setFieldOfView(fov);

	string ra_format = (pointingRA_str.find(":") == string::npos) ? "(degrees)               " : "(hours:minuntes:seconds)";
	string dec_format = (pointingDec_str.find(":") == string::npos) ? "(degrees)                      " : "(degrees:arcminutes:arcseconds)";
	cout << "Pointing direction RA "+ra_format+"                 :" << pointingRA_str << endl;
	cout << "Pointing direction declination "+dec_format+" :" << pointingDec_str << endl;

	if (m_params->GetAsBool("redundantHardware")) m_data->setRedundantHardware();
	m_data->setSerialReadRate(timeConf.getExposureTimeAsDouble() <= 2.326 ? 230. : 100.); //Will be overwritten by DarkCurrentGenerator if run
	cout << endl << "Readout via " << (m_params->GetAsBool("redundantHardware") ? "right" : "left") << " amplifier at " << m_data->getSerialReadRate() << " kHz" << endl << endl;

	//Set the nominal values for each of the bias voltages, extracting them from the REF_APP_GainCorrection file
	string gainCorrectionFilename = m_params->GetAsString("gainCorrectionFilename");
	if (gainCorrectionFilename[gainCorrectionFilename.length()-10] != 'V') {
		throw runtime_error("name of REF_APP_GainCorrection filename must end with Vxxxx.fits");
	}
	if (m_data->redundantHardware()) {
		boost::replace_last(gainCorrectionFilename,"0.fits","2.fits");
		boost::replace_last(gainCorrectionFilename,"2018","2050");
	}
	m_data->setNominalVoltages(gainCorrectionFilename);

	m_data->setPhotometryParams(m_params->GetAsDouble("radius_barycentre"),
								m_params->GetAsDouble("radius_bkgInner"),
								m_params->GetAsDouble("radius_bkgOuter"),
								m_params->GetAsDouble("radius_psf"),
								m_params->GetAsBool("subtractBackground"));

	//Get the wavelength dependence configuration
	if (m_params->GetAsString("modulesToRun") != "DataReduction") { //No need to do this for case of standalone photometric extraction

		bool applyThroughput = m_params->GetAsBool("applyThroughput");
		string throughputFilename = m_params->GetAsString("throughputFilename");
		bool applyQE = m_params->GetAsBool("applyQE");
		string QEFilename = m_params->GetAsString("QEFilename");
		string SEDFilename = m_params->GetAsString("SEDFilename");
		string targetSEDFilename = m_params->GetAsString("targetSEDFilename");

		//Initialize the wavelength dependence
		m_data->setWavelengthDependence(new WavelengthDependence(applyThroughput,throughputFilename,applyQE,QEFilename,SEDFilename,targetSEDFilename));

	}

	//Open the output MPS_PRE_Visits file and set the header keywords
	UTC startTime_utc = timeConf.getVisitStartTimeUTC();
	UTC endTime_utc = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());
	string dirname = m_data->getOutputDirectory()+"/data";
	boost::filesystem::create_directory(dirname);
	MpsPreVisits * mpsVisits_out = new MpsPreVisits(buildFileName(dirname, startTime_utc,
			VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visit.m_visitCtr),
			PassId(), std::string(), MpsPreVisits::getExtName()),"CREATE");
	mpsVisits_out->setKeyProcNum(visit.m_versionNum);
	mpsVisits_out->setKeyArchRev(0);
	mpsVisits_out->setKeyVStrtU(startTime_utc);
	mpsVisits_out->setKeyVStrtM(startTime_utc.getMjd());
	mpsVisits_out->setKeyVStopU(endTime_utc);
	mpsVisits_out->setKeyVStopM(endTime_utc.getMjd());

	//Set the column values for the MPS visit
	string email = "";
	ifstream email_in("email.txt", ios::in);
	if (email_in.good()) {
	    string line;
        getline(email_in, line);
        stringstream(line) >> email;
	}
	email_in.close();
	mpsVisits_out->setCellProgrammeType(visit.m_progType);
	mpsVisits_out->setCellProgrammeId(visit.m_progId);
	mpsVisits_out->setCellRequestId(visit.m_reqId);
	mpsVisits_out->setCellObsid(visit.m_obsId);
	mpsVisits_out->setCellUtcTimeStart(startTime_utc);
	mpsVisits_out->setCellUtcTimeStop(endTime_utc);
	mpsVisits_out->setCellPiName(visit.m_piName);
	mpsVisits_out->setCellPiUid(visit.m_piUid);
	mpsVisits_out->setCellEMail(email);
	mpsVisits_out->setCellPropFirstVisit(visit.m_prpFirst);
	mpsVisits_out->setCellPropLastVisit(visit.m_prpLast);
	mpsVisits_out->setCellTargetName(visit.m_targName);
	mpsVisits_out->setCellObsCategory(visit.m_obsCat);
	mpsVisits_out->setCellReadoutMode(visit.m_readMode);
	mpsVisits_out->setCellExptime(timeConf.getExposureTimeAsDouble());
	mpsVisits_out->setCellNexp(timeConf.getExposuresPerStack());
	mpsVisits_out->setCellNFullframeExp(m_doFullFrame);
	mpsVisits_out->setCellNWindowframeExp(timeConf.getNumberOfStackedImages() * timeConf.getExposuresPerStack());
	mpsVisits_out->setCellRaTarg(m_data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	mpsVisits_out->setCellDecTarg(m_data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
	mpsVisits_out->setCellStrayLightThreshold(visit.m_strayLightThreshold);
	m_data->setMpsPreVisits(mpsVisits_out);

//	UTC utc = UTC(to_iso_extended_string(startTime));
//	MJD mjd = utc.getMjd();
//	cout << "MJD " << setprecision(15) << mjd  << " " << mjd+2400000.5<< endl;
//	BarycentricOffset offset(pointingRA/3600.,pointingDec/3600.);
//	BJD bjd = mjd+offset;
//	cout << "BJD " << setprecision(15) << bjd << " " << bjd+2400000.5 << endl;
//	MJD mjd2 = bjd-offset;
//	cout << "MJD2 " << setprecision(15) << mjd2  << " " << mjd2+2400000.5<< endl;

//	string orbitFile = "CH_PR990001_TG000001_TU2018-08-01T15-30-00_MPS_PRE_Visits_V0000.fits";
//	OrbitInterpolation orbitInterpolation(orbitFile,string(getenv("CHEOPS_SW"))+"/resources/CH_TU2020-03-11T00-00-00_EXT_APP_SAAMap-700km_V0102.fits");
//	OrbitInterpolation::OrbitData orbitData = orbitInterpolation.interpolate(UTC("2018-08-01T15:34:00.000000"));
//	cout << orbitData.m_moonAngle << " " << orbitData.m_sunAngle << " " << orbitData.m_earthLimbAngle << " " << orbitData.m_earthOccultationFlag << " " << orbitData.m_SAAFlag << endl;
//	orbitData = orbitInterpolation.interpolate(UTC("2018-08-01T15:34:10.000000"));
//	cout << orbitData.m_moonAngle << " " << orbitData.m_sunAngle << " " << orbitData.m_earthLimbAngle << " " << orbitData.m_earthOccultationFlag << " " << orbitData.m_SAAFlag << endl;

}

void Simulator::connectToDatabase(string databaseServer) {

	//get the job ID from the directory name
	string jobdir = string(boost::filesystem::basename(getenv("PWD")));
	string jobdir_prefix = "CHEOPSim_job";
	size_t found = jobdir.find(jobdir_prefix);
	if (found!=string::npos) {
		m_jobid = jobdir.substr(found+jobdir_prefix.size());

		//check that the extracted job ID is an integer
		bool jobid_isValid = true;
		for (unsigned i=0; i<m_jobid.length(); i++) {
			if (!isdigit(m_jobid[i])) jobid_isValid = false;
		}

		//Open the database connection
		if (jobid_isValid) {
			try {
				if (databaseServer.substr(0,1)=="p") {
					m_database = new pqxx::connection("host="+databaseServer+".isdc.unige.ch dbname=p_cheopsim user=p_cheopsim password=ch30ps");
				} else if (databaseServer.substr(0,1)=="d") {
					m_database = new pqxx::connection("host="+databaseServer+".isdc.unige.ch dbname=d_cheopsim user=d_cheopsim password=ch30ps");
				} else {
					m_database = new pqxx::connection("host=localhost dbname=cheopsim user=futyan password=ch30ps");
				}
			} catch(std::exception& e) {
				throw runtime_error(e.what());
			}
		}
	}

}

void Simulator::process() {

	//initialize all modules added to the simulator
	if (!m_modulesInitialized) initializeModules();

	//Generate bias voltages and CCD temperature with drift and random fluctuations at 1.2s cadence, to be used for image metadata and HK data
	generateVoltagesAndTemperatures();

	//Loop over modules which should run before the time loop
	processBegin();
	//Loop over time, looping over modules which should run for each time step
	processTimeLoop();
	//Loop over modules which should run after the time loop
	processEnd();

}

void Simulator::initializeModules() {

    cout << endl << "=================== Module Initialization ========================" << endl << endl;

	// Instantiate required modules via ModuleRegistry and add them to m_modules
    ModuleRegistry moduleRegistry;
    moduleRegistry.instanciateModules(m_params);
    m_modules = moduleRegistry.getModules();

    // Call the Module::initialize() method for all instantiated modules
    for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {

    	string moduleName = (*module)->getName();
		cout << endl << "Initializing module " << moduleName << endl;
		const ModuleParams & moduleParams = m_params->Module(moduleName);
    	(*module)->initialize(moduleParams);
    	if (moduleName=="FrameTransferSmearer") m_data->setFrameTransferSmearing();
    	if (moduleName=="ChargeTransferSimulator") {
    		if (moduleParams.GetAsBool("endOfLife")) m_data->setChargeTransferEoL();
    	}
    	if (moduleName=="ImageWriter") {
    		if (moduleParams.GetAsBool("omitSAA")) m_data->setOmitSAA();
    		if (moduleParams.GetAsBool("omitEarthOccultation")) m_data->setOmitEarthOccultation();
    		if (moduleParams.GetAsBool("omitStrayLight")) m_data->setOmitStrayLight();
    		m_data->setNumberOfCEsToOmitAtStart(moduleParams.GetAsInt("numberOfCEsToOmitAtStart"));
    		m_data->setOnlyWriteEveryNthCE(moduleParams.GetAsInt("onlyWriteEveryNthCE"));
    	}

    	//Set the margin mode for the visit if ImageWriter is to be run (unless already set via MPS_PRE_Visits)
    	if (moduleName=="ImageWriter") {
    		Data::Visit visit = m_data->getVisit();
    		if (visit.m_marginMode=="undefined") {
    			visit.m_marginMode = moduleParams.GetAsString("marginMode");
    			m_data->setVisit(visit);
    		}
		}

    }
    m_modulesInitialized = true;

}


void Simulator::generateVoltagesAndTemperatures() {

	//Set up the Gaussian random number generator
	boost::mt19937 randomNumberEngine(12345);
    boost::normal_distribution<double> normalDistribution(SatelliteData::kDefaultCcdTemperature,SatelliteData::kCcdTemperatureSigma[SatelliteData::main]/1.E3);
    RANDOM_GAUSS * gaussianNoiseGenerator = new RANDOM_GAUSS(randomNumberEngine,normalDistribution);

	SatelliteData::READOUT_HARDWARE readout_hardware = m_data->redundantHardware() ? SatelliteData::redundant :SatelliteData::main;

	//Get the total duration of the simulation in seconds
	TimeConfiguration timeConf = m_data->getTimeConfiguration();
	UTC startTime_utc = timeConf.getVisitStartTimeUTC();
	UTC endTime_utc = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());
	double totalDuration = (endTime_utc-startTime_utc).getSeconds();

	//Loop over the total duration (plus kAveragedHKCadence because Data::calculateHKAverages calculates the average for the previous kAveragedHKCadence seconds,
	//plus a further kAveragedHKCadence because in HkWriter there is a point after the end of the visit to ensure the visit is covered and the SDS counter includes all omitted CEs,
	//plus Data::kRawHKCadence to avoid a problem for certain exposure durations)
	//in steps of 1.2s (the measurement cadence on board for voltages and temperatures)
	for (double timeSinceVisitStart = 0; timeSinceVisitStart <= totalDuration+2*Data::kAveragedHKCadence+Data::kRawHKCadence; timeSinceVisitStart += Data::kRawHKCadence) {

		gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(SatelliteData::kDefaultCcdTemperature,SatelliteData::kCcdTemperatureSigma[readout_hardware]/1.E3));
		double ccdTemp = (*gaussianNoiseGenerator)();
		//cout << "fluctuated CCD temperature: " << ccdTemp << " " << SatelliteData::kDefaultCcdTemperature << " " << SatelliteData::kCcdTemperatureSigma[readout_hardware] << " ";

		gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(SatelliteData::kDefaultFeeBiasTemperature,SatelliteData::kFeeBiasTemperatureSigma[readout_hardware]/1.E3));
		double biasTemp = (*gaussianNoiseGenerator)();

		gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(SatelliteData::kDefaultFeeAdcTemperature,SatelliteData::kFeeAdcTemperatureSigma[readout_hardware]/1.E3));
		double adcTemp = (*gaussianNoiseGenerator)();

		double voltages[4];
		gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(m_data->getVoltageWithDrift(timeSinceVisitStart,SatelliteData::VSS),SatelliteData::kVoltageSigma[readout_hardware][SatelliteData::VSS]/1.E6));
		voltages[SatelliteData::VSS] = (*gaussianNoiseGenerator)();

		for (unsigned i=0; i<3; i++) {
			SatelliteData::VOLTAGE_TYPE type = static_cast<SatelliteData::VOLTAGE_TYPE>(i);
			//cout << "fluctuated voltage[" << type << "]: " << " " << m_data->getVoltageWithDrift(timeSinceVisitStart,type) << " " << SatelliteData::kVoltageSigma[readout_hardware][type] << " ";
			//rms from payload calibration document is for Voltage-VSS, so use voltage-vss for the gaussian mean
			gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(m_data->getVoltageWithDrift(timeSinceVisitStart,type)-voltages[SatelliteData::VSS],SatelliteData::kVoltageSigma[readout_hardware][type]/1.E6));
			//add vss subtracted above to assign the absolute voltage
			voltages[type] = voltages[SatelliteData::VSS] + (*gaussianNoiseGenerator)();
		}

		//cout << setprecision(10) << timeSinceVisitStart << " " << ccdTemp << " " << biasTemp << " " << adcTemp << " " << voltages[SatelliteData::VOD] << " " << voltages[SatelliteData::VRD] << " " << voltages[SatelliteData::VOG] << " " << voltages[SatelliteData::VSS] << endl;
		m_data->appendRawHKData(Data::HKData(ccdTemp,biasTemp,adcTemp,voltages));

	}

	//Calculate the average HK data at 20s intervals (the cadence for HK data sent to ground), taking the mean over the previous 16 raw measurements (HK cadence / onboard cadence = 20/1.2 = 16.67)
	m_data->calculateHKAverages();

}

void Simulator::processBegin() {

	if (!m_modulesInitialized) initializeModules();

    cout << endl << "=================== Module pre-processing ========================" << endl << endl;

    Module * strayLightModule = nullptr;
    for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
    	if ((*module)->getName() =="StrayLightGenerator") strayLightModule = *module;
    }

	//Loop over modules. doBegin must be called for StrayLightGenerator before OrbitSimulator
    //(in order to set the stray light flux in MPS_PRE_VisitConstraints) and before FocalPlaneGenerator
    //(in order to define the number of images to omit due to SL before defining the image cube size)
    for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
    	if (strayLightModule!=nullptr) {
    		if ((*module)->getName() == "OrbitSimulator" || ((*module)->getName() == "FocalPlaneGenerator" && !m_data->moduleIsActive("OrbitSimulator"))) {
    			cout << "Processing doBegin for module StrayLightGenerator" << endl;
    			strayLightModule->doBegin(m_data);
    			cout << endl;
    		}
    	}
    	if ((*module)->getName() != "StrayLightGenerator") {
    		cout << "Processing doBegin for module " << (*module)->getName() << endl;
    		(*module)->doBegin(m_data);
    		cout << endl;
    	}
    }

}

void Simulator::processTimeLoop() {

	if (!m_modulesInitialized) initializeModules();

	//Check that there is at least one timeLoop module
	bool modulesInTimeLoop = false;
	for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
		if ((*module)->schedule()==Module::timeLoop) {
			modulesInTimeLoop = true;
			break;
		}
	}
	if (!modulesInTimeLoop) return;

	if (m_doFullFrame) {
		//Generate full frame image for first exposure
	    cout << "=================== Generating full frame image for first exposure ========================" << endl;
		for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
			if ((*module)->getName()=="StarProducer") {
				cout << endl << "Processing StarProducer for full frame image" << endl;
				(*module)->doBegin(m_data,true);
			} else if ((*module)->getName()=="BiasGenerator") {
				cout << endl << "Re-initializing bias frame and CCD non-linearity for full frame readout frequency" << endl;
				(*module)->doBegin(m_data,true);
			}
			cout << endl;
			if ((*module)->schedule()==Module::timeLoop && (*module)->getName()!="IdealLightCurveGenerator") {
				cout << "Processing timestep module " << (*module)->getName() << " for full frame image" << endl;
				int fullFrameTimestep = m_data->getTimeConfiguration().getFullFrameTimeStep();
				if ((*module)->getName().find("TransitFluxModulator_star") != string::npos ||
						(*module)->getName()=="StellarNoiseFluxModulator_star0" ||
						(*module)->getName()=="StellarVariationFluxModulator_star0" ||
						(*module)->getName()=="UserFluxModifier_star0" ||
						(*module)->getName()=="StrayLightGenerator" ||
						(*module)->getName()=="HaloGenerator" ||
						(*module)->getName()=="PSFGenerator" ||
						(!m_data->moduleIsActive("JitterProducer")&&!m_data->moduleIsActive("OrbitSimulator"))) {
					fullFrameTimestep = 0;
				}
				(*module)->process(m_data,fullFrameTimestep,true);
			}
		}
		cout << endl << "======================== Full frame image generation complete =============================" << endl;

		//re-run doBegin for StarProducer to generate star list for sub-array only
		//re-run doBegin for BiasGenerator to read non-linearity and bias frame for potentially different readout frequency
		for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
			if ((*module)->getName()=="StarProducer") {
				cout << endl << "Processing StarProducer for sub-array images" << endl;
				(*module)->doBegin(m_data);
			} else if ((*module)->getName()=="BiasGenerator")  {
				cout << endl << "Re-initializing bias frame and CCD non-linearity for sub-array readout frequency" << endl;
				(*module)->doBegin(m_data);
			}
		}
	}

	//Loop over time steps
    cout << endl << "================================ Processing time loop =====================================" << endl << endl;
	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
		cout << "exposure " << iTimeStep+1 << " / " << nTimesteps << endl;

		//update the progress in the database
		if (m_database != nullptr) {
			double fraction_done = (boost::lexical_cast<double>(iTimeStep)/boost::lexical_cast<double>(nTimesteps));
			string percent_done = boost::lexical_cast<std::string>(floor(fraction_done*1000.)/10.);
			try {
				pqxx::work transaction(*m_database);
				transaction.exec("UPDATE cheopsimjobs SET percent_done = '" + percent_done + "' WHERE id = " + m_jobid + ";");
				transaction.commit();
			} catch(std::exception& e) {
				cerr << "unable to connect to database" << endl;
			}
		}

		//Loop over modules
		//clock_t startTime = clock();
		for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
	    	//Check that this module should run within the time loop
			if ((*module)->schedule()==Module::timeLoop) {
				//Run the module to process the data
				//cout << "Processing timestep module " << (*module)->getName() << endl;
				//clock_t moduleStartTime = clock();
    			(*module)->process(m_data,iTimeStep);
    			//cout << (*module)->getName() << ": " << double( clock() - moduleStartTime ) / (double)CLOCKS_PER_SEC << " seconds" << endl;
    		}
    	}
		//cout << "Total processing time: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds" << endl;

    }

}

void Simulator::processEnd() {


	if (!m_modulesInitialized) initializeModules();

	//Loop over modules
    for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) {
    	if ((*module)->getName()=="DataReduction") {
    	    cout << endl << "=================================== Data Reduction ========================================" << endl << endl;
    	}
		(*module)->doEnd(m_data);
    }

}

Simulator::~Simulator() {

	//Set MPS_PRE_Visits column values which could not be set in the constructor since relevant modules had not run, and write the row
	MpsPreVisits * mpsVisits_out = m_data->getMpsPreVisits();
	if (m_data->getFieldOfView()->getStars().size()>0) {
		Star * targetStar = m_data->getFieldOfView()->getStars()[0];
		mpsVisits_out->setCellSpectralType(targetStar->getSpectralTypeString());
		mpsVisits_out->setCellMagTarg(targetStar->getGaiaMagnitude());
		mpsVisits_out->setCellTransitTime(targetStar->getTransitTime());
		double transitPeriod = targetStar->getTransitPeriod();
		if (transitPeriod>0.) {
			mpsVisits_out->setCellTransitPeriod(transitPeriod);
		} else {
			mpsVisits_out->setNullTransitPeriod();
		}
	} else {
		mpsVisits_out->setNullSpectralType();
		mpsVisits_out->setNullMagTarg();
		mpsVisits_out->setNullTransitTime();
		mpsVisits_out->setNullTransitPeriod();
	}
	mpsVisits_out->setCellTargetLocationX(m_data->getSubarrayDimensions().m_targetLocationX);
	mpsVisits_out->setCellTargetLocationY(m_data->getSubarrayDimensions().m_targetLocationY);
	mpsVisits_out->setCellMarginMode(m_data->getVisit().m_marginMode);
	mpsVisits_out->WriteRow();

	//Update the database now that the job has completed
	delete m_data;
	if (m_database != nullptr) {
		string execution_end_time = boost::posix_time::to_iso_extended_string(boost::posix_time::second_clock::local_time());
		cout << "Job " << m_jobid << " completed at " << execution_end_time << endl;
		pqxx::work transaction(*m_database);
		transaction.exec("UPDATE cheopsimjobs SET percent_done='100', end_time='"+execution_end_time+"' WHERE id="+m_jobid+";");
		transaction.commit();
		m_database->disconnect();
	}
	delete m_database;

	//Clean up
    for (vector<Module*>::const_iterator module = m_modules.begin(); module!=m_modules.end(); ++module) delete (*module);
	boost::filesystem::remove_all(string(getenv("PWD"))+"/lc*");
}
