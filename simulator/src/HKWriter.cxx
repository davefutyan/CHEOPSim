/*
 * HKWriter.cxx
 *
 *  Created on: 6 Sep 2017
 *      Author: futyand
 */

#include <fstream>

#include "CreateFitsFile.hxx"
#include "Hk2RawProcessing.hxx"
#include "SCI_RAW_HkCe.hxx"

#include "HKWriter.hxx"

void HKWriter::doEnd(Data* data) const {

	//Create the output directory
	string dirname = data->getOutputDirectory()+"/hk";
	boost::filesystem::create_directory(dirname);

	//Instantiate the class to generate the HK data
	string hkEnumConversionFilename = "CH_TU2015-01-01T00-00-00_REF_APP_HkEnumConversion_V0103.fits";
	Data::Visit visit = data->getVisit();
	Hk2RawProcessing hk2RawProcessing(
			string(getenv("CHEOPS_SW"))+"/resources/"+hkEnumConversionFilename,
			dirname, VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visit.m_visitCtr));

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << hkEnumConversionFilename << endl;
	referenceFilesList.close();

	//Write out the housekeeping data with cadence kAveragedHKCadence, with the
	//timings for each type of HK data offset w.r.t. one another
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double startTimeOffset = (timeConf.getStartTimeUTC() - timeConf.getVisitStartTimeUTC()).getSeconds();
	for (double time = -startTimeOffset; time <= timeConf.getDuration()+Data::kAveragedHKCadence; time+=Data::kAveragedHKCadence) {

        int timeStep = static_cast<int>(floor(time/timeConf.getRepetitionPeriod()));
        if (timeStep == int(timeConf.getNumberOfTimeSteps())) timeStep = timeConf.getNumberOfTimeSteps()-1;
        if (timeStep<-1*int(timeConf.getNumberOfPreSubArrayTimeSteps())) timeStep=-1*int(timeConf.getNumberOfPreSubArrayTimeSteps());
        UTC currentTime = timeConf.getStartTimeUTC() + DeltaTime(time);
        //cout << currentTime.getUtc() << " " << timeStep << " " << timeConf.getNumberOfTimeSteps() << " " << time << " " << timeConf.getDuration() << endl;

		unsigned numberOfDiscardedCEs = m_sdsOffset;
		for (int i=0; i<=timeStep;  i++) {
			if ((i+1)%timeConf.getExposuresPerStack() == 0 && i<int(timeConf.getNumberOfTimeSteps())) { //i corresponds to the last exposure in a stack
				if (data->hasSatelliteData()) {
					if (data->getSatelliteData(i)->getDiscardFlag()) numberOfDiscardedCEs++;
				}
			}
		}

		double telescopeTemperature = SatelliteData::kDefaultTelescopeTemperature;
		double dpuTemperature = SatelliteData::kDefaultDpuTemperature;
		double dpuVoltage = SatelliteData::kDefaultDpuVoltage;
		if (data->hasSatelliteData() && data->moduleIsActive("OrbitSimulator") && timeStep<int(timeConf.getNumberOfTimeSteps())) {
			SatelliteData * satData = data->getSatelliteData(timeStep);
			dpuTemperature = satData->getDpuTemperature();
			telescopeTemperature = satData->getTelescopeTemperature();
		}

		Data::HKData hkData = data->getClosestAveragedHKData(time+startTimeOffset);

		hk2RawProcessing.writeRow(currentTime, "SCI_RAW_HkDefault");
		hk2RawProcessing.updateParameter("TEMP_FEE_CCD",hkData.m_ccdTemp-273.15);
		hk2RawProcessing.updateParameter("TEMP_FEE_ADC",hkData.m_adcTemp-273.15);
		hk2RawProcessing.updateParameter("TEMP_FEE_BIAS",hkData.m_biasTemp-273.15);
		hk2RawProcessing.updateParameter("VOLT_FEE_VOD",hkData.m_vod);
		hk2RawProcessing.updateParameter("VOLT_FEE_VRD",hkData.m_vrd);
		hk2RawProcessing.updateParameter("VOLT_FEE_VOG",hkData.m_vog);
		hk2RawProcessing.updateParameter("VOLT_FEE_VSS",hkData.m_vss);
		hk2RawProcessing.updateParameter("ADC_TEMP1_U",dpuTemperature-273.15);
		hk2RawProcessing.updateParameter("ADC_N5V_L",dpuVoltage);
		hk2RawProcessing.updateParameter("sdsCounter",numberOfDiscardedCEs);
		hk2RawProcessing.updateParameter("observationId",visit.m_obsId);
		hk2RawProcessing.writeRow(currentTime, "SCI_RAW_HkExtended");
		hk2RawProcessing.updateParameter("ttc1AvTempAft",telescopeTemperature-273.15+0.1+0.5);
		hk2RawProcessing.updateParameter("ttc1AvTempFrt",telescopeTemperature-273.15+0.5);
		if (time+2 <= timeConf.getDuration()) hk2RawProcessing.writeRow(currentTime+DeltaTime(2), "SCI_RAW_HkIaswDg");
		//correspondence with SCI_RAW_ImageMetadata thermfornt and thermaft values according to https://redmine.astro.unige.ch/issues/17499
		hk2RawProcessing.updateParameter("ADC_TEMPOH4B",telescopeTemperature-273.15+0.4+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH3B",telescopeTemperature-273.15+0.3+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH2B",telescopeTemperature-273.15+0.2+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH1B",telescopeTemperature-273.15+0.1+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH1A",telescopeTemperature-273.15+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH2A",telescopeTemperature-273.15-0.1+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH3A",telescopeTemperature-273.15-0.2+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMPOH4A",telescopeTemperature-273.15-0.3+0.5);
		hk2RawProcessing.updateParameter("ADC_TEMP1",dpuTemperature-273.15);
		hk2RawProcessing.updateParameter("ADC_N5V",dpuVoltage);
		if (time-2 >= -startTimeOffset) hk2RawProcessing.writeRow(currentTime+DeltaTime(-2), "SCI_RAW_HkIfsw");
		if (time-1 >= -startTimeOffset) hk2RawProcessing.writeRow(currentTime+DeltaTime(-1), "SCI_RAW_HkIaswPar");
		if (time+3 <= timeConf.getDuration()) hk2RawProcessing.writeRow(currentTime+DeltaTime(3), "SCI_RAW_HkIbswDg");
		if (time+1 <= timeConf.getDuration()) hk2RawProcessing.writeRow(currentTime+DeltaTime(1), "SCI_RAW_HkIbswPar");

	}

	//Set the HK header keywords
	string specType = "undefined";
	double Gmag = 0.;
	double cheopsMag = 0.;
	double GmagErr = 0.;
	double cheopsMagErr = 0.;
	double TEff = 0;
	if (data->getFieldOfView()->getStars().size()>0) {
		Star * targetStar = data->getFieldOfView()->getStars()[0];
		specType = targetStar->getSpectralTypeString();
		TEff = targetStar->getEffectiveTemperature();
		Gmag = targetStar->getGaiaMagnitude();
		cheopsMag = targetStar->getCheopsMagnitude();
		GmagErr = targetStar->getGaiaMagnitudeError();
		cheopsMagErr = targetStar->getCheopsMagnitudeError();
	}
	hk2RawProcessing.setKeyProcChn("CHEOPSim");
	hk2RawProcessing.setKeyCheopsid(0);
	hk2RawProcessing.setKeyTargname(visit.m_targName);
	hk2RawProcessing.setKeySpectype(specType);
	//hk2RawProcessing.setKeyTEff(TEff); //Uncomment after common_sw 3.1.5 or later has been released
	hk2RawProcessing.setKeyMagG(Gmag);
	hk2RawProcessing.setKeyMagChps(cheopsMag);
	hk2RawProcessing.setKeyMagGerr(GmagErr);
	hk2RawProcessing.setKeyMagCerr(cheopsMagErr);
	hk2RawProcessing.setKeyPiName(visit.m_piName);
	hk2RawProcessing.setKeyPiUid(visit.m_piUid);
	hk2RawProcessing.setKeyObsCat(visit.m_obsCat);
	hk2RawProcessing.setKeyObsid(visit.m_obsId);
	hk2RawProcessing.setKeyPrpVst1(visit.m_prpFirst);
	hk2RawProcessing.setKeyPrpVstn(visit.m_prpLast);
	hk2RawProcessing.setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	hk2RawProcessing.setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);

	//Get pointers to the stacked and unstacked image metadata
	SciRawImagemetadata * stackedMetaData = data->getFitsImageMetaData(Data::stacked);
	SciRawImagemetadata * unstackedMetaData = data->getFitsImageMetaData(Data::unstacked);

	if (unstackedMetaData != nullptr && stackedMetaData != nullptr) {

		//Open the SCI_RAW_HkCe output file and define the header keywords
		UTC startTime = stackedMetaData->getKeyVStrtU();
		UTC endTime = stackedMetaData->getKeyVStopU();
		SciRawHkce * sciRawHkce = createFitsFile<SciRawHkce>(dirname, startTime,
				VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visit.m_visitCtr), PassId());
		sciRawHkce->setKeyProcChn("CHEOPSim");
		sciRawHkce->setKeyProcNum(visit.m_versionNum);
		sciRawHkce->setKeyArchRev(0);
		sciRawHkce->setKeyVStrtU(startTime);
		sciRawHkce->setKeyVStrtM(startTime.getMjd());
		sciRawHkce->setKeyVStopU(endTime);
		sciRawHkce->setKeyVStopM((endTime).getMjd());
		sciRawHkce->setKeyTargname(visit.m_targName);
		sciRawHkce->setKeySpectype(specType);
		sciRawHkce->setKeyTEff(TEff);
		sciRawHkce->setKeyMagG(Gmag);
		sciRawHkce->setKeyMagChps(cheopsMag);
		sciRawHkce->setKeyMagGerr(GmagErr);
		sciRawHkce->setKeyMagCerr(cheopsMagErr);
		sciRawHkce->setKeyPiName(visit.m_piName);
		sciRawHkce->setKeyPiUid(visit.m_piUid);
		sciRawHkce->setKeyProgtype(visit.m_progType);
		sciRawHkce->setKeyProgId(visit.m_progId);
		sciRawHkce->setKeyReqId(visit.m_reqId);
		sciRawHkce->setKeyVisitctr(visit.m_visitCtr);
		sciRawHkce->setKeyObsCat(visit.m_obsCat);
		sciRawHkce->setKeyObsid(visit.m_obsId);
		sciRawHkce->setKeyPrpVst1(visit.m_prpFirst);
		sciRawHkce->setKeyPrpVstn(visit.m_prpLast);
		sciRawHkce->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
		sciRawHkce->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);

		//Assign the SCI_RAW_HkCe data by taking the average over each stacked image of the temperatures and voltages assigned to the unstacked metadata (closest raw HK values)
		unsigned nexp = timeConf.getExposuresPerStack();
		unsigned istack = 0;
		double ccdTemp = 0.;
		double biasTemp = 0.;
		double adcTemp = 0.;
		double dpuTemp = 0.;
		double telescopeTemp = 0.;
		double vod = 0.;
		double vrd = 0.;
		double vog = 0.;
		double vss = 0.;
		while (unstackedMetaData->ReadRow()) {
			istack++;
			ccdTemp += unstackedMetaData->getCellHkTempFeeCcd();
			biasTemp += unstackedMetaData->getCellHkTempFeeBias();
			adcTemp += unstackedMetaData->getCellHkTempFeeAdc();
			dpuTemp += unstackedMetaData->getCellAdcTemp1();
			telescopeTemp += unstackedMetaData->getCellThermfront1();
			vod += unstackedMetaData->getCellHkVoltFeeVod();
			vrd += unstackedMetaData->getCellHkVoltFeeVrd();
			vog += unstackedMetaData->getCellHkVoltFeeVog();
			vss += unstackedMetaData->getCellHkVoltFeeVss();
			if (istack%nexp == 0) {
				stackedMetaData->ReadRow();
				sciRawHkce->setCellObtTime(stackedMetaData->getCellObtTime());
				sciRawHkce->setCellUtcTime(stackedMetaData->getCellUtcTime());
				sciRawHkce->setCellMjdTime(stackedMetaData->getCellMjdTime());
				sciRawHkce->setCellCeTempFeeCcd(ccdTemp/nexp);
				sciRawHkce->setCellCeTempFeeBias(biasTemp/nexp);
				sciRawHkce->setCellCeTempFeeAdc(adcTemp/nexp);
				sciRawHkce->setCellCeAdcTemp1(dpuTemp/nexp);
				sciRawHkce->setCellCeAdcN5v(SatelliteData::kDefaultDpuVoltage);
				sciRawHkce->setCellCeThermaft4(telescopeTemp/nexp+0.4);
				sciRawHkce->setCellCeThermaft3(telescopeTemp/nexp+0.3);
				sciRawHkce->setCellCeThermaft2(telescopeTemp/nexp+0.2);
				sciRawHkce->setCellCeThermaft1(telescopeTemp/nexp+0.1);
				sciRawHkce->setCellCeThermfront1(telescopeTemp/nexp);
				sciRawHkce->setCellCeThermfront2(telescopeTemp/nexp-0.1);
				sciRawHkce->setCellCeThermfront3(telescopeTemp/nexp-0.2);
				sciRawHkce->setCellCeThermfront4(telescopeTemp/nexp-0.3);
				sciRawHkce->setCellCeVoltFeeVod(vod/nexp);
				sciRawHkce->setCellCeVoltFeeVrd(vrd/nexp);
				sciRawHkce->setCellCeVoltFeeVog(vog/nexp);
				sciRawHkce->setCellCeVoltFeeVss(vss/nexp);
				sciRawHkce->WriteRow();
				istack = 0;
				ccdTemp = 0.;
				biasTemp = 0.;
				adcTemp = 0.;
				dpuTemp = 0.;
				telescopeTemp = 0.;
				vod = 0.;
				vrd = 0.;
				vog = 0.;
				vss = 0.;
			}
		}

		delete sciRawHkce;

	}

}
