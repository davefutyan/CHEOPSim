/*
 * HaloGenerator.cxx
 *
 *  Created on: 14 Jul 2016
 *      Author: futyand
 */

#include "HaloGenerator.hxx"
#include "PSFGenerator.hxx"
#include "detector/include/CosmicRayGenerator.hxx"

void HaloGenerator::initialize(const ModuleParams& params) {
	//m_oversampleJitter = params.GetAsBool("oversampleJitter"); //For now not configurable via web interface: jitter oversampling is off
	m_doGhosts = params.GetAsBool("doGhosts");
}

void HaloGenerator::doBegin(Data* data, bool fullFrame) {
}

void HaloGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	unsigned nJitterPerExposure = data->getTimeConfiguration().getJitterSamplingsPerExposure();

	Image * image = data->getImages().back();

	//Get the image dimensions and offsets w.r.t. CCD edges
	int xOffset = image->getXOffset();
	int yOffset = image->getYOffset();
	int xDim = image->getXDim();
	int yDim = image->getYDim();
	//If end of life charge transfer is to be simulated, include also stars below the sub-array
	if (data->doChargeTransferEoL()) {
		yDim += yOffset;
		yOffset = 0.;
	}

	//Loop over all stars in the field of view
	for (unsigned istar=0; istar<data->getFieldOfView()->getStars().size(); istar++) {

		if (fullFrame) cout << "processing star " << istar << endl;

		Star * star = data->getFieldOfView()->getStars()[istar];
		StarData * starTimeData = star->getTimeSeries()[timeStep];

		//Get the the number of PSF positions over the exposure (due to FOV rotation and/or jitter)
		unsigned nJitter = star->getFocalPlaneX().size();
		if (nJitter != nJitterPerExposure && nJitter != 1) {
			throw runtime_error("Error in HaloGenerator::process: FOV rotation vector size must equal exposure time or 1");
		}

		//Calculate the mean position of the star over the exposure, taking into account jitter
		double starPositionX = 0.;
		double starPositionY = 0.;
		for (unsigned ijitter=0; ijitter<nJitter; ijitter++) {
			starPositionX += star->getFocalPlaneX()[ijitter];
			starPositionY += star->getFocalPlaneY()[ijitter];
		}
		starPositionX /= nJitter;
		starPositionY /= nJitter;

		//Require that the star position is within the sub-array plus a margin
		if (starPositionX<xOffset-PSFGenerator::kPsfMargin || starPositionX>=xOffset+xDim+PSFGenerator::kPsfMargin ||
			starPositionY<yOffset-PSFGenerator::kPsfMargin || starPositionY>=yOffset+yDim+PSFGenerator::kPsfMargin) {
			continue;
		}

		//Obtain the flux of the star in photons/cm2
		double flux = star->getMeanPhotonFlux();
		double fluxFactor = starTimeData->getCombinedFluxFactor();
		flux = flux * fluxFactor * data->getTimeConfiguration().getExposureTimeAsDouble();

		//Width of one pixel in mm
		double pixelWidth = CosmicRayGenerator::kPixelWidth/1000.;

		//clock_t startTime = clock();
		//Loop over jitter within the exposure and increment the halo for each jitter position
		double integratedScatter = 0.;
		nJitter = (m_oversampleJitter ? nJitter : 1);
		for (unsigned ijitter=0; ijitter<nJitter; ijitter++) {

			double jitteredPositionX = starPositionX;
			double jitteredPositionY = starPositionY;
			if (m_oversampleJitter) {
				//Get the PSF centre position taking into account jitter and FOV rotation
				//(Note: it is already required that the mean position is within the sub-array plus margin)
				jitteredPositionX = star->getFocalPlaneX()[ijitter];
				jitteredPositionY = star->getFocalPlaneY()[ijitter];
			}

			//Increment the halo flux for each pixel
			for (int ix=xOffset; ix<xOffset+xDim; ix++) {
				for (int iy=yOffset; iy<yOffset+yDim; iy++) {
					double deltaX = static_cast<double>(ix)+0.5 - jitteredPositionX; //add 0.5 to shift to pixel centre
					double deltaY = static_cast<double>(iy)+0.5 - jitteredPositionY; //add 0.5 to shift to pixel centre
					double radius = pixelWidth*sqrt(deltaX*deltaX + deltaY*deltaY);

					//Simplified Peterson model taken from CHEOPS-INAF-INST-RP-001 Issue 3 Revision 6 page 57
					//Multiplying factor is 395 instead of 315 because the latter value assumes
					//throughput*QE already applied with value 0.85x(1-0.008)^8
					//This gives the scatter flux for an incident flux of 1 photon/mm2
					double simplifiedPeterson = 395.*pow(1+19*radius*radius,-1.135);

					//factor 100 to convert flux from photons/cm2 to photons/mm2
					double pixel_value = pixelWidth*pixelWidth * flux/100. * simplifiedPeterson;

					//Increment the pixel value for the current jitter
					image->incrementPixelValue(ix,iy,pixel_value/double(nJitter));
					integratedScatter += pixel_value/double(nJitter);
				}
			}

		}
		//cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds" << endl;

		//double peakScatter = pixelWidth*pixelWidth * flux/100. * 395.*pow(1,-1.135);
		//cout << istar << " " << integratedScatter << " " << peakScatter << endl;

		star->setHaloFlux(integratedScatter);

		if (m_doGhosts) generateGhostFlux(image,flux,starPositionX,starPositionY,data->doChargeTransferEoL());

	}

}

void HaloGenerator::generateGhostFlux(Image* image, double incidentFlux, double starPositionX, double starPositionY, bool doChargeTransferEoL) const {

	//Get the image dimensions and offsets w.r.t. CCD edges
	int xOffset = image->getXOffset();
	int yOffset = image->getYOffset();
	int xDim = image->getXDim();
	int yDim = image->getYDim();
	//If end of life charge transfer is to be simulated, include also stars below the sub-array
	if (doChargeTransferEoL) {
		yDim += yOffset;
		yOffset = 0.;
	}

	//The total ghost flux as a fraction (in %) of the total PSF flux is defined as a function of the
	//angle w.r.t. the pointing direction in Figure 36 of CHEOPS-INAF-INST-RP-001 Issue 3 Revision 6
	double fieldAngles[10] = {0.,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18};
	double ghostFluxFractions[10] = {0.0318818,0.0377916,0.04934603,0.05605599,0.05867185,0.05479005,0.04593002,0.02427683,0.00726283,0.00719129};

	double distanceFromFieldCentre = sqrt(pow(starPositionX-512,2) + pow(starPositionY-512,2));

	double fieldAngle = distanceFromFieldCentre/3600.; //Convert from distance in pixels to degrees, using 1 arcsec/pixel plate scale

	//Identify the first point on the ghostFluxFraction vs field angle curve with field angle greater than the current field angle
	unsigned index = 0;
	double mu = 0.;
	for (unsigned i=1; i<10; i++) {
		if (fieldAngle<fieldAngles[i]) {
			index = i;
			mu = (fieldAngle-fieldAngles[i-1])/(fieldAngles[i]-fieldAngles[i-1]);
			break;
		}
	}

	//Determine the ghost flux fraction for the current field angle by linear interpolation between the two nearest
	//points on the ghost flux fraction vs field angle curve
	double ghostFluxFraction;
	if (index != 0) {
		ghostFluxFraction = ghostFluxFractions[index-1]*(1.-mu)+ghostFluxFractions[index]*mu;
	} else { //Field angle is greater than maximum value in fieldAngles array (0.18)
		ghostFluxFraction =  ghostFluxFractions[9];
	}
	ghostFluxFraction /= 100.; //Convert from percent

	//Total ghost flux is the total incident flux multiplied by the ghost flux fraction (this is how ghost flux fraction is defined)
	double ghostFlux = incidentFlux * PSFGenerator::kTelescopeArea * ghostFluxFraction;

	//Spread the ghost flux uniformly over the image
	double ghostFluxPerPixel = ghostFlux/(Image::kXDim*Image::kYDim);
	for (int ix=xOffset; ix<xOffset+xDim; ix++) {
		for (int iy=yOffset; iy<yOffset+yDim; iy++) {
			image->incrementPixelValue(ix,iy,ghostFluxPerPixel);
		}
	}

	//cout << starPositionX << " " << starPositionY << " " << distanceFromFieldCentre << " " << fieldAngle << " " << index << " " << mu << " " <<  ghostFluxFraction << endl;
	//cout << ghostFlux << " " << ghostFluxPerPixel << " " << ghostFlux/(PSFGenerator::kTelescopeArea*incidentFlux) << endl;

}
