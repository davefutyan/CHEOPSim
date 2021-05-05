/*
 * FrameTransferSmearer.cxx
 *
 *  Created on: 22 Jul 2014
 *      Author: futyand
 */

#include "FrameTransferSmearer.hxx"

void FrameTransferSmearer::initialize(const ModuleParams& params) {

	m_transferClockPeriod = params.GetAsDouble("transferClockPeriod")/1.E6;

}

void FrameTransferSmearer::process(Data* data, int timeStep, bool fullFrame) const {

	//Three images are output by PSFGenerator: data->getImages().back() is an image with
	//brightness corresponding to the exposure duration, containing jittered PSF(s).
	//data->getImageToSmearUp() is an image with brightness corresponding to 1 second,
	//with PSF position(s) corresponding to the last second (i.e. last jitter sampling) of the exposure.
	//data->getImageToSmearDown() is an image with brightness corresponding to 1 second,
	//with PSF position(s) corresponding to the first second (i.e. first jitter sampling) of the exposure.

	Image * image = data->getImages().back();

	if (data->getImagesToSmearUp().size()==0 || data->getImagesToSmearDown().size()==0 ) {
		throw runtime_error("Error: FrameTransferSmearer cannot be run without running PSFGenerator or DarkCurrentGenerator first");
	}

	Image * imageToSmearUp = data->getImagesToSmearUp().back();
	Image * imageToSmearDown = data->getImagesToSmearDown().back();

	int xOffset_smear = imageToSmearUp->getXOffset();
	int yOffset_smear = imageToSmearUp->getYOffset();
	int xDim_smear = imageToSmearUp->getXDim();
	int yDim_smear = imageToSmearUp->getYDim();

	//Scale the brightness of the image to be smeared by the transfer clock period
	(*imageToSmearUp) *= m_transferClockPeriod;
	(*imageToSmearDown) *= m_transferClockPeriod;

	//Make a copy of the image to be smeared to be subtracted from the smear trail (the initial image is part of the exposure, not the trail)
	Image * imageToSmearUp_copy = new Image(xDim_smear,yDim_smear,xOffset_smear,yOffset_smear);
	*imageToSmearUp_copy += *imageToSmearUp;
	*imageToSmearUp_copy *= -1.;

	//Loop over the image to smear, adding to each row the values in the row below/above in order to generate the smear trails
	for (int ix=xOffset_smear; ix<xOffset_smear+xDim_smear; ix++) {
		for (int iy=max(1,yOffset_smear); iy<yOffset_smear+yDim_smear+Image::kNDarkRows+Image::kNOverscanRows; iy++) {
			imageToSmearUp->incrementPixelValue(ix,iy,imageToSmearUp->getPixelValue(ix,iy-1));
		}
		for (int iy=min(Image::kYDim-1,yOffset_smear+yDim_smear); iy>=yOffset_smear; iy--) {
			imageToSmearDown->incrementPixelValue(ix,iy,imageToSmearDown->getPixelValue(ix,iy+1));
		}
	}

	//Subtract the original image to be smeared from the smear trail (the initial image is part of the exposure, not the trail)
	*imageToSmearUp += *imageToSmearUp_copy;
	delete imageToSmearUp_copy;

	//Store the truth information for the smear trails
	vector<double> smearTrailTruthRow;
	int top = Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows-1;
	for (int ix=xOffset_smear; ix<xOffset_smear+xDim_smear; ix++) {
		smearTrailTruthRow.push_back(imageToSmearUp->getPixelValue(ix,top));
	}
	image->getTruthData()->setSmearTrailRow(smearTrailTruthRow);

	//Add the smear trails to the original image
	for (int ix=0; ix<Image::kXDim; ix++) {
		for (int iy=0; iy<Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows; iy++) {
			image->incrementPixelValue(ix,iy,imageToSmearUp->getPixelValue(ix,iy));
			image->incrementPixelValue(ix,iy,imageToSmearDown->getPixelValue(ix,iy));
		}
	}

}
