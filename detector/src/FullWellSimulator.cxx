/*
 * FullWellSimulator.cxx
 *
 *  Created on: 5 Aug 2014
 *      Author: futyand
 */

#include <fstream>

#include "FullWellSimulator.hxx"

void FullWellSimulator::initialize(const ModuleParams& params) {

	m_fullWellCapacity = params.GetAsInt("fullWellCapacity");

}

void FullWellSimulator::process(Data* data, int timeStep, bool fullFrame) const {

	Image* image = data->getImages().back();

	double delta = 1.E-8;
	double fullWellCapacity = static_cast<double>(m_fullWellCapacity);
	bool saturated = false;

	//Uncomment the following line to append the time to saturation to an ascii file.
	//This can be used together with the timeToSaturation.sh script to generate an array of
	//values for the time to saturation for different magnitudes and spectral types
	//calculateTimeToSaturation(data);

	//Bleeding occurs only vertically, so treat each column independently
	for (int ix=0; ix<Image::kXDim; ix++) {

		//Iterate until there are no pixels in the column with content exceeding the full well capacity
		double nElectronsMax = fullWellCapacity+1.;
		while (nElectronsMax > fullWellCapacity+delta) {

			//Find the pixel within the column that has the maximum no. of electrons
			nElectronsMax = 0.;
			int iymax = 0;
			for (int iy=0; iy<Image::kYDim; iy++) {
				double nElectrons = image->getPixelValue(ix,iy);
				if (nElectrons > nElectronsMax) {
					iymax = iy;
					nElectronsMax = nElectrons;
				}
			}
			//cout << iymax << " " << nElectronsMax << endl;

			//If the pixel is saturated, distribute the excess to the fill the nearest unsaturated pixels
			//in the vertical direction, with half being sent in each direction
			if (nElectronsMax > fullWellCapacity+delta) {
				saturated = true;
				double nOverflow = nElectronsMax - fullWellCapacity;
				image->setPixelValue(ix,iymax,fullWellCapacity);
				overflow(image,ix,iymax,nOverflow/2.,true); //half the excess is sent upwards
				overflow(image,ix,iymax,nOverflow/2.,false); //half the excess is sent downwards
			}

		}

	}

	if (saturated) image->getTruthData()->flagFullWellSaturation();

}

void FullWellSimulator::overflow(Image * image, int ix, int iy, double nOverflow, bool up) const {

	double fullWellCapacity = static_cast<double>(m_fullWellCapacity);
	double delta = 1.E-8;

	//Keep shifting vertically until all the excess electrons have been added to an unsaturated pixel,
	//or the edge of the CCD is reached
	while (nOverflow>0. && (up?iy<Image::kYTotal:iy>=0)) {
		up ? iy++ : iy--;
		double nAvailable = fullWellCapacity - image->getPixelValue(ix,iy);
		if (nAvailable > delta) { //pixel has space for nAvailable electrons
			image->incrementPixelValue(ix,iy,min(nOverflow,nAvailable)); //add the excess electrons to the pixel
			nOverflow -= nAvailable; //remaining excess to be sent on to the next pixel
		}
	}

}

void FullWellSimulator::calculateTimeToSaturation(const Data * data) const {

	double nElectrons_max = 0.;
	for (int ix=0; ix<Image::kXDim; ix++) {
		for (int iy=0; iy<Image::kYDim; iy++) {
			double nElectrons = data->getImages().back()->getPixelValue(ix,iy);
			if (nElectrons > nElectrons_max) nElectrons_max = nElectrons;
		}
	}

	cout << "max nElectrons in image = " << round(nElectrons_max) << endl;

	double exposureTime = data->getTimeConfiguration().getExposureTimeAsDouble();
	double timeToSaturation = static_cast<double>(m_fullWellCapacity)*exposureTime/nElectrons_max;
	cout << "time to saturation: " << timeToSaturation << " seconds" << endl;

	double magnitude = data->getFieldOfView()->getStars()[0]->getGaiaMagnitude();
	Star::spectralType spectralType = data->getFieldOfView()->getStars()[0]->getSpectralType();

	ofstream saturation_file("saturation.txt", ios::app);
	if (spectralType==Star::O7) saturation_file << magnitude;
	saturation_file << "	" << timeToSaturation;
	if (spectralType==Star::M9) saturation_file << endl;
	saturation_file.close();

	ofstream saturation_array_file("saturation_array.txt", ios::app);
	if (spectralType==Star::O7) saturation_array_file << "      [" << magnitude;
	saturation_array_file << "," << timeToSaturation;
	if (spectralType==Star::M9) saturation_array_file << "]," << endl;
	saturation_array_file.close();

}
