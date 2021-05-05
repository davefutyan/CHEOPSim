/*
 * CosmicRayGenerator.cxx
 *
 *  Created on: 12 May 2014
 *      Author: futyand
 */

#include <cmath>
#include <iostream>
#include <fstream>

#include "CosmicRayGenerator.hxx"

void CosmicRayGenerator::initialize(const ModuleParams& params) {

	//For HST (Instrument Science Report ACS 2007-07), rate is 2.2% of pixels affected per 1000s.
	//From this module, mean number of pixels per cosmic ray is 4.9 for pixel depth 15 microns.
	//Implies mean no. of cosmic rays per 1000s should be no. of pixels * (0.022/4.9) = 180 for 200x200
	//Implies mean no. of cosmic rays per minute (stacked image time) = 10.8
	unsigned seed = params.GetAsInt("cosmicSeed");
	m_cosmicsPerMinute = params.GetAsDouble("meanCosmicsPerMinute");
	m_cosmicsPerMinute *= static_cast<double>(Image::kXDim * Image::kYDim) / static_cast<double>(200*200);
	m_SAAFluxFactor = params.GetAsDouble("SAAFluxFactor");
	m_cosmicEnergyScaleFactor = params.GetAsDouble("cosmicEnergyScaleFactor");

	//Initialize random number generators. In principle only one Poisson RNG is needed, but keep
	//two so that the distribution of cosmics is the same as in earlier releases in which
	//a Gaussian RNG was used for the number of electrons per pixel in the cosmic trails
	boost::mt19937 randomNumberEngine(seed);
	boost::uniform_real<double> uniform_random(0.,1.);
	m_uniformRandomGenerator = new RANDOM_UNIFORM_REAL(randomNumberEngine,uniform_random);
    boost::poisson_distribution<int,double> poissonDistribution(1.);
    m_poissonRandomGenerator = new RANDOM_POISSON(randomNumberEngine,poissonDistribution);
    m_poissonRandomGenerator2 = new RANDOM_POISSON(randomNumberEngine,poissonDistribution);

    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    m_gslRandomGenerator = gsl_rng_alloc(T);

    //hist_file = new TFile("cosmics.root");
    //hist_logNElectrons = new TH1F("hist_logNElectrons","log of number of electrons per cosmic ray",180,1.4,5.);
}

CosmicRayGenerator::~CosmicRayGenerator() {
	//hist_file->Write();
	//hist_file->Close();
}

void CosmicRayGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	if (m_cosmicsPerMinute <= 0.) {
		cout << "mean rate for cosmics is set to zero or negative - no cosmics will be generated" << endl;
		return;
	}

	//Get the pointer to the cosmic ray truth image, or create it if no image is currently in memory (first exposure of a new stacked image)
	Image* image = data->getImages().back();
	Image* image_cr = data->getCosmicRayImage();
	if (image_cr == nullptr) {
		image_cr = new Image(image->getXDim(),image->getYDim(),image->getXOffset(),image->getYOffset());
		data->setCosmicRayImage(image_cr);
	}

	//double mean_length = 0.;

	//Get the mean number of cosmics to generate as the cosmics rate per second multiplied by the time for which the image is exposed to cosmics
	//This time is the sum of the exposure time and half of the readout time (averaging over the fact that the top of the image waits longer to be read out than the bottom)
	double meanNumberOfCosmics = (m_cosmicsPerMinute/60.) * (data->getTimeConfiguration().getExposureTimeAsDouble() + data->getReadoutTime()/2.);

	//Increase the flux if the satellite is over the SAA
	if (data->getSatelliteData(timeStep)->getSAAFlag()) meanNumberOfCosmics *= m_SAAFluxFactor;

	//generate the number of cosmics for the current image by drawing randomly from a Poisson
	m_poissonRandomGenerator->distribution().param(boost::poisson_distribution<int,double>::param_type(meanNumberOfCosmics));
	int numberOfCosmics = (*m_poissonRandomGenerator)();

	//ofstream hist_file;
	//hist_file.open("cosmics.dat");

	//cout << "numberOfCosmics = " << numberOfCosmics << endl;

	for (int i=0; i<numberOfCosmics; i++) {

		//Generate the cosmic ray path
		double cosmic_x = ((*m_uniformRandomGenerator)()*(Image::kXDim+2*Image::kNDarkCols))-Image::kNDarkCols; //x-coordinate of impact position of cosmic on CCD surface
		double cosmic_y = (*m_uniformRandomGenerator)()*(Image::kYDim+Image::kNDarkRows); //y-coordinate of impact position of cosmic on CCD surface
		double cosmic_azimuth = (*m_uniformRandomGenerator)()*2.*M_PI; //azimuthal angle of cosmic direction of travel
		double cosmic_polar = (*m_uniformRandomGenerator)()*M_PI/2.; //polar angle of cosmic direction of travel
		CosmicPath cosmicPath(cosmic_x,cosmic_y,cosmic_polar,cosmic_azimuth);
		//cout << endl << "New Cosmic, length = " << cosmicPath.m_length << " " << atan2(cosmicPath.m_length,kPixelDepth/kPixelWidth) << endl;

		//Skip the cosmic if no part of it is inside the sub-array
		if ((cosmicPath.m_xstart<image->getXOffset() && cosmicPath.m_xend<image->getXOffset()) ||
			(cosmicPath.m_xstart>image->getXOffset()+image->getXDim() && cosmicPath.m_xend>image->getXOffset()+image->getXDim()) ||
			(cosmicPath.m_ystart<image->getYOffset() && cosmicPath.m_yend<image->getYOffset()) ||
			(cosmicPath.m_ystart>image->getYOffset()+image->getYDim() && cosmicPath.m_yend>image->getYOffset()+image->getYDim())) continue;

		//Draw the total number of electrons randomly from a Landau distribution
		double nElectrons_tot = cosmicPath.m_length3D *kPixelWidth * 650./kPixelDepth;
		nElectrons_tot += 150.*(gsl_ran_landau(m_gslRandomGenerator));
		if (nElectrons_tot > 50000.) {i--; continue;} //cut off at 50000 electrons to avoid extreme tail

		//Scale the energy according to the configuration (default factor 1)
		nElectrons_tot *= m_cosmicEnergyScaleFactor;

		//Propagate along the line in steps of kPropagationStep pixels and count the number of steps within
		//each pixel crossed as a measure of the length traversed within each pixel
		double x = cosmicPath.m_xstart;
		double y = cosmicPath.m_ystart;
		double z = 0.; //define z=0 at upper surface of CCD. Polar angle is in range 0,pi/2 so is always from above
		vector<CosmicRayPixel*> cosmicRay;
		int ix_previous = -999;
		int iy_previous = -999;
		double propagationStep = kPropagationStep;
		if (data->getSatelliteData(timeStep)->getSAAFlag()) propagationStep*=10.; //reduce the number of steps when in SAA to improve execution time
		while (propagateCosmic(cosmicPath,propagationStep,x,y,z)) {
			//cout << cosmicPath.m_polar << " " << cosmicPath.m_azimuth << " " << cosmicPath.m_xstart << " " << cosmicPath.m_xend << " " << cosmicPath.m_ystart << " " << cosmicPath.m_yend << " " << x << " " << y << " " << z << endl;
			int ix_pixel = static_cast<int>(floor(x));
			int iy_pixel = static_cast<int>(floor(y));
			if (ix_pixel != ix_previous || iy_pixel != iy_previous) {
				if (ix_pixel > Image::kXDim+Image::kNDarkCols-1 || iy_pixel > Image::kYDim+Image::kNDarkRows-1 || ix_pixel < -Image::kNDarkCols || iy_pixel < 0) break;
				CosmicRayPixel * pixel = new CosmicRayPixel(ix_pixel,iy_pixel,propagationStep);
				cosmicRay.push_back(pixel);
				ix_previous = ix_pixel;
				iy_previous = iy_pixel;
			} else {
				cosmicRay.back()->incrementValue(propagationStep);
			}
		}

		//mean_length += cosmicRay.size();

		for (vector<CosmicRayPixel*>::const_iterator pixel = cosmicRay.begin(); pixel!=cosmicRay.end(); ++pixel) {
			double nElectrons = nElectrons_tot * ((*pixel)->getValue()/cosmicPath.m_length3D);
			m_poissonRandomGenerator2->distribution().param(boost::poisson_distribution<int,double>::param_type(nElectrons));
			nElectrons = (*m_poissonRandomGenerator2)();
			int ix =(*pixel)->getXIndex();
			int iy =(*pixel)->getYIndex();
			if ((ix>=image->getXOffset() && ix<image->getXOffset()+image->getXDim() &&
					((iy>=image->getYOffset() && iy<image->getYOffset()+image->getYDim()) ||
					 (iy>=Image::kYDim && iy<Image::kYDim+Image::kNDarkRows))) ||
				((iy>=image->getYOffset() && iy<image->getYOffset()+image->getYDim() &&
					((ix>=-Image::kNDarkCols && ix<0) ||
					 (ix>=Image::kXDim && ix<Image::kXDim+Image::kNDarkCols))))) {
				image->incrementPixelValue(ix,iy,round(nElectrons));
				image_cr->incrementPixelValue(ix,iy,round(nElectrons));
				image->getTruthData()->addCosmicPixel(ix,iy,round(nElectrons));
			}
			delete (*pixel);
		}

		//hist_file << (nElectrons_tot>0 ? log10(nElectrons_tot) : 0.) << " " << cosmicRay.size() << endl;
		//hist_logNElectrons->Fill(log10(nElectrons_tot));

		cosmicRay.clear();
	}

	//hist_file.close();

	//if (numberOfCosmics > 0) mean_length /= numberOfCosmics;
	//cout << mean_length << endl;

}

bool CosmicRayGenerator::propagateCosmic(const CosmicPath & cosmicPath, double step, double & x, double & y, double & z) const {

	x += step*sin(cosmicPath.m_polar)*cos(cosmicPath.m_azimuth);
	y += step*sin(cosmicPath.m_polar)*sin(cosmicPath.m_azimuth);
	z -= step*cos(cosmicPath.m_polar);

	if ((cosmicPath.m_xend > cosmicPath.m_xstart && x > cosmicPath.m_xend) || (cosmicPath.m_xend < cosmicPath.m_xstart && x < cosmicPath.m_xend)) return false;
	if ((cosmicPath.m_yend > cosmicPath.m_ystart && y > cosmicPath.m_yend) || (cosmicPath.m_yend < cosmicPath.m_ystart && y < cosmicPath.m_yend)) return false;
	if (z < -kPixelDepth/kPixelWidth) return false;

	return true;

}
