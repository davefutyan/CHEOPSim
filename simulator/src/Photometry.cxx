/*
 * Photometry.cxx
 *
 *  Created on: 21 Nov 2014
 *      Author: futyand
 */

#include "CommonTools.hxx"
#include "Photometry.hxx"

double Photometry::extractFlux(const Image* image, std::pair<double,double> barycentre) const {

	//Evaluate the background as the median within an annulus centred on the barycentre of the PSF
	double background_psfCentred = medianBackground(image,barycentre.first,barycentre.second);

	//Extract total flux within a circle of radius m_radius_psf centred on the barycentre of the PSF, subtracting the background
	double flux = 0.;
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=image->getYOffset(); iy<image->getYOffset()+image->getYDim(); iy++) {
			double x = double(ix) - floor(m_targetLocationX + barycentre.first);
			double y = double(iy) - floor(m_targetLocationY + barycentre.second);
			double radius = sqrt(x*x + y*y);
			if (radius < m_radius_psf) {
				double pixelValue = image->getPixelValue(ix,iy);
				if (m_subtractBackground) pixelValue -= background_psfCentred;
				if (m_flatField != nullptr) pixelValue /= m_flatField->getPixelValue(ix,iy);
				flux += pixelValue;
			}
		}
	}

	return flux;

}

std::pair<double,double> Photometry::getBarycentre(const Image* image) const {

	//Evaluate the background as the median within an annulus around the intended target location
	double background_ccdCentred = medianBackground(image,0.,0.);

	//Determine the barycentre of the PSF, subtracting the background
	double xBarycentre = 0;
	double yBarycentre = 0;
	double flux_tot = 0.;
	for (int ix=static_cast<int>(m_targetLocationX-(m_radius_barycentre+1.)); ix<static_cast<int>(m_targetLocationX+(m_radius_barycentre+1.)); ix++) {
		for (int iy=static_cast<int>(m_targetLocationY-(m_radius_barycentre+1.)); iy<static_cast<int>(m_targetLocationY+(m_radius_barycentre+1.)); iy++) {
			//-0.5 to take into account offset of 0.5 pixels between centre of pixel (512,512)
			//and the centre of the image, which is at the boundary between pixels 511 and 512 (i.e. 511.5):
			double x = double(ix)+0.5 - m_targetLocationX;
			double y = double(iy)+0.5 - m_targetLocationY;
			double radius = sqrt(x*x + y*y);
			if (radius < m_radius_barycentre) {
				double pixelflux = image->getPixelValue(ix,iy) - background_ccdCentred;
				//Weight by square of pixel value, for consistency with IntensityWeightedCenterOfGravity2D calculation in IFSW
				xBarycentre += (x * pixelflux*pixelflux);
				yBarycentre += (y * pixelflux*pixelflux);
				flux_tot += pixelflux*pixelflux;
			}
		}
	}
	xBarycentre /= flux_tot;
	yBarycentre /= flux_tot;

	return std::make_pair(xBarycentre,yBarycentre);

}

double Photometry::medianBackground(const Image * image, double xOffset, double yOffset) const {

	std::vector<double> pixelValues;
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=image->getYOffset(); iy<image->getYOffset()+image->getYDim(); iy++) {
			double x = double(ix) - floor(m_targetLocationX + xOffset);
			double y = double(iy) - floor(m_targetLocationY + yOffset);
			double radius = sqrt(x*x + y*y);
			if (radius > m_radius_bkgInner && radius < m_radius_bkgOuter) {
				pixelValues.push_back(image->getPixelValue(ix,iy));
			}
		}
	}

	return CommonTools::median(pixelValues);

}
