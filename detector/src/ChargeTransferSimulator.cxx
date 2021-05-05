/*
 * ChargeTransferSimulator.cxx
 *
 *  Created on: 12 Mar 2014
 *      Author: futyand
 */

#include <boost/math/distributions/exponential.hpp>

#include "ChargeTransferSimulator.hxx"

void ChargeTransferSimulator::initialize(const ModuleParams& params) {

	m_endOfLife = params.GetAsBool("endOfLife");
	m_cte_vertical = params.GetAsDouble("cte_vertical")/100.;
	m_cte_horizontal = params.GetAsDouble("cte_horizontal")/100.;
	m_cti_trailFraction = params.GetAsDouble("cti_trailFraction");
	m_cti_intensityScaling = params.GetAsDouble("cti_intensityScaling");

	//Fill an array to hold the value of an exponential PDF with lambda=cti_trailLength in unit (1 pixel) steps
	boost::math::exponential_distribution<double> exponentialPDF = boost::math::exponential_distribution<double>(1./params.GetAsDouble("cti_trailLength"));
	for (int i=0;i<Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows; i++) {
		m_exponentialDistribution[i] = boost::math::pdf(exponentialPDF,i);
	}

}

void ChargeTransferSimulator::process(Data* data, int timeStep, bool fullFrame) const {

	Image * image = data->getImages().back();

	m_endOfLife ? processEndOfLife(image) : processBeginOfLife(image);

}

void ChargeTransferSimulator::processBeginOfLife(Image * image) const {

	int xOffset = image->getXOffset();
	int yOffset = image->getYOffset();
	int xDim = image->getXDim();
	int yDim = image->getYDim();

	//ccdArray is an array representing the entire CCD including the storage area, margins and pseudo-pixels
	//Initially ccdArray has the storage region empty and the exposed region contains the input image
	double ** ccdArray = new double*[Image::kXTotal];
	for (int i = 0; i < Image::kXTotal; ++i) ccdArray[i] = new double[Image::kYDim+Image::kYTotal];
	for (int ix=0; ix<Image::kXTotal; ix++) {
		for (int iy=0; iy<Image::kYDim; iy++) { //storage region
			ccdArray[ix][iy] = 0.;
		}
		for (int iy=Image::kYDim; iy<Image::kYDim+Image::kYTotal; iy++) { //exposed region
			ccdArray[ix][iy] = image->getPixelValue(ix-Image::kLeftMargin,iy-Image::kYDim);
		}
	}

	//First step is to simulate the frame transfer
	//Loop over all rows
	for (int iVert = 0; iVert<Image::kYTotal; iVert++) {
		if ((iVert >= yOffset && iVert<yOffset+yDim) || iVert >= Image::kYDim) { //Sub-array plus top margin

			//Transfer the current row down by Image::kYDim rows in successive shifts, taking into account CTE
			for (int iv = Image::kYDim+iVert; iv>iVert; iv--) {
				for (int iHoriz = 0; iHoriz<Image::kXTotal; iHoriz++) {
					//Only process pixels within the sub-array or left and right margins
					if ((iHoriz >= Image::kLeftMargin+xOffset && iHoriz<Image::kLeftMargin+xOffset+xDim) ||
							iHoriz < Image::kLeftMargin || iHoriz >= Image::kLeftMargin+Image::kXDim) {
						ccdArray[iHoriz][iv-1] = ccdArray[iHoriz][iv]*m_cte_vertical + ccdArray[iHoriz][iv-1]*(1.-m_cte_vertical);
					}
				}
			}

		}
	}

	//Second step is to simulate the readout
	//Loop over all rows
	for (int iVert = 0; iVert<Image::kYTotal; iVert++) {
		if ((iVert >= yOffset && iVert<yOffset+yDim) || iVert >= Image::kYDim) { //Sub-array plus top margin

			//Transfer the current row down to the bottom of the frame in successive shifts, taking into account CTE
			for (int iv = iVert; iv>0; iv--) {
				for (int iHoriz = 0; iHoriz<Image::kXTotal; iHoriz++) {
					//Only process pixels within the sub-array or left and right margins
					if ((iHoriz >= Image::kLeftMargin+xOffset && iHoriz<Image::kLeftMargin+xOffset+xDim) ||
							iHoriz < Image::kLeftMargin || iHoriz >= Image::kLeftMargin+Image::kXDim) {
						ccdArray[iHoriz][iv-1] = ccdArray[iHoriz][iv]*m_cte_vertical + ccdArray[iHoriz][iv-1]*(1.-m_cte_vertical);
					}
				}
			}

			//After transferring the current row to the bottom, read out the bottom row
			for (int iHoriz = 0; iHoriz<Image::kXTotal; iHoriz++) {
				//Only process pixels within the sub-array and left and right margins
				if ((iHoriz >= Image::kLeftMargin+xOffset && iHoriz<Image::kLeftMargin+xOffset+xDim) ||
						iHoriz < Image::kLeftMargin || iHoriz >= Image::kLeftMargin+Image::kXDim) {

					//Transfer the contents of the pixel at location iHoriz to the left in successive shifts, taking into account CTE
					for (int ih = iHoriz; ih>0; ih--) {
						ccdArray[ih-1][0] = ccdArray[ih][0]*m_cte_horizontal + ccdArray[ih-1][0]*(1.-m_cte_horizontal);
					}

					//After transferring the contents of pixel iHoriz to the leftmost pixel, read out that pixel, taking into account CTE
					double readoutValue = m_cte_horizontal*ccdArray[0][0];

					//Set the value for pixel (iHoriz,iVert) to the readout value
					image->setPixelValue(iHoriz-Image::kLeftMargin,iVert,readoutValue);

				}
			}

		}
	}

	for(int i = 0; i < Image::kXTotal; ++i) delete [] ccdArray[i];
	delete [] ccdArray;

}

void ChargeTransferSimulator::processEndOfLife(Image * image) const {

	int xOffset = image->getXOffset();
	int yOffset = image->getYOffset();
	int xDim = image->getXDim();
	int yDim = image->getYDim();

	Image * image_trails = new Image(xDim,yDim,xOffset,yOffset);

	for (int ix=xOffset; ix<xOffset+xDim; ix++) {

		//double tail_integral = 0.;
		//double PSF_integral = 0.;
		//double origPSF_integral = 0.;
		for (int iy=0; iy<Image::kYDim+Image::kNDarkRows; iy++) {
			//double integral = 0.;

			//Get the value of the current pixel
			double pixel_value = image->getPixelValue(ix,iy);
			//if (ix==512) origPSF_integral += pixel_value;

			//Calculate the fraction of the pixel intensity which will be transfered to the tail
			double tail_fraction = m_cti_trailFraction*pow(pixel_value/10000.,m_cti_intensityScaling);
			tail_fraction = min(tail_fraction,1.); //Limit the trail fraction to 100% of the original pixel intensity

			//Reduce the pixel value by the fraction which will be transfered to the tail
			image->setPixelValue(ix,iy,pixel_value*(1.-tail_fraction));
			//if (ix==512) PSF_integral += pixel_value*(1.-tail_fraction);

			//Generate trails for the main image by looping from the pixel above the current pixel to the top of the image,
			//adding a number of electrons to each pixel according to an exponential distribution
			for (int j=iy+1; j<yOffset+yDim; j++) {
				double value = tail_fraction*pixel_value*m_exponentialDistribution[j-iy];
				//integral += value;
				//if (ix==512) tail_integral += value;
				//if (ix==512 && iy==512) cout << j-iy << " " << pixel_value << " " << tail_fraction*pixel_value << " " << pixel_value*(1.-tail_fraction) << " "
				//							   << m_exponentialDistribution[j-iy] << " " << integral << " " << (pixel_value*(1.-tail_fraction)+integral)/pixel_value << endl;
				image_trails->incrementPixelValue(ix,j,value);
			}

			//Extend tails into the CCD top margins
			for (int j=Image::kYDim; j<Image::kYTotal; j++) {
				if (j>iy) image_trails->incrementPixelValue(ix,j,tail_fraction*pixel_value*m_exponentialDistribution[j-iy]);
			}

		}
		//if (ix==512) cout << origPSF_integral << " " << PSF_integral << " " << tail_integral << " " << (PSF_integral+tail_integral)/origPSF_integral << endl;
	}

	//Add the CTI trails to the original image
	for (int ix=xOffset; ix<xOffset+xDim; ix++) {
		for (int iy=0; iy<Image::kYTotal; iy++) {
			image->incrementPixelValue(ix,iy,image_trails->getPixelValue(ix,iy));
		}
	}

	delete image_trails;

}

