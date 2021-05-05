/*
 * Image.cxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#include <cmath>

#include "Image.hxx"

Image::Image(int xDim,int yDim,int xOffset,int yOffset) :
							m_xDim(xDim),m_yDim(yDim),
							m_xOffset(xOffset),m_yOffset(yOffset) {

	for (int ix=0; ix<kXTotal; ix++) {
		for (int iy=0; iy<kYTotal; iy++) {
			m_image[ix][iy] = 0.;
		}
	}

	m_truthData = new TruthData();

}

void Image::saturate(unsigned nBits) {

	for (int ix=0; ix<kXTotal; ix++) {
		for (int iy=0; iy<kYTotal; iy++) {
			if (m_image[ix][iy] >= pow(2.,nBits)) {
				m_image[ix][iy] = pow(2.,nBits)-1.;
				m_truthData->flagAdcSaturation();
			}
		}
	}

}

void Image::setTruthData(TruthData* truthData) {
	*(m_truthData) = *(truthData);
}

void Image::roundToIntegers() {

	for (int ix=0; ix<kXTotal; ix++) {
		for (int iy=0; iy<kYTotal; iy++) {
			m_image[ix][iy] = nearbyint(m_image[ix][iy]);
		}
	}

}

Image & Image::operator +=(const Image & image) {

	for (int ix=0; ix<kXTotal; ix++) {
		for (int iy=0; iy<kYTotal; iy++) {
			this->m_image[ix][iy] += image.m_image[ix][iy];
		}
	}

	*(this->m_truthData) += *(image.m_truthData);

	return *this;

}

Image & Image::operator *=(double scaleFactor) {

	for (int ix=0; ix<kXTotal; ix++) {
		for (int iy=0; iy<kYTotal; iy++) {
			this->m_image[ix][iy] *= scaleFactor;
		}
	}

	return *this;

}
