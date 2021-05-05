/*
 * TruthData.hxx
 *
 *  Created on: 3 Dec 2014
 *      Author: futyand
 */

#ifndef TRUTHDATA_HXX_
#define TRUTHDATA_HXX_

#include <vector>
#include <numeric>
using namespace std;

#include "PSF.hxx"

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container class for truth information
///
/// An instance of this class is contained in the Data class
///
////////////////////////////////////////////////////////////////////////

class TruthData {
public:

	/// @brief Type of hot pixel
	enum HOTPIXEL_TYPE{hot,warm,telegraphicActive,telegraphicInactive};

	/// @brief Cosmic ray pixel
	struct CosmicPixel {

		/** *************************************************************************
		 *  @brief CosmicPixel constructor
		 *
		 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
		 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
		 *  @param [in] val  Number of electrons generated in the pixel by the CR
		 */
		CosmicPixel(int ix, int iy, int val) : m_xPixel(ix), m_yPixel(iy), m_nElectrons(val) {};

		int m_xPixel; ///< x position of the pixel (exposed full frame coordinates)
		int m_yPixel; ///< y position of the pixel (exposed full frame coordinates)
		int m_nElectrons; ///< Number of electrons generated in the pixel by the CR
	};

	/// @brief Hot pixel
	struct HotPixel {

		/** *************************************************************************
		 *  @brief HotPixel constructor
		 *
		 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
		 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
		 *  @param [in] val  Number of electrons in the hot pixel for the current exposure
		 *  @param [in] type  Type of hot pixel (hot, warm, telegraphicActive or telegraphicInactive)
		 */
		HotPixel(int ix, int iy, int val, HOTPIXEL_TYPE type) : m_xPixel(ix), m_yPixel(iy), m_nElectrons(val), m_type(type) {};
		int m_xPixel; ///< x position of the pixel (exposed full frame coordinates)
		int m_yPixel; ///< y position of the pixel (exposed full frame coordinates)
		int m_nElectrons; ///< Number of electrons in the hot pixel for the current exposure
		HOTPIXEL_TYPE m_type; ///< Type of hot pixel (hot, warm, telegraphicActive or telegraphicInactive)

	};

	/// @brief Dead pixel
	struct DeadPixel {

		/** *************************************************************************
		 *  @brief DeadPixel constructor
		 *
		 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
		 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
		 *  @param [in] val  Quantum efficiency of the dead pixel
		 */
		DeadPixel(int ix, int iy, double val) : m_xPixel(ix), m_yPixel(iy), m_quantumEfficiency(val) {};
		int m_xPixel; ///< x position of the pixel (exposed full frame coordinates)
		int m_yPixel; ///< y position of the pixel (exposed full frame coordinates)
		double m_quantumEfficiency; ///< Quantum efficiency of the dead pixel
	};

	TruthData() : m_fullWellSaturated(false),m_adcSaturated(false),m_zodiacalLight(0.),m_strayLight(0.) {};
	virtual ~TruthData();

	/** *************************************************************************
	 *  @brief Addition operator: combines the information contained in two
	 *  	   TruthData objects: For PSFs, fluxes are summed and list of
	 *  	   jittered positions concatenated; flux values for
	 *  	   hot/warm/telegraphic pixels are incremented; the list of cosmics
	 *  	   is concatenated; the smear trail truth row data is incremented.
	 */
	TruthData & operator +=(const TruthData & truthData);

	/// @brief Assignment operator: Sets the contents of a TruthData object to be equal to those of another
	TruthData & operator =(const TruthData & truthData);

	/** *************************************************************************
	 *  @brief Adds a PSF to the list
	 *
	 *  @param [in] x  Vector of truth x positions of the PSF centre (exposed full frame coordinates, one value for each second of the exposure)
	 *  @param [in] y  Vector of truth y positions of the PSF centre (exposed full frame coordinates, one value for each second of the exposure)
	 *  @param [in] flux  Truth PSF flux
	 */
	void addPSF(vector<double> x, vector<double> y, double flux) {m_psfs.push_back(PSF(x,y,flux));}

	/** *************************************************************************
	 *  @brief Adds a CosmicPixel to the list
	 *
	 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
	 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
	 *  @param [in] nElectrons  Number of electrons generated in the pixel by the CR
	 */
	void addCosmicPixel(int ix, int iy, int nElectrons) {m_cosmicPixels.push_back(CosmicPixel(ix,iy,nElectrons));}

	/** *************************************************************************
	 *  @brief Adds a HotPixel to the list
	 *
	 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
	 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
	 *  @param [in] nElectrons  Number of electrons in the hot pixel for the
	 *  						current exposure
	 *  @param [in] type  Type of hot pixel (hot, warm, telegraphicActive or telegraphicInactive)
	 */
	void addHotPixel(int ix, int iy, int nElectrons, HOTPIXEL_TYPE type) {m_hotPixels.push_back(HotPixel(ix,iy,nElectrons,type));}

	/** *************************************************************************
	 *  @brief Adds a DeadPixel to the list
	 *
	 *  @param [in] ix  x position of the pixel (exposed full frame coordinates)
	 *  @param [in] iy  y position of the pixel (exposed full frame coordinates)
	 *  @param [in] QE  Quantum efficiency of the dead pixel
	 */
	void addDeadPixel(int ix, int iy, double QE) {m_deadPixels.push_back(DeadPixel(ix,iy,QE));}

	/** *************************************************************************
	 *  @brief Sets the smear trail row truth data
	 *
	 *  @param [in] smearTrailRow  Array of pixels making up the smear trail
	 *  						   truth row (top overscan row)
	 */
	void setSmearTrailRow(vector<double> smearTrailRow) {m_smearTrailRow = smearTrailRow;}

	/// @brief Sets the full well saturated image flag to true
	void flagFullWellSaturation() {m_fullWellSaturated = true;}

	/// @brief Sets the ADC saturated image flag to true
	void flagAdcSaturation() {m_adcSaturated = true;}

	/// @brief Sets the global throughput
	void setGlobalThroughput(double globalThroughput) {m_globalThroughput.push_back(globalThroughput);}

	/// @brief Sets the ADC gain
	void setGain(double gain) {m_gain.push_back(gain);}

	/// @brief Returns the list of PSFs corresponding to the exposure (1 per second)
	vector<PSF> getPSFs() {return m_psfs;}

	/// @brief Returns the list of CosmicPixels corresponding to the exposure
	vector<CosmicPixel> getCosmicPixels() {return m_cosmicPixels;}

	/// @brief Returns the list of HotPixels
	vector<HotPixel> getHotPixels() {return m_hotPixels;}

	/// @brief Returns the list of DeadPixels
	vector<DeadPixel> getDeadPixels() {return m_deadPixels;}

	/// @brief Returns the smear trail row truth data
	vector<double> getSmearTrailRow() {return m_smearTrailRow;}

	/// @brief Returns whether or not an image contains full well saturated pixels
	bool isFullWellSaturated() {return m_fullWellSaturated;}

	/// @brief Returns whether or not an image contains ADC saturated pixels
	bool isAdcSaturated() {return m_adcSaturated;}

	/// @brief Returns the wavelength integral of Blackbody(target star) * Optical throughput * QE (incident photons -> electrons conversion factor).
	///		   The value is averaged in the case of image stacking
	double getGlobalThroughput() {return m_globalThroughput.size()>0 ? (std::accumulate(m_globalThroughput.begin(),m_globalThroughput.end(),0.0))/m_globalThroughput.size() : 1.;}

	/// @brief Sets the ADC gain value at the CCD temperature corresponding to the image (electrons->ADU conversion factor)
	///		   The value is averaged in the case of image stacking
	double getGain() {return m_globalThroughput.size()>0 ? (std::accumulate(m_gain.begin(),m_gain.end(),0.0))/m_gain.size() : 1.;}

	/// @brief Sets the zodiacal light photons per pixel
	void setZodiacalLight(double zodiacalLight) {m_zodiacalLight = zodiacalLight;}

	/// @brief Returns the zodiacal light photons per pixel
	double getZodiacalLight() {return m_zodiacalLight;}

	/// @brief Sets the stray light photons per pixel
	void setStrayLight(double strayLight) {m_strayLight = strayLight;}

	/// @brief Returns the stray light photons per pixel
	double getStrayLight() {return m_strayLight;}

private:
	vector<PSF> m_psfs; ///< list of PSFs corresponding to the exposure (1 per second)
	vector<CosmicPixel> m_cosmicPixels; ///< list of CosmicPixels corresponding to the exposure
	vector<HotPixel> m_hotPixels; ///< List of HotPixels
	vector<DeadPixel> m_deadPixels; ///< List of DeadPixels
	vector<double> m_smearTrailRow; ///< Smear trail row truth data
	bool m_fullWellSaturated; ///< Boolean to flag images with full well saturated pixels
	bool m_adcSaturated; ///< Boolean to flag images with ADC saturated pixels
	vector<double> m_globalThroughput; ///< Wavelength integral of Blackbody(target star) * Optical throughput * QE (incident photons -> electrons conversion factor)
	vector<double> m_gain; ///< ADC gain value at the CCD temperature corresponding to the image (electrons->ADU conversion factor)
	double m_zodiacalLight; ///< Zodiacal light photons per pixel
	double m_strayLight; ///< Stray light photons per pixel
};

inline TruthData operator +(TruthData truth1, const TruthData & truth2) {
	truth1+=truth2;
	return truth1;
}

#endif /* TRUTHDATA_HXX_ */
