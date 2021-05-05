/*
 * Data.hxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Data container for CHEOPSim
///
/// Container for all data which can be output by CHEOPSim
/// A pointer to the data is passed as input argument to  Module::process().
/// In this way all processing modules have access to all the data
/// available at the point where Module::process() is called.
///
/// The data consists of the following four components, for which accessor
/// methods are provided:
///		1. TimeConfiguration: Container class for simulation start time,
///						   duration, time resolution and exposure time
///						   (initialized in the constructor)
///		2. SkyFieldOfView: The field of view (null pointer until initialized  by
///						   StarProducer), which in turn contains a vector of Stars
///		3. vector<const SatelliteData*>: Container for information relating to the
///						   satellite conditions at a given time (temperature, roll
///                        angle, jitter etc).  Each vector element corresponds to
///                        a different time step (the vector is initially empty).
///		4. vector<const Image*>: 2D image array.  Each vector element corresponds to
///                        a different time step (the vector is initially empty).
///
////////////////////////////////////////////////////////////////////////

#ifndef _DATA_HXX_
#define _DATA_HXX_

#include "boost/filesystem.hpp"

#include "ProgramParams.hxx"
#include "MPS_PRE_Visits.hxx"
#include "SCI_RAW_SubArray.hxx"
#include "SIM_RAW_UnstackedSubArray.hxx"
#include "SCI_RAW_ImageMetadata.hxx"
#include "SCI_RAW_UnstackedImageMetadata.hxx"
#include "SCI_RAW_Imagette.hxx"
#include "SCI_RAW_ImagetteMetadata.hxx"
#include "SIM_RAW_DoublePrecisionSubArray.hxx"
#include "SCI_RAW_Centroid.hxx"
#include "SIM_TRU_SubArray.hxx"
#include "SIM_TRU_FlatField.hxx"
#include "REF_APP_BadPixelMap.hxx"
#include "REF_APP_BadPixelMapLeft.hxx"
#include "REF_APP_BadPixelMapRight.hxx"
#include "REF_APP_BadPixelMapTop.hxx"
#include "REF_APP_PhotPixelMap.hxx"
#include "REF_APP_PhotPixelMapLeft.hxx"
#include "REF_APP_PhotPixelMapRight.hxx"
#include "REF_APP_PhotPixelMapTop.hxx"
#include "SCI_RAW_DarkLeft.hxx"
#include "SCI_RAW_DarkRight.hxx"
#include "SCI_RAW_DarkTop.hxx"
#include "SCI_RAW_BlankLeft.hxx"
#include "SCI_RAW_BlankRight.hxx"
#include "SCI_RAW_OverscanLeft.hxx"
#include "SCI_RAW_OverscanRight.hxx"
#include "SCI_RAW_OverscanTop.hxx"
#include "SIM_RAW_UnstackedDarkLeftImage.hxx"
#include "SIM_RAW_UnstackedDarkRightImage.hxx"
#include "SIM_RAW_UnstackedDarkTopImage.hxx"
#include "SIM_RAW_UnstackedBlankLeftImage.hxx"
#include "SIM_RAW_UnstackedBlankRightImage.hxx"
#include "SIM_RAW_UnstackedOverscanLeftImage.hxx"
#include "SIM_RAW_UnstackedOverscanRightImage.hxx"
#include "SIM_RAW_UnstackedOverscanTopImage.hxx"

#include "TimeConfiguration.hxx"
#include "SkyFieldOfView.hxx"
#include "Star.hxx"
#include "StarData.hxx"
#include "SatelliteData.hxx"
#include "Image.hxx"
#include "WavelengthDependence.hxx"

class Data {
public:

	static const int kProgramType; ///< Program type used in output data structure keywords
	static constexpr double kRawHKCadence = 1.2; ///< Cadence for on board HK measurements
	static constexpr double kAveragedHKCadence = 20.; ///< Cadence for writing HK data to ground

	/// @brief Telegraphic pixel
	struct TelegraphicPixel {
		/** *************************************************************************
		 *  @brief Constructor for telegraphic pixel
		 *
		 *  @param [in] ix  x position of the pixel
		 *  @param [in] iy  y position of the pixel
		 *  @param [in] meanActiveRate  Mean dark current at 233K for active state
		 *  @param [in] active  Boolean to indicate if the pixel is initially active
		 *  @param [in] transitions  List of times, in seconds since the start of the simulation, at which transitions occur
		 */
		TelegraphicPixel(int ix, int iy, double meanActiveRate, bool active, vector<float> transitions) :
						 m_xPixel(ix), m_yPixel(iy), m_meanActiveRate(meanActiveRate), m_active(active), m_transitions(transitions) {};
		void flip() {m_active = (m_active ? false : true);} ///< Switch the state between active/inactive
		int m_xPixel; ///< x position of the pixel
		int m_yPixel; ///< y position of the pixel
		double m_meanActiveRate; ///< Mean dark current at 233K for active state
		bool m_active; ///< boolean to indicate if the pixel is currently active
		vector<float> m_transitions; ///< list of times, in seconds since the start of the simulation, at which transitions occur
	};

	/// @brief Visit data
	struct Visit {
		/** *************************************************************************
		 *  @brief Constructor for visit data
		 *
		 *  @param [in] progType  Program Type (default kProgramType)
		 *  @param [in] progId  Program ID (default 1)
		 *  @param [in] reqId  Request ID (default 0)
		 *  @param [in] obsId  Observation ID (default -1)
		 *  @param [in] visitCtr  Visit Counter (default 1)
		 *  @param [in] obsCat  Observation Category (default "time critical")
		 */
		Visit(int progType=kProgramType, int progId=1, int reqId=0, int obsId=-1, int visitCtr=1, string obsCat="time critical") :
			m_progType(progType), m_progId(progId), m_reqId(reqId), m_obsId(obsId), m_visitCtr(visitCtr), m_obsCat(obsCat) {

			// Set default values for simulation
			m_versionNum = 0;
			m_piName = "CHEOPSim";
			m_piUid = 0;
			m_prpFirst = 365;
			m_prpLast = 548;
			m_targName = "simulation";
			m_readMode = "undefined";
			m_marginMode = "undefined";
			m_gaiaMagError = 0.;
			m_strayLightThreshold = 1.E9;

			// If no observation ID is provided, set it equal to the job ID obtained from the directory name
			if (m_obsId==-1) {
				m_obsId = 0;
				string jobdir = string(boost::filesystem::basename(getenv("PWD")));
				string jobdir_prefix = "CHEOPSim_job";
				size_t found = jobdir.find(jobdir_prefix);
				if (found!=string::npos) {
					string jobid = jobdir.substr(found+jobdir_prefix.size());
					//check that the extracted job ID is an integer
					bool jobid_isValid = true;
					for (unsigned i=0; i<jobid.length(); i++) {
						if (!isdigit(jobid[i])) jobid_isValid = false;
					}
					if (jobid_isValid) m_obsId = boost::lexical_cast<uint32_t>(jobid);
				}
			}

		};

		int m_progType; ///< Program Type
		int m_progId; ///< Program ID
		int m_reqId; ///< Request ID
		int m_obsId; ///< Observation ID
		int m_visitCtr; ///< Visit counters
		int m_versionNum; ///< Version number
		string m_obsCat; ///< Observation category
		string m_piName; ///< Name of the PI of the observing program
		int m_piUid; ///< Account ID of the PI at UGE
		int m_prpFirst; ///< Proprietary period, depending on first visit
		int m_prpLast; ///< Proprietary period, depending on last visit
		string m_targName; ///< Name of the target as provided by the proposal
		string m_readMode; ///< readout mode of the CCD:  "faint", "faint fast", "bright" or "ultrabright"
		string m_marginMode; ///< margin mode of the CCD: "image", "reduced" or "total collapsed"
		double m_gaiaMagError; ///< Error on the Gaia Magnitude, as specified in the MPS_PRE_Visits file
		double m_strayLightThreshold; ///< Stray light threshold in photons/s/pixel

	};

	/// @brief Struct to contain subarray dimensions
	struct SubarrayDimensions {
		/** *************************************************************************
		 *  @brief Constructor for SubarrayDimensions struct
		 *
		 *  @param [in] xDim  Number of pixels for x dimension (default 200)
		 *  @param [in] yDim  Number of pixels for y dimension (default 200)
		 *  @param [in] xOffset  Offset in x direction of the left edge of the image
		 *  					 w.r.t. the left edge of the exposed part of the CCD
		 *  					 (default 412)
		 *  @param [in] yOffset  Offset in y direction of the bottom edge of the image
		 *  					 w.r.t. the bottom edge of the exposed part of the CCD
		 *  					 (default 412)
		 *  @param [in] targetLocationX  Location of centre of FOV rotation in x direction
		 *  					 w.r.t. the bottom edge of the exposed part of the CCD
		 *  					 (default 512)
		 *  @param [in] targetLocationY  Location of centre of FOV rotation in y direction
		 *  					 w.r.t. the bottom edge of the exposed part of the CCD
		 *  					 (default 512)
		 */
		SubarrayDimensions(int xDim=200, int yDim=200, int xOffset=412, int yOffset=412,
						   double targetLocationX=512, double targetLocationY=512) :
									m_xDim(xDim),m_yDim(yDim),
									m_xOffset(xOffset),m_yOffset(yOffset),
									m_targetLocationX(targetLocationX), m_targetLocationY(targetLocationY){};
		int m_xDim; ///< Number of pixels in X dimension
		int m_yDim; ///< Number of pixels in Y dimension
		int m_xOffset; ///< Offset of first pixel in X dimension w.r.t. CCD edge
		int m_yOffset; ///< Offset of first pixel in Y dimension w.r.t. CCD edge
		double m_targetLocationX; ///< Location of centre of FOV rotation in X dimension w.r.t. CCD edge
		double m_targetLocationY; ///< Location of centre of FOV rotation in Y dimension w.r.t. CCD edge
	};

	/// @brief Struct to contain parameters for generating the ideal light curve
	struct IdealLightCurveParams {
		/// @brief Constructor for IdealLightCurveParams struct
		IdealLightCurveParams() : m_photonNoise(false), m_biasOffset(0.), m_readNoise(0.) {};
		bool m_photonNoise; ///< Boolean to indicate whether or not PhotonNoiseGenerator is run
		double m_biasOffset; ///< Value of biasMean parameter in BiasGenerator
		double m_readNoise; ///< Value of biasWidth parameter in BiasGenerator
	};

	/// @brief Struct to contain parameters for photometric extraction
	struct PhotometryParams {
		/** *************************************************************************
		 *  @brief Constructor for PhotometryParams struct
		 *
		 *  @param [in] radius_barycentre  Radius of circle centred on image centre used to determine PSF barycentre
		 *  @param [in] radius_bkgInner  Inner radius of annulus centred on image centre used to evaluate the background
		 *  @param [in] radius_bkgOuter  Outer radius of annulus centred on image centre used to evaluate the background
		 *  @param [in] radius_psf  Radius of circle centred on PSF barycentre used to extract the signal flux
		 *  @param [in] subtractBackground  Flag to indicate whether or not to perform background subtraction
		 */
		PhotometryParams(double radius_barycentre=0., double radius_bkgInner=0., double radius_bkgOuter=0., double radius_psf=0, bool subtractBackground=true) :
			m_radius_barycentre(radius_barycentre), m_radius_bkgInner(radius_bkgInner),
			m_radius_bkgOuter(radius_bkgOuter), m_radius_psf(radius_psf),
			m_subtractBackground(subtractBackground) {};
		double m_radius_barycentre; ///< Radius of circle centred on image centre used to determine PSF barycentre
		double m_radius_bkgInner; ///< Inner radius of annulus centred on image centre used to evaluate the background
		double m_radius_bkgOuter; ///< Inner radius of annulus centred on image centre used to evaluate the background
		double m_radius_psf; ///< Radius of circle centred on PSF barycentre used to extract the signal flux
		bool m_subtractBackground; ///< Flag to indicate whether or not to perform background subtraction
	};

	/// @brief Struct to contain voltages and temperatures for HK data
	struct HKData {
		/// @brief Default constructor
		HKData() : m_ccdTemp(SatelliteData::kDefaultCcdTemperature), m_biasTemp(SatelliteData::kDefaultFeeBiasTemperature), m_adcTemp(SatelliteData::kDefaultFeeAdcTemperature),
				m_vod(0.), m_vrd(0.), m_vog(0.), m_vss(0.) {};
		/// @brief Constructor for HKData struct
		HKData(double ccdTemp, double biasTemp, double adcTemp, double * voltages) : m_ccdTemp(ccdTemp), m_biasTemp(biasTemp), m_adcTemp(adcTemp),
				m_vod(voltages[SatelliteData::VOD]), m_vrd(voltages[SatelliteData::VRD]), m_vog(voltages[SatelliteData::VOG]), m_vss(voltages[SatelliteData::VSS]) {};
		double m_ccdTemp; ///< CCD temperature in Kelvin
		double m_biasTemp; ///< Bias temperature in Kelvin
		double m_adcTemp; ///< ADC temperature in Kelvin
		double m_vod; ///< Output drain level voltage
		double m_vrd; ///< Reset drain level voltage
		double m_vog; ///< Output gate level voltage
		double m_vss; ///< Substrate level voltage
	};

	/// @brief Type of image: stacked, unstacked, imagette or full frame
	enum IMAGETYPE{stacked,unstacked,imagette,fullframe};

	/** *************************************************************************
	 *  @brief Constructor. Initializes all pointers to nullptr
	 *
	 *  @param [in] timeConf  Time configuration for the simulation
	 *  @param [in] modulesToRun  Comma separated list of modules to be run
	 */
	Data(const TimeConfiguration & timeConf, string modulesToRun);

	virtual ~Data();

	/** *************************************************************************
	 *  @brief Returns whether or not the specified module is being run in the configuration
	 */
	bool moduleIsActive(string moduleName);

	/** *************************************************************************
	 *  @brief Resets the vectors of stacked and unstacked images held in memory
	 */
	void clearImages();

	/** *************************************************************************
	 *  @brief Returns the name of the output directory
	 */
	string getOutputDirectory() const;

	/** *************************************************************************
	 *  @brief Set the time configuration for the simulation
	 */
	void setTimeConfiguration(const TimeConfiguration & timeConf) {m_timeConf = timeConf;}

	/** *************************************************************************
	 *  @brief Sets the number of time steps (repetition periods) between the visit start time
	 *  and the start of the first sub-array image in the time configuration.
	 */
	void setPreSubArrayTimeSteps() {m_timeConf.setPreSubArrayTimeSteps();}

	/** *************************************************************************
	 *  @brief Set the field of view
	 */
	void setFieldOfView(SkyFieldOfView * fov) {if (m_fov!=nullptr) delete m_fov; m_fov = fov;}

	/** *************************************************************************
	 *  @brief Set the wavelength dependence information
	 */
	void setWavelengthDependence(WavelengthDependence * wavelengthDependence) {m_wavelengthDependence = wavelengthDependence;}

	/** *************************************************************************
	 *  @brief Set the temperature range corresponding to breathing of the PSF
	 *
	 *  @param [in] min  lowest value for telescope temperature,
	 *  				 corresponding to the hot1 thermal map for the PSF
	 *  @param [in] max  highest value for telescope temperature,
	 *  				 corresponding to the hot2 thermal map for the PSF
	 */
	void setBreathingTemperatureRange(double min, double max) {m_breathingTemperatureRange = make_pair(min,max);}

	/** *************************************************************************
	 *  @brief Add a telegraphic pixel to the list
	 */
	void addTelegraphicPixel(TelegraphicPixel telegraphicPixel) {m_telegraphicPixels->push_back(telegraphicPixel);}

	/** *************************************************************************
	 *  @brief Set the FITS table containing mission planning visit information (output)
	 */
	void setMpsPreVisits(MpsPreVisits * mpsVisits) {m_MPSVisits = mpsVisits;}

	/** *************************************************************************
	 *  @brief Set the FITS table containing mission planning visit constraints (input)
	 */
	void setMpsVisitConstraints(MpsPreVisitconstraints * mpsVisitConstraints) {m_mpsVisitConstraints = mpsVisitConstraints;}

	/** *************************************************************************
	 *  @brief Add an unstacked image to the list held in memory
	 */
	void addImage(Image * image) {m_images.push_back(image);}

	/** *************************************************************************
	 *  @brief Add an stacked image to the list held in memory
	 */
	void addStackedImage(Image * image) {m_stackedImages.push_back(image);}

	/** *************************************************************************
	 *  @brief Add an unstacked imagette to the list held in memory
	 */
	void addImagette(Image * image) {m_imagettes.push_back(image);}

	/** *************************************************************************
	 *  @brief Add an stacked imagette to the list held in memory
	 */
	void addStackedImagette(Image * image) {m_stackedImagettes.push_back(image);}

	/** *************************************************************************
	 *  @brief Add an image corresponding to the last second of the current exposure,
	 *  	   used to define frame transfer smear trails in the upward direction
	 */
	void addImageToSmearUp(Image * image) {m_imagesToSmearUp.push_back(image);}

	/** *************************************************************************
	 *  @brief Add an image corresponding to the first second of the current exposure,
	 *  	   used to define frame transfer smear trails in the downward direction
	 */
	void addImageToSmearDown(Image * image) {m_imagesToSmearDown.push_back(image);}

	/** *************************************************************************
	 *  @brief Set the cosmic ray truth image
	 */
	void setCosmicRayImage(Image * image_cr) {m_cosmicRayImage = image_cr;}

	/** *************************************************************************
	 *  @brief Set the flat field image, combined over wavelength
	 */
	void setFlatField(Image * flatField) {m_flatField = flatField;}

	/** *************************************************************************
	 *  @brief Set the (optionally imperfect) flat field image
	 *  	   to be used for flat field correction
	 */
	void setFlatFieldToSubtract(Image * flatField) {m_flatFieldToSubtract = flatField;}

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cube to contain imagettes
	 *  	   and to the corresponding metadata
	 *
	 *  @param [in] imageCube  pointer to the imagette image cube
	 *  @param [in] metaData  pointer to the imagette metadata
	 */
	void setFitsImageCube(SciRawImagette* imageCube, SciRawImagettemetadata* metaData);

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cube to contain stacked images
	 *  	   (32-bit integers) and to the corresponding metadata
	 *
	 *  @param [in] imageCube  pointer to the image cube for stacked images
	 *  @param [in] metaData  pointer to the metadata for stacked images
	 *  @param [in] metaData_unstacked  pointer to the unstacked metadata for stacked images
	 */
	void setFitsImageCube(SciRawSubarray* imageCube, SciRawImagemetadata* metaData, SciRawUnstackedimagemetadata* metaData_unstacked);

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cube to contain unstacked images
	 *  	   (16-bit integers) and to the corresponding metadata
	 *
	 *  @param [in] imageCube  pointer to the image cube for unstacked images
	 *  @param [in] metaData  pointer to the metadata for unstacked images
	 */
	void setFitsImageCube(SimRawUnstackedsubarray* imageCube, SciRawImagemetadata* metaData);

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cube to contain
	 *  	   double precision stacked images and to the corresponding metadata
	 *
	 *  @param [in] imageCube  pointer to the image cube for double precision stacked images
	 *  @param [in] metaData  pointer to the metadata for double precision stacked images
	 *  @param [in] metaData_unstacked  pointer to the unstacked metadata for double precision stacked images
	 */
	void setFitsImageCube(SimRawDoubleprecisionsubarray* imageCube, SciRawImagemetadata* metaData, SciRawUnstackedimagemetadata* metaData_unstacked);

	/** *************************************************************************
	 *  @brief Set the pointer to the FITS image cube to contain
	 *  	   stacked images with smear trails only
	 *
	 *  @param [in] imageCube  pointer to the image cube for stacked images with smear trails only
	 */
	void setFitsImageCube_smear(SimRawDoubleprecisionsubarray* imageCube);

	/** *************************************************************************
	 *  @brief Set the pointer to the FITS image cube to contain
	 *  	   stacked images with cosmic rays only
	 *
	 *  @param [in] imageCube  pointer to the image cube for stacked images with cosmic rays only
	 */
	void setFitsImageCube_cosmics(SimRawDoubleprecisionsubarray* imageCube);

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cubes
	 *  	   for each type of CCD margin, for stacked images (32-bit integers)
	 */
	void setFitsCcdMargins(SciRawBlankleft * blankLeftMarginCube,
						   SciRawBlankright * blankRightMarginCube,
						   SciRawDarkleft * darkLeftMarginCube,
						   SciRawDarkright * darkRightMarginCube,
						   SciRawDarktop * darkTopMarginCube,
						   SciRawOverscanleft * overscanLeftMarginCube,
						   SciRawOverscanright * overscanRightMarginCube,
						   SciRawOverscantop * overscanTopMarginCube);

	/** *************************************************************************
	 *  @brief Set the pointers to the FITS image cubes
	 *  	   for each type of CCD margin, for unstacked images (16-bit integers)
	 */
	void setFitsCcdMarginImages(SimRawUnstackedblankleftimage * blankLeftImageCube,
							    SimRawUnstackedblankrightimage * blankRightImageCube,
							    SimRawUnstackeddarkleftimage * darkLeftImageCube,
							    SimRawUnstackeddarkrightimage * darkRightImageCube,
							    SimRawUnstackeddarktopimage * darkTopImageCube,
							    SimRawUnstackedoverscanleftimage * overscanLeftImageCube,
							    SimRawUnstackedoverscanrightimage * overscanRightImageCube,
							    SimRawUnstackedoverscantopimage * overscanTopImageCube);

	/** *************************************************************************
	 *  @brief Set the pointer to the FITS table storing the centroid positions
	 */
	void setFitsCentroid(SciRawCentroid * centroid) {m_centroid = centroid;}

	/** *************************************************************************
	 *  @brief Set the pointer to the FITS table storing truth data
	 */
	void setTruthMetaData(SimTruSubarray * truthMetaData, IMAGETYPE type);

	/** *************************************************************************
	 *  @brief Get the visit data
	 */
	Visit getVisit() const {return m_visit;}

	/** *************************************************************************
	 *  @brief Set the visit data
	 */
	void setVisit(const Visit visit) {m_visit = visit;}

	/** *************************************************************************
	 *  @brief Return the time configuration for the simulation
	 */
	TimeConfiguration getTimeConfiguration() const {return m_timeConf;}

	/** *************************************************************************
	 *  @brief Return the field of view
	 */
	SkyFieldOfView * getFieldOfView() const;

	/** *************************************************************************
	 *  @brief Return the wavelength dependence information
	 */
	WavelengthDependence * getWavelengthDependence() {return m_wavelengthDependence;}

	/** *************************************************************************
	 *  @brief Set the temperature range corresponding to breathing of the PSF
	 *
	 *  @return the lowest and highest values for telescope temperature,
	 *  		corresponding to the hot1 and hot2 thermal maps for the PSF, respectively
	 */
	pair<double,double> getBreathingTemperatureRange() {return m_breathingTemperatureRange;}

	/** *************************************************************************
	 *  @brief Get the satellite data for a specified time step
	 *
	 *  @param [in] timeStep  time step of the simulation
	 *  @param [in] create 	 boolean to indicate whether or not to create a new Satellite data instance
	 *  @param [in] checkTimeStep 	 boolean to indicate whether or not compare the timeStep against the size of the SatelliteData vector
	 */
	SatelliteData* getSatelliteData(int timeStep, bool create=false, bool checkTimeStep=true);

	/** *************************************************************************
	 *  @brief Get the list of telegraphic pixels
	 */
	vector<TelegraphicPixel> * getTelegraphicPixels() {return m_telegraphicPixels;}

	/** *************************************************************************
	 *  @brief Get the FITS table containing mission planning visit information (output)
	 */
	MpsPreVisits* getMpsPreVisits() const {return m_MPSVisits;}

	/** *************************************************************************
	 *  @brief Get the FITS table containing mission planning visit constraints (input)
	 */
	MpsPreVisitconstraints* getMpsVisitConstraints() const {return m_mpsVisitConstraints;}

	/** *************************************************************************
	 *  @brief Get the list of unstacked images held in memory
	 */
	vector<Image*> getImages() const {return m_images;}

	/** *************************************************************************
	 *  @brief Get the list of stacked images held in memory
	 */
	vector<Image*> getStackedImages() const {return m_stackedImages;}

	/** *************************************************************************
	 *  @brief Get the list of unstacked imagettes held in memory
	 */
	vector<Image*> getImagettes() const {return m_imagettes;}

	/** *************************************************************************
	 *  @brief Get the list of stacked imagettes held in memory
	 */
	vector<Image*> getStackedImagettes() const {return m_stackedImagettes;}

	/** *************************************************************************
	 *  @brief Get the set of images corresponding to the last second of the exposure,
	 *  	   used to define frame transfer smear trails in the upward direction
	 */
	vector<Image*> getImagesToSmearUp() const {return m_imagesToSmearUp;}

	/** *************************************************************************
	 *  @brief Get the set of images corresponding to the first second of each exposure,
	 *  	   used to define frame transfer smear trails in the downward direction
	 */
	vector<Image*> getImagesToSmearDown() const {return m_imagesToSmearDown;}

	/** *************************************************************************
	 *  @brief Get the cosmic ray truth image
	 */
	Image* getCosmicRayImage() const {return m_cosmicRayImage;}

	/** *************************************************************************
	 *  @brief Get the flat field image, combined over wavelength
	 */
	Image* getFlatField() const {return m_flatField;}

	/** *************************************************************************
	 *  @brief Get the (optionally imperfect) flat field image
	 *  	   to be used for flat field correction
	 */
	Image* getFlatFieldToSubtract() const {return m_flatFieldToSubtract;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing imagettes
	 */
	void getFitsImageCube(SciRawImagette ** imageCube) const {*imageCube = m_imagetteCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing stacked images
	 *  	   (32-bit integers)
	 */
	void getFitsImageCube(SciRawSubarray ** imageCube) const {*imageCube = m_stackedImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing unstacked images
	 *  	   (16-bit integers)
	 */
	void getFitsImageCube(SimRawUnstackedsubarray ** imageCube) const {*imageCube = m_unstackedImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing double precision
	 *  	   stacked images
	 */
	void getFitsImageCube(SimRawDoubleprecisionsubarray ** imageCube) const {*imageCube = m_stackedImageCube_DP;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing stacked images
	 *  	   with smear trails only
	 */
	void getFitsImageCube_smear(SimRawDoubleprecisionsubarray ** imageCube) const {*imageCube = m_stackedImageCube_smear;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube containing stacked images
	 *  	   with cosmic rays only
	 */
	void getFitsImageCube_cosmics(SimRawDoubleprecisionsubarray ** imageCube) const {*imageCube = m_stackedImageCube_cosmics;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left blank columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawBlankleft ** blankLeftMarginCube) {*blankLeftMarginCube = m_stackedBlankLeftMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right blank columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawBlankright ** blankRightMarginCube) {*blankRightMarginCube = m_stackedBlankRightMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left dark columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawDarkleft ** darkLeftMarginCube) {*darkLeftMarginCube = m_stackedDarkLeftMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right dark columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawDarkright ** darkRightMarginCube) {*darkRightMarginCube = m_stackedDarkRightMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the top dark rows,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawDarktop ** darkTopMarginCube) {*darkTopMarginCube = m_stackedDarkTopMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left overscan columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawOverscanleft ** overscanLeftMarginCube) {*overscanLeftMarginCube = m_stackedOverscanLeftMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right overscan columns,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawOverscanright ** overscanRightMarginCube) {*overscanRightMarginCube = m_stackedOverscanRightMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the top overscan rows,
	 *  	   for stacked images (32-bit integers)
	 */
	void getFitsCcdMargin(SciRawOverscantop ** overscanTopMarginCube) {*overscanTopMarginCube = m_stackedOverscanTopMarginCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left blank columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackedblankleftimage ** blankLeftImageCube) {*blankLeftImageCube = m_unstackedBlankLeftImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right blank columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackedblankrightimage ** blankRightImageCube) {*blankRightImageCube = m_unstackedBlankRightImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left dark columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackeddarkleftimage ** darkLeftImageCube) {*darkLeftImageCube = m_unstackedDarkLeftImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right dark columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackeddarkrightimage ** darkRightImageCube) {*darkRightImageCube = m_unstackedDarkRightImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the top dark rows,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackeddarktopimage ** darkTopImageCube) {*darkTopImageCube = m_unstackedDarkTopImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the left overscan columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackedoverscanleftimage ** overscanLeftImageCube) {*overscanLeftImageCube = m_unstackedOverscanLeftImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the right overscan columns,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackedoverscanrightimage ** overscanRightImageCube) {*overscanRightImageCube = m_unstackedOverscanRightImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS image cube for the top overscan rows,
	 *  	   for unstacked images (16-bit integers)
	 */
	void getFitsCcdMargin(SimRawUnstackedoverscantopimage ** overscanTopImageCube) {*overscanTopImageCube = m_unstackedOverscanTopImageCube;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS table containing metadata
	 *  	   for the specified image type (stacked, unstacked or fullframe)
	 */
	SciRawImagemetadata * getFitsImageMetaData(IMAGETYPE type);

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS table containing unstacked metadata
	 *  	   for the specified image type (stacked or fullframe)
	 */
	SciRawUnstackedimagemetadata * getFitsUnstackedImageMetaData(IMAGETYPE type);

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS table containing metadata for imagettes
	 */
	SciRawImagettemetadata * getFitsImagetteMetaData() {return m_imagetteMetaData;};

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS table storing the centroid positions
	 */
	SciRawCentroid * getFitsCentroid() {return m_centroid;}

	/** *************************************************************************
	 *  @brief Get the pointer to the FITS table storing truth data
	 */
	SimTruSubarray * getFitsTruthMetaData(IMAGETYPE type);

	/** *************************************************************************
	 *  @brief Initializes the REF_APP_BadPixelMap data structure and its
	 *  	   margin extensions.
	 *
	 *	@param [in] reference_file   name of flat field or dark current
	 *								 reference file (depending on the calling
	 *								 module) used to set the corresponding
	 *								 header keyword
	 */
	void initializeBadPixelMap(string reference_file);

	/** *************************************************************************
	 *  @brief Assigns the value to a pixel in the appropriate extension of the
	 *  	   REF_APP_BadPixelMap data structure given the pixel location in
	 *  	   global CCD coordinates (x=0 at left edge of left margin)
	 *
	 *  @param [in] i  x-index of pixel position in global CCD coordinates
	 *  @param [in] j  y-index of pixel position in global CCD coordinates
	 *  @param [in] value  bad pixel value: 0=good, 1=hot (e-/s>5), 3=telegraphic,
	 *  				   -1=partially dead (flat field < 0.8)
	 */
	void setBadPixel(int i, int j, int value);

	/** *************************************************************************
	 *  @brief Delete all pointers to FITS image cubes and tables, setting their
	 *  	   values to nullptr
	 */
	void resetFitsImageCube(IMAGETYPE type);

	/** *************************************************************************
	 *  @brief Get the number of stacked images per image cube
	 */
	unsigned getNumberOfImagesPerCube(IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Set the number of stacked images for which the payload is in the loop
	 */
	void setNumberOfStackedImagesPITL(unsigned numberOfStackedImages_PITL) {m_timeConf.setNumberOfStackedImagesPITL(numberOfStackedImages_PITL);}

	/** *************************************************************************
	 *  @brief Determine and set the imagette stacking number
	 */
	void setImagetteStackingNumber() {m_timeConf.setImagetteStackingNumber();}

	/** *************************************************************************
	 *  @brief Returns a boolean to indicate whether or not any satellite data
	 *  	   has been defined
	 */
	bool hasSatelliteData() const {return (m_satelliteData.size()!=0);}

	/** *************************************************************************
	 *  @brief Increments the counter for the number of stacked images
	 *  	   that have been written out
	 */
	void incrementStackedImageCount() {m_stackedImageCount++;}

	/** *************************************************************************
	 *  @brief Returns the value of the counter for the number of stacked images
	 *  	   that have been written out
	 */
	unsigned stackedImageCount() const {return m_stackedImageCount;}

	/** *************************************************************************
	 *  @brief Sets whether or not frame transfer smearing is switched on
	 */
	void setFrameTransferSmearing() {m_frameTransferSmear = true;}

	/** *************************************************************************
	 *  @brief Returns whether or not frame transfer smearing is switched on
	 */
	bool doFrameTransferSmearing() const {return m_frameTransferSmear;}

	/** *************************************************************************
	 *  @brief Sets whether or not frame transfer smearing is switched on
	 */
	void setChargeTransferEoL() {m_chargeTransferEol = true;}

	/** *************************************************************************
	 *  @brief Returns whether or not frame transfer smearing is switched on
	 */
	bool doChargeTransferEoL() const {return m_chargeTransferEol;}

	/** *************************************************************************
	 *  @brief Returns whether or not the sub-array images are preceded by a full frame image
	 */
	bool doFullFrame() const {return m_doFullFrame;}

	/** *************************************************************************
	 *  @brief Sets the readout time
	 */
	void setReadoutTime(double readoutTime) {m_readoutTime = readoutTime;}

	/** *************************************************************************
	 *  @brief Returns the readout time
	 */
	double getReadoutTime() const {return m_readoutTime;}

	/** *************************************************************************
	 *  @brief Sets the bias offset of the electronic readout
	 */
	void setBiasOffset(double offset) {m_biasOffset = offset;}

	/** *************************************************************************
	 *  @brief Returns the bias offset of the electronic readout
	 */
	double getBiasOffset() const {return m_biasOffset;}

	/** *************************************************************************
	 *  @brief Sets the name of the file containing the electronic gain formula
	 */
	void setGainFilename(string filename) {m_gainFilename = filename;}

	/** *************************************************************************
	 *  @brief Returns the name of the file containing the electronic gain formula
	 */
	string getGainFilename() const {return m_gainFilename;}

	/** *************************************************************************
	 *  @brief Extracts the nominal values for each of the four VOLTAGE_TYPE bias
	 *  	   voltages and the nominal CCD temperature from the REF_APP_GainCorrection file
	 */
	void setNominalVoltages(string gainFilename);

	/** *************************************************************************
	 *  @brief Returns the nominal voltage for the specified bias voltage type
	 */
	double getNominalVoltage(SatelliteData::VOLTAGE_TYPE voltageType) const {return m_nominalVoltage[voltageType];}

	/** *************************************************************************
	 *  @brief Return the specified bias voltage at the specified temperature
	 *  	   using temperature-voltage relationship defined in
	 *  	   CHEOPS-DLR-INST-TR-038 Issue 1.0, Table 13
	 *
	 *  @param [in] temperature  FEE bias temperature in Kelvin
	 *  @param [in] type  Bias voltage type
	 */
	double getTemperatureDependentVoltage(double temperature, SatelliteData::VOLTAGE_TYPE type) const;

	/** *************************************************************************
	 *  @brief Return the specified bias voltage at the specified time,
	 *  	   without Gaussian fluctuation, but taking into account the
	 *  	   drift of the voltage with time defined in Tables 5-3 and 5-4 of
	 *  	   CHEOPS-UGE-SYS-TR-020 Issue 1.0 (payload calibration report)
	 *
	 *  @param [in] timeSinceVisitStart  time since the start of the visit in seconds
	 *  @param [in] type  Bias voltage type
	 */
	double getVoltageWithDrift(double timeSinceVisitStart, SatelliteData::VOLTAGE_TYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the flag to omit images when the spacecraft is within the SAA
	 *  	   to true
	 */
	void setOmitSAA() {m_omitSAA = true;}

	/** *************************************************************************
	 *  @brief Returns true if images should be omitted when the spacecraft is
	 *  	   within the SAA
	 */
	bool omitSAA() const {return m_omitSAA;}

	/** *************************************************************************
	 *  @brief Sets the flag to omit images when the target is occulted by the
	 *  	   Earth to true
	 */
	void setOmitEarthOccultation() {m_omitEarthOccultation = true;}

	/** *************************************************************************
	 *  @brief Returns true if images should be omitted when the target is
	 *  	   occulted by the Earth
	 */
	bool omitEarthOccultation() const {return m_omitEarthOccultation;}

	/** *************************************************************************
	 *  @brief Sets the flag to omit images when the stray light is above
	 *  	   threshold to true
	 */
	void setOmitStrayLight() {m_omitStrayLight = true;}

	/** *************************************************************************
	 *  @brief Returns true if images should be omitted when the stray light is
	 *  	   above threshold
	 */
	bool omitStrayLight() const {return m_omitStrayLight;}

	/** *************************************************************************
	 *  @brief Sets the number of stacked images to be discarded at the start of the simulation
	 */
	void setNumberOfCEsToOmitAtStart(unsigned numberOfCEsToOmitAtStart) {m_numberOfCEsToOmitAtStart = numberOfCEsToOmitAtStart;}

	/** *************************************************************************
	 *  @brief Returns the number of stacked images to be discarded at the start of the simulation
	 */
	unsigned getNumberOfCEsToOmitAtStart() const {return m_numberOfCEsToOmitAtStart;}

	/** *************************************************************************
	 *  @brief sets the value is N, in the case that only every Nth stacked image
	 *  is to be written out. All images are written out if N is 0 or 1.
	 */
	void setOnlyWriteEveryNthCE(unsigned onlyWriteEveryNthCE) {m_onlyWriteEveryNthCE = onlyWriteEveryNthCE==0?1:onlyWriteEveryNthCE;}

	/** *************************************************************************
	 *  @brief Returns the value is N, in the case that only every Nth stacked image
	 *  is to be written out. All images are written out if N is 0 or 1.
	 */
	unsigned getOnlyWriteEveryNthCE() const {return m_onlyWriteEveryNthCE;}

	/** *************************************************************************
	 *  @brief Increments the number of images which have been discarded due to
	 *  	   Earth occultation or SAA
	 */
	void incrementDiscardedImages() {m_discardedImages++;}

	/** *************************************************************************
	 *  @brief Returns the number of images which have been discarded due to
	 *  	   Earth occultation or SAA
	 */
	unsigned getNumberOfDiscardedImages() const {return m_discardedImages;}

	/** *************************************************************************
	 *  @brief Sets the image sub-array dimensions
	 */
	void setSubarrayDimensions(int xDim, int yDim, int xOffset, int yOffset, double targetLocationX, double targetLocationY);

	/** *************************************************************************
	 *  @brief Sets the target location
	 */
	void setTargetLocation(double targetLocationX, double targetLocationY);

	/** *************************************************************************
	 *  @brief Returns the image sub-array dimensions
	 */
	SubarrayDimensions getSubarrayDimensions() const {return m_subarrayDimensions;}

	/** *************************************************************************
	 *  @brief Sets the parameters for photometric extraction
	 */
	void setPhotometryParams(double radius_barycentre, double radius_bkgInner, double radius_bkgOuter, double radius_psf, bool subtractBackground);

	/** *************************************************************************
	 *  @brief Returns the parameters for photometric extraction
	 */
	PhotometryParams getPhotometryParams() const {return m_photometryParams;}

	/** *************************************************************************
	 *  @brief Returns the parameters for for generating the ideal light curve
	 */
	IdealLightCurveParams * getIdealLightCurveParams() {return m_idealLightCurveParams;}

	/** *************************************************************************
	 *  @brief Add a point to the ideal light curve
	 */
	void appendIdealLightCurve(double flux) {m_idealLightCurve.push_back(flux);}

	/** *************************************************************************
	 *  @brief Returns the ideal light curve
	 */
	vector<double> getIdealLightCurve() const {return m_idealLightCurve;}

	/** *************************************************************************
	 *  @brief Sets the flag to indicate whether or not to use the redundant
	 *  	   right amplifier rather than left amplifier for CCD readout
	 */
	void setRedundantHardware() {m_redundantHardware = true;}

	/** *************************************************************************
	 *  @brief Returns true if the redundant right amplifier is used rather than
	 *  	   the left amplifier for CCD readout
	 */
	bool redundantHardware() const {return m_redundantHardware;}

	/** *************************************************************************
	 *  @brief Sets the flag to indicate whether or not gain non-linearity has been simulated
	 */
	void setGainNonLinearity() {m_gainNonLinearity = true;}

	/** *************************************************************************
	 *  @brief Returns true if gain non-linearity has been simulated
	 */
	bool gainNonLinearity() const {return m_gainNonLinearity;}

	/** *************************************************************************
	 *  @brief Sets the serial readout rate
	 *
	 *  @param [in] readRate  serial read rate in kHz
	 */
	void setSerialReadRate(double readRate) {m_serialReadRate = readRate;}

	/** *************************************************************************
	 *  @brief Returns the serial readout rate in kHz
	 */
	double getSerialReadRate() const {return m_serialReadRate;}

	/** *************************************************************************
	 *  @brief Appends the vector of raw HK data recorded at 1.2s intervals
	 *
	 *  @param [in] hkData  instance of the HKData struct
	 */
	void appendRawHKData(HKData hkData) {m_rawHKData.push_back(hkData);}

	/** *************************************************************************
	 *  @brief Calculates the average HK data at 20s intervals, taking the mean
	 *  	   over the previous 16 raw measurements (raw cadence is 1.2s)
	 */
	void calculateHKAverages();

	/** *************************************************************************
	 *  @brief Returns the raw HK data closest in time to the specified time since
	 *  	   the start of the visit
	 */
	HKData getClosestRawHKData(double timeSinceVisitStart) const;

	/** *************************************************************************
	 *  @brief Returns the HK data averaged over the previous 20s, for the time
	 *         closest to the specified time since the start of the visit
	 */
	HKData getClosestAveragedHKData(double timeSinceVisitStart) const;

private:

	//Voltage temperature coefficients taken from CHEOPS-DLR-INST-TR-038 Issue 1.0 table 13
	static constexpr double kVODTempCoeff =55.; ///< VOD temperature coefficient [microvolts/K]
	static constexpr double kVRDTempCoeff =-12.; ///< VRD temperature coefficient [microvolts/K]
	static constexpr double kVOGTempCoeff =-0.8; ///< VOG temperature coefficient [microvolts/K] //slope at -10C (nominal FEE bias temperature) estimated from Figure 20 rather than Table 13
	static constexpr double kVSSTempCoeff =-13.; ///< VSS temperature coefficient [microvolts/K]

	/** *************************************************************************
	 *  @brief Set the values of the pointers to FITS cubes containing
	 *  	   CCD margin data for stacked images to nullptr
	 */
	void initializeStackedCcdMargins();

	/** *************************************************************************
	 *  @brief Set the values of the pointers to FITS cubes containing
	 *  	   CCD margin data for unstacked images to nullptr
	 */
	void initializeUnstackedCcdMargins();

	/** *************************************************************************
	 *  @brief Delete the pointers to FITS cubes containing
	 *  	   CCD margin data for stacked images
	 */
	void deleteStackedCcdMargins();

	/** *************************************************************************
	 *  @brief Delete the pointers to FITS cubes containing
	 *  	   CCD margin data for unstacked images
	 */
	void deleteUnstackedCcdMargins();

	string m_modulesToRun; ///< List of modules being run as a comma separated string
	Visit m_visit; ///< Visit data for the simulation
	TimeConfiguration m_timeConf;  ///< Time configuration for the simulation
	SkyFieldOfView * m_fov;  ///< Field of view
	WavelengthDependence * m_wavelengthDependence;  ///< Wavelength dependence information
	SubarrayDimensions m_subarrayDimensions;  ///< Dimensions of the image sub-array
	IdealLightCurveParams * m_idealLightCurveParams;  ///< Parameters for photometric extraction
	PhotometryParams m_photometryParams;  ///< Parameters for photometric extraction
	vector<double> m_idealLightCurve;  ///< Ideal light curve flux time series vector
	pair<double,double> m_breathingTemperatureRange;  ///< minimum and maximum telescope temperatures, corresponding to the hot1 and hot2 thermal maps for the PSF
	vector<SatelliteData*> m_satelliteData;  ///< Satellite data
	vector<TelegraphicPixel> * m_telegraphicPixels;  ///< List of telegraphic pixels
	vector<Image*> m_images;  ///< List of stacked images held in memory
	vector<Image*> m_stackedImages;  ///< List of unstacked images held in memory
	unsigned m_stackedImageCount;  ///< Counter for the number of stacked images written out
	bool m_frameTransferSmear; ///< Boolean to indicate whether or not frame transfer smearing is switched on
	bool m_chargeTransferEol; ///< Boolean to indicate whether or not charge transfer simulation in end of life mode is switched on
	bool m_doFullFrame; ///< Boolean to indicate whether or not the sub-array images are preceded by a full frame image
	double m_readoutTime; ///< Time to readout one exposure
	double m_biasOffset; ///< Bias Offset for the electronic readout readout, defined in BiasGenerator and used to set PIX_DATA_OFFSET in the image metadata
	string m_gainFilename; ///< Name of the file containing the electronic gain formula
	double m_nominalVoltage[5]; ///< nominal values for each of the four VOLTAGE_TYPE bias voltages, and the nominal CCD temperature
	bool m_omitSAA; ///< Flag to indicate whether or not to omit images when the spacecraft is within the SAA
	bool m_omitEarthOccultation; ///< Flag to indicate whether or not to omit images when the target is occulted by the Earth
	bool m_omitStrayLight; ///< Flag to indicate whether or not to omit images when the stray light is above threshold
	unsigned m_numberOfCEsToOmitAtStart; ///< Number of stacked images to be discarded at the start of the simulation
	unsigned m_onlyWriteEveryNthCE; ///< If the value is N, only every Nth stacked image will be written out. All images are written out if the value is 0 or 1.
	unsigned m_discardedImages; ///< Number of stacked images which have been discarded due to SAA, Earth occultation, option to omit first N images, or option to only write out every Nth image
	bool m_redundantHardware; ///< Flag to indicate whether or not to use the redundant right amplifier rather than left amplifier for CCD readout
	bool m_gainNonLinearity; ///< Flag to indicate whether or not gain non-linearity has been simulated
	double m_serialReadRate; ///< Serial read rate in kHz

	MpsPreVisits * m_MPSVisits; ///< Output FITS table containing mission planning visit information (output)
	MpsPreVisitconstraints * m_mpsVisitConstraints; ///< Output FITS table containing mission planning visit constraints (input)

	vector<HKData> m_rawHKData; ///< Vector of raw HK data recorded at 1.2s intervals, beginning at the start of the visit
	vector<HKData> m_averagedHKData; ///< Vector of HK data averaged over the previous 16 raw intervals (19.2s), at 20s intervals, beginning at the start of the visit

	SciRawSubarray * m_stackedImageCube;  ///< FITS Image cube containing stacked images
	SciRawImagemetadata * m_stackedMetaData;  ///< FITS table containing meta data for stacked images
	SciRawUnstackedimagemetadata * m_stackedMetaData_unstacked;  ///< FITS table containing unstacked meta data for stacked images
	SimRawUnstackedsubarray * m_unstackedImageCube;  ///< FITS Image cube containing unstacked images
	SciRawImagemetadata * m_unstackedMetaData;  ///< FITS table containing meta data for unstacked images
	SciRawImagette * m_imagetteCube;  ///< FITS Image cube containing imagettes
	SciRawImagettemetadata * m_imagetteMetaData;  ///< FITS table containing meta data for imagettes
	SimRawDoubleprecisionsubarray * m_stackedImageCube_DP;  ///< FITS Image cube containing double precision stacked images
	SimRawDoubleprecisionsubarray * m_stackedImageCube_smear;  ///< FITS Image cube containing stacked images with smear trails only
	SimRawDoubleprecisionsubarray * m_stackedImageCube_cosmics;  ///< FITS Image cube containing stacked images with cosmic rays only

	SciRawBlankleft * m_stackedBlankLeftMarginCube;  ///< FITS image cube containing left blank reference columns
	SciRawBlankright * m_stackedBlankRightMarginCube;  ///< FITS image cube containing right blank reference columns
	SciRawDarkleft * m_stackedDarkLeftMarginCube;  ///< FITS image cube containing left dark reference columns
	SciRawDarkright * m_stackedDarkRightMarginCube;  ///< FITS image cube containing right dark reference columns
	SciRawDarktop * m_stackedDarkTopMarginCube;  ///< FITS image cube containing top dark reference rows
	SciRawOverscanleft * m_stackedOverscanLeftMarginCube;  ///< FITS image cube containing left overscan reference columns
	SciRawOverscanright * m_stackedOverscanRightMarginCube;  ///< FITS image cube containing right overscan reference columns
	SciRawOverscantop * m_stackedOverscanTopMarginCube;  ///< FITS image cube containing top overscan reference rows

	SimRawUnstackedblankleftimage * m_unstackedBlankLeftImageCube;  ///< FITS image cube containing left blank reference columns for unstacked images
	SimRawUnstackedblankrightimage * m_unstackedBlankRightImageCube;  ///< FITS image cube containing right blank reference columns for unstacked images
	SimRawUnstackeddarkleftimage * m_unstackedDarkLeftImageCube;  ///< FITS image cube containing left dark reference columns for unstacked images
	SimRawUnstackeddarkrightimage * m_unstackedDarkRightImageCube;  ///< FITS image cube containing left dark reference columns for unstacked images
	SimRawUnstackeddarktopimage * m_unstackedDarkTopImageCube;  ///< FITS image cube containing right dark reference rows for unstacked images
	SimRawUnstackedoverscanleftimage * m_unstackedOverscanLeftImageCube;  ///< FITS image cube containing left overscan reference columns for unstacked images
	SimRawUnstackedoverscanrightimage * m_unstackedOverscanRightImageCube;  ///< FITS image cube containing right overscan reference columns for unstacked images
	SimRawUnstackedoverscantopimage * m_unstackedOverscanTopImageCube;  ///< FITS image cube containing top overscan reference rows for unstacked images

	SciRawCentroid * m_centroid;  ///< FITS table containing centroid positions
	SimTruSubarray * m_stackedTruthMetaData;  ///< FITS table containing stacked image truth data
	SimTruSubarray * m_unstackedTruthMetaData;  ///< FITS table containing unstacked image truth data
	SimTruSubarray * m_imagetteTruthMetaData;  ///< FITS table containing imagette truth data
	RefAppBadpixelmap * m_badPixels;  ///< FITS image containing set of bad pixels
	RefAppBadpixelmapleft * m_badPixelsLeft;  ///< FITS image containing set of bad pixels (left margin)
	RefAppBadpixelmapright * m_badPixelsRight;  ///< FITS image containing set of bad pixels (right margin)
	RefAppBadpixelmaptop * m_badPixelsTop;  ///< FITS image containing set of bad pixels (top margin)
	RefAppPhotpixelmap * m_photPixels;  ///< FITS image containing set of bad pixels for photometry
	RefAppPhotpixelmapleft * m_photPixelsLeft;  ///< FITS image containing set of bad pixels for photometry (left margin)
	RefAppPhotpixelmapright * m_photPixelsRight;  ///< FITS image containing set of bad pixels for photometry (right margin)
	RefAppPhotpixelmaptop * m_photPixelsTop;  ///< FITS image containing set of bad pixels for photometry (top margin)

	vector <Image *> m_imagettes; ///< Vector of unstacked imagettes
	vector <Image *> m_stackedImagettes; ///< Vector of stacked imagettes
	vector <Image *> m_imagesToSmearUp;  ///< Set of images corresponding to the last second of each exposure
	vector <Image *> m_imagesToSmearDown;  ///< Set of images corresponding to the first second of each exposure
	Image * m_cosmicRayImage;  ///< Image containing cosmic rays only, corresponding the current stacked image
	Image * m_flatField;  ///< Flat field image, combined over wavelength
	Image * m_flatFieldToSubtract;  ///< Flat field image to be used for flat field correction
};

#endif /* _DATA_HXX_ */
