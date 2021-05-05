/*
 * ImageWriter.hxx
 *
 *  Created on: 17 Nov 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This Module is used to write CHEOSPim output from the Data in
 *  	   memory to FITS output files on disk.
 *
 *  Information written out by this Module includes all types of image
 *  (unstacked, stacked, imagettes, double precision stacked images, full
 *  frame images), including CCD margin data and image metadata as FITS
 *  extensions. With the exception of the full frame
 *  image, which is written once at the start of the simulation, all image
 *  data is written out each time there is a new stacked image, rather than
 *  for every time step. The class performs initialization of all output
 *  FITS data structures, including the setting of header keywords. The truth
 *  data associated to each image is also output. Houskeeping data is also
 *  output with user defined cadence.
 */

#ifndef IMAGEWRITER_HXX_
#define IMAGEWRITER_HXX_

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/multi_array.hpp"

#include "SCI_RAW_FullArray.hxx"

#include "Photometry.hxx"
#include "Module.hxx"

class ImageWriter: public Module {
public:
	static const int32_t kNullInt = -2147483648; ///< NULL value for truth data
	static constexpr float kNullFloat = -9999.; ///< NULL value for truth data

	/// @brief Margin mode: image, reduced (3 values per row) or total collapsed (4 values in total)
	enum MARGIN_MODE{image,reduced,total_collapsed};

	/// @brief Struct to hold mean, median, standard deviation and MAD values for CCD margins
	struct Averages {
		/** *************************************************************************
		 *  @brief Constructor
		 *
		 *  @param [in] meanVal  Mean value of pixels within a dark, blank or
		 *  					 overscan region of the CCD margins
		 *	@param [in] medianVal  Median value of pixels within a dark, blank or
		 *  					   overscan region of the CCD margins
		 *	@param [in] stddevVal  Standard deviation of pixel values within a dark,
		 *						   blank or overscan region of the CCD margins
		 *  @param [in] madVal  Median absolute deviation of pixels values within a
		 *  					dark, blank or overscan region of the CCD margins
		 */
		Averages(double meanVal, double medianVal, double stddevVal, double madVal) :
			mean(meanVal),median(medianVal),stddev(stddevVal),mad(madVal) {};
		Averages() :
			mean(0),median(0),stddev(0),mad(0) {};
		double mean; ///< Mean value of pixels within a dark, blank or overscan region of the CCD margins
		double median; ///< Median value of pixels within a dark, blank or overscan region of the CCD margins
		double stddev; ///< Standard deviation of pixel values within a dark, blank or overscan region of the CCD margins
		double mad; ///< Median absolute deviation of pixels values within a dark, blank or overscan region of the CCD margins
	};

	typedef boost::multi_array<double, 2> array2D; ///< 2D array type used for CCD margin sub-images

	ImageWriter() : Module("ImageWriter",timeLoop), m_writeUnstackedImages(true),
												    m_writeStackedImages(true),
												    m_writeImagettes(true),
												    m_doublePrecisionStackedImages(false),
												    m_marginMode(image),
												    m_writeCentroid(true),
												    m_writeTruthData(true),
													m_redundantHardware(false) {};
	virtual ~ImageWriter() {delete m_photometry;}

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;
	void doEnd(Data * data) const;

private:
	static const unsigned kJitterTruthSize = 600; ///< Size of array to hold jitter truth information for sub-array images
	static const unsigned kPsfTruthSize = 500; ///< Size of array to hold PSF truth information for sub-array images
	static const unsigned kCosmicTruthSize = 300; ///< Size of array to hold cosmic ray truth information for sub-array images
	static const unsigned kHotTruthSize = 500; ///< Size of array to hold hot pixel truth information for sub-array images
	static const unsigned kDeadTruthSize = 200; ///< Size of array to hold dead pixel truth information for sub-array images
	static const unsigned kCosmicTruthSize_fullFrame = 2000; ///< Size of array to hold cosmic ray truth information for the full frame image
	static const unsigned kHotTruthSize_fullFrame = 12500; ///< Size of array to hold hot pixel truth information for the full frame image
	static const unsigned kDeadTruthSize_fullFrame = 5000; ///< Size of array to hold dead pixel truth information for the full frame image

	/** *************************************************************************
	 *  @brief Creates the output FITS file for the templated data type,
	 *  	   sets the header keywords and adds FITS extensions to hold the
	 *  	   CCD margin images. Also creates the output FITS file to contain
	 *  	   the  corresponding truth metadata and sets its header keywords.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] numberOfImages  Number of images to be stored in the current
	 *  						    image cube
	 */
	template <class T> void initializeFitsOutput(Data * data, Data::IMAGETYPE type,
			boost::posix_time::ptime startTime, unsigned numberOfImages) const;

	/** *************************************************************************
	 *  @brief Adds FITS extensions to the image cube to hold CCD margin images
	 *  	   (dark, blank and overscan), and sets the corresponding header keywords
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] xdim  Number of pixels for x dimension of the sub-array image
	 *  @param [in] ydim  Number of pixels for y dimension of the sub-array image
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] numberOfImages  Number of images to be stored in the current
	 *  						    image cube
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void initializeCcdMargins(T * imageCube, Data * data, int xdim, int ydim,
			boost::posix_time::ptime startTime, unsigned numberOfImages, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Template specialization of initializeCcdMarginImages<class T> for
	 *  T = SimRawUnstackedsubarray (required because FITS image types differ)
	 */
	void initializeCcdMargins(SimRawUnstackedsubarray * imageCube, Data * data, int xdim, int ydim,
			boost::posix_time::ptime startTime, unsigned numberOfImages, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Template specialization of initializeCcdMarginImages<class T> for
	 *  T = SciRawImagette: performs no action since Imagette cube should not
	 *  contain CCD margin extensions
	 */
	void initializeCcdMargins(SciRawImagette * imageCube, Data * data, int xdim, int ydim,
			boost::posix_time::ptime startTime, unsigned numberOfImages, Data::IMAGETYPE type) const {};

	/** *************************************************************************
	 *  @brief Creates the output FITS file for the onboard centroid and sets the
	 *  	   header keywords
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void initializeCentroid(Data *data) const;

	/** *************************************************************************
	 *  @brief Called from initializeFitsOutput. Sets header keywords for the
	 *  	   image cube
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void setImageKeywords(T * imageCube, Data * data,
			boost::posix_time::ptime startTime, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords specifying the readout setup
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template<class T> void setReadoutKeywords(T * imageCube, Data * data, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Returns the readout timing script ID
	 *
	 *  @param [in] readMode  readout mode (ultrabright, bright, faint fast or faint)
	 *  @param [in] exposureTime  Exposure duration in seconds
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	int readoutScript(string readMode, double exposureTime, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords corresponding to the end time of the image cube.
	 *  	   Called from writeImages after each write, so that the values
	 *  	   are updated until the output file is closed.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] timeStep  Current time step of the simulation
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void setEndTimeKeywords(Data * data, int timeStep, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords corresponding to the end time of the image cube
	 *  	   for the CCD margin images. Called from setEndTimeKeywords.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] endTime_utc  UTC time corresponding to the end of the last
	 *  						 image in the image cube
	 *  @param [in] lastImageMidTime_utc  UTC time corresponding to the midpoint
	 *  						 		  of the last image in the image cube
	 */
	template <class T> void setCcdMarginImageEndTime(Data *data, UTC endTime_utc, UTC lastImageMidTime_utc) const;

	/** *************************************************************************
	 *  @brief Sets the keywords for CCD margin (dark, blank and overscan) rows
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void setCcdMarginRowKeywords(T * imageCube, Data * data,
			boost::posix_time::ptime startTime, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords for CCD margin (dark, blank and overscan) columns
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void setCcdMarginColumnKeywords(T * imageCube, Data * data,
			boost::posix_time::ptime startTime, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords to describe the format of the CCD margin columns
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] fullFrame  Flag to indicate whether or not the image is full frame
	 *  @param [in] nAverages  Number of values for averages:
	 *  					   can be 4 (mean, median, stddev, mad),
	 *  					   3 (mean, median, stddev),
	 *  					   or the default value 0 (averges not calculated)
	 */
	template <class T> void setCcdMarginDescriptionKeywords(T * imageCube, Data * data, bool fullFrame, int nAverages=0) const;

	/** *************************************************************************
	 *  @brief Sets the keywords that are common to all types of iamge cube
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] startTime  Time corresponding to the start of the first image
	 *  					   in the current image cube
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void setCommonKeywords(T * imageCube, Data * data,
			boost::posix_time::ptime startTime, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Sets the keywords for the x and y position of the image origin
	 *  	   (defined as the offset of the image sub-array w.r.t. the edge of
	 *  	   the exposed part of the CCD)
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 */
	template<class T> void setImageOriginKeywords(T * imageCube, Data * data) const;

	/** *************************************************************************
	 *  @brief Template specialization of setImageOriginKeywords<class T> for
	 *  T = SciRawFullarray: performs no action since full frame image cubes
	 *  do not use image origin header keywords
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] data  Pointer to the Data
	 */
	void setImageOriginKeywords(SciRawFullarray * imageCube, Data * data) const {}

	/** *************************************************************************
	 *  @brief Sets header keywords relating to the visit
	 *
	 *  @param [in] fitsFile  Pointer to the FITS file
	 */
	template <class T> void setVisitKeywords(T * fitsFile) const;

	/** *************************************************************************
	 *  @brief Sets header keywords relating to the target star
	 *
	 *  @param [in] fitsFile  Pointer to the FITS file
	 */
	template <class T> void setTargetKeywords(T * fitsFile) const;

	/** *************************************************************************
	 *  @brief Appends the metadata as a FITS extension to the image cube,
	 *  	   sets its header keywords, and adds the image cube and metadata
	 *  	   to the Data
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] imageCube  Pointer to the image cube
	 */
	template <class T> void setFitsImageCubeWithMetaData(Data * data, Data::IMAGETYPE type, T * imageCube) const;

	/** *************************************************************************
	 *  @brief Performs ADU saturation of images by calling Image::saturate for
	 *  	   each image
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] nBits  Number of bits for the ADC conversion
	 */
	void saturateImages(Data* data, unsigned nBits) const;

	/** *************************************************************************
	 *  @brief Performs integer rounding followed by image stacking
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void stackImages(Data* data) const;

	/** *************************************************************************
	 *  @brief Extracts an imagette from each unstacked image and stores it in the Data
	 *  	   Imagette is centred on the centre of the sub-array, or on the PSF
	 *  	   barycentre, depending on the value of the dymanaicImagettes parameter
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void extractImagettes(Data* data) const;

	/** *************************************************************************
	 *  @brief Performs stacking of imagettes
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void stackImagettes(Data* data) const;

	/** *************************************************************************
	 *  @brief Performs writing of image output, called each time there a set of
	 *  	   unstacked images corresponding to one stacked image in memory.
	 *  	   Starts new image cubes on the first call or whenever an image cube
	 *  	   is full. Calls saturateImages, followed by stackImages, followed
	 *  	   by writeImageData, followed by setEndTimeKeywords.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] timeStep  Current time step of the simulation
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void writeImages(Data * data, int timeStep, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Called by writeImages. Writes image data to the image cube,
	 *  	   including CCD margins and truth metadata
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] imageCount  Number of images since the start of the simulation
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void writeImageData(Data * data, int imageCount, int cubeLayer, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Writes the ADU value of the current pixel to the FITS image cube,
	 *  	   casting to uint_16 or uint_32 according to the type of image
	 *
	 *  @param [in] imageCube  Pointer to the image cube
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 *  @param [in] i  Horizontal index of pixel (full CCD coordinates including margins)
	 *  @param [in] j  Vertical index of pixel (full CCD coordinates including margins)
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] pixelValue  ADU value for the pixel
	 *  @param [in] round  Flag to indicate whether or not to round to integer
	 *  				   (set to false for CCD margins)
	 */
	template <class T> void writePixelData(T *imageCube, int cubeLayer, int i, int j, Data::IMAGETYPE type, double pixelValue, bool round=true) const;

	/** *************************************************************************
	 *  @brief Called by writeImageData. Writes CCD margin data according to the
	 *  	   margin mode to the corresponding FITS extension of the image cube
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] image  Pointer to the CCD margin image array
	 *  @param [in] averages  Mean, median, standard deviation and MAD values for the margin
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] marginMode  Margin processing mode (image, reduced or total_collapsed)
	 */
	template <class T> void writeCcdMargin(Data * data, array2D * image, ImageWriter::Averages averages, int cubeLayer, Data::IMAGETYPE type, MARGIN_MODE marginMode) const;

	/** *************************************************************************
	 *  @brief Called by writeImageData. Writes image data for CCD margins to the
	 *		   corresponding FITS extension of the image cube
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] image  Pointer to the CCD margin image array
	 *  @param [in] image  Pointer to the CCD margin image array
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 */
	template <class T> void writeCcdMarginImage(Data * data, array2D * image, int cubeLayer, Data::IMAGETYPE type) const;

	/** *************************************************************************
	 *  @brief Called by writeFullFrameImage and writeImageData.
	 *  	   Returns the mean, median, standard deviation and MAD values for a
	 *  	   CCD margin region collapsed to a single row/column (reduced mode)
	 *  	   or to a single value (total collapsed mode).
	 *
	 *  @param [in] image  Pointer to the CCD margin image array
	 *  @param [in] rowcol_index  Reduced mode: Index of column or row for which
	 *  				  the average is to be calculated (should be called for
	 *  				  each (column/row) of the margin).
	 *  				  Total collapsed mode: For default value -1, calculates the
	 *  				  average over all rows/columns of the margin
	 *  @return mean, median, standard deviation and MAD values for the CCD margin region
	 */
	Averages ccdMarginAverages(array2D * image, int rowcol_index=-1) const;

	/** *************************************************************************
	 *  @brief Returns the PSF barycentre either from truth or from photometric
	 *  	   extraction depending on the truthBarycentre input parameter
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @param [in] subarray  Dimensions of the sub-array
	 *  @return PSF barycentre either from truth or from photometric extraction
	 */
	pair<double,double> getBarycentre(const Image * image, const Data::SubarrayDimensions & subarray) const;

	/** *************************************************************************
	 *  @brief Writes the image metadata for the specified image, including
	 *  	   CCD margins collapsed to single values
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] imageCount  Number of images since the start of the simulation
	 */
	void writeMetaData(Data * data, Data::IMAGETYPE type, int imageCount) const;

	/** *************************************************************************
	 *  @brief Writes the centroid the current image to a SCI_RAW_Centroid FITS
	 *  	   file. The centroid is calculated using the getBarycentre method.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] timeStep  Current time step of the simulation
	 */
	void writeCentroid(Data * data, int timeStep) const;

	/** *************************************************************************
	 *  @brief Writes the truth metadata for the current image to a
	 *  	   SIM_TRU_Subarray or SIM_TRU_FullArray FITS file.
	 *  	   Calls fillTruthData.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] truthData  Pointer to the TruthData object associated to the
	 *  					   current Image
	 *  @param [in] imageXOffset  Offset of the left edge of the image sub-array
	 *  						  w.r.t. the left edge of the exposed part of the CCD
	 *  @param [in] imageYOffset  Offset of the bottom edge of the image sub-array
	 *  						  w.r.t. the bottom edge of the exposed part of the CCD
	 *  @param [in] imageXSize  Number of pixels for x dimension of the sub-array image
	 *  @param [in] imageYSize  Number of pixels for y dimension of the sub-array image
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] imageCount  Number of images since the start of the simulation
	 */
	void writeTruthMetaData(Data * data, TruthData * truthData, int imageXOffset, int imageYOffset,
							int imageXSize, int imageYSize, Data::IMAGETYPE type, int imageCount) const;

	/** *************************************************************************
	 *  @brief Copies the truth information from the TruthData object associated
	 *  	   to the Image in memory to the SIM_TRU_Subarray or
	 *  	   SIM_TRU_FullArray FITS file. Called by writeTruthMetaData.
	 *
	 *  @param [in] truthMetaData  Pointer to the truth metadata FITS table
	 *  						   SIM_TRU_SubArray or SIM_TRU_FullArray
	 *  @param [in] truthData  Pointer to the TruthData object associated to the
	 *  					   current Image
	 *  @param [in] satData  Pointer to the SateliteData
	 *  @param [in] imageXOffset  Offset of the left edge of the image sub-array
	 *  						  w.r.t. the left edge of the exposed part of the CCD
	 *  @param [in] imageYOffset  Offset of the bottom edge of the image sub-array
	 *  						  w.r.t. the bottom edge of the exposed part of the CCD
	 *  @param [in] imageXSize  Number of pixels for x dimension of the sub-array image
	 *  @param [in] imageYSize  Number of pixels for y dimension of the sub-array image
	 *  @param [in] fullFrame  Flag to indicate whether or not the image is full frame
	 *  @param [in] imagette  Flag to indicate whether or not the image is an imagette
	 */
	template <class T> void fillTruthData(T * truthMetaData, TruthData * truthData, vector<SatelliteData*> satData,
										  int imageXOffset, int imageYOffset, unsigned imageXSize, unsigned imageYSize,
										  bool fullFrame, bool imagette) const;

	/** *************************************************************************
	 *  @brief Prints the truth information to the log file
	 *
	 *  @param [in] truthData  Pointer to the TruthData object associated to the
	 *  					   current Image
	 */
	void printTruthData(TruthData* truthData) const;

	/** *************************************************************************
	 *  @brief Called by writeImages. Writes images for frame transfer smear trails
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] imageCount  Number of images since the start of the simulation
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 */
	void writeTruthSmearTrails(Data * data, int imageCount, int cubeLayer) const;

	/** *************************************************************************
	 *  @brief Called by writeImages. Writes images containing cosmic rays only
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] imageCount  Number of images since the start of the simulation
	 *  @param [in] cubeLayer  Layer of the current image cube to write to
	 */
	void writeTruthCosmicImage(Data * data, int imageCount, int cubeLayer) const;

	/** *************************************************************************
	 *  @brief Creates a SCI_RAW_FullArray FITS file, defines its metadata, and
	 *  	   writes the full frame image data to it
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] timeStep  Time step of the simulation (this method is called for
	 *  					  images in the time loop if the sub-array dimension
	 *  					  corresponds to the full frame)
	 *  @param [in] initialFullFrame  Boolean to indicate whether or not the image
	 *  							  to be written is the full frame image before
	 *  							  the time loop
	 */
	void writeFullFrameImage(Data * data, int timeStep, bool initialFullFrame) const;


	/** *************************************************************************
	 *  @brief Returns the duration of an image in seconds. For stacked images,
	 *  	   the duration is from the start of the first exposure to the end of
	 *  	   the last exposure. It does not include the gap after the end of
	 *  	   the last exposure in the case where the repetition period is
	 *  	   longer than the exposure duration
	 *
	 *  @param [in] type  Type of image (stacked, unstacked, imagette or fullframe)
	 *  @param [in] timeConf  The time configuration
	 */
	double getImageDuration(TimeConfiguration timeConf, Data::IMAGETYPE type) const;

	bool m_writeUnstackedImages; ///< Flag to indicate whether or not to output unstacked images
	bool m_writeStackedImages; ///< Flag to indicate whether or not to output stacked images
	bool m_writeImagettes; ///< Flag to indicate whether or not to output imagettes
	bool m_dynamicImagettes; ///< Flag to indicate if imagette extraction should be dynamic (the centre of the imagette follows the centroid) rather than static (fixed position on CCD - required if m_stackImagettes is true)
	bool m_stackImagettes; ///< Flag to indicate whether or not to stack the imagettes (requires m_dynamicImagettes to be false)
	string m_imagetteShape; ///< Shape of the imagettes: either CIRCULAR or RECTANGULAR
	string m_subarrayShape; ///< Shape of the sub-array images: either CIRCULAR or RECTANGULAR
	bool m_doublePrecisionStackedImages; ///< Flag to indicate whether or not stacked images should be output in double precision
	MARGIN_MODE m_marginMode; ///< Margin mode for the visit (image, reduced or total collapsed)
	bool m_writeCentroid; ///< Flag to indicate whether or not to output onboard centroid data
	bool m_writeTruthData; ///< Flag to indicate whether or not to output truth data
	unsigned m_imagetteSize; ///< Size of imagette in number of pixels along each side of a square
	string m_stackingMethod; ///< Stacking method: coadd (coaddition) or mean (mean value for each pixel)
	string m_marginStackingMethod; ///< Stacking method: coadd (coaddition) or mean (mean value for each pixel) to be applied to the CCD margins
	bool m_truthBarycentre; ///< Flag to indicate whether or not to use to the truth barycentre rather than the barycentre from photometry
	bool m_targetLocationFromJitter; ///< Flag to indicate whether or not to define the target location using the jitter offsets in order to have a different target location for each image
	unsigned m_imagettesPerStackedSubarray; ///< Number of imagettes corresponding to each stacked sub-array
	bool m_redundantHardware; ///< Flag to indicate whether or not to use the redundant right amplifier rather than left amplifier for CCD readout
	unsigned m_ceCounterOffset; ///< CE counter initial offset, intended for Dark Frame M&C data, for which CE counter increments through consecutive visits

	Data::Visit m_visit; ///< Visit ID data
	string m_specType; ///< Spectral type to be used for header keywords
	double m_TEff; ///< Effective temperature to be used for header keywords
	double m_Gmag; ///< Gaia magnitude of the star to be used for header keywords
	double m_cheopsMag; ///< CHEOPS magnitude of the star to be used for header keywords
	double m_GmagErr; ///< Error on Gaia magnitude of the star to be used for header keywords
	double m_cheopsMagErr; ///< Error on CHEOPS magnitude of the star to be used for header keywords
	double m_nominalGain; ///< Nominal gain
	unsigned m_mpsVisit; ///< MPS visit ID
	string m_dirname; ///< Name of output directory

	Photometry * m_photometry; ///< Pointer to an instance of the Photometry class, used to extract the barycentre

};

#endif /* IMAGEWRITER_HXX_ */
