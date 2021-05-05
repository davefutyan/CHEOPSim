#!/usr/bin/python

##imports:

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import nnls
from astropy.io import fits

import matplotlib.pyplot as plt
from skimage import exposure

## Global Variables

# Input files for code
FILE_QE =         "CH_TU2020-01-29T00-00-00_REF_APP_QE_V0101.fits"
FILE_TP =         "CH_TU2018-01-01T00-00-00_REF_APP_Throughput-BeginOfLife_V0100.fits"
FILE_SEDS =       "CH_TU2015-01-01T00-00-00_REF_APP_SEDTeff_V0101.fits"
FILE_FF_FILTER  = "CH_TU2018-01-01T00-00-00_REF_APP_FlatFieldFilter-PointSource_V0100.fits"
FILE_SED_FILTER = "CH_TU2018-01-01T00-00-00_REF_APP_SEDFilter_V0100.fits"
FILE_FF_TEFF_IN = "CH_TU2018-01-01T00-00-00_REF_APP_FlatFieldTeff-PointSource_V0100_empty.fits"
FILE_FF_TEFF_OUT= "CH_TU2020-01-29T00-00-00_REF_APP_FlatFieldTeff-PointSource_V0102.fits"
CUT_MONO_FF_NM  = 740

## My functions:


def get_global_throuput(wl):
    """
    Returns global throughput (optical x QE) sampled for wl (e.g. broadband SED)
    """ 
    #Not considering dependence on temperature for QE
    qe = fits.getdata(FILE_QE)
    tp = fits.getdata(FILE_TP)
    qe_wl = interp1d(qe["WAVELENGTH"],qe["QE"],kind='linear')(wl)
    tp_wl = interp1d(tp["WAVELENGTH"],tp["THROUGHPUT"],kind='linear')(wl)
    gt = qe_wl * tp_wl
    return gt


def get_sed_blackbody(teff):
    """
    Returns the flux of Planck black body with the specified
    effective temperature, at the specified wavelength.
    """
    h = 6.6260693e-27   #; ///< Planck constant, for energy in ergs
    k = 1.380658e-16    #; ///< Boltzmann constant, for energy in ergs
    c = 299792458.0     #; ///< Speed of light

    wavelength = np.arange(200,1150,1)
    planck = (2.*h*c*c/(wavelength*1.e-9)**5) / (np.exp(h*c/(wavelength*1.e-9*k*teff))-1.)
    planck_photons = planck/(h*c/(wavelength*1.e-9))
    return {"wavelength":wavelength, "flux_energy": planck, "flux_photons": planck_photons}

def get_sed(teff_star, sedProduct):
    """
    Get the sed for a stellar temperature
    """
    h = 6.6260693e-27   #; ///< Planck constant, for energy in ergs
    k = 1.380658e-16    #; ///< Boltzmann constant, for energy in ergs
    c = 299792458.0     #; ///< Speed of light

    teff_sed = []
    sed_data = sedProduct[1].data
    i=0 
    while (sed_data.shape[0]>i):
        teff_sed.append(sed_data['TEMPERATUR'][i])
        i+=1
    teff_sed = np.array(teff_sed)
    index = np.argmin(abs(teff_sed-teff_star))
    
    wavelength = np.array(sed_data['WAVELENGTH'][index])
    flux = np.array(sed_data['FLUX'][index])
    flux_photons = flux/(h*c/(wavelength*1.e-9))
    return {"wavelength":wavelength, "flux_energy": flux, "flux_photons": flux_photons}


def get_sed_ff(teff):
    """
    Interface to get the correct SED for a given temperature
    """
    sedProduct = fits.open(FILE_SEDS)
    if teff < sedProduct[1].data['TEMPERATUR'][0] or teff > sedProduct[1].data['TEMPERATUR'][-1]:
        sed_out = get_sed_blackbody(teff)
    else:
        sed_out = get_sed(teff, sedProduct)
    sedProduct.close()
    ### NEED TO CONFIRM IF PHOTONS OR ENERGY (For now is passing energy CHECK)
    sed_tuple_arr = (sed_out["wavelength"], sed_out["flux_energy"])
    return sed_tuple_arr



def read_FFs_data():
    """
    Read all the Flat Field data from the input data products
    """

    ff_filter_product = fits.open(FILE_FF_FILTER)

    ff_images = ff_filter_product[1].data
    ff_metadata = ff_filter_product[2].data

    mono_wl    = []
    mono_fwhm  = []
    mono_ff    = []
    mono_fferr = []

    for i in range(len(ff_metadata)):

        ## Get the monochromatic FFs only bellow 740 nm CUT

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i].replace(".","").isdigit() and float(ff_metadata["FILTER"][i]) < CUT_MONO_FF_NM:
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            mono_wl.append(float(ff_metadata["FILTER"][i]))
            mono_fwhm.append(ff_metadata["BANDWIDTH"][i])
            mono_ff.append(ff_images[i,:,:])
            #assuming that there is half of FF and half FFERROR as it should:
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            mono_fferr.append(ff_images[i_err,:,:])


        ## Get the broadband FF

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i]=="U":
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            u_ff = ff_images[i,:,:]
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            u_fferr = ff_images[i_err,:,:]          

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i]=="B":
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            b_ff = ff_images[i,:,:]
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            b_fferr = ff_images[i_err,:,:]          

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i]=="V":
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            v_ff = ff_images[i,:,:]
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            v_fferr = ff_images[i_err,:,:]          

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i]=="R":
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            r_ff = ff_images[i,:,:]
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            r_fferr = ff_images[i_err,:,:]          

        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD" and ff_metadata["FILTER"][i]=="I":
            print(ff_metadata["DATA_TYPE"][i],ff_metadata["FILTER"][i] )
            i_ff = ff_images[i,:,:]
            i_err = int(i+len(ff_metadata)/2)
            print(ff_metadata["DATA_TYPE"][i_err],ff_metadata["FILTER"][i_err] )
            i_fferr = ff_images[i_err,:,:]          
    ff_filter_product.close()


    #Get broadband SEDs:

    bb_sed_filter_product = fits.open(FILE_SED_FILTER)
    sed_filter_data = bb_sed_filter_product[1].data
    # assuming that the wavelenght array is the same for all bb seds. Taking furs from U
    wl = sed_filter_data["WAVELENGTH"][0]

    for i in range(len(sed_filter_data)):
        if sed_filter_data["FILTER"][i] == "U":
            u_sed = sed_filter_data["FLUX"][i]
        if sed_filter_data["FILTER"][i] == "B":
            b_sed = sed_filter_data["FLUX"][i]
        if sed_filter_data["FILTER"][i] == "V":
            v_sed = sed_filter_data["FLUX"][i]
        if sed_filter_data["FILTER"][i] == "R":
            r_sed = sed_filter_data["FLUX"][i]
        if sed_filter_data["FILTER"][i] == "I":
            i_sed = sed_filter_data["FLUX"][i]
    bb_sed_filter_product.close()

    broad_seds = u_sed,b_sed,v_sed,r_sed,i_sed

    mono_wl    = np.array(mono_wl)
    mono_fwhm  = np.array(mono_fwhm)
    mono_ff    = np.array(mono_ff)
    mono_fferr = np.array(mono_fferr)
    mono_flats = mono_wl, mono_fwhm, mono_ff, mono_fferr

    broad_ff = u_ff,b_ff,v_ff,r_ff,i_ff 
    broad_ffer = u_fferr,b_fferr,v_fferr,r_fferr,i_fferr
    broad_flats = wl, broad_seds, broad_ff, broad_ffer

    return mono_flats, broad_flats


# Normalised cumulative SED for monochromatic flat fields
def mono_cumspd(wl0, width, wl, gt):
    """
    Creating gaussian sed for monochromatic flats
    """
    sigma = width/(2*np.sqrt(2*np.log(2)))  # width is interpreted as the FWHM
    gauss_sed = np.exp(-(wl-wl0)**2/(2.*sigma**2))
    gauss_spd = gauss_sed*gt*wl
    cumspd = np.cumsum(gauss_spd)
    cumspd /= cumspd.max()
    return cumspd



def ff_synthesis(sed_in, mono_flats, broad_flats, gt):
    """
    The code to generate a FF for a given SED shared by Adrien passed into this function:
    The code converts all Spectral Energy Distributions (SED) into Spectral Photon
    Distributions (SPD) because the intrinsic flat field units are photons.
    Then, the fit is done with the cumulative SPD (may also work with basic SPD).
    """

# Input parameters
# Parameters extracted from reference FITS files
# star_wl = "wavelengths of the stellar SED points"
# star_sed = "stellar SED (units must be energy (e.g. erg), not photons!)"
    star_wl, star_sed = sed_in

#mono_wl = "wavelenghts of the monochromatic FF"
#mono_fwhm = "widths of the monochromatic FF"    # nearly 30nm for all
#mono_ff = "monochromatic FF data cube"
#mono_fferr = "monochromatic FF error data cube"
    mono_wl, mono_fwhm, mono_ff, mono_fferr = mono_flats

#    gt = "global throughput (optical x QE) sampled as broadband SED" (Comes as input)

#wl = "wavelengths of the broadband SED points"
#u_sed,b_sed,v_sed,r_sed,i_sed = "(U,B,V,R,I) broadband SED"
#u_ff,b_ff,v_ff,r_ff,i_ff = "(U,B,V,R,I) broadband FF"
#u_fferr,b_fferr,v_fferr,r_fferr,i_fferr = "(U,B,V,R,I) broadband FF errors"
    wl, broad_seds, broad_ff, broad_ffer = broad_flats
    u_sed,b_sed,v_sed,r_sed,i_sed = broad_seds
    u_ff,b_ff,v_ff,r_ff,i_ff = broad_ff
    u_fferr,b_fferr,v_fferr,r_fferr,i_fferr = broad_ffer

    # Normalised cumulative SPD for broadband flat fields
    u_spd = u_sed*gt*wl
    u_cumspd = np.cumsum(u_spd)
    u_cumspd /= u_cumspd.max()
    b_spd = b_sed*gt*wl
    b_cumspd = np.cumsum(b_spd)
    b_cumspd /= b_cumspd.max()
    v_spd = v_sed*gt*wl
    v_cumspd = np.cumsum(v_spd)
    v_cumspd /= v_cumspd.max()
    r_spd = r_sed*gt*wl
    r_cumspd = np.cumsum(r_spd)
    r_cumspd /= r_cumspd.max()
    i_spd = i_sed*gt*wl
    i_cumspd = np.cumsum(i_spd)
    i_cumspd /= i_cumspd.max()


    # Normalised cumulative SED for the stellar spectra
    star_sed = interp1d(star_wl,star_sed,kind='linear')(wl) # interpolating stellar SED to match broadband SED sampling
    star_spd = star_sed*gt*wl
    star_cumspd = np.cumsum(star_spd)
    star_cumspd /= star_cumspd.max()


    # Extracting the coefficients (non-negative least-square solution)
    mono_matrix = np.array([mono_cumspd(mono_wl[i],mono_fwhm[i], wl, gt) for i in range(len(mono_wl))])
    broad_matrix = np.array([u_cumspd,b_cumspd,v_cumspd,r_cumspd,i_cumspd])
    matrix = np.concatenate((mono_matrix,broad_matrix))

    coeff,res = nnls(matrix.T,star_cumspd)
    mono_coeff = coeff[:len(mono_wl)]
    broad_coeff = coeff[len(mono_wl):]


    # Generating the synthetic flat field
    ff_syn = ((mono_ff.T*mono_coeff).T.sum(axis=0)+
              (np.array([u_ff,b_ff,v_ff,r_ff,i_ff]).T*broad_coeff).T.sum(axis=0))
    ff_syn /= ff_syn.mean()


    # Estimating error on the synthetic flat field
    mono_fferr *= mono_ff
    u_fferr *= u_ff
    b_fferr *= b_ff
    v_fferr *= v_ff
    r_fferr *= r_ff
    i_fferr *= i_ff
    ff_syn_err = np.sqrt(((mono_fferr**2).T*mono_coeff).T.sum(axis=0)+
                         ((np.array([u_fferr,b_fferr,v_fferr,r_fferr,i_fferr])**2).T*broad_coeff).T.sum(axis=0))

    return ff_syn, ff_syn_err




def get_teff_array():
    """
    Create the temperature sampling for the Flat Field library
    Based on the Teff from the SED files + extra temperatures to cover stars listed for CHEOPSim
    """

## Aditional table from table 7 and 8 Cheopsim user manual
    teff_lowm = np.array([2450,2500, 2650, 2800])
    teff_hots = np.array([7440, 7500,  7800,  8000,  8080,  8270,  8550,  8840,  9200,  9700, 10700, 
                          12500, 14000, 14500,15700, 16700, 17000, 20600, 26000, 31500, 32500, 34500, 36500])

    sedProduct = fits.open(FILE_SEDS)
    teff_sed = sedProduct[1].data['TEMPERATUR']
    sedProduct.close()
    #teff_vec = np.concatenate((teff_lowm, teff_sed, teff_hots)) #No need for teff_lowm now that seds go down to 2300K
    teff_vec = np.concatenate((teff_sed, teff_hots))
    return teff_vec


def create_all_flats():
    """
    It creates the data cube for the flat images and flat_error images, returns also the temperature array
    """

    teffarr = get_teff_array()
    mono_flats, broad_flats = read_FFs_data()
    wl, broad_seds, broad_ff, broad_ffer = broad_flats
    gt = get_global_throuput(wl)

    flat_list = []
    flat_list_er = []
    for teff in teffarr:
        print("Calculating flat for TEFF: ", teff)
        sed_in = get_sed_ff(teff)
        flat_test, flat_test_er = ff_synthesis(sed_in, mono_flats, broad_flats, gt)
        flat_list.append(flat_test)
        flat_list_er.append(flat_test_er)
    return teffarr, np.array(flat_list), np.array(flat_list_er)


def fill_ref_ff(filename_out=FILE_FF_TEFF_OUT, template=FILE_FF_TEFF_IN):
    """
    Main function that puts everything together and writes the file following the data structure scheme
    """

    ffteffprod = fits.open(template)
    teffarr, ffs,ffs_err=create_all_flats()

    number_flats = ffs.shape[0]
    ff_image_cube = np.concatenate((ffs,ffs_err))
    ffteffprod[1].data = ff_image_cube

    cols = ffteffprod[2].data.columns
    header = ffteffprod[2].header
    ff_table = fits.BinTableHDU.from_columns(cols,nrows=number_flats*2)

    for iline in range(ff_image_cube.shape[0]):
        if iline < number_flats:
            data_type = "FLAT FIELD"
            teffin = teffarr[iline]
        else:
            data_type = "FLAT FIELD ERROR"
            teffin = teffarr[iline-number_flats]
        ff_table.data['DATA_TYPE'][iline] = data_type
        ff_table.data['T_EFF'][iline] = teffin
        ff_table.data['STATUS'][iline] = 0
    ff_table.header = header
    ffteffprod[2] = ff_table

    ffteffprod[1].header["FF_RF"]    = FILE_FF_FILTER
    ffteffprod[1].header["SED_T_RF"] = FILE_SEDS
    ffteffprod[1].header["SED_F_RF"] = FILE_SED_FILTER
    ffteffprod[1].header["THRGH_RF"] = FILE_TP
    ffteffprod[1].header["QE_RF"]    = FILE_QE

    ffteffprod.writeto(filename_out, overwrite=True)



def test_sed_reading():
    """
    test function
    """
    teff=6000
    out = get_sed_ff(teff)
    sed1_wl, sed1_fl = out
    sedbb = get_sed_blackbody(teff)
    plt.plot(sed1_wl, sed1_fl/np.max(sed1_fl))
    plt.plot(sedbb['wavelength'], sedbb['flux_energy']/np.max(sedbb['flux_energy']))
    plt.show()


def plot_image(image, title):
    """
    auxiliary function to plot ff image
    """
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(111)
    test = exposure.equalize_hist(image)
    cax = ax1.imshow(test, cmap=plt.cm.gray, origin='lower')
    ax1.set_title(title)
    cbar = fig.colorbar(cax)

def plot_all_ff():
    """
    auxiliary function to plot ref FF images
    """
    ff_filter_product = fits.open(FILE_FF_FILTER)

    ff_images = ff_filter_product[1].data
    ff_metadata = ff_filter_product[2].data

    for i in range(len(ff_metadata)):
        if ff_metadata["DATA_TYPE"][i] == "FLAT FIELD":
            plot_image(ff_images[i,:,:],ff_metadata["FILTER"][i])
            plt.show()


### Main program:
def main():
    fill_ref_ff()


if __name__ == "__main__":
    main()
