#
#                                                               SOAR Observatory
#
#                                                                        Dragons
#                                                               soar_instruments
#                                                                soar.adclass.py
# ------------------------------------------------------------------------------
import re
import math
import datetime
import dateutil.parser

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle

from astrodata import AstroData
from astrodata import astro_data_tag
from astrodata import astro_data_descriptor
from astrodata import TagSet, Section

from ..utilities import section_to_tuple

#I copy these from gemini.py

soar_keyword_names = dict( #Those with the '#' aren't implemented yet
    adc_status = 'ADCSTAT',#
    adc_position = 'ADCPOS',#
    airmass = 'AIRMASS', 
    azimuth = 'MOUNT_AZ', 
    binning = 'CCDSUM',
    camera = 'CAMERA', 
    camera_angle_in_deg = 'CAM_ANG',#
    camera_target_in_deg = 'CAM_TARG',#
    camera_temp_in_C = 'CAM_TEMP',
    camera_focus = 'CAM_FOC', 
    collimator_focus = 'COLL_FOC',#
    collimator_temp_in_C = 'COL_TEMP',#
    dec = 'DEC',
    detector_size = 'DETSIZE', #
    dispersion_axis = 'DISPAXIS', #
    dome_shutter_azimuth = 'DOME_AZ',#
    elevation = 'MOUNT_EL', 
    epoch = 'EQUINOX', 
    exposure_time = 'EXPTIME', 
    filter_1 = 'FILTER1', ##
    filter_2 = 'FILTER2', ##
    focus = 'FOCUS',
    gain = 'GAIN', 
    grating = 'GRATING', #
    grating_angle = 'GRT_ANG', #
    grating_target = 'GRT_TARG', #
    hour_angle = 'HA', 
    instrument = 'INSTRUME', 
    local_sidereal_time = 'LST', 
    environment_wind_at_start = 'ENVWIN',
    environment_atm_pressure_at_start = 'ENVPRE', 
    environment_wind_direction_at_start = 'ENVDIR', 
    environment_outside_temperature_in_C_at_start = 'ENVTEM', 
    environment_relative_humidity_at_start = 'ENVHUM', 
    obs_ra = 'OBSRA',
    obs_dec = 'OBSDEC',
    observation_id = 'OBSID',
    observation_type = 'OBSTYPE',
    pos_angle = 'POSANGLE',
    ra = 'RA', 
    read_noise = 'RDNOISE', 
    rotator = 'ROTATOR',
    seeing = 'SEEING',
    slit = 'SLIT',
    trim_section = 'TRIMSEC',
    ut_time = 'UT', 
    wavelength_config = 'WAVMODE'
)

def use_keyword_if_prepared(fn):
    """
    A decorator for descriptors. If decorated, the descriptor will bypass its
    main code on "PREPARED" data in favour of simply returning the value of
    the associated header keyword (as defined by the "_keyword_for" method)
    if this exists in all the headers (if the keyword is missing, it will
    execute the code in the descriptor method).
    """
    def gn(self):
        if "PREPARED" in self.tags:
            try:
                return self.hdr[self._keyword_for(fn.__name__)]
            except (KeyError, AttributeError):
                pass
        return fn(self)
    return gn


# ------------------------------------------------------------------------------
class AstroDataSOAR(AstroData):
    __keyword_dict = soar_keyword_names

    @staticmethod
    def _matches_data(source):
        obs = source[0].header.get('OBSERVAT', '').upper()
        tel = source[0].header.get('TELESCOP', '').upper()

        isSOAR = (obs == 'SOAR' or tel == 'SOAR 4.1M')
        return isSOAR

    @astro_data_tag
    def _type_observatory(self):
        return TagSet(['SOAR'])

    # @astro_data_tag #Pendiente, Ver como clasificar entre imagen y espectro
    # def _type_mode(self):
    #     mode = self.phu.get(self._keyword_for('observation_mode'), '').upper()

    #     if mode:
    #         tags = [mode]
    #         if mode != 'IMAGE':
    #             # assume SPECT
    #             tags.append('SPECT')
    #         return TagSet(tags)

    @astro_data_tag
    def _status_prepared(self):
        if any(('PREPAR' in kw) for kw in self.phu):
            return TagSet(['PREPARED'])
        else:
            return TagSet(['UNPREPARED'])

    @astro_data_tag #No se que es GEM-TLM, pendiente
    def _status_raw(self):
        if 'GEM-TLM' not in self.phu:
            return TagSet(['RAW'])

    @astro_data_tag #Pendiente ver como se añaden estos tags
    def _status_overscan(self):
        found = []
        for pattern, tag in (('TRIMOVER', 'OVERSCAN_TRIMMED'), ('SUBOVER', 'OVERSCAN_SUBTRACTED')):
            if any((pattern in kw) for kw in self.phu):
                found.append(tag)
        if found:
            return TagSet(found)

    @astro_data_tag #Entiedno que aqui se definene los procesos intermedios tras los recipes
    def _status_processed_cals(self):
        kwords = {'PROCARC', 'GBIAS', 'PROCBIAS', 'PROCDARK',
                      'GIFLAT', 'PROCFLAT', 'GIFRINGE', 'PROCFRNG', 'PROCSTND', 'PROCILLM'}

        if set(self.phu.keys()) & kwords:
            return TagSet(['PROCESSED'])

    @astro_data_tag #Pendiente 
    def _status_processed_science(self):
        kwords = {'GMOSAIC', 'PROCSCI'}

        if self.phu['OBSTYPE'] == 'OBJECT' and set(self.phu.keys()) & kwords:
            return TagSet(['PROCESSED_SCIENCE', 'PROCESSED'], blocks=['RAW'])

    @astro_data_tag #Pendiente
    def _type_extracted(self):
        if 'EXTRACT' in self.phu:
            return TagSet(['EXTRACTED'])

    # # GCALFLAT and the LAMPON/LAMPOFF are kept separated because the        #Entiendo que esta parte tiene que ver con las lamparas
    # # PROCESSED status will cancel the tags for lamp status, but the
    # # GCALFLAT is still needed
    # @astro_data_tag
    # def _type_gcalflat(self):
    #     gcallamp = self.phu.get('GCALLAMP')
    #     if gcallamp == 'IRhigh' or (gcallamp is not None and gcallamp.startswith('QH')):
    #         return TagSet(['GCALFLAT', 'FLAT', 'CAL'])

    # @astro_data_tag
    # def _type_gcal_lamp(self):
    #     if self.phu['INSTRUME'].startswith('GMOS'):
    #         return

    #     gcallamp = self.phu.get('GCALLAMP')
    #     if gcallamp == 'IRhigh' or (gcallamp is not None and gcallamp.startswith('QH')):
    #         shut = self.phu.get('GCALSHUT')
    #         if shut == 'OPEN':
    #             return TagSet(['GCAL_IR_ON', 'LAMPON'], blocked_by=['PROCESSED'])
    #         elif shut == 'CLOSED':
    #             return TagSet(['GCAL_IR_OFF', 'LAMPOFF'], blocked_by=['PROCESSED'])
    #     elif self.phu.get('GCALLAMP') == 'No Value' and \
    #          self.phu.get('GCALSHUT') == 'CLOSED':
    #         return TagSet(['GCAL_IR_OFF', 'LAMPOFF'], blocked_by=['PROCESSED'])

    def _ra(self):
        """
        Parse RA from header.

        Utility method to pull the right ascension from the header, parsing text if appropriate.

        Returns
        -------
        float : right ascension in degrees, or None
        """
        ra = self.phu.get(self._keyword_for('ra'), None)
        if type(ra) == str:
            # maybe it's just a float
            try:
                return float(ra)
            except:
                try:
                    if not ra.endswith('hours') and not ra.endswith('degrees'):
                        rastr = f'{ra} hours'
                    else:
                        rastr = ra
                    return Angle(rastr).degree
                except:
                    self._logger.warning(f"Unable to parse RA from {ra}")
                    return None
        return ra

    def _dec(self):
        """
        Parse DEC from header.

        Utility method to pull the declination from the header, parsing text if appropriate.

        Returns
        -------
        float : declination in degrees, or None
        """
        dec = self.phu.get(self._keyword_for('dec'), None)
        if type(dec) == str:
            # maybe it's just a float
            try:
                return float(dec)
            except:
                try:
                    if not dec.endswith('degrees'):
                        decstr = f'{dec} degrees'
                    else:
                        decstr = dec
                    return Angle(decstr).degree
                except:
                    self._logger.warning(f"Unable to parse dec from {dec}")
                    return None
        return dec

    # Utilities
    def _parse_section(self, keyword, pretty):
        try:
            value_filter = (str if pretty else section_to_tuple)
            process_fn = lambda x: (None if x is None else value_filter(x))
            # Dummy keyword FULLFRAME returns shape of full data array
            if keyword == 'FULLFRAME':
                try:
                    sections = '[1:{1},1:{0}]'.format(*self.data.shape)
                except AttributeError:
                    sections = ['[1:{1},1:{0}]'.format(*ext.shape)
                                for ext in self.data]
            else:
                sections = self.hdr.get(keyword)
            if self.is_single:
                return process_fn(sections)
            else:
                return [process_fn(raw) for raw in sections]
        except (KeyError, TypeError):
            return None

    @astro_data_descriptor
    def airmass(self):
        """
        Returns the airmass of the observation.

        Returns
        -------
        float
            Airmass value.
        """
        return self.phu.get(self._keyword_for('airmass')) #hdr devuelve una lista, [] y phu el valor?

    @astro_data_descriptor
    def azimuth(self):
        """
        Returns the azimuth of the telescope, in degrees

        Returns
        -------
        float
            azimuth

        """
        return self.phu.get(self._keyword_for('azimuth'))

    @astro_data_descriptor
    def binning(self):
        """
        Returns the binning of the observation.

        Returns
        -------
        float
            Binning value.
        """
        return self.hdr.get(self._keyword_for('binning'))

    @astro_data_descriptor
    def camera(self, stripID=False, pretty=False): #Esto esta mal
        """
        Returns the name of the camera.  The component ID can be removed
        with either 'stripID' or 'pretty' set to True.

        Parameters
        ----------
        stripID : bool
            If True, removes the component ID and returns only the name of
            the camera.
        pretty : bool
            Same as for stripID.  Pretty here does not do anything more.

        Returns
        -------
        str
            The name of the camera with or without the component ID.

        """
        return self._may_remove_component(self._keyword_for('camera'),
                                          stripID, pretty)

    @astro_data_descriptor
    def camera_focus(self, stripID=False, pretty=False):
        """
        Returns the value of the camera focus.  

        Returns
        -------
        float
            camera focus

        """
        return self.phu.get(self._keyword_for('camera_focus'))

    @astro_data_descriptor
    def dec(self): 
        """
        Returns the Declination of the center of the field, in degrees.

        Returns
        -------
        float
            declination in degrees
        """
        dec = self.wcs_dec()
        if dec is None:
            dec = self._dec()
        return dec
   
    @astro_data_descriptor
    def elevation(self):
        """
        Returns the elevation of the telescope, in degrees

        Returns
        -------
        float
            elevation
        """
        return self.phu.get(self._keyword_for('elevation'))

    @astro_data_descriptor
    def environment_wind_at_start(self):
        """
        Returns the environment wind at the start of the exposure, in km/h

        Returns
        -------
        float
            environmental wind
        """
        return self.phu.get(self._keyword_for('environment_wind_at_start'))

    @astro_data_descriptor
    def environment_wind_direction_at_start(self): 
        """
        Returns the environment wind direction at the start of the exposure, 

        Returns
        -------
        float
            environmental wind direction
        """
        return self.phu.get(self._keyword_for('environment_wind_direction_at_start'))

    @astro_data_descriptor
    def environment_atm_pressure_at_start(self):
        """
        Returns the environment atmospheric pressure at the start of the exposure, pending units

        Returns
        -------
        float
            environmental atmospheric pressure
        """
        return self.phu.get(self._keyword_for('environment_atm_pressure_at_start'))  

    @astro_data_descriptor
    def environment_outside_temperature_in_C_at_start(self):
        """
        Returns the temperature outside the observaratory at the start of the exposure, in C

        Returns
        -------
        float
            outside temperature
        """
        return self.phu.get(self._keyword_for('environment_outside_temperature_in_C_at_start'))   
     
    @astro_data_descriptor
    def environment_relative_humidity_at_start(self):
        """
        Returns the relative humidity outside the observaratory at the start of the exposure, pending units

        Returns
        -------
        float
            outside humidity
        """
        return self.phu.get(self._keyword_for('environment_relative_humidity_at_start'))   

    @astro_data_descriptor
    def epoch(self):
        """
        Returns the current epoch for celestial coordinates

        Returns
        -------
        float
            epoch
        """
        return self.phu.get(self._keyword_for('epoch'))

    @astro_data_descriptor
    def focus(self):
        """
        Returns the focus of the telescope

        Returns
        -------
        float
            focus used for the observation
        """
        return self.phu.get(self._keyword_for('focus'))

    @astro_data_descriptor
    def gain(self):
        """
        Returns the gain (electrons/ADU) for each extension

        Returns
        -------
        list of floats/float
            Gains used for the observation
        """
        return self.phu.get(self._keyword_for('gain'))    

    @astro_data_descriptor
    def hour_angle(self):
        """
        Returns the hour angle of the target

        Returns
        -------
        
            hour angle
        """
        return self.phu.get(self._keyword_for('hour_angle')) 

    @astro_data_descriptor
    def instrument(self):
        """
        Returns the current instrument used in SOAR

        Returns
        -------
        
            str
        """
        return self.phu.get(self._keyword_for('instrument')) 

    @astro_data_descriptor
    def local_sidereal_time(self): 
        """
        Returns the local time stored at the time of the observation.

        Returns
        -------
        datetime.datetime.time()
            Local time of the observation.
        """
        try:
            local_time = self.phu[self._keyword_for('local_sidereal_time')]
            return dateutil.parser.parse(local_time).time()
        except (ValueError, TypeError, KeyError):
            return None

    @astro_data_descriptor
    def observation_type(self):
        """
        Returns the type of an observation, e.g., 'OBJECT', 'FLAT', 'ARC'.

        Returns
        -------
        str
            the observation type
        """
        return self.phu.get('OBSTYPE')

    @astro_data_descriptor 
    def ra(self):
        """
        Returns the Right Ascension of the center of the field, in degrees.

        Returns
        -------
        float
            right ascension in degrees
        """
        ra = self.wcs_ra()
        if ra is None:
            ra = self._ra()
        return ra   

    @astro_data_descriptor
    def read_noise(self):
        """
        Returns the read noise in electrons for each extension. A list is
        returned unless called on a single-extension slice, when a float

        Returns
        -------
        float/list of floats
            the read noise
        """
        return self.phu.get(self._keyword_for('read_noise'))

    @astro_data_descriptor
    def rotator(self):
        """
        Returns the rotator. Pending docstring

        Returns
        -------
        float/list of floats
            the read noise
        """
        return self.phu.get(self._keyword_for('rotator'))

    @astro_data_descriptor
    def seeing(self):
        """
        Returns the value of the seeing, in arcsec

        Returns
        -------
        float
            seeing
        """
        return self.phu.get(self._keyword_for('seeing'))

    @astro_data_descriptor
    def slit(self):
        """
        Returns the name of the entrance slit used for the observation

        Returns
        -------
        str
            the slit name
        """
        return self.phu.get(self._keyword_for('slit'))

    @astro_data_descriptor
    def trim_section(self):
        """
        Returns the list of values of the image trim section

        Returns
        -------
        list, int
            trim section
        """
        return self.phu.get(self._keyword_for('trim_section'))

    @astro_data_descriptor
    def ut_time(self):
        """
        Returns the date of observation in UT

        Returns. Pending docstring
        -------
        list, float
            trim section
        """
        return self.phu.get(self._keyword_for('ut_time'))

    @astro_data_descriptor
    def wavelength_config(self):
        """
        Returns the mode of observation in SOAR. For goodman spectropraph it
         will return the current configurations modes for the grating. Typical
          modes includes 400_M1, 400_M2.

        Returns string?
        -------
        str
            wavelength mode
        """
        return self.phu.get(self._keyword_for('wavelength'))

    def _get_wcs_coords(self):
        """
        Returns the RA and dec of the location around which the celestial
        sphere is being projected (CRVALi in a FITS representation)

        Returns
        -------
        dict
            {'lon': right ascension, 'lat': declination} plus other coords
        """
        wcs = self.wcs if self.is_single else self[0].wcs
        if wcs is None:
            return None

        ra = dec = None
        coords = {name: None for name in wcs.output_frame.axes_names}
        for m in wcs.forward_transform:
            try:
                ra = m.lon.value
                dec = m.lat.value
            except AttributeError:
                pass

        if 'NON_SIDEREAL' in self.tags and ra is not None and dec is not None:
            ra, dec = gmu.toicrs('APPT', ra, dec, ut_datetime=self.ut_datetime())

        coords["lon"] = ra
        coords["lat"] = dec
        return coords

    @astro_data_descriptor
    def wcs_ra(self):
        """
        Returns the Right Ascension of the center of the field based on the
        WCS rather than the RA header keyword.

        Returns
        -------
        float
            right ascension in degrees
        """
        # Return None if the WCS isn't sky coordinates
        try:
            return self._get_wcs_coords()['lon']
        except (KeyError, TypeError):
            return None

    @astro_data_descriptor
    def wcs_dec(self):
        """
        Returns the Declination of the center of the field based on the
        WCS rather than the DEC header keyword.

        Returns
        -------
        float
            declination in degrees
        """
        # Return None if the WCS isn't sky coordinates
        try:
            return self._get_wcs_coords()['lat']
        except (KeyError, TypeError):
            return None