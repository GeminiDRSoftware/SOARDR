from astrodata import (astro_data_tag, astro_data_descriptor, returns_list,
                       TagSet, Section)
from ..soar import AstroDataSOAR
from . import lookup

class AstroDataGoodman(AstroDataSOAR):

    @staticmethod
    def _matches_data(source):
        return source[0].header.get('INSTRUME', '') in ('ghts_red', 'ghts_blue', 'ghts_red_imager', 'ghts_blue_imager')

    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['GOODMAN'])
    
    @astro_data_tag
    def _tag_arc(self):
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['ARC', 'CAL'])

    @astro_data_tag
    def _tag_bias(self):
        if self.phu.get('OBSTYPE') == 'BIAS':
            return TagSet(['BIAS', 'CAL'], blocks=['IMAGE','SPECT'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('OBSTYPE') == 'LAMPFLAT':
            return TagSet(['FLAT', 'CAL'])
    
    @astro_data_tag
    def _tag_domeflat(self):
        if self.phu.get('OBJECT', '').upper() == 'DFLAT':
            return TagSet(['DOMEFLAT', 'CAL', 'FLAT'])

    @astro_data_tag
    def _tag_camera(self):
        if self.phu.get('INSTCONF') == 'Red':
            return TagSet(['RED'])
        elif self.phu.get('INSTCONF') == 'Blue':
            return TagSet(['BLUE'])
        else: #No se si esta parte es necesaria, podria usar un raise
            try:
                return TagSet(['NEITHER_BLUE_OR_RED'])
            except TypeError:  # either is None
                return None

    @astro_data_tag
    def _tag_grating(self):
        grating_map = {
            '400_SYZY': 'GRATING_400',
            '600_SYZY': 'GRATING_600',
            '930_SYZY': 'GRATING_930',
            '1200_SYZY': 'GRATING_1200',
            '2100_SYZY': 'GRATING_2100'
        }

        grating_value = self.phu.get('GRATING')
        if grating_value in grating_map:
            return TagSet([grating_map[grating_value]])
        else:
            # Si 'GRATING' no coincide con ninguno de los valores conocidos, devuelve 'GRATING_CUSTOM'
            return TagSet(['GRATING_CUSTOM'])

    # @astro_data_tag
    # def _tag_blue(self):
    #     if self.phu.get('INSTCONF') == 'Blue':
    #         return TagSet(['BLUE'])

    @astro_data_tag
    def _tag_mode(self):
        if self.phu.get('WAVMODE') == 'IMAGING':
            return TagSet(['IMAGE'])
        elif self.phu.get('WAVMODE') != 'IMAGING':
            return TagSet(['SPECT'])    
        else: #No se si esta parte es necesaria, podria usar un raise
            return TagSet(['Error_WAVMODE'])

    # @astro_data_tag
    # def _tag_spect(self):
    #     if self.phu.get('CAM_ANG') != 0 and self.phu.get('WAVMOD') != 'Imaging':
    #         return TagSet(['SPECT'])

    # @astro_data_descriptor
    # def dec(self):
    #     """
    #     Returns the Declination of the center of field in degrees.  Since a
    #     fiber is used it coincides with the position of the target. For code
    #     re-used, use target_dec() if you really want the position of the target
    #     rather than the center of the field.

    #     Returns
    #     -------
    #     float
    #         declination in degrees
    #     """
    #     return self.target_dec()

    @astro_data_descriptor
    def filter_1(self):
        """
        Returns the name of the primary filter in the wheel.  

        Returns
        -------
        str
            primary filter

        """
        return self.phu.get('FILTER')

    @astro_data_descriptor
    def filter_2(self):
        """
        Returns the name of the secondary filter in the wheel.  

        Returns
        -------
        str
            secondary filter

        """
        return self.phu.get('FILTER2')

    @astro_data_descriptor
    def grating(self):
        """
        Returns the name of the VPH grating + manufacturer name.  

        Returns
        -------
        str
            grating

        """
        return self.phu.get('GRATING')

    @astro_data_descriptor
    def instrument(self, generic=False):
        """
        Returns the name of the instrument making the observation

        Parameters
        ----------
        generic: boolean
            If set, don't specify the specific instrument if there are clones
            (e.g., return "ghst" rather than "ghst_red" or "ghst_blue")

        Returns
        -------
        str
            instrument name
        """
        return 'ghst' if generic else self.phu.get('INSTRUME')

    @astro_data_descriptor
    def wavelength_config(self):
        """
        Returns the name of the VPH grating + configuration mode.  

        Returns
        -------
        str
            grating mode

        """
        return self.phu.get('WAVMODE')

    # @astro_data_descriptor
    # def ra(self):
    #     """
    #     Returns the Right Ascension of the center of field in degrees.  Since a
    #     fiber is used it coincides with the position of the target. For code
    #     re-used, use target_ra() if you really want the position of the target
    #     rather than the center of the field.

    #     Returns
    #     -------
    #     float
    #         right ascension in degrees
    #     """
    #     return self.target_ra()