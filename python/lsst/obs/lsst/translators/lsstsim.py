# This file is part of obs_lsst.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Metadata translation code for LSST simulations"""

__all__ = ("LsstSimTranslator", )

import warnings
import logging

import astropy.utils.exceptions
from astropy.time import Time
from astropy.coordinates import AltAz

from astro_metadata_translator import cache_translation, StubTranslator

from .lsst import LSST_LOCATION, read_detector_ids, compute_detector_exposure_id

log = logging.getLogger(__name__)


class LsstSimTranslator(StubTranslator):
    """Shared routines for LSST Simulated Data.
    """
    cameraPolicyFile = None
    """Path to policy file relative to obs_lsst root."""

    detectorMapping = None
    """Mapping of detector name to detector number."""

    @cache_translation
    def to_telescope(self):
        # Docstring will be inherited. Property defined in properties.py
        telescope = None
        if "OUTFILE" in self._header and self._header["OUTFILE"].startswith("lsst"):
            telescope = "LSST"
            self._used_these_cards("OUTFILE")
        elif "LSST_NUM" in self._header:
            telescope = "LSST"
            self._used_these_cards("LSST_NUM")
        return telescope

    @cache_translation
    def to_location(self):
        # Docstring will be inherited. Property defined in properties.py
        location = None
        tel = self.to_telescope()
        if tel == "LSST":
            location = LSST_LOCATION
        return location

    @cache_translation
    def to_datetime_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        self._used_these_cards("MJD-OBS")
        return Time(self._header["MJD-OBS"], scale="tai", format="mjd")

    @cache_translation
    def to_datetime_end(self):
        # Docstring will be inherited. Property defined in properties.py
        return self.to_datetime_begin() + self.to_exposure_time()

    @cache_translation
    def to_detector_num(self):
        # Docstring will be inherited. Property defined in properties.py
        raft = self.to_detector_group()
        detector = self.to_detector_name()
        fullname = f"{raft}_{detector}"

        num = None
        if self.cameraPolicyFile is not None:
            if self.detectorMapping is None:
                self.__class__.detectorMapping = read_detector_ids(self.cameraPolicyFile)
            if fullname in self.detectorMapping:
                num = self.detectorMapping[fullname]
            else:
                log.warning("Unable to determine detector number from detector name {fullname}")

        return num

    @cache_translation
    def to_detector_exposure_id(self):
        exposure_id = self.to_exposure_id()
        num = self.to_detector_num()
        return compute_detector_exposure_id(exposure_id, num, max_num=1000, mode="concat")

    @cache_translation
    def to_observation_type(self):
        # Docstring will be inherited. Property defined in properties.py
        obstype = self._header["IMGTYPE"]
        self._used_these_cards("IMGTYPE")
        obstype = obstype.lower()
        if obstype == "skyexp":
            obstype = "science"
        return obstype

    @cache_translation
    def to_altaz_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        if self.to_observation_type() == "science":
            # Derive from RADec in absence of any other information
            radec = self.to_tracking_radec()
            if radec is not None:
                # This triggers warnings because of the future dates
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=astropy.utils.exceptions.AstropyWarning)
                    altaz = radec.transform_to(AltAz)
                return altaz
        return None
