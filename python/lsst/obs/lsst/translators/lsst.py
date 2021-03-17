# This file is currently part of obs_lsst but is written to allow it
# to be migrated to the astro_metadata_translator package at a later date.
#
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file in this directory for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

"""Metadata translation support code for LSST headers"""

__all__ = ("TZERO", "SIMONYI_LOCATION", "read_detector_ids",
           "compute_detector_exposure_id_generic", "LsstBaseTranslator",
           "SIMONYI_TELESCOPE")

import os.path
import yaml
import logging
import re
import datetime
import hashlib

import astropy.coordinates
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation

from lsst.utils import getPackageDir

from astro_metadata_translator import cache_translation, FitsTranslator
from astro_metadata_translator.translators.helpers import tracking_from_degree_headers, \
    altaz_from_degree_headers


TZERO = Time("2015-01-01T00:00", format="isot", scale="utc")
TZERO_DATETIME = TZERO.to_datetime()

# Delimiter to use for multiple filters/gratings
FILTER_DELIMITER = "~"

# Regex to use for parsing a GROUPID string
GROUP_RE = re.compile(r"^(\d\d\d\d\-\d\d\-\d\dT\d\d:\d\d:\d\d)\.(\d\d\d)(?:[\+#](\d+))?$")

# LSST Default location in the absence of headers
SIMONYI_LOCATION = EarthLocation.from_geodetic(-70.749417, -30.244639, 2663.0)

# Name of the main survey telescope
SIMONYI_TELESCOPE = "Simonyi Survey Telescope"

obs_lsst_packageDir = getPackageDir("obs_lsst")

log = logging.getLogger(__name__)


def read_detector_ids(policyFile):
    """Read a camera policy file and retrieve the mapping from CCD name
    to ID.

    Parameters
    ----------
    policyFile : `str`
        Name of YAML policy file to read, relative to the obs_lsst
        package.

    Returns
    -------
    mapping : `dict` of `str` to (`int`, `str`)
        A `dict` with keys being the full names of the detectors, and the
        value is a `tuple` containing the integer detector number and the
        detector serial number.

    Notes
    -----
    Reads the camera YAML definition file directly and extracts just the
    IDs and serials.  This routine does not use the standard
    `~lsst.obs.base.yamlCamera.YAMLCamera` infrastructure or
    `lsst.afw.cameraGeom`.  This is because the translators are intended to
    have minimal dependencies on LSST infrastructure.
    """

    file = os.path.join(obs_lsst_packageDir, policyFile)
    try:
        with open(file) as fh:
            # Use the fast parser since these files are large
            camera = yaml.load(fh, Loader=yaml.CSafeLoader)
    except OSError as e:
        raise ValueError(f"Could not load camera policy file {file}") from e

    mapping = {}
    for ccd, value in camera["CCDs"].items():
        mapping[ccd] = (int(value["id"]), value["serial"])

    return mapping


def compute_detector_exposure_id_generic(exposure_id, detector_num, max_num):
    """Compute the detector_exposure_id from the exposure id and the
    detector number.

    Parameters
    ----------
    exposure_id : `int`
        The exposure ID.
    detector_num : `int`
        The detector number.
    max_num : `int`
        Maximum number of detectors to make space for.

    Returns
    -------
    detector_exposure_id : `int`
        Computed ID.

    Raises
    ------
    ValueError
        The detector number is out of range.
    """

    if detector_num is None:
        raise ValueError("Detector number must be defined.")
    if detector_num >= max_num or detector_num < 0:
        raise ValueError(f"Detector number out of range 0 <= {detector_num} < {max_num}")

    return max_num*exposure_id + detector_num


class LsstBaseTranslator(FitsTranslator):
    """Translation methods useful for all LSST-style headers."""

    _const_map = {}
    _trivial_map = {}

    # Do not specify a name for this translator
    cameraPolicyFile = None
    """Path to policy file relative to obs_lsst root."""

    detectorMapping = None
    """Mapping of detector name to detector number and serial."""

    detectorSerials = None
    """Mapping of detector serial number to raft, number, and name."""

    DETECTOR_MAX = 1000
    """Maximum number of detectors to use when calculating the
    detector_exposure_id.

    Note that because this is the maximum number *of* detectors, for
    zero-based ``detector_num`` values this is one greater than the maximum
    ``detector_num``.  It is also often rounded up to the nearest power of
    10 anyway, to allow ``detector_exposure_id`` values to be easily decoded by
    humans.
    """

    _DEFAULT_LOCATION = SIMONYI_LOCATION
    """Default telescope location in absence of relevant FITS headers."""

    _ROLLOVER_TIME = TimeDelta(12*60*60, scale="tai", format="sec")
    """Time delta for the definition of a Rubin Observatory start of day.
    Used when the header is missing. See LSE-400 for details."""

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """Ensure that subclasses clear their own detector mapping entries
        such that subclasses of translators that use detector mappings
        do not pick up the incorrect values from a parent."""

        cls.detectorMapping = None
        cls.detectorSerials = None

        super().__init_subclass__(**kwargs)

    def search_paths(self):
        """Search paths to use for LSST data when looking for header correction
        files.

        Returns
        -------
        path : `list`
            List with a single element containing the full path to the
            ``corrections`` directory within the ``obs_lsst`` package.
        """
        return [os.path.join(obs_lsst_packageDir, "corrections")]

    @classmethod
    def compute_detector_exposure_id(cls, exposure_id, detector_num):
        """Compute the detector exposure ID from detector number and
        exposure ID.

        This is a helper method to allow code working outside the translator
        infrastructure to use the same algorithm.

        Parameters
        ----------
        exposure_id : `int`
            Unique exposure ID.
        detector_num : `int`
            Detector number.

        Returns
        -------
        detector_exposure_id : `int`
            The calculated ID.
        """
        return compute_detector_exposure_id_generic(exposure_id, detector_num, max_num=cls.DETECTOR_MAX)

    @classmethod
    def max_detector_exposure_id(cls):
        """The maximum detector exposure ID expected to be generated by
        this instrument.

        Returns
        -------
        max_id : `int`
            The maximum value.
        """
        max_exposure_id = cls.max_exposure_id()
        # We subtract 1 from DETECTOR_MAX because LSST detector_num values are
        # zero-based, and detector_max is the maximum number *of* detectors,
        # while this returns the (inclusive) maximum ID value.
        return cls.compute_detector_exposure_id(max_exposure_id, cls.DETECTOR_MAX - 1)

    @classmethod
    def max_exposure_id(cls):
        """The maximum exposure ID expected from this instrument.

        Returns
        -------
        max_exposure_id : `int`
            The maximum value.
        """
        max_date = "2050-12-31T23:59.999"
        max_seqnum = 99_999
        max_controller = "C"  # This controller triggers the largest numbers
        return cls.compute_exposure_id(max_date, max_seqnum, max_controller)

    @classmethod
    def detector_mapping(cls):
        """Returns the mapping of full name to detector ID and serial.

        Returns
        -------
        mapping : `dict` of `str`:`tuple`
            Returns the mapping of full detector name (group+detector)
            to detector number and serial.

        Raises
        ------
        ValueError
            Raised if no camera policy file has been registered with this
            translation class.

        Notes
        -----
        Will construct the mapping if none has previously been constructed.
        """
        if cls.cameraPolicyFile is not None:
            if cls.detectorMapping is None:
                cls.detectorMapping = read_detector_ids(cls.cameraPolicyFile)
        else:
            raise ValueError(f"Translation class '{cls.__name__}' has no registered camera policy file")

        return cls.detectorMapping

    @classmethod
    def detector_serials(cls):
        """Obtain the mapping of detector serial to detector group, name,
        and number.

        Returns
        -------
        info : `dict` of `tuple` of (`str`, `str`, `int`)
            A `dict` with the serial numbers as keys and values of detector
            group, name, and number.
        """
        if cls.detectorSerials is None:
            detector_mapping = cls.detector_mapping()

            if detector_mapping is not None:
                # Form mapping to go from serial number to names/numbers
                serials = {}
                for fullname, (id, serial) in cls.detectorMapping.items():
                    raft, detector_name = fullname.split("_")
                    if serial in serials:
                        raise RuntimeError(f"Serial {serial} is defined in multiple places")
                    serials[serial] = (raft, detector_name, id)
                cls.detectorSerials = serials
            else:
                raise RuntimeError("Unable to obtain detector mapping information")

        return cls.detectorSerials

    @classmethod
    def compute_detector_num_from_name(cls, detector_group, detector_name):
        """Helper method to return the detector number from the name.

        Parameters
        ----------
        detector_group : `str`
            Name of the detector grouping.  This is generally the raft name.
        detector_name : `str`
            Detector name.

        Returns
        -------
        num : `int`
            Detector number.
        """
        fullname = f"{detector_group}_{detector_name}"

        num = None
        detector_mapping = cls.detector_mapping()
        if detector_mapping is None:
            raise RuntimeError("Unable to obtain detector mapping information")

        if fullname in detector_mapping:
            num = detector_mapping[fullname]
        else:
            log.warning(f"Unable to determine detector number from detector name {fullname}")
            return None

        return num[0]

    @classmethod
    def compute_detector_info_from_serial(cls, detector_serial):
        """Helper method to return the detector information from the serial.

        Parameters
        ----------
        detector_serial : `str`
            Detector serial ID.

        Returns
        -------
        info : `tuple` of (`str`, `str`, `int`)
            Detector group, name, and number.
        """
        serial_mapping = cls.detector_serials()
        if serial_mapping is None:
            raise RuntimeError("Unable to obtain serial mapping information")

        if detector_serial in serial_mapping:
            info = serial_mapping[detector_serial]
        else:
            raise RuntimeError("Unable to determine detector information from detector serial"
                               f" {detector_serial}")

        return info

    @staticmethod
    def compute_exposure_id(dayobs, seqnum, controller=None):
        """Helper method to calculate the exposure_id.

        Parameters
        ----------
        dayobs : `str`
            Day of observation in either YYYYMMDD or YYYY-MM-DD format.
            If the string looks like ISO format it will be truncated before the
            ``T`` before being handled.
        seqnum : `int` or `str`
            Sequence number.
        controller : `str`, optional
            Controller to use. If this is "O", no change is made to the
            exposure ID. If it is "C" a 1000 is added to the year component
            of the exposure ID. If it is "H" a 2000 is added to the year
            component.
            `None` indicates that the controller is not relevant to the
            exposure ID calculation (generally this is the case for test
            stand data).

        Returns
        -------
        exposure_id : `int`
            Exposure ID in form YYYYMMDDnnnnn form.
        """
        if "T" in dayobs:
            dayobs = dayobs[:dayobs.find("T")]

        dayobs = dayobs.replace("-", "")

        if len(dayobs) != 8:
            raise ValueError(f"Malformed dayobs: {dayobs}")

        # Expect no more than 99,999 exposures in a day
        maxdigits = 5
        if seqnum >= 10**maxdigits:
            raise ValueError(f"Sequence number ({seqnum}) exceeds limit")

        # Camera control changes the exposure ID
        if controller is not None:
            if controller == "O":
                pass
            elif controller == "C":
                # Add 1000 to the year component
                dayobs = int(dayobs)
                dayobs += 1000_00_00
            elif controller == "H":
                # Add 2000 to the year component for pHosim
                dayobs = int(dayobs)
                dayobs += 2000_00_00
            else:
                raise ValueError(f"Supplied controller, '{controller}' is neither 'O' nor 'C' nor 'H'")

        # Form the number as a string zero padding the sequence number
        idstr = f"{dayobs}{seqnum:0{maxdigits}d}"

        # Exposure ID has to be an integer
        return int(idstr)

    def _is_on_mountain(self):
        """Indicate whether these data are coming from the instrument
        installed on the mountain.

        Returns
        -------
        is : `bool`
            `True` if instrument is on the mountain.
        """
        if "TSTAND" in self._header:
            return False
        return True

    def is_on_sky(self):
        """Determine if this is an on-sky observation.

        Returns
        -------
        is_on_sky : `bool`
            Returns True if this is a observation on sky on the
            summit.
        """
        # For LSST we think on sky unless tracksys is local
        if self.is_key_ok("TRACKSYS"):
            if self._header["TRACKSYS"].lower() == "local":
                # not on sky
                return False

        # These are obviously not on sky
        if self.to_observation_type() in ("bias", "dark", "flat"):
            return False

        return self._is_on_mountain()

    @cache_translation
    def to_location(self):
        # Docstring will be inherited. Property defined in properties.py
        if not self._is_on_mountain():
            return None
        try:
            # Try standard FITS headers
            return super().to_location()
        except KeyError:
            return self._DEFAULT_LOCATION

    @cache_translation
    def to_datetime_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        # We prefer to use MJD-OBS but it seems that this can be missing
        # and so we fall back to DATE-OBS
        if self.is_key_ok("MJD-OBS"):
            self._used_these_cards("MJD-OBS")
            return Time(self._header["MJD-OBS"], scale="tai", format="mjd")
        return super().to_datetime_begin()

    @cache_translation
    def to_datetime_end(self):
        # Docstring will be inherited. Property defined in properties.py
        if self.is_key_ok("DATE-END"):
            return super().to_datetime_end()

        return self.to_datetime_begin() + self.to_exposure_time()

    @cache_translation
    def to_detector_num(self):
        # Docstring will be inherited. Property defined in properties.py
        raft = self.to_detector_group()
        detector = self.to_detector_name()
        return self.compute_detector_num_from_name(raft, detector)

    @cache_translation
    def to_detector_exposure_id(self):
        # Docstring will be inherited. Property defined in properties.py
        exposure_id = self.to_exposure_id()
        num = self.to_detector_num()
        return self.compute_detector_exposure_id(exposure_id, num)

    @cache_translation
    def to_observation_type(self):
        # Docstring will be inherited. Property defined in properties.py
        obstype = self._header["IMGTYPE"]
        self._used_these_cards("IMGTYPE")
        obstype = obstype.lower()
        if obstype in ("skyexp", "object"):
            obstype = "science"
        return obstype

    @cache_translation
    def to_observation_reason(self):
        # Docstring will be inherited. Property defined in properties.py
        if self.is_key_ok("TESTTYPE"):
            reason = self._header["TESTTYPE"]
            self._used_these_cards("TESTTYPE")
            return reason.lower()
        # no specific header present so use the default translation
        return super().to_observation_reason()

    @cache_translation
    def to_dark_time(self):
        """Calculate the dark time.

        If a DARKTIME header is not found, the value is assumed to be
        identical to the exposure time.

        Returns
        -------
        dark : `astropy.units.Quantity`
            The dark time in seconds.
        """
        if self.is_key_ok("DARKTIME"):
            darktime = self._header["DARKTIME"]*u.s
            self._used_these_cards("DARKTIME")
        else:
            log.warning("%s: Unable to determine dark time. Setting from exposure time.",
                        self._log_prefix)
            darktime = self.to_exposure_time()
        return darktime

    @cache_translation
    def to_exposure_id(self):
        """Generate a unique exposure ID number

        This is a combination of DAYOBS and SEQNUM, and optionally
        CONTRLLR.

        Returns
        -------
        exposure_id : `int`
            Unique exposure number.
        """
        if "CALIB_ID" in self._header:
            self._used_these_cards("CALIB_ID")
            return None

        dayobs = self._header["DAYOBS"]
        seqnum = self._header["SEQNUM"]
        self._used_these_cards("DAYOBS", "SEQNUM")

        if self.is_key_ok("CONTRLLR"):
            controller = self._header["CONTRLLR"]
            self._used_these_cards("CONTRLLR")
        else:
            controller = None

        return self.compute_exposure_id(dayobs, seqnum, controller=controller)

    @cache_translation
    def to_visit_id(self):
        """Calculate the visit associated with this exposure.

        Notes
        -----
        For LATISS and LSSTCam the default visit is derived from the
        exposure group.  For other instruments we return the exposure_id.
        """

        exposure_group = self.to_exposure_group()
        # If the group is an int we return it
        try:
            visit_id = int(exposure_group)
            return visit_id
        except ValueError:
            pass

        # A Group is defined as ISO date with an extension
        # The integer must be the same for a given group so we can never
        # use datetime_begin.
        # Nominally a GROUPID looks like "ISODATE+N" where the +N is
        # optional.  This can be converted to seconds since epoch with
        # an adjustment for N.
        # For early data lacking that form we hash the group and return
        # the int.
        matches_date = GROUP_RE.match(exposure_group)
        if matches_date:
            iso_str = matches_date.group(1)
            fraction = matches_date.group(2)
            n = matches_date.group(3)
            if n is not None:
                n = int(n)
            else:
                n = 0
            iso = datetime.datetime.strptime(iso_str, "%Y-%m-%dT%H:%M:%S")

            tdelta = iso - TZERO_DATETIME
            epoch = int(tdelta.total_seconds())

            # Form the integer from EPOCH + 3 DIGIT FRAC + 0-pad N
            visit_id = int(f"{epoch}{fraction}{n:04d}")
        else:
            # Non-standard string so convert to numbers
            # using a hash function. Use the first N hex digits
            group_bytes = exposure_group.encode("us-ascii")
            hasher = hashlib.blake2b(group_bytes)
            # Need to be big enough it does not possibly clash with the
            # date-based version above
            digest = hasher.hexdigest()[:14]
            visit_id = int(digest, base=16)

            # To help with hash collision, append the string length
            visit_id = int(f"{visit_id}{len(exposure_group):02d}")

        return visit_id

    @cache_translation
    def to_physical_filter(self):
        """Calculate the physical filter name.

        Returns
        -------
        filter : `str`
            Name of filter. Can be a combination of FILTER, FILTER1 and FILTER2
            headers joined by a "~". Returns "unknown" if no filter is declared
        """
        joined = self._join_keyword_values(["FILTER", "FILTER1", "FILTER2"], delim=FILTER_DELIMITER)
        if not joined:
            joined = "unknown"

        return joined

    @cache_translation
    def to_tracking_radec(self):
        if not self.is_on_sky():
            return None

        # RA/DEC are *derived* headers and for the case where the DATE-BEG
        # is 1970 they are garbage and should not be used.
        if self._header["DATE-OBS"] == self._header["DATE"]:
            # A fixed up date -- use AZEL as source of truth
            altaz = self.to_altaz_begin()
            radec = astropy.coordinates.SkyCoord(altaz.transform_to(astropy.coordinates.ICRS()),
                                                 obstime=altaz.obstime,
                                                 location=altaz.location)
        else:
            radecsys = ("RADESYS",)
            radecpairs = (("RASTART", "DECSTART"), ("RA", "DEC"))
            radec = tracking_from_degree_headers(self, radecsys, radecpairs)

        return radec

    @cache_translation
    def to_altaz_begin(self):
        if not self._is_on_mountain():
            return None

        # ALTAZ always relevant unless bias or dark
        if self.to_observation_type() in ("bias", "dark"):
            return None

        return altaz_from_degree_headers(self, (("ELSTART", "AZSTART"),),
                                         self.to_datetime_begin(), is_zd=False)

    @cache_translation
    def to_exposure_group(self):
        """Calculate the exposure group string.

        For LSSTCam and LATISS this is read from the ``GROUPID`` header.
        If that header is missing the exposure_id is returned instead as
        a string.
        """
        if self.is_key_ok("GROUPID"):
            exposure_group = self._header["GROUPID"]
            self._used_these_cards("GROUPID")
            return exposure_group
        return super().to_exposure_group()

    @staticmethod
    def _is_filter_empty(filter):
        """Return true if the supplied filter indicates an empty filter slot

        Parameters
        ----------
        filter : `str`
            The filter string to check.

        Returns
        -------
        is_empty : `bool`
            `True` if the filter string looks like it is referring to an
            empty filter slot. For example this can be if the filter is
            "empty" or "empty_2".
        """
        return bool(re.match(r"empty_?\d*$", filter.lower()))

    def _determine_primary_filter(self):
        """Determine the primary filter from the ``FILTER`` header.

        Returns
        -------
        filter : `str`
            The contents of the ``FILTER`` header with some appropriate
            defaulting.
        """

        if self.is_key_ok("FILTER"):
            physical_filter = self._header["FILTER"]
            self._used_these_cards("FILTER")

            if self._is_filter_empty(physical_filter):
                physical_filter = "empty"
        else:
            # Be explicit about having no knowledge of the filter
            # by setting it to "unknown". It should always have a value.
            physical_filter = "unknown"

            # Warn if the filter being unknown is important
            obstype = self.to_observation_type()
            if obstype not in ("bias", "dark"):
                log.warning("%s: Unable to determine the filter",
                            self._log_prefix)

        return physical_filter

    @cache_translation
    def to_observing_day(self):
        """Return the day of observation as YYYYMMDD integer.

        For LSSTCam and other compliant instruments this is the value
        of the DAYOBS header.

        Returns
        -------
        obs_day : `int`
            The day of observation.
        """
        if self.is_key_ok("DAYOBS"):
            self._used_these_cards("DAYOBS")
            return int(self._header["DAYOBS"])

        # Calculate it ourselves correcting for the Rubin offset
        date = self.to_datetime_begin().tai
        date -= self._ROLLOVER_TIME
        return int(date.strftime("%Y%m%d"))

    @cache_translation
    def to_observation_counter(self):
        """Return the sequence number within the observing day.

        Returns
        -------
        counter : `int`
            The sequence number for this day.
        """
        if self.is_key_ok("SEQNUM"):
            # Some older LATISS data may not have the header
            # but this is corrected in fix_header for LATISS.
            self._used_these_cards("SEQNUM")
            return int(self._header["SEQNUM"])

        # This indicates a problem so we warn and return a 0
        log.warning("%s: Unable to determine the observation counter so returning 0",
                    self._log_prefix)
        return 0

    @cache_translation
    def to_boresight_rotation_coord(self):
        """Boresight rotation angle.

        Only relevant for science observations.
        """
        unknown = "unknown"
        if not self.is_on_sky():
            return unknown

        self._used_these_cards("ROTCOORD")
        coord = self._header.get("ROTCOORD", unknown)
        if coord is None:
            coord = unknown
        return coord

    @cache_translation
    def to_boresight_airmass(self):
        """Calculate airmass at boresight at start of observation.

        Notes
        -----
        Early data are missing AMSTART header so we fall back to calculating
        it from ELSTART.
        """
        if not self.is_on_sky():
            return None

        # This observation should have AMSTART
        amkey = "AMSTART"
        if self.is_key_ok(amkey):
            self._used_these_cards(amkey)
            return self._header[amkey]

        # Instead we need to look at azel
        altaz = self.to_altaz_begin()
        if altaz is not None:
            return altaz.secz.to_value()

        log.warning("%s: Unable to determine airmass of a science observation, returning 1.",
                    self._log_prefix)
        return 1.0
