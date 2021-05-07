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

"""Gen3 Butler Formatters for LSST raw data.
"""

__all__ = (
    "LsstCamRawFormatter",
    "LatissRawFormatter",
    "LsstCamImSimRawFormatter",
    "LsstCamPhoSimRawFormatter",
    "LsstTS8RawFormatter",
    "LsstTS3RawFormatter",
    "LsstComCamRawFormatter",
    "LsstUCDCamRawFormatter",
)

from astro_metadata_translator import fix_header, merge_headers

import lsst.afw.fits
from lsst.obs.base import FitsRawFormatterBase

from ._instrument import LsstCam, Latiss, \
    LsstCamImSim, LsstCamPhoSim, LsstTS8, \
    LsstTS3, LsstUCDCam, LsstComCam
from .translators import LatissTranslator, LsstCamTranslator, \
    LsstUCDCamTranslator, LsstTS3Translator, LsstComCamTranslator, \
    LsstCamPhoSimTranslator, LsstTS8Translator, LsstCamImSimTranslator
from .assembly import fixAmpsAndAssemble, readRawAmps


class LsstCamRawFormatter(FitsRawFormatterBase):
    translatorClass = LsstCamTranslator
    filterDefinitions = LsstCam.filterDefinitions
    _instrument = LsstCam

    def readMetadata(self):
        """Read all header metadata directly into a PropertyList.

        Specialist version since some of our data does not
        set INHERIT=T so we have to merge the headers manually.

        Returns
        -------
        metadata : `~lsst.daf.base.PropertyList`
            Header metadata.
        """
        file = self.fileDescriptor.location.path
        phdu = lsst.afw.fits.readMetadata(file, 0)
        if "INHERIT" in phdu:
            # Trust the inheritance flag
            return super().readMetadata()

        # Merge ourselves
        md = merge_headers([phdu, lsst.afw.fits.readMetadata(file)],
                           mode="overwrite")
        fix_header(md)
        return md

    def getDetector(self, id):
        # Docstring inherited.
        return self._instrument.getCamera()[id]

    def readImage(self):
        # Docstring inherited.
        return self.readFull().getImage()

    def readFull(self):
        # Docstring inherited.
        rawFile = self.fileDescriptor.location.path
        ccd = self.getDetector(self.observationInfo.detector_num)
        ampExps = readRawAmps(rawFile, ccd)
        exposure = fixAmpsAndAssemble(ampExps, rawFile)
        self.attachComponentsFromMetadata(exposure)
        return exposure


class LatissRawFormatter(LsstCamRawFormatter):
    translatorClass = LatissTranslator
    _instrument = Latiss
    filterDefinitions = Latiss.filterDefinitions
    wcsFlipX = True


class LsstCamImSimRawFormatter(LsstCamRawFormatter):
    translatorClass = LsstCamImSimTranslator
    _instrument = LsstCamImSim
    filterDefinitions = LsstCamImSim.filterDefinitions


class LsstCamPhoSimRawFormatter(LsstCamRawFormatter):
    translatorClass = LsstCamPhoSimTranslator
    _instrument = LsstCamPhoSim
    filterDefinitions = LsstCamPhoSim.filterDefinitions


class LsstTS8RawFormatter(LsstCamRawFormatter):
    translatorClass = LsstTS8Translator
    _instrument = LsstTS8
    filterDefinitions = LsstTS8.filterDefinitions


class LsstTS3RawFormatter(LsstCamRawFormatter):
    translatorClass = LsstTS3Translator
    _instrument = LsstTS3
    filterDefinitions = LsstTS3.filterDefinitions


class LsstComCamRawFormatter(LsstCamRawFormatter):
    translatorClass = LsstComCamTranslator
    _instrument = LsstComCam
    filterDefinitions = LsstComCam.filterDefinitions


class LsstUCDCamRawFormatter(LsstCamRawFormatter):
    translatorClass = LsstUCDCamTranslator
    _instrument = LsstUCDCam
    filterDefinitions = LsstUCDCam.filterDefinitions
