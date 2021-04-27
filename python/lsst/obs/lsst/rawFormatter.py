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
from lsst.afw.cameraGeom import makeUpdatedDetector

from ._instrument import LsstCam, Latiss, \
    LsstCamImSim, LsstCamPhoSim, LsstTS8, \
    LsstTS3, LsstUCDCam, LsstComCam
from .translators import LatissTranslator, LsstCamTranslator, \
    LsstUCDCamTranslator, LsstTS3Translator, LsstComCamTranslator, \
    LsstCamPhoSimTranslator, LsstTS8Translator, LsstCamImSimTranslator
from .assembly import fixAmpsAndAssemble, fixAmpGeometry, readRawAmps, warn_once


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
        in_detector = self._instrument.getCamera()[id]
        # The detectors attached to the Camera object represent the on-disk
        # amplifier geometry, not the assembled raw.  But Butler users
        # shouldn't know or care about what's on disk; they want the Detector
        # that's equivalent to `butler.get("raw", ...).getDetector()`, so we
        # adjust it accordingly.  This parallels the logic in
        # fixAmpsAndAssemble, but that function and the ISR AssembleCcdTask it
        # calls aren't set up to handle bare bounding boxes with no pixel.  We
        # also can't remove those without API breakage.  So this is fairly
        # duplicative, and hence fragile.
        # We start by fixing amp bounding boxes based on the size of the amp
        # images themselves, because the camera doesn't have the right overscan
        # regions for all images.
        filename = self.fileDescriptor.location.path
        temp_detector = in_detector.rebuild()
        temp_detector.clear()
        with warn_once(filename) as logCmd:
            for n, in_amp in enumerate(in_detector):
                reader = lsst.afw.image.ImageFitsReader(filename, hdu=(n + 1))
                out_amp, _ = fixAmpGeometry(in_amp,
                                            bbox=reader.readBBox(),
                                            metadata=reader.readMetadata(),
                                            logCmd=logCmd)
                temp_detector.append(out_amp)
        adjusted_detector = temp_detector.finish()
        # Now we need to apply flips and offsets to reflect assembly.  The
        # function call that does this in fixAmpsAndAssemble is down inside
        # ip.isr.AssembleCcdTask.
        return makeUpdatedDetector(adjusted_detector)

    def readImage(self):
        """Read just the image component of the Exposure.

        Returns
        -------
        image : `~lsst.afw.image.Image`
            In-memory image component.
        """
        return self.readFull().getImage()

    def readFull(self):
        """Read the complete exposure.

        This correctly fixes amplifier bounding box deviations from
        the camera definitions, and so should provide the safest
        interface to the data.

        Returns
        -------
        exposure : `~lsst.afw.image.Exposure`
            Complete in-memory exposure representation.

        """
        rawFile = self.fileDescriptor.location.path
        ccd = self._instrument.getCamera()[self.observationInfo.detector_num]
        ampExps = readRawAmps(rawFile, ccd)
        exposure = fixAmpsAndAssemble(ampExps, rawFile)

        mask = self.readMask()
        if mask is not None:
            exposure.setMask(mask)
        variance = self.readVariance()
        if variance is not None:
            exposure.setVariance(variance)

        info = exposure.getInfo()
        info.setFilterLabel(self.makeFilterLabel())
        info.setVisitInfo(self.makeVisitInfo())
        info.setWcs(self.makeWcs(info.getVisitInfo(), info.getDetector()))

        exposure.setMetadata(self.metadata)

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
