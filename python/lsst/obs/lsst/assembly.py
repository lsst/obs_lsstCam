# This file is part of obs_lsst
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

__all__ = ("attachRawWcsFromBoresight", "fixAmpGeometry", "assembleUntrimmedCcd",
           "fixAmpsAndAssemble", "readRawAmps")

import lsst.log
import lsst.afw.image as afwImage
from lsst.obs.base import bboxFromIraf, MakeRawVisitInfoViaObsInfo, createInitialSkyWcs
from lsst.geom import Box2I, Extent2I
from lsst.ip.isr import AssembleCcdTask
from astro_metadata_translator import ObservationInfo

logger = lsst.log.Log.getLogger("obs.lsst.assembly")


def attachRawWcsFromBoresight(exposure, dataIdForErrMsg=None):
    """Attach a WCS by extracting boresight, rotation, and camera geometry from
    an Exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Image object with attached metadata and detector components.

    Return
    ------
    attached : `bool`
        If True, a WCS component was successfully created and attached to
        ``exposure``.
    """
    md = exposure.getMetadata()
    # Use the generic version since we do not have a mapper available to
    # tell us a specific translator to use.
    obsInfo = ObservationInfo(md)
    visitInfo = MakeRawVisitInfoViaObsInfo.observationInfo2visitInfo(obsInfo, log=logger)
    exposure.getInfo().setVisitInfo(visitInfo)

    # LATISS (and likely others) need flipping, DC2 etc do not
    flipX = False
    if obsInfo.instrument in ("LATISS",):
        flipX = True

    if visitInfo.getBoresightRaDec().isFinite():
        exposure.setWcs(createInitialSkyWcs(visitInfo, exposure.getDetector(), flipX=flipX))
        return True

    if obsInfo.observation_type == "science":
        logger.warn("Unable to set WCS from header as RA/Dec/Angle are unavailable%s",
                    ("" if dataIdForErrMsg is None else " for dataId %s" % dataIdForErrMsg))
    return False


def fixAmpGeometry(inAmp, bbox, metadata, logCmd=None):
    """Make sure a camera geometry amplifier matches an on-disk bounding box.

    Bounding box differences that are consistent with differences in overscan
    regions are assumed to be overscan regions, which gives us enough
    information to correct the camera geometry.

    Parameters
    ----------
    inAmp : `lsst.afw.cameraGeom.Amplifier`
        Amplifier description from camera geometry.
    bbox : `lsst.geom.Box2I`
        The on-disk bounding box of the amplifer image.
    metadata : `lsst.daf.base.PropertyList`
        FITS header metadata from the amplifier HDU.
    logCmd : `function`, optional
        Call back to use to issue log messages about patching.  Arguments to
        this function should match arguments to be accepted by normal logging
        functions.  Warnings about bad EXTNAMES are always sent directly to
        the module-level logger.

    Return
    ------
    outAmp : `~lsst.afw.cameraGeom.Amplifier.Builder`
    modified : `bool`
        `True` if ``amp`` was modified; `False` otherwise.

    Raises
    ------
    RuntimeError
        Raised if the bounding boxes differ in a way that is not consistent
        with just a change in overscan.
    """
    if logCmd is None:
        # Define a null log command
        def logCmd(*args):
            return

    # check that the book-keeping worked and we got the correct EXTNAME
    extname = metadata.get("EXTNAME")
    predictedExtname = f"Segment{inAmp.getName()[1:]}"
    if extname is not None and predictedExtname != extname:
        logger.warning('expected to see EXTNAME == "%s", but saw "%s"', predictedExtname, extname)

    modified = False

    outAmp = inAmp.rebuild()
    if outAmp.getRawBBox() != bbox:  # Oh dear. cameraGeom is wrong -- probably overscan
        if outAmp.getRawDataBBox().getDimensions() != outAmp.getBBox().getDimensions():
            raise RuntimeError("Active area is the wrong size: %s v. %s" %
                               (outAmp.getRawDataBBox().getDimensions(), outAmp.getBBox().getDimensions()))

        logCmd("outAmp.getRawBBox() != data.getBBox(); patching. (%s v. %s)", outAmp.getRawBBox(), bbox)

        w, h = bbox.getDimensions()
        ow, oh = outAmp.getRawBBox().getDimensions()  # "old" (cameraGeom) dimensions
        #
        # We could trust the BIASSEC keyword, or we can just assume that
        # they've changed the number of overscan pixels (serial and/or
        # parallel).  As Jim Chiang points out, the latter is safer
        #
        fromCamGeom = outAmp.getRawHorizontalOverscanBBox()
        hOverscanBBox = Box2I(fromCamGeom.getBegin(),
                              Extent2I(w - fromCamGeom.getBeginX(), fromCamGeom.getHeight()))
        fromCamGeom = outAmp.getRawVerticalOverscanBBox()
        vOverscanBBox = Box2I(fromCamGeom.getBegin(),
                              Extent2I(fromCamGeom.getWidth(), h - fromCamGeom.getBeginY()))
        outAmp.setRawBBox(bbox)
        outAmp.setRawHorizontalOverscanBBox(hOverscanBBox)
        outAmp.setRawVerticalOverscanBBox(vOverscanBBox)
        #
        # This gets all the geometry right for the amplifier, but the size
        # of the untrimmed image will be wrong and we'll put the amp sections
        # in the wrong places, i.e.
        #   outAmp.getRawXYOffset()
        # will be wrong.  So we need to recalculate the offsets.
        #
        xRawExtent, yRawExtent = outAmp.getRawBBox().getDimensions()

        x0, y0 = outAmp.getRawXYOffset()
        ix, iy = x0//ow, y0/oh
        x0, y0 = ix*xRawExtent, iy*yRawExtent
        outAmp.setRawXYOffset(Extent2I(ix*xRawExtent, iy*yRawExtent))

        modified = True

    #
    # Check the "IRAF" keywords, but don't abort if they're wrong
    #
    # Only warn about the first amp, use debug for the others
    #
    d = metadata.toDict()
    detsec = bboxFromIraf(d["DETSEC"]) if "DETSEC" in d else None
    datasec = bboxFromIraf(d["DATASEC"]) if "DATASEC" in d else None
    biassec = bboxFromIraf(d["BIASSEC"]) if "BIASSEC" in d else None

    if detsec and outAmp.getBBox() != detsec:
        logCmd("DETSEC doesn't match (%s != %s)", outAmp.getBBox(), detsec)
    if datasec and outAmp.getRawDataBBox() != datasec:
        logCmd("DATASEC doesn't match for (%s != %s)", outAmp.getRawDataBBox(), detsec)
    if biassec and outAmp.getRawHorizontalOverscanBBox() != biassec:
        logCmd("BIASSEC doesn't match for (%s != %s)", outAmp.getRawHorizontalOverscanBBox(), detsec)

    return outAmp, modified


def assembleUntrimmedCcd(ccd, exposures):
    """Assemble an untrimmmed CCD from per-amp Exposure objects.

    Parameters
    ----------
    ccd : `~lsst.afw.cameraGeom.Detector`
        The detector geometry for this ccd that will be used as the
        framework for the assembly of the input amplifier exposures.
    exposures : sequence of `lsst.afw.image.Exposure`
        Per-amplifier images, in the same order as ``amps``.

    Returns
    -------
    ccd : `lsst.afw.image.Exposure`
        Assembled CCD image.
    """
    ampDict = {}
    for amp, exposure in zip(ccd, exposures):
        ampDict[amp.getName()] = exposure
    config = AssembleCcdTask.ConfigClass()
    config.doTrim = False
    assembleTask = AssembleCcdTask(config=config)
    return assembleTask.assembleCcd(ampDict)


def fixAmpsAndAssemble(ampExps, msg):
    """Fix amp geometry and assemble into exposure.

    Parameters
    ----------
    ampExps : sequence of `lsst.afw.image.Exposure`
        Per-amplifier images.
    msg : `str`
        Message to add to log and exception output.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        Exposure with the amps combined into a single image.

    Notes
    -----
    The returned exposure does not have any metadata or WCS attached.

    """
    if not len(ampExps):
        raise RuntimeError(f"Unable to read raw_amps for {msg}")

    ccd = ampExps[0].getDetector()      # the same (full, CCD-level) Detector is attached to all ampExps
    #
    # Check that the geometry in the metadata matches cameraGeom
    #
    warned = False

    def logCmd(s, *args):
        nonlocal warned
        if warned:
            logger.debug(f"{msg}: {s}", *args)
        else:
            logger.warn(f"{msg}: {s}", *args)
            warned = True

    # Rebuild the detector and the amplifiers to use their corrected geometry.
    tempCcd = ccd.rebuild()
    tempCcd.clear()
    for amp, ampExp in zip(ccd, ampExps):
        outAmp, _ = fixAmpGeometry(amp,
                                   bbox=ampExp.getBBox(),
                                   metadata=ampExp.getMetadata(),
                                   logCmd=logCmd)
        tempCcd.append(outAmp)

    ccd = tempCcd.finish()

    # Update the data to be combined to point to the newly rebuilt detector.
    for ampExp in ampExps:
        ampExp.setDetector(ccd)

    exposure = assembleUntrimmedCcd(ccd, ampExps)
    return exposure


def readRawAmps(fileName, detector):
    """Given a file name read the amps and attach the detector.

    Parameters
    ----------
    fileName : `str`
        The full path to a file containing data from a single CCD.
    detector : `lsst.afw.cameraGeom.Detector`
        Detector to associate with the amps.

    Returns
    -------
    ampExps : `list` of `lsst.afw.image.Exposure`
       All the individual amps read from the file.
    """
    amps = []
    for hdu in range(1, len(detector)+1):
        exp = afwImage.makeExposure(afwImage.makeMaskedImage(afwImage.ImageF(fileName, hdu=hdu)))
        exp.setDetector(detector)
        amps.append(exp)
    return amps
