# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Utilities and a class to help the camera team use LSST DM camera utilities

Mostly support mapping between amp/ccd/focal plane coordinates, but
also the odd utility function
"""
import lsst.geom as geom
import lsst.afw.cameraGeom as cameraGeom

__all__ = ["LsstCameraTransforms", "getAmpImage", "channelToAmp"]


class LsstCameraTransforms():
    """Package up DM's camera geometry in ways useful to the camera team"""
    def __init__(self, camera, detectorName=None):
        r"""Make a new LsstCameraTransforms

        Parameters
        ----------
        camera : `lsst.afw.cameraGeom.Camera`
           The object describing a/the LSST Camera
        detectorName : `str`
           A default detector name (e.g. "R12_S01"); used when you don't
           pass a detectorName (None to disable the default detectorName)
        """
        self.camera = camera
        self.__detectorName = detectorName

    def setDetectorName(self, detectorName):
        r"""Set the default detector

        Parameters
        ----------
        detectorName : `str`
           The detector (e.g. "R34_S00") to use when none is passed to a mathod
           Set to None to disable the use of a default
        """

        self.__detectorName = detectorName

    def getDetector(self, detectorName=None):
        r"""Return a specified Detector, or the default Detector

        Parameters
        ----------
        detectorName : `str`
           Name of detector (or default from setDetectorName() if None)

        Returns
        -------
        detector : `lsst.afw.cameraGeom.Detector`
           The requested detector
        """

        if detectorName is None:
            detectorName = self.__detectorName

        if detectorName is None:
            raise RuntimeError("Please specify a detector name, "
                               "or set a default with LsstCameraTransforms.setDetectorName()")

        return self.camera[detectorName]

    def ampPixelToCcdPixel(self, ampX, ampY, channel, detectorName=None):
        r"""Given raw amplifier position return position on full detector

        Parameters
        ----------
        ampX : `int`
           column on amp segment
        ampY : `int`
           row on amp segment
        channel: `int`
           Channel number of amplifier (1-indexed; identical to HDU)
        detectorName : `str`
           Name of detector (or default from setDetectorName() if None)

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0)

        Returns
        -------
        ccdX : `int`
           The column pixel position relative to the corner of the detector
        ccdY : `int`
           The row pixel position relative to the corner of the detector
        """

        return ampPixelToCcdPixel(ampX, ampY, self.getDetector(detectorName), channel)

    def ccdPixelToAmpPixel(self, ccdX, ccdY, detectorName=None):
        r"""Given position within a detector, return the amplifier position

        Parameters
        ----------
        ccdX : `int`
           column pixel position within detector
        ccdY : `int`
           row pixel position within detector
        detectorName : `str`
           Name of detector (or default from setDetectorName() if None)

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0)

        Returns
        -------
        channel: `int`
           Channel number of amplifier (1-indexed; identical to HDU)
        ampX : `int`
           The column coordinate relative to the corner of the single-amp image
        ampY : `int`
           The row coordinate relative to the corner of the single-amp image

        Raises
        ------
        RuntimeError
            If the requested pixel doesn't lie on the detector
        """

        amp, ampXY = ccdPixelToAmpPixel(geom.PointD(ccdX, ccdY), self.getDetector(detectorName))

        ampX, ampY = ampXY
        return ampToChannel(amp), ampX, ampY

    def ccdPixelToFocalMm(self, ccdX, ccdY, detectorName=None):
        r"""Given position within a detector return the focal plane position

        Parameters
        ----------
        ccdX : `int`
           column pixel position within detector
        ccdY : `int`
           row pixel position within detector
        detectorName : `str`
           Name of detector (or default from setDetectorName() if None)

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0), and CCD columns are taken to be vertical, with "R00"
        in the bottom left of the image

        Returns
        -------
        focalPlaneX : `float`
           The x-focal plane position (mm) relative to the centre of the camera
        focalPlaneY : `float`
           The y-focal plane position (mm) relative to the centre of the camera
        """
        detector = self.getDetector(detectorName)

        return detector.transform(geom.Point2D(ccdX, ccdY), cameraGeom.PIXELS, cameraGeom.FOCAL_PLANE)

    def ampPixelToFocalMm(self, ampX, ampY, channel, detectorName=None):
        r"""Given position within an amplifier return the focal plane position

        ampX : `int`
           column on amp segment
        ampY : `int`
           row on amp segment
        channel: `int`
           Channel number of amplifier (1-indexed; identical to HDU)
        detectorName : `str`
           Name of detector (or default from setDetectorName() if None)

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0), and CCD columns are taken to be vertical, with "R00"
        in the bottom left of the image

        Returns
        -------
        focalPlaneX : `float`
           The x-focal plane position (mm) relative to the centre of the camera
        focalPlaneY : `float`
           The y-focal plane position (mm) relative to the centre of the camera

        Raises
        ------
        RuntimeError
            If the requested pixel doesn't lie on the detector
        """

        detector = self.getDetector(detectorName)

        ccdX, ccdY = ampPixelToCcdPixel(ampX, ampY, detector, channel)

        return detector.transform(geom.Point2D(ccdX, ccdY), cameraGeom.PIXELS, cameraGeom.FOCAL_PLANE)

    def focalMmToCcdPixel(self, focalPlaneX, focalPlaneY):
        r"""Given focal plane position return the detector position

        Parameters
        ----------
        focalPlaneX : `float`
           The x-focal plane position (mm) relative to the centre of the camera
        focalPlaneY : `float`
           The y-focal plane position (mm) relative to the centre of the camera

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0) and CCD columns are taken to be vertical, with "R00" in
        the bottom left of the image

        Returns
        -------
        detectorName : `str`
           The name of the detector
        ccdX : `int`
           The column pixel position relative to the corner of the detector
        ccdY : `int`
           The row pixel position relative to the corner of the detector

        Raises
        ------
        RuntimeError
            If the requested position doesn't lie on a detector
        """
        detector, ccdXY = focalMmToCcdPixel(self.camera, geom.PointD(focalPlaneX, focalPlaneY))
        ccdX, ccdY = ccdXY

        return detector.getName(), ccdX, ccdY

    def focalMmToAmpPixel(self, focalPlaneX, focalPlaneY):
        r"""Given a focal plane position plane return the amplifier position

        Parameters
        ----------
        focalPlaneX : `float`
           The x-focal plane position (mm) relative to the centre of the camera
        focalPlaneY : `float`
           The y-focal plane position (mm) relative to the centre of the camera

        N.b. all pixel coordinates have the centre of the bottom-left pixel
        at (0.0, 0.0) and CCD columns are taken to be vertical, with "R00" in
        the bottom left of the image

        Returns
        -------
        detectorName : `str`
           The name of the detector
        channel: `int`
           Channel number of amplifier (1-indexed; identical to HDU)
        ampX : `int`
           The column coordinate relative to the corner of the single-amp image
        ampY : `int`
           The row coordinate relative to the corner of the single-amp image

        Raises
        ------
        RuntimeError
            If the requested position doesn't lie on a detector
        """

        detector, ccdXY = focalMmToCcdPixel(self.camera, geom.PointD(focalPlaneX, focalPlaneY))
        amp, ampXY = ccdPixelToAmpPixel(ccdXY, detector)

        ampX, ampY = ampXY
        return detector.getName(), ampToChannel(amp), ampX, ampY

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def getAmpImage(butler, dataId, channel):
    r"""Return a amplifier image

    Parameters
    ----------
    butler : `lsst.daf.persistence.butler.Butler`
       The butler which will perform the I/O
    dataId : `dict` or `dafPersist.butler.DataId`
       The dataId specifying the desired detector
    channel : `int`
       The 1-indexed channel ID for the desired amplifier

    Returns
    -------
    detectorName : `lsst.afw.image.ExposureF`
       The desired image

    Notes
    -----
    The Gen2 butler can't quite do this natively as it doesn't
    know how to only lookup things that it doesn't know, so we
    so that for it, poor lamb

    """
    keys = ["run", "detector"]
    did = dataId.copy()
    did.update(dict(zip(keys, butler.queryMetadata("raw", keys, did)[0])))
    # the only snap is snap==0
    did["snap"] = 0

    return butler.get('raw_amp', did, channel=channel)


def ampToChannel(amp):
    r"""Given an Amplifier, return the channel

    Parameters
    ----------
    amp : `lsst.afw.table.AmpInfoRecord`
       The amplifier in question

    Returns
    -------
    channel : `int`
       The 1-indexed channel ID for the desired amplifier
    """
    return amp.get("hdu")


def channelToAmp(detector, channel):
    r"""Given a Detector and channel, return the Amplifier

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`
        The requested detector
    channel : `int`
       The 1-indexed channel ID for the desired amplifier

    Returns
    -------
    amp : `lsst.afw.table.AmpInfoRecord`
       The amplifier in question
    """
    return [amp for amp in detector if amp.get("hdu") == channel][0]


def ampPixelToCcdPixel(x, y, detector, channel):
    r"""Given a position in a raw amplifier return position on full detector

    Parameters
    ----------
    x : `int`
       column on amp segment
    y : `int`
       row on amp segment
    detector : `lsst.afw.cameraGeom.Detector`
        The requested detector
    channel: `int`
       Channel number of amplifier (1-indexed; identical to HDU)

    Returns
    -------
    ccdX : `int`
       The column pixel position relative to the corner of the detector
    ccdY : `int`
       The row pixel position relative to the corner of the detector
    """

    amp = channelToAmp(detector, channel)
    bbox = amp.getRawDataBBox()

    # Allow for flips (due e.g. to physical location of the amplifiers)
    x, y = geom.PointI(x, y)            # definitely ints
    w, h = bbox.getDimensions()
    if amp.getRawFlipX():
        x = w - x - 1

    if amp.getRawFlipY():
        y = h - y - 1

    return amp.getBBox().getBegin() + geom.ExtentI(x, y)


def ccdPixelToAmpPixel(xy, detector):
    r"""Given an position within a detector return position within an amplifier

    Parameters
    ----------
    xy : `lsst.geom.PointD`
       pixel position within detector
    detector : `lsst.afw.cameraGeom.Detector`
        The requested detector

    N.b. all pixel coordinates have the centre of the bottom-left pixel
    at (0.0, 0.0)

    Returns
    -------
    amp : `lsst.afw.table.AmpInfoRecord`
       The amplifier that the pixel lies in
    ampXY : `lsst.geom.PointI`
       The pixel coordinate relative to the corner of the single-amp image

    Raises
    ------
    RuntimeError
        If the requested pixel doesn't lie on the detector
    """
    xy = geom.PointI(xy)                # use pixel coordinates

    found = False
    for amp in detector:
        if amp.getBBox().contains(xy):
            x, y = xy - amp.getBBox().getBegin()  # coordinates within data segment
            found = True
            break

    if not found:
        raise RuntimeError("Point (%g, %g) does not lie on detector %s" % (xy[0], xy[1], detector.getName()))

    # Allow for flips (due e.g. to physical location of the amplifiers)
    w, h = amp.getBBox().getDimensions()

    if amp.getRawFlipX():
        x = w - x - 1

    if amp.getRawFlipY():
        y = h - y - 1

    return amp, geom.PointI(x, y)


def focalMmToCcdPixel(camera, focalPlaneXY):
    r"""Given an position in the focal plane, return the position on a detector

    Parameters
    ----------
    camera : `lsst.afw.cameraGeom.Camera`
       The object describing a/the LSST Camera
    focalPlaneXY : `lsst.geom.PointD`
       The focal plane position (mm) relative to the centre of the camera

    N.b. all pixel coordinates have the centre of the bottom-left pixel
    at (0.0, 0.0) and CCD columns are taken to be vertical, with "R00"
    in the bottom left of the image

    Returns
    -------
    detector : `lsst.afw.cameraGeom.Detector`
       The requested detector
    ccdPos : `lsst.geom.Point2D`
       The pixel position relative to the corner of the detector

    Raises
    ------
    RuntimeError
        If the requested position doesn't lie on a detector
    """

    for detector in camera:
        ccdXY = detector.transform(focalPlaneXY, cameraGeom.FOCAL_PLANE, cameraGeom.PIXELS)
        if geom.BoxD(detector.getBBox()).contains(ccdXY):
            return detector, ccdXY

    raise RuntimeError("Failed to map focal plane position (%.3f, %.3f) to a detector" % (focalPlaneXY))
