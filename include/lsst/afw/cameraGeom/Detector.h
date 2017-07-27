/*
 * LSST Data Management System
 * Copyright 2014 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#if !defined(LSST_AFW_CAMERAGEOM_DETECTOR_H)
#define LSST_AFW_CAMERAGEOM_DETECTOR_H

#include <string>
#include <vector>
#include <unordered_map>
#include "lsst/base.h"
#include "lsst/afw/geom/TransformMap.h"
#include "lsst/afw/table/AmpInfo.h"
#include "lsst/afw/cameraGeom/CameraSys.h"
#include "lsst/afw/cameraGeom/CameraPoint.h"
#include "lsst/afw/cameraGeom/Orientation.h"

namespace lsst {
namespace afw {
namespace cameraGeom {

/**
 * Type of imaging detector
 */
enum DetectorType {
    SCIENCE,
    FOCUS,
    GUIDER,
    WAVEFRONT,
};

/**
 * Information about a CCD or other imaging detector
 *
 * Supports conversion of CameraPoint between FOCAL_PLANE and pixel-based coordinate systems.
 * Also an iterator over amplifiers (in C++ use begin(), end(), in Python use "for amplifier in detector").
 *
 * @todo: this would probably be a bit more robust if it used a ConstAmpInfoCatalog
 * (a catalog with const records) but I don't think const catalogs really work yet;
 * for instance it is not possible to construct one from a non-const catalog,
 * so I don't know how to construct one.
 *
 * @note code definitions use lsst::afw::table:: instead of table:: because the latter confused swig
 * when I tried it. This is a known issue: ticket #2461.
 */
class Detector {
public:
    typedef ndarray::Array<float const, 2> CrosstalkMatrix;

    /**
     * Make a Detector
     *
     * @warning
     * * The keys for the detector-specific coordinate systems in the transform registry
     *   must include the detector name (even though this is redundant).
     *
     * @throws lsst::pex::exceptions::InvalidParameterError if:
     * - any amplifier names are not unique
     * - any CamerSys in transformMap has a detector name other than "" or this detector's name
     */
    explicit Detector(
            std::string const &name,    ///< name of detector's location in the camera
            int id,                     ///< detector integer ID; used as keys in some tables
            DetectorType type,          ///< type of detector
            std::string const &serial,  ///< serial "number" that identifies the physical detector
            geom::Box2I const &bbox,    ///< bounding box
            lsst::afw::table::AmpInfoCatalog const &ampInfoCatalog,  ///< catalog of amplifier information
            Orientation const &orientation,   ///< detector position and orientation in focal plane
            geom::Extent2D const &pixelSize,  ///< pixel size (mm)
            CameraTransformMap::Transforms const &transforms,  ///< map of
            ///< CameraSys: lsst::afw::geom::XYTransform, where each transform's "forwardTransform" method
            ///< transforms from PIXELS to the specified camera system
            CrosstalkMatrix const &crosstalk=CrosstalkMatrix() ///< matrix of crosstalk coefficients
            );

    ~Detector() {}

    /** Get the detector name */
    std::string getName() const { return _name; }

    /** Get the detector ID */
    int getId() const { return _id; }

    DetectorType getType() const { return _type; }

    /** Get the detector serial "number" */
    std::string getSerial() const { return _serial; }

    /** Get the bounding box */
    lsst::afw::geom::Box2I getBBox() const { return _bbox; }

    /** Get the corners of the detector in the specified camera coordinate system */
    std::vector<geom::Point2D> getCorners(CameraSys const &cameraSys) const;

    /** Get the corners of the detector in the specified camera coordinate system prefix */
    std::vector<geom::Point2D> getCorners(CameraSysPrefix const &cameraSysPrefix) const;

    /** Get the center of the detector in the specified camera coordinate system */
    CameraPoint getCenter(CameraSys const &cameraSys) const;

    /** Get the center of the detector in the specified camera coordinate system prefix */
    CameraPoint getCenter(CameraSysPrefix const &cameraSysPrefix) const;

    /** Get the amplifier information catalog */
    lsst::afw::table::AmpInfoCatalog const getAmpInfoCatalog() const { return _ampInfoCatalog; }

    /** Get detector's orientation in the focal plane */
    Orientation const getOrientation() const { return _orientation; }

    /** Get size of pixel along (mm) */
    geom::Extent2D getPixelSize() const { return _pixelSize; }

    /** Get the transform registry */
    CameraTransformMap const getTransformMap() const { return _transformMap; }

    /** Have we got crosstalk coefficients? */
    bool hasCrosstalk() const {
        return !(_crosstalk.isEmpty() || _crosstalk.getShape() == ndarray::makeVector(0, 0));
    }

    /** Get the crosstalk coefficients */
    CrosstalkMatrix const getCrosstalk() const { return _crosstalk; }

    /** Get iterator to beginning of amplifier list */
    lsst::afw::table::AmpInfoCatalog::const_iterator begin() const { return _ampInfoCatalog.begin(); }

    /** Get iterator to end of amplifier list */
    lsst::afw::table::AmpInfoCatalog::const_iterator end() const { return _ampInfoCatalog.end(); }

    /**
     * Get the amplifier specified by index
     *
     * @throws std::out_of_range if index is out of range
     */
    lsst::afw::table::AmpInfoRecord const &operator[](size_t i) const { return _ampInfoCatalog.at(i); }

    /**
     * Get the amplifier specified by name
     *
     * @throws lsst::pex::exceptions::InvalidParameterError if no such amplifier
     */
    lsst::afw::table::AmpInfoRecord const &operator[](std::string const &name) const;

    /**
     * Get the amplifier specified by index, returning a shared pointer to an AmpInfo record
     *
     * @warning Intended only for use by pybind11. This exists because
     * Operator[] returns a data type that is difficult for pybind11 to use.
     * Since we have it, we also take advantage of the fact that it is only for pybind11
     * to support negative indices in the python style.
     *
     * @param[in] i  Ampifier index; if < 0 then treat as an offset from the end (the Python convention)
     *
     * @throws std::out_of_range if index is out of range
     */
    std::shared_ptr<lsst::afw::table::AmpInfoRecord const> _get(int i) const;

    /**
     * Get the amplifier specified by name, returning a shared pointer to an AmpInfo record
     *
     * @warning Intended only for internal and pybind11 use. This exists because
     * Operator[] returns a data type that is difficult for pybind11 to use.
     *
     * @param[in] name  Amplifier name
     *
     * @throws std::out_of_range if index is out of range
     */
    std::shared_ptr<lsst::afw::table::AmpInfoRecord const> _get(std::string const &name) const;

    /**
     * Get the number of amplifiers. Renamed to `__len__` in Python.
     */
    size_t size() const { return _ampInfoCatalog.size(); }

    /** Does the specified camera coordinate system exist in the transform registry? */
    bool hasTransform(CameraSys const &cameraSys) const;

    /** Does the specified camera coordiante system prefix exist in the transform registry? */
    bool hasTransform(CameraSysPrefix const &cameraSysPrefix) const;

    /**
     * Get an XYTransform that transforms from `cameraSys` to the native system in the forward direction
     *
     * @param[in] cameraSys  camera coordinate system
     * @returns a shared_ptr to an lsst::afw::XYTransform
     *
     * @throws pex::exceptions::InvalidParameterError if `cameraSys` is unknown
     */
    std::shared_ptr<afw::geom::XYTransform const> getTransform(CameraSys const &cameraSys) const;

    /**
     * Get an XYTransform that transforms from cameraSysPrefix to the native system in the forward direction
     *
     * @param[in] cameraSysPrefix  camera coordinate system prefix
     * @returns a shared_ptr to an lsst::afw::geom::XYTransform
     *
     * @throws pex::exceptions::InvalidParameterError if coordSys is unknown
     */
    std::shared_ptr<afw::geom::XYTransform const> getTransform(CameraSysPrefix const &cameraSysPrefix) const;

    /**
     * Make a CameraPoint from a point and a camera system
     *
     * @param[in] point  2D point
     * @param[in] cameraSys  Camera coordinate system
     *
     * @note the CameraSysPrefix version needs the detector name, which is why this is not static.
     */
    CameraPoint makeCameraPoint(geom::Point2D const &point,  ///< 2-d point
                                CameraSys const &cameraSys   ///< coordinate system
                                ) const {
        return CameraPoint(point, cameraSys);
    }

    /**
     * Make a CameraPoint from a point and a camera system prefix
     *
     * @param[in] point  2D point
     * @param[in] cameraSysPrefix  Camera coordinate system prefix
     */
    CameraPoint makeCameraPoint(geom::Point2D const &point, CameraSysPrefix const &cameraSysPrefix) const {
        return CameraPoint(point, makeCameraSys(cameraSysPrefix));
    }

    /**
     * Get a coordinate system from a coordinate system (return input unchanged and untested)
     *
     * @param[in] cameraSys  Camera coordinate system
     * @return `cameraSys` unchanged
     *
     * @note the CameraSysPrefix version needs the detector name, which is why this is not static.
     */
    CameraSys const makeCameraSys(CameraSys const &cameraSys) const { return cameraSys; }

    /**
     * Get a coordinate system from a detector system prefix (add detector name)
     *
     * @param[in] cameraSysPrefix  Camera coordinate system prefix
     * @return `cameraSysPrefix` with the detector name added
     */
    CameraSys const makeCameraSys(CameraSysPrefix const &cameraSysPrefix) const {
        return CameraSys(cameraSysPrefix.getSysName(), _name);
    }

    /**
     * Convert a CameraPoint from one coordinate system to another
     *
     * @param[in] fromCameraPoint  Camera point to transform
     * @param[in] toSys  To camera coordinate system
     * @return The transformed camera point
     *
     * @throws pex::exceptions::InvalidParameterError if from or to coordinate system is unknown
     */
    CameraPoint transform(CameraPoint const &fromCameraPoint, CameraSys const &toSys) const {
        return CameraPoint(
                _transformMap.transform(fromCameraPoint.getPoint(), fromCameraPoint.getCameraSys(), toSys),
                toSys);
    }

    /**
     * Convert a CameraPoint from one coordinate system to a coordinate system prefix
     *
     * @param[in] fromCameraPoint  Camera point to transform
     * @param[in] toSys  To camera coordinate system prefix (this detector is used)
     * @return The transformed camera point; the detector will be this detector
     *
     * @throws pex::exceptions::InvalidParameterError if from or to coordinate system is unknown
     */
    CameraPoint transform(CameraPoint const &fromCameraPoint, CameraSysPrefix const &toSys) const {
        return transform(fromCameraPoint, makeCameraSys(toSys));
    }

private:
    typedef std::unordered_map<std::string, table::AmpInfoCatalog::const_iterator> _AmpInfoMap;
    /**
     * Finish constructing this object
     *
     * Set _ampNameIterMap from _ampInfoCatalog
     * Check detector name in the CoordSys in the transform registry
     *
     * @throws lsst::pex::exceptions::InvalidParameterError if:
     * - any amplifier names are not unique
     * - any CamerSys in transformMap has a detector name other than "" or this detector's name
     */
    void _init();

    std::string _name;                      ///< name of detector's location in the camera
    int _id;                                ///< detector numeric ID
    DetectorType _type;                     ///< type of detectorsize_t
    std::string _serial;                    ///< serial "number" that identifies the physical detector
    geom::Box2I _bbox;                      ///< bounding box
    table::AmpInfoCatalog _ampInfoCatalog;  ///< list of amplifier data
    _AmpInfoMap _ampNameIterMap;            ///< map of amplifier name: catalog iterator
    Orientation _orientation;               ///< position and orientation of detector in focal plane
    geom::Extent2D _pixelSize;              ///< pixel size (mm)
    CameraTransformMap _transformMap;       ///< registry of coordinate transforms
    CrosstalkMatrix _crosstalk;             ///< crosstalk coefficients
};
}  // namespace cameraGeom
}  // namespace afw
}  // namespace lsst

#endif
