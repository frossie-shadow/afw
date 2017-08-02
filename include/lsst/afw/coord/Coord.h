// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#if !defined(LSST_AFW_COORD_COORD_H)
#define LSST_AFW_COORD_COORD_H
/*
 * Functions to handle coordinates
 */
#include <iostream>
#include <limits>
#include <map>
#include <memory>

#include "lsst/base.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/SpherePoint.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/coord/Observatory.h"
#include "lsst/daf/base.h"

namespace lsst {
namespace afw {
namespace coord {

/*
 * Information about the coordinate system we support
 */
enum CoordSystem { UNKNOWN = -1, FK5, ICRS, GALACTIC, ECLIPTIC, TOPOCENTRIC };
/**
 * A utility function to get the enum value of a coordinate system from a string name.
 */
CoordSystem makeCoordEnum(std::string const system);

/**
 * Return average of a list of coordinates
 *
 * @param[in] coords  list of coords to average
 */
std::shared_ptr<geom::SpherePoint> averageCoord(
        std::vector<std::shared_ptr<geom::SpherePoint>> const& coords);

/**
 * Get the offset on the tangent plane from one SpherePoint to another
 *
 * This is suitable only for small angles.
 *
 * @param from  Coordinate from which to compute offset
 * @param to  Coordinate to which to compute offset
 * @returns pair of Angles: Longitude and Latitude offsets
 */
std::pair<geom::Angle, geom::Angle> getTangentPlaneOffset(geom::SpherePoint const &from,
                                                          geom::SpherePoint const &to);

/**
 * Convert a dd:mm:ss string to Angle
 *
 * @param dms Coord as a string in dd:mm:ss format
 */
lsst::afw::geom::Angle dmsStringToAngle(std::string const dms);
/// Convert a hh:mm:ss string to Angle
lsst::afw::geom::Angle hmsStringToAngle(std::string const hms);
/**
 * a Function to convert a coordinate in decimal degrees to a string with form dd:mm:ss
 *
 * @todo allow a user specified format
 */
std::string angleToDmsString(lsst::afw::geom::Angle const deg);
/// a function to convert decimal degrees to a string with form hh:mm:ss.s
std::string angleToHmsString(lsst::afw::geom::Angle const deg);

}
}
}

#endif
