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

/*
 * Provide functions to handle coordinates
 *
 * Most (nearly all) algorithms adapted from Astronomical Algorithms, 2nd ed. (J. Meeus)
 *
 */
#include <cmath>
#include <limits>
#include <cstdio>
#include <iomanip>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "lsst/pex/exceptions.h"
#include "boost/algorithm/string.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/format.hpp"

#include "lsst/afw/coord/Coord.h"
#include "lsst/daf/base/DateTime.h"

namespace ex = lsst::pex::exceptions;
namespace dafBase = lsst::daf::base;

namespace lsst {
namespace afw {
namespace coord {

namespace {

typedef std::map<std::string, CoordSystem> CoordSystemMap;

CoordSystemMap const getCoordSystemMap() {
    CoordSystemMap idMap;
    idMap["FK5"] = FK5;
    idMap["ICRS"] = ICRS;
    idMap["ECLIPTIC"] = ECLIPTIC;
    idMap["GALACTIC"] = GALACTIC;
    idMap["ELON"] = ECLIPTIC;
    idMap["GLON"] = GALACTIC;
    idMap["TOPOCENTRIC"] = TOPOCENTRIC;
    return idMap;
}

}  // end anonymous namespace

CoordSystem makeCoordEnum(std::string const system) {
    static CoordSystemMap idmap = getCoordSystemMap();
    if (idmap.find(system) != idmap.end()) {
        return idmap[system];
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterError, "System " + system + " not defined.");
    }
}

namespace {

/*
 * A local class to handle dd:mm:ss coordinates
 *
 * This class allows a decimal or dd:mm:ss coordinate to be
 * disassembed into d, m, and s components.
 * It's in an anonymous namespace, but there are public functions
 * which perform the transformations directly:
 *
 * --> std::string dmsStr = degreesToDmsString(double deg);
 * --> double deg = dmsStringToDegrees(std::string dms);
 */
class Dms {
public:
    Dms(){};

    // note that isSouth is needed to specify coords between dec = 0, and dec = -1
    // otherwise, d = -0 gets carried as d = 0 ... need a way to specify it explicitly
    Dms(int const d, int const m, double const s, bool const isSouth = false) {
        sign = (d < 0 || isSouth) ? -1 : 1;
        deg = std::abs(d);
        min = m;
        sec = s;
    };
    // unit could be "degrees" or "hours"
    Dms(geom::Angle const deg00, geom::AngleUnit const unit = geom::degrees) {
        double deg0 = deg00.asAngularUnits(unit);
        double const absVal = std::fabs(deg0);
        sign = (deg0 >= 0) ? 1 : -1;
        deg = static_cast<int>(std::floor(absVal));
        min = static_cast<int>(std::floor((absVal - deg) * 60.0));
        sec = ((absVal - deg) * 60.0 - min) * 60.0;
    }

    int deg;
    int min;
    double sec;
    int sign;
};


}  // end anonymous namespace

static std::string angleToXmsString(geom::Angle const a, geom::AngleUnit const unit) {
    Dms dms(a, unit);

    // make sure rounding won't give 60.00 for sec or min
    if ((60.00 - dms.sec) < 0.005) {
        dms.sec = 0.0;
        dms.min += 1;
        if (dms.min == 60) {
            dms.min = 0;
            dms.deg += 1;
            if (dms.deg == 360) {
                dms.deg = 0;
            }
        }
    }

    std::string fmt("%02d:%02d:%05.2f");
    std::string s = (boost::format(fmt) % dms.deg % dms.min % dms.sec).str();
    if (dms.sign < 0) {
        s = "-" + s;
    }
    return s;
}

std::string angleToDmsString(geom::Angle const a) { return angleToXmsString(a, geom::degrees); }

std::string angleToHmsString(geom::Angle const a) { return angleToXmsString(a, geom::hours); }

/**
 * @internal Convert a XX:mm:ss string to Angle
 *
 * @param dms Coord as a string in dd:mm:ss format
 * @param unit the units assumed for the first part of `dms`. The second and third
 *             parts shall be defined to be 1/60 and 1/3600 of `unit`, respectively.
 */
static geom::Angle xmsStringToAngle(std::string const dms, geom::AngleUnit unit) {
    if (dms.find(":") == std::string::npos) {
        throw LSST_EXCEPT(ex::InvalidParameterError,
                          (boost::format("String is not in xx:mm:ss format: %s") % dms).str());
    }
    std::vector<std::string> elements;
    boost::split(elements, dms, boost::is_any_of(":"));
    if (elements.size() != 3) {
        throw LSST_EXCEPT(ex::InvalidParameterError,
                          (boost::format("Could not parse string as xx:mm:ss format: %s") % dms).str());
    }
    int const deg = abs(atoi(elements[0].c_str()));
    int const min = atoi(elements[1].c_str());
    double const sec = atof(elements[2].c_str());

    geom::Angle ang = (deg + min / 60.0 + sec / 3600.0) * unit;
    if ((elements[0].c_str())[0] == '-') {
        ang *= -1.0;
    }
    return ang;
}

geom::Angle hmsStringToAngle(std::string const hms) { return xmsStringToAngle(hms, geom::hours); }

geom::Angle dmsStringToAngle(std::string const dms) { return xmsStringToAngle(dms, geom::degrees); }

std::shared_ptr<geom::SpherePoint> averageCoord(
        std::vector<std::shared_ptr<geom::SpherePoint>> const &coords) {
    if (coords.size() == 0) {
        throw LSST_EXCEPT(ex::LengthError, "No coordinates provided to average");
    }

    geom::Point3D sum(0, 0, 0);
    geom::Point3D corr(0, 0, 0);  // Kahan summation correction
    for (auto &&cc : coords) {
        geom::Point3D const point = cc->getVector();
        // Kahan summation
        geom::Extent3D const add = point - corr;
        geom::Point3D const temp = sum + add;
        corr = (temp - geom::Extent3D(sum)) - add;
        sum = temp;
    }
    sum.scale(1.0 / coords.size());
    return std::make_shared<geom::SpherePoint>(sum);
}

std::pair<geom::Angle, geom::Angle> getTangentPlaneOffset(geom::SpherePoint const &from,
                                                          geom::SpherePoint const &to) {
    geom::Angle const alpha1 = to.getRa();
    geom::Angle const delta1 = to.getDec();
    geom::Angle const alpha2 = from.getRa();
    geom::Angle const delta2 = from.getDec();

    double const sinDelta1 = std::sin(delta1);
    double const cosDelta1 = std::cos(delta1);
    double const sinDelta2 = std::sin(delta2);
    double const cosDelta2 = std::cos(delta2);
    double const cosAlphaDiff = std::cos(alpha2 - alpha1);
    double const sinAlphaDiff = std::sin(alpha2 - alpha1);

    double const div = cosDelta1 * cosAlphaDiff * cosDelta2 + sinDelta1 * sinDelta2;
    double const xi = cosDelta1 * sinAlphaDiff / div;
    double const eta = (cosDelta1 * cosAlphaDiff * sinDelta2 - sinDelta1 * cosDelta2) / div;

    return std::make_pair(xi * geom::radians, eta * geom::radians);
}

}
}
}  // end lsst::afw::coord
