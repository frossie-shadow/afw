/*
 * LSST Data Management System
 * Copyright 2008-2016  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "lsst/afw/coord/Coord.h"

namespace py = pybind11;
using namespace py::literals;

using namespace lsst::afw::coord;

PYBIND11_DECLARE_HOLDER_TYPE(MyCoordType, std::shared_ptr<MyCoordType>);

PYBIND11_PLUGIN(_coord) {
    py::module mod("_coord", "Python wrapper for afw _coord library");

    /* Types and enums */
    py::enum_<CoordSystem>(mod, "CoordSystem")
            .value("UNKNOWN", CoordSystem::UNKNOWN)
            .value("FK5", CoordSystem::FK5)
            .value("ICRS", CoordSystem::ICRS)
            .value("GALACTIC", CoordSystem::GALACTIC)
            .value("ECLIPTIC", CoordSystem::ECLIPTIC)
            .value("TOPOCENTRIC", CoordSystem::TOPOCENTRIC)
            .export_values();

    mod.def("makeCoordEnum", makeCoordEnum);

    mod.def("averageCoord", averageCoord, py::arg("coords"));

    mod.def("getTangentPlaneOffset", getTangentPlaneOffset, "from"_a, "to"_a);

    return mod.ptr();
}
