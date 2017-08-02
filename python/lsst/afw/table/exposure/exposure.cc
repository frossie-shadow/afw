/*
 * LSST Data Management System
 * Copyright 2008-2017  AURA/LSST.
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

#include "pybind11/pybind11.h"

#include <memory>

#include "lsst/afw/detection/Psf.h"  // forward-declared by Exposure.h
#include "lsst/afw/fits.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/SpherePoint.h"
#include "lsst/afw/geom/polygon/Polygon.h"  // forward-declared by Exposure.h
#include "lsst/afw/image/ApCorrMap.h"       // forward-declared by Exposure.h
#include "lsst/afw/image/Calib.h"           // forward-declared by Exposure.h
#include "lsst/afw/image/VisitInfo.h"       // forward-declared by Exposure.h
#include "lsst/afw/image/Wcs.h"             // forward-declared by Exposure.h
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/python/catalog.h"
#include "lsst/afw/table/python/columnView.h"
#include "lsst/afw/table/python/sortedCatalog.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace afw {
namespace table {

namespace {

using PyExposureRecord = py::class_<ExposureRecord, std::shared_ptr<ExposureRecord>, BaseRecord>;
using PyExposureTable = py::class_<ExposureTable, std::shared_ptr<ExposureTable>, BaseTable>;
using PyExposureCatalog =
        py::class_<ExposureCatalogT<ExposureRecord>, std::shared_ptr<ExposureCatalogT<ExposureRecord>>,
                   SortedCatalogT<ExposureRecord>>;

PyExposureRecord declareExposureRecord(py::module &mod) {
    PyExposureRecord cls(mod, "ExposureRecord");
    cls.def("getId", &ExposureRecord::getId);
    cls.def("setId", &ExposureRecord::setId, "id"_a);
    cls.def("getBBox", &ExposureRecord::getBBox);
    cls.def("setBBox", &ExposureRecord::setBBox, "bbox"_a);
    cls.def("getTable", &ExposureRecord::getTable);
    cls.def_property_readonly("table", &ExposureRecord::getTable);
    cls.def("contains",
            (bool (ExposureRecord::*)(geom::SpherePoint const &, bool) const) & ExposureRecord::contains,
            "coord"_a, "includeValidPolygon"_a = false);
    cls.def("contains", (bool (ExposureRecord::*)(geom::Point2D const &, image::Wcs const &, bool) const) &
                                ExposureRecord::contains,
            "point"_a, "wcs"_a, "includeValidPolygon"_a = false);
    cls.def("getWcs", &ExposureRecord::getWcs);
    cls.def("setWcs", &ExposureRecord::setWcs, "wcs"_a);
    cls.def("getPsf", &ExposureRecord::getPsf);
    cls.def("setPsf", &ExposureRecord::setPsf, "psf"_a);
    cls.def("getCalib", &ExposureRecord::getCalib);
    cls.def("setCalib", &ExposureRecord::setCalib, "calib"_a);
    cls.def("getApCorrMap", &ExposureRecord::getApCorrMap);
    cls.def("setApCorrMap", &ExposureRecord::setApCorrMap, "appCorrMap"_a);
    cls.def("getValidPolygon", &ExposureRecord::getValidPolygon);

    // Workaround for DM-10289.
    cls.def("setValidPolygon",
            [](ExposureRecord &self, py::object polygon) {
                if (polygon == py::none()) {
                    self.setValidPolygon(nullptr);
                } else {
                    self.setValidPolygon(py::cast<std::shared_ptr<afw::geom::polygon::Polygon>>(polygon));
                }
            },
            "polygon"_a);

    cls.def("getVisitInfo", &ExposureRecord::getVisitInfo);
    cls.def("setVisitInfo", &ExposureRecord::setVisitInfo, "visitInfo"_a);
    return cls;
}

PyExposureTable declareExposureTable(py::module &mod) {
    PyExposureTable cls(mod, "ExposureTable");
    cls.def_static("make", &ExposureTable::make);
    cls.def_static("makeMinimalSchema", &ExposureTable::makeMinimalSchema);
    cls.def_static("checkSchema", &ExposureTable::checkSchema, "schema"_a);

    cls.def_static("getIdKey", &ExposureTable::getIdKey);
    cls.def_static("getBBoxMinKey", &ExposureTable::getBBoxMinKey);
    cls.def_static("getBBoxMaxKey", &ExposureTable::getBBoxMaxKey);

    cls.def("clone", &ExposureTable::clone);
    cls.def("makeRecord", &ExposureTable::makeRecord);
    cls.def("copyRecord", (std::shared_ptr<ExposureRecord> (ExposureTable::*)(BaseRecord const &)) &
                                  ExposureTable::copyRecord);
    cls.def("copyRecord",
            (std::shared_ptr<ExposureRecord> (ExposureTable::*)(BaseRecord const &, SchemaMapper const &)) &
                    ExposureTable::copyRecord);
    return cls;
}

PyExposureCatalog declareExposureCatalog(py::module &mod) {
    using Catalog = ExposureCatalogT<ExposureRecord>;
    table::python::declareSortedCatalog<ExposureRecord>(mod, "Exposure", true);

    // We need py::dynamic_attr() below to support our Python-side caching of the associated ColumnView.
    PyExposureCatalog cls(mod, "ExposureCatalog", py::dynamic_attr());

    cls.def(py::init<Schema const &>(), "schema"_a);
    cls.def(py::init<std::shared_ptr<ExposureTable> const &>(), "table"_a = std::shared_ptr<ExposureTable>());
    cls.def(py::init<Catalog const &>(), "other"_a);
    // Constructor taking C++ iterators not wrapped; we recommend .extend() (defined in pure Python) instead.
    cls.def_static("readFits", (Catalog(*)(std::string const &, int, int)) & Catalog::readFits, "filename"_a,
                   "hdu"_a = INT_MIN, "flags"_a = 0);
    cls.def_static("readFits", (Catalog(*)(fits::MemFileManager &, int, int)) & Catalog::readFits,
                   "manager"_a, "hdu"_a = INT_MIN, "flags"_a = 0);
    // readFits taking Fits objects not wrapped, because Fits objects are not wrapped.

    cls.def("subset", (Catalog (Catalog::*)(ndarray::Array<bool const, 1> const &) const) & Catalog::subset,
            "mask"_a);
    cls.def("subset",
            (Catalog (Catalog::*)(std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t) const) & Catalog::subset,
            "startd"_a, "stopd"_a, "step"_a);
    cls.def("subsetContaining",
            (Catalog (Catalog::*)(geom::SpherePoint const &, bool) const) & Catalog::subsetContaining, "coord"_a,
            "includeValidPolygon"_a = false);
    cls.def("subsetContaining",
            (Catalog (Catalog::*)(geom::Point2D const &, image::Wcs const &, bool) const) &
                    Catalog::subsetContaining,
            "point"_a, "wcs"_a, "includeValidPolygon"_a = false);
    return cls;
};

}  // anonymous namespace

PYBIND11_PLUGIN(exposure) {
    py::module mod("exposure");
    py::module::import("lsst.afw.table.simple");
    py::module::import("lsst.afw.geom");
    // afw.image and afw.detection cannot be imported due to circular dependencies

    auto clsExposureRecord = declareExposureRecord(mod);
    auto clsExposureTable = declareExposureTable(mod);
    auto clsExposureColumnView = table::python::declareColumnView<ExposureRecord>(mod, "Exposure");
    auto clsExposureCatalog = declareExposureCatalog(mod);

    clsExposureRecord.attr("Table") = clsExposureTable;
    clsExposureRecord.attr("ColumnView") = clsExposureColumnView;
    clsExposureRecord.attr("Catalog") = clsExposureCatalog;
    clsExposureTable.attr("Record") = clsExposureRecord;
    clsExposureTable.attr("ColumnView") = clsExposureColumnView;
    clsExposureTable.attr("Catalog") = clsExposureCatalog;
    clsExposureCatalog.attr("Record") = clsExposureRecord;
    clsExposureCatalog.attr("Table") = clsExposureTable;
    clsExposureCatalog.attr("ColumnView") = clsExposureColumnView;

    return mod.ptr();
}
}
}
}  // namespace lsst::afw::table
