#ifndef AFW_TABLE_PYBIND11_SORTEDCATALOG_H_INCLUDED
#define AFW_TABLE_PYBIND11_SORTEDCATALOG_H_INCLUDED
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

#include "lsst/afw/table/SortedCatalog.h"
#include "lsst/afw/table/python/catalog.h"

namespace lsst {
namespace afw {
namespace table {
namespace python {

template <typename Record>
using PySortedCatalog =
        pybind11::class_<SortedCatalogT<Record>, std::shared_ptr<SortedCatalogT<Record>>, CatalogT<Record>>;

/**
Wrap an instantiation of lsst::afw::table::SortedCatalogT<Record>.

In addition to calling this method (which also instantiates and wraps the CatalogT base class),
you must call addCatalogMethods on the class object in Python.

@tparam Record  Record type, e.g. BaseRecord or SimpleRecord.

@param[in] mod    Module object class will be added to.
@param[in] name   Name prefix of the record type, e.g. "Base" or "Simple".
@param[in] isBase Whether this instantiation is only being used as a base class (used to set the class name).

*/
template <typename Record>
PySortedCatalog<Record> declareSortedCatalog(pybind11::module &mod, std::string const &name,
                                             bool isBase = false) {
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Catalog = SortedCatalogT<Record>;
    using Table = typename Record::Table;

    auto clsBase = declareCatalog<Record>(mod, name, true);

    std::string fullName;
    if (isBase) {
        fullName = "_" + name + "SortedCatalogBase";
    } else {
        fullName = name + "Catalog";
    }

    // We need py::dynamic_attr() below to support our Python-side caching of the associated ColumnView.
    PySortedCatalog<Record> cls(mod, fullName.c_str(), py::dynamic_attr());

    /* Constructors */
    cls.def(pybind11::init<Schema const &>());
    cls.def(pybind11::init<std::shared_ptr<Table> const &>(), "table"_a = std::shared_ptr<Table>());
    cls.def(pybind11::init<Catalog const &>());

    /* Overridden and Variant Methods */
    cls.def_static("readFits", (Catalog(*)(std::string const &, int, int)) & Catalog::readFits, "filename"_a,
                   "hdu"_a = INT_MIN, "flags"_a = 0);
    cls.def_static("readFits", (Catalog(*)(fits::MemFileManager &, int, int)) & Catalog::readFits,
                   "manager"_a, "hdu"_a = INT_MIN, "flags"_a = 0);
    // readFits taking Fits objects not wrapped, because Fits objects are not wrapped.

    cls.def("subset", (Catalog (Catalog::*)(ndarray::Array<bool const, 1> const &) const) & Catalog::subset);
    cls.def("subset",
            (Catalog (Catalog::*)(std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t) const) & Catalog::subset);

    // The following three methods shadow those in the base class in C++ (unlike the base class versions,
    // they do not require a ley argument because we assume it's the ID key).  In Python, we make that appear
    // as though the key argument is available but has a default value.  If that key is not None, we delegate
    // to the base class.
    cls.def("isSorted",
            [clsBase](py::object const &self, py::object key) -> py::object {
                if (key == py::none()) {
                    key = self.attr("table").attr("getIdKey")();
                }
                return clsBase.attr("isSorted")(self, key);
            },
            "key"_a = py::none());
    cls.def("sort",
            [clsBase](py::object const &self, py::object key) -> py::object {
                if (key == py::none()) {
                    key = self.attr("table").attr("getIdKey")();
                }
                return clsBase.attr("sort")(self, key);
            },
            "key"_a = py::none());
    cls.def("find",
            [clsBase](py::object const &self, py::object const &value, py::object key) -> py::object {
                if (key == py::none()) {
                    key = self.attr("table").attr("getIdKey")();
                }
                return clsBase.attr("find")(self, value, key);
            },
            "value"_a, "key"_a = py::none());

    return cls;
};
}
}
}
}  // lsst::afw::table::python

#endif  // !AFW_TABLE_PYBIND11_CATALOG_H_INCLUDED
