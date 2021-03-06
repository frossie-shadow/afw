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
#include "pybind11/stl.h"

#include "lsst/afw/table/Field.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/afw/table/detail/SchemaImpl.h"
#include "lsst/afw/table/SchemaMapper.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace afw {
namespace table {
namespace {

using PySchemaMapper = py::class_<SchemaMapper, std::shared_ptr<SchemaMapper>>;

template <typename T>
void declareSchemaMapperOverloads(PySchemaMapper &cls, std::string const &suffix) {
    cls.def("getMapping", (Key<T> (SchemaMapper::*)(Key<T> const &) const) & SchemaMapper::getMapping);
};

PYBIND11_PLUGIN(schemaMapper) {
    py::module mod("schemaMapper");
    py::module::import("lsst.afw.table.schema");

    PySchemaMapper cls(mod, "SchemaMapper");

    cls.def(py::init<>());
    cls.def(py::init<Schema const &, Schema const &>());
    cls.def(py::init<Schema const &, bool>(), "input"_a, "shareAliasMap"_a = false);

    cls.def("getInputSchema", &SchemaMapper::getInputSchema);
    cls.def("getOutputSchema", &SchemaMapper::getOutputSchema);
    cls.def("editOutputSchema", &SchemaMapper::editOutputSchema, py::return_value_policy::reference_internal);
    cls.def("addMinimalSchema", &SchemaMapper::addMinimalSchema, "minimal"_a, "doMap"_a = true);
    cls.def_static("removeMinimalSchema", &SchemaMapper::removeMinimalSchema);
    cls.def_static("join", &SchemaMapper::join, "inputs"_a, "prefixes"_a = std::vector<std::string>());
    declareSchemaMapperOverloads<std::uint16_t>(cls, "U");
    declareSchemaMapperOverloads<std::int32_t>(cls, "I");
    declareSchemaMapperOverloads<std::int64_t>(cls, "L");
    declareSchemaMapperOverloads<float>(cls, "F");
    declareSchemaMapperOverloads<double>(cls, "D");
    declareSchemaMapperOverloads<std::string>(cls, "String");
    declareSchemaMapperOverloads<lsst::afw::geom::Angle>(cls, "Angle");
    declareSchemaMapperOverloads<lsst::afw::table::Flag>(cls, "Flag");
    declareSchemaMapperOverloads<lsst::afw::table::Array<std::uint16_t>>(cls, "ArrayU");
    declareSchemaMapperOverloads<lsst::afw::table::Array<int>>(cls, "ArrayI");
    declareSchemaMapperOverloads<lsst::afw::table::Array<float>>(cls, "ArrayF");
    declareSchemaMapperOverloads<lsst::afw::table::Array<double>>(cls, "ArrayD");

    return mod.ptr();
}
}
}
}
}  // namespace lsst::afw::table::<anonymous>
