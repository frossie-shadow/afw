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
//#include <pybind11/operators.h>
//#include <pybind11/stl.h>

#include "numpy/arrayobject.h"

#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"
#include "lsst/daf/base.h"

#include "lsst/afw/fits.h"

namespace py = pybind11;

using namespace lsst::afw::fits;
using namespace pybind11::literals;


void defineImageCompression(py::module & mod) {
    py::class_<ImageCompressionOptions> cls(mod, "ImageCompressionOptions");

    py::enum_<ImageCompressionOptions::CompressionScheme> scheme(cls, "CompressionScheme");
    scheme.value("NONE", ImageCompressionOptions::CompressionScheme::NONE);
    scheme.value("GZIP", ImageCompressionOptions::CompressionScheme::GZIP);
    scheme.value("GZIP", ImageCompressionOptions::CompressionScheme::GZIP_SHUFFLE);
    scheme.value("RICE", ImageCompressionOptions::CompressionScheme::RICE);
    scheme.value("PLIO", ImageCompressionOptions::CompressionScheme::PLIO);
    scheme.export_values();

    cls.def(py::init<ImageCompressionOptions::CompressionScheme, ImageCompressionOptions::Tiles, int>(),
            "scheme"_a, "tiles"_a, "quantizeLevel"_a=0.0);

    cls.def(py::init<lsst::afw::image::Image<unsigned char> const&>());
    cls.def(py::init<lsst::afw::image::Image<unsigned short> const&>());
    cls.def(py::init<lsst::afw::image::Image<short> const&>());
    cls.def(py::init<lsst::afw::image::Image<int> const&>());
    cls.def(py::init<lsst::afw::image::Image<unsigned int> const&>());
    cls.def(py::init<lsst::afw::image::Image<float> const&>());
    cls.def(py::init<lsst::afw::image::Image<double> const&>());
    cls.def(py::init<lsst::afw::image::Image<std::uint64_t> const&>());

    cls.def(py::init<lsst::afw::image::Mask<unsigned char> const&>());
    cls.def(py::init<lsst::afw::image::Mask<unsigned short> const&>());
    cls.def(py::init<lsst::afw::image::Mask<short> const&>());
    cls.def(py::init<lsst::afw::image::Mask<std::int32_t> const&>());

    cls.def_readonly("scheme", &ImageCompressionOptions::scheme);
    cls.def_readonly("tiles", &ImageCompressionOptions::tiles);
    cls.def_readonly("quantizeLevel", &ImageCompressionOptions::quantizeLevel);
}


// Wrapping for lsst::afw::fits::Fits
//
// Not every feature is wrapped, only those that we guess might be useful.
// In particular, the header keyword read/write and table read/write are not wrapped.
void defineFits(py::module & mod) {
    py::class_<Fits> cls(mod, "Fits");

    cls.def(py::init<std::string const&, std::string const&, int>(), "filename"_a, "mode"_a,
            "behavior"_a=Fits::AUTO_CLOSE | Fits::AUTO_CHECK);
    cls.def(py::init<MemFileManager &, std::string const&, int>(), "manager"_a, "mode"_a,
            "behavior"_a=Fits::AUTO_CLOSE | Fits::AUTO_CHECK);

    cls.def("closeFile", &Fits::closeFile);
    cls.def("getFileName", &Fits::getFileName);
    cls.def("getHdu", &Fits::getHdu);
    cls.def("setHdu", &Fits::setHdu, "hdu"_a, "relative"_a=false);
    cls.def("countHdus", &Fits::countHdus);

    cls.def("writeMetadata", &Fits::writeMetadata);
    cls.def("readMetadata", [](Fits & self, bool strip=false) { return readMetadata(self, strip); },
            "strip"_a=false);
    cls.def("createEmpty", &Fits::createEmpty);

    cls.def("setImageCompression", &Fits::setImageCompression);
    cls.def("getImageCompression", &Fits::getImageCompression);
    cls.def("checkCompressedImagePhu", &Fits::checkCompressedImagePhu);

    cls.def_readonly("status", &Fits::status);
}


PYBIND11_PLUGIN(_fits) {
    py::module mod("_fits", "Python wrapper for afw _fits library");

    // Need to import numpy for ndarray and eigen conversions
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    py::class_<MemFileManager> clsMemFileManager(mod, "MemFileManager");

    lsst::pex::exceptions::python::declareException<FitsError, lsst::pex::exceptions::IoError>(
            mod, "FitsError", "IoError");
    //    lsst::pex::exceptions::python::declareException<FitsTypeError, FitsError>(mod, "FitsTypeError",
    //    "FitsError");

    clsMemFileManager.def(py::init<>());
    clsMemFileManager.def(py::init<size_t>());

    /* TODO: We should really revisit persistence and pickling as this is quite ugly.
     * But it is what Swig did (sort of, it used the cdata.i extension), so I reckon this
     * is cleaner because it does not expose casting to the Python side. */
    clsMemFileManager.def("getLength", &MemFileManager::getLength);
    clsMemFileManager.def("getData", [](MemFileManager &m) {
        return py::bytes(static_cast<char *>(m.getData()), m.getLength());
    });
    clsMemFileManager.def("setData", [](MemFileManager &m, py::bytes const &d, size_t size) {
        memcpy(m.getData(), PyBytes_AsString(d.ptr()), size);
    });
    clsMemFileManager.def("readMetadata",
                          [](MemFileManager & self, int hdu=INT_MIN, bool strip=false) {
                              return readMetadata(self, hdu, strip);
                          }, "hdu"_a=INT_MIN, "strip"_a=false);

    defineImageCompression(mod);
    defineFits(mod);

    mod.def("readMetadata",
            [](std::string const& filename, int hdu=INT_MIN, bool strip=false) {
                return readMetadata(filename, hdu, strip);
            }, "fileName"_a, "hdu"_a=INT_MIN, "strip"_a=false);

    return mod.ptr();
}