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
#include <pybind11/stl.h>

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/exceptions/python/Exception.h"
#include "lsst/daf/base.h"
#include "lsst/afw/image/Image.h"

#include "lsst/afw/fits.h"

namespace py = pybind11;

using namespace lsst::afw::fits;
using namespace lsst::afw::fits::detail;
using namespace pybind11::literals;


void definePixelArray(py::module & mod) {
    // Very minimal because this carries internal details and shouldn't be needed from python
    py::class_<detail::PixelArrayBase, std::shared_ptr<detail::PixelArrayBase>> cls(mod, "PixelArrayBase");
    cls.def("getNum", &PixelArrayBase::getNum);
};


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


template <typename T>
void defineImageScalingOptionsTemplates(py::class_<ImageScalingOptions> & cls) {
    cls.def("determine", &ImageScalingOptions::determine<T>);
    cls.def("apply", &ImageScalingOptions::apply<T>);
}

void defineImageScalingOptions(py::module & mod) {
    py::class_<ImageScalingOptions> cls(mod, "ImageScalingOptions");
    py::enum_<ImageScalingOptions::ScalingScheme> scheme(cls, "CompressionScheme");
    scheme.value("NONE", ImageScalingOptions::ScalingScheme::NONE);
    scheme.value("RANGE", ImageScalingOptions::ScalingScheme::RANGE);
    scheme.value("STDEV_POSITIVE", ImageScalingOptions::ScalingScheme::STDEV_POSITIVE);
    scheme.value("STDEV_NEGATIVE", ImageScalingOptions::ScalingScheme::STDEV_NEGATIVE);
    scheme.value("STDEV_BOTH", ImageScalingOptions::ScalingScheme::STDEV_BOTH);
    scheme.value("MANUAL", ImageScalingOptions::ScalingScheme::MANUAL);
    scheme.export_values();

    cls.def(py::init<>());
    cls.def(py::init<ImageScalingOptions::ScalingScheme, int, std::vector<std::string> const&, unsigned long,
                     float, float, bool, double, double>(),
            "scheme"_a, "bitpix"_a, "maskPlanes"_a, "seed"_a=1, "quantizeLevel"_a=4.0, "quantizePad"_a=5.0,
            "fuzz"_a=true, "bscale"_a=1.0, "bzero"_a=0.0);

    cls.def_readonly("scheme", &ImageScalingOptions::scheme);
    cls.def_readonly("bitpix", &ImageScalingOptions::bitpix);
    cls.def_readonly("maskPlanes", &ImageScalingOptions::maskPlanes);
    cls.def_readonly("seed", &ImageScalingOptions::seed);
    cls.def_readonly("quantizeLevel", &ImageScalingOptions::quantizeLevel);
    cls.def_readonly("quantizePad", &ImageScalingOptions::quantizePad);
    cls.def_readonly("fuzz", &ImageScalingOptions::fuzz);
    cls.def_readonly("bscale", &ImageScalingOptions::bscale);
    cls.def_readonly("bzero", &ImageScalingOptions::bzero);

    defineImageScalingOptionsTemplates<float>(cls);
    defineImageScalingOptionsTemplates<double>(cls);
}

template <typename T>
void defineImageScaleTemplates(py::class_<ImageScale> & cls, std::string const& suffix) {
    cls.def("toDisk", &ImageScale::toDisk<T>, "image"_a, "fuzz"_a=true, "seed"_a=1);
    cls.def("fromDisk", &ImageScale::fromDisk<T>);
}

void defineImageScale(py::module & mod) {
    py::class_<ImageScale> cls(mod, "ImageScale");
    cls.def(py::init<int, double, double>(), "bitpix"_a, "bscale"_a, "bzero"_a);
    cls.def_readonly("bitpix", &ImageScale::bitpix);
    cls.def_readonly("bscale", &ImageScale::bscale);
    cls.def_readonly("bzero", &ImageScale::bzero);

    defineImageScaleTemplates<float>(cls, "F");
    defineImageScaleTemplates<double>(cls, "D");
}

void defineImageWriteOptions(py::module & mod) {
    py::class_<ImageWriteOptions> cls(mod, "ImageWriteOptions");

    cls.def(py::init<lsst::afw::image::Image<std::uint16_t>>());
    cls.def(py::init<lsst::afw::image::Image<std::int32_t>>());
    cls.def(py::init<lsst::afw::image::Image<std::uint64_t>>());
    cls.def(py::init<lsst::afw::image::Image<float>>());
    cls.def(py::init<lsst::afw::image::Image<double>>());

    cls.def(py::init<lsst::afw::image::Mask<lsst::afw::image::MaskPixel>>());

    cls.def(py::init<ImageCompressionOptions const&, ImageScalingOptions const&>(),
            "compression"_a, "scaling"_a=ImageScalingOptions());
    cls.def(py::init<ImageScalingOptions const&>());
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

    cls.def("gotoFirstHdu", [](Fits & self) { self.setHdu(INT_MIN); });

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
    defineImageScalingOptions(mod);
    defineImageScale(mod);
    defineImageWriteOptions(mod);
    defineFits(mod);

    mod.def("readMetadata",
            [](std::string const& filename, int hdu=INT_MIN, bool strip=false) {
                return readMetadata(filename, hdu, strip);
            }, "fileName"_a, "hdu"_a=INT_MIN, "strip"_a=false);
    mod.def("setAllowImageCompression", &setAllowImageCompression, "allow"_a);
    mod.def("getAllowImageCompression", &getAllowImageCompression);

    return mod.ptr();
}