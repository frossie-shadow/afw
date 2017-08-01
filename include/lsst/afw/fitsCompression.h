// -*- lsst-c++ -*-
#ifndef LSST_AFW_fitsCompression_h_INCLUDED
#define LSST_AFW_fitsCompression_h_INCLUDED

#include <string>

#include "lsst/pex/exceptions.h"
#include "lsst/daf/base.h"
#include "ndarray.h"
#include "ndarray/eigen.h"

#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Mask.h"

namespace lsst {
namespace afw {
namespace fits {

// Forward declarations
class Fits;

/// Options for tile compression of image pixels
struct ImageCompressionOptions {
    // cfitsio compression schemes are #define-ed; these have a namespace.
    //
    // We deliberately don't support HCOMPRESS: it doesn't appear to be useful to us (e.g., lossy)
    // and it requires extra configuration.
    enum CompressionScheme {
        COMPRESS_NONE,              ///< No compression
        COMPRESS_GZIP,              ///< Standard GZIP compression
        COMPRESS_GZIP_SORTED,       ///< GZIP compression with sorting
        COMPRESS_RICE,              ///< RICE compression
        COMPRESS_PLIO,              ///< PLIO compression
    };
    typedef ndarray::Array<long, 1> Tiles;

    CompressionScheme scheme;
    Tiles tiles;           ///< tile size
    float quantizeLevel;   ///< quantization level: 0.0 = none

    explicit ImageCompressionOptions(
        CompressionScheme scheme_,
        Tiles tiles_,
        int quantizeLevel_=0.0
    ) : scheme(scheme_), tiles(tiles_), quantizeLevel(quantizeLevel_) {}

    /// Default compression for a particular style of image
    template <typename T>
    explicit ImageCompressionOptions(image::Image<T> const& image)
      : scheme(image.getBBox().getArea() > 0 ? COMPRESS_GZIP_SORTED : COMPRESS_NONE),
        quantizeLevel(0.0) {
        tiles = ndarray::allocate(2);
        tiles.asEigen() = image.getDimensions().asEigen().template cast<long>();
    }
    template <typename T>
    explicit ImageCompressionOptions(image::Mask<T> const& mask)
      : scheme(mask.getBBox().getArea() > 0 ? COMPRESS_GZIP : COMPRESS_NONE), quantizeLevel(0.0) {
        tiles = ndarray::allocate(2);
        tiles.asEigen() = mask.getDimensions().asEigen().template cast<long>();
    }
};

ImageCompressionOptions::CompressionScheme compressionSchemeFromString(std::string const& name);
std::string compressionSchemeToString(ImageCompressionOptions::CompressionScheme scheme);
ImageCompressionOptions::CompressionScheme compressionSchemeFromCfitsio(int cfitsio);
int compressionSchemeToCfitsio(ImageCompressionOptions::CompressionScheme scheme);

}}} // namespace lsst::afw::fits

#endif // ifndef LSST_AFW_fitsCompression_h_INCLUDED