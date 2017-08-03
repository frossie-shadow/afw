// -*- lsst-c++ -*-
#ifndef LSST_AFW_fitsCompression_h_INCLUDED
#define LSST_AFW_fitsCompression_h_INCLUDED

#include <string>
#include <limits>

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
        NONE,              ///< No compression
        GZIP,              ///< Standard GZIP compression
        GZIP_SHUFFLE,      ///< GZIP compression with shuffle (most-significant byte first)
        RICE,              ///< RICE compression
        PLIO,              ///< PLIO compression
    };
    typedef ndarray::Array<long, 1> Tiles;

    CompressionScheme scheme;
    Tiles tiles;           ///< tile size
    float quantizeLevel;   ///< quantization level: 0.0 = none

    /// Custom compression
    explicit ImageCompressionOptions(
        CompressionScheme scheme_,
        Tiles tiles_,
        int quantizeLevel_=0.0
    ) : scheme(scheme_), tiles(tiles_), quantizeLevel(quantizeLevel_) {}

    /// Compression by rows
    explicit ImageCompressionOptions(
        CompressionScheme scheme_,
        int quantizeLevel_=0.0
    );

    /// Full image compression
    ImageCompressionOptions(
        CompressionScheme scheme_,
        geom::Extent2I const& dims,
        int quantizeLevel_=0.0
    ) : scheme(scheme_), quantizeLevel(quantizeLevel_) {
        tiles = ndarray::allocate(2);
        tiles.asEigen() = dims.asEigen().template cast<long>();
    }

    /// Default compression for a particular style of image
    template <typename T>
    explicit ImageCompressionOptions(image::Image<T> const& image) :
        ImageCompressionOptions(image.getBBox().getArea() > 0 ? GZIP : NONE, 0.0) {}
    template <typename T>
    explicit ImageCompressionOptions(image::Mask<T> const& mask) :
        ImageCompressionOptions(mask.getBBox().getArea() > 0 ? GZIP : NONE, 0.0) {}

    /// Disable compression for int64: cfitsio won't compress them
    explicit ImageCompressionOptions(image::Image<std::int64_t> const& image) :
        ImageCompressionOptions(NONE, 0.0) {}
    explicit ImageCompressionOptions(image::Mask<std::int64_t> const& mask) :
        ImageCompressionOptions(NONE, 0.0) {}
    explicit ImageCompressionOptions(image::Image<std::uint64_t> const& image) :
        ImageCompressionOptions(NONE, 0.0) {}
    explicit ImageCompressionOptions(image::Mask<std::uint64_t> const& mask) :
        ImageCompressionOptions(NONE, 0.0) {}

};

ImageCompressionOptions::CompressionScheme compressionSchemeFromString(std::string const& name);
std::string compressionSchemeToString(ImageCompressionOptions::CompressionScheme scheme);
ImageCompressionOptions::CompressionScheme compressionSchemeFromCfitsio(int cfitsio);
int compressionSchemeToCfitsio(ImageCompressionOptions::CompressionScheme scheme);

}}} // namespace lsst::afw::fits

#endif // ifndef LSST_AFW_fitsCompression_h_INCLUDED