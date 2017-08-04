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


class ImageScale {
  public:
    int bitpix;
    double bscale;
    double bzero;

    ImageScale(int bitpix_, double bscale_, double bzero_) :
      bitpix(bitpix_), bscale(bscale_), bzero(bzero_) {}

    template <typename DiskT, typename MemT>
    image::Image<DiskT> toDisk(image::Image<MemT> const& image, bool fuzz, unsigned long seed);

    template <typename MemT, typename DiskT>
    image::Image<MemT> fromDisk(image::Image<DiskT> const& image);
};


class ImageScalingOptions {
  public:
    enum ScalingScheme {
        NONE,
        RANGE,
        STDEV_POSITIVE,
        STDEV_NEGATIVE,
        STDEV_BOTH,
        MANUAL,
    };
    ScalingScheme scheme;
    int bitpix;
    bool fuzz;                          ///< Fuzz the values when quantising floating-point values?
    unsigned long seed;
    std::vector<std::string> maskPlanes;
    float quantizeLevel;
    float quantizePad;                     ///< Number of standard deviations to pad off the edge
    double bscale, bzero;               ///< Manually specified BSCALE and BZERO (for SCALE_MANUAL)

    ImageScalingOptions(
        ScalingScheme scheme_,
        int bitpix_,
        std::vector<std::string> const& maskPlanes_,
        unsigned long seed_,
        float quantizeLevel_=4.0,
        float quantizePad_=5.0,
        bool fuzz_=true,
        double bscale_=1.0,
        double bzero_=0.0
    ) : scheme(scheme_), bitpix(bitpix_), fuzz(fuzz_), seed(seed_), maskPlanes(maskPlanes_),
        quantizeLevel(quantizeLevel_), quantizePad(quantizePad_), bscale(bscale_), bzero(bzero_) {}

    template <typename T>
    ImageScale determine(
        image::Image<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel>> mask=std::shared_ptr<image::Mask<image::MaskPixel>>()
    );

    template <typename DiskT, typename MemT>
    image::Image<DiskT> apply(
        image::Image<MemT> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel>> mask=std::shared_ptr<image::Mask<image::MaskPixel>>()
    ) {
        return determine(image, mask).template toDisk<DiskT>(image, fuzz, seed);
    }

  private:
    template <typename T>
    ImageScale determineFromRange(
        image::Image<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel>> mask
    );
    template <typename T>
    ImageScale determineFromStdev(
        image::Image<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel>> mask
    );

};

ImageScalingOptions::ScalingScheme scalingSchemeFromString(std::string const& name);
std::string scalingSchemeToString(ImageScalingOptions::ScalingScheme scheme);

}}} // namespace lsst::afw::fits

#endif // ifndef LSST_AFW_fitsCompression_h_INCLUDED