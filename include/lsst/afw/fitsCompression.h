// -*- lsst-c++ -*-
#ifndef LSST_AFW_fitsCompression_h_INCLUDED
#define LSST_AFW_fitsCompression_h_INCLUDED

#include <string>
#include <limits>

#include "boost/cstdfloat.hpp"

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

namespace detail {

template <typename T>
struct Bitpix;

template <>
struct Bitpix<std::uint8_t> {
    static int const value = 8;
};
template <>
struct Bitpix<std::int16_t> {
    static int const value = 16;
};
template <>
struct Bitpix<std::int32_t> {
    static int const value = 32;
};
template <>
struct Bitpix<std::int64_t> {
    static int const value = 64;
};
template <>
struct Bitpix<std::uint16_t> {
    static int const value = 16;
};
template <>
struct Bitpix<std::uint32_t> {
    static int const value = 32;
};
template <>
struct Bitpix<std::uint64_t> {
    static int const value = 64;
};
template <>
struct Bitpix<float> {
    static int const value = -32;
};
template <>
struct Bitpix<double> {
    static int const value = -64;
};


template <int bitpix>
struct BitpixType;

template <>
struct BitpixType<16> {
    typedef std::int16_t type;
};
template <>
struct BitpixType<32> {
    typedef std::int32_t type;
};
template <>
struct BitpixType<64> {
    typedef std::int64_t type;
};
template <>
struct BitpixType<-32> {
    typedef boost::float32_t type;
};
template <>
struct BitpixType<-64> {
    typedef boost::float64_t type;
};

class PixelArrayBase {
  public:
    PixelArrayBase(std::size_t num) : _num(num) {}
    virtual ~PixelArrayBase() {}
    virtual void * getData() const = 0;
    std::size_t getNum() const { return _num; }

  private:
    std::size_t _num;
};

template <typename T>
class PixelArray : public PixelArrayBase {
  public:
    PixelArray(std::size_t num)
      : PixelArrayBase(num),
        _pixels(reinterpret_cast<PixelT *>(new typename std::remove_const<T>::type [num])) {}

    template <typename U>
    PixelArray(ndarray::Array<U, 1, 1> const& array)
      : PixelArrayBase(array.getNumElements()),
        _pixels(reinterpret_cast<PixelT *>(new typename std::remove_const<T>::type [getNum()])) {
        std::copy(array.begin(), array.end(),
                  reinterpret_cast<typename std::remove_const<T>::type *>(_pixels));
    }

    virtual ~PixelArray() {
        delete [] _pixels;
    }

    void * getData() const override { return reinterpret_cast<void *>(_pixels); }

  private:
    typedef std::uint8_t PixelT;  // bytes, because 'void' doesn't work, and bytes are next most generic
    PixelT * _pixels;
};

inline std::shared_ptr<PixelArrayBase> makePixelArray(int bitpix, int num) {
    switch (bitpix) {
      case 8: return std::make_shared<PixelArray<std::uint8_t>>(num);
      case 16: return std::make_shared<PixelArray<std::int16_t>>(num);
      case 32: return std::make_shared<PixelArray<std::int32_t>>(num);
      case 64: return std::make_shared<PixelArray<std::int64_t>>(num);
      case -32: return std::make_shared<PixelArray<boost::float32_t>>(num);
      case -64: return std::make_shared<PixelArray<boost::float64_t>>(num);
      default:
        std::ostringstream os;
        os << "Unrecognized bitpix: " << bitpix;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}

template <typename T>
std::shared_ptr<PixelArrayBase> makePixelArray(
    int bitpix,
    ndarray::Array<T, 1, 1> const& array
) {
    switch (bitpix) {
      case 0: return std::make_shared<PixelArray<T>>(array);
      case 8: return std::make_shared<PixelArray<std::uint8_t>>(array);
      case 16: return std::make_shared<PixelArray<std::int16_t>>(array);
      case 32: return std::make_shared<PixelArray<std::int32_t>>(array);
      case 64: return std::make_shared<PixelArray<std::int64_t>>(array);
      case -32: return std::make_shared<PixelArray<boost::float32_t>>(array);
      case -64: return std::make_shared<PixelArray<boost::float64_t>>(array);
      default:
        std::ostringstream os;
        os << "Unrecognized bitpix: " << bitpix;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}



} // namespace detail

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

    template <typename T>
    std::shared_ptr<detail::PixelArrayBase> toDisk(
        ndarray::Array<T const, 2, 2> const& image,
        bool fuzz,
        unsigned long seed=1
    ) const;

    template <typename T>
    ndarray::Array<T, 2, 2> fromDisk(ndarray::Array<T, 2, 2> const& image) const;
};


// value in memory = BZERO + BSCALE * value on disk
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

    // Disable scaling by default
    explicit ImageScalingOptions()
      : scheme(NONE), bitpix(0), fuzz(false), seed(1), quantizeLevel(4.0), quantizePad(5.0),
        bscale(std::numeric_limits<double>::quiet_NaN()), bzero(std::numeric_limits<double>::quiet_NaN()) {}

    ImageScalingOptions(
        ScalingScheme scheme_,
        int bitpix_,
        std::vector<std::string> const& maskPlanes_={},
        unsigned long seed_=1,
        float quantizeLevel_=4.0,
        float quantizePad_=5.0,
        bool fuzz_=true,
        double bscale_=1.0,
        double bzero_=0.0
    ) : scheme(scheme_), bitpix(bitpix_), fuzz(fuzz_), seed(seed_), maskPlanes(maskPlanes_),
        quantizeLevel(quantizeLevel_), quantizePad(quantizePad_), bscale(bscale_), bzero(bzero_) {}

    ImageScalingOptions(int bitpix_, double bscale_=1.0, double bzero_=0.0)
      : scheme(MANUAL), bitpix(bitpix_), bscale(bscale_), bzero(bzero_) {}

    template <typename T>
    ImageScale determine(
        image::ImageBase<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel> const> mask=
            std::shared_ptr<image::Mask<image::MaskPixel>>()
    ) const {
        auto const arrays = _toArray(image, mask);
        return determine(arrays.first, arrays.second);
    }

    template <typename T, int N>
    ImageScale determine(
        ndarray::Array<T const, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask
    ) const;

#if 0
    template <typename T>
    std::shared_ptr<detail::PixelArrayBase> apply(
        image::ImageBase<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel>const> mask=
            std::shared_ptr<image::Mask<image::MaskPixel>>()
    ) const {
        auto const arrays = _toArray(image, mask);
        return apply(arrays.first, arrays.second);
    }

    template <typename T, int N>
    std::shared_ptr<detail::PixelArrayBase> apply(
        ndarray::Array<T, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask
    ) const {
        return determine(image, mask).toDisk(image, fuzz, seed);
    }
#endif

  private:
    template <typename T>
    std::pair<ndarray::Array<T const, 2, 2>, ndarray::Array<bool, 2, 2>> _toArray(
        image::ImageBase<T> const& image,
        std::shared_ptr<image::Mask<image::MaskPixel> const> mask=
            std::shared_ptr<image::Mask<image::MaskPixel>>()
    ) const {
        if (mask && image.getDimensions() != mask->getDimensions()) {
            std::ostringstream os;
            os << "Size mismatch between image and mask: ";
            os << image.getWidth() << "x" << image.getHeight();
            os << " vs ";
            os << mask->getWidth() << "x" << mask->getHeight();
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
        }
        ndarray::Array<T const, 2, 2> imageArray = ndarray::dynamic_dimension_cast<2>(image.getArray());
        if (imageArray.empty()) imageArray = ndarray::copy(image.getArray());
        ndarray::Array<bool, 2, 2> maskArray = ndarray::allocate(imageArray.getShape());
        if (mask) {
            maskArray.deep() = (mask->getArray() & mask->getPlaneBitMask(maskPlanes));
        } else {
            maskArray.deep() = false;
        }
        return std::make_pair(imageArray, maskArray);
    }

    template <typename T, int N>
    ImageScale determineFromRange(
        ndarray::Array<T const, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask,
        bool isUnsigned=false
    ) const;
    template <typename T, int N>
    ImageScale determineFromStdev(
        ndarray::Array<T const, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask
    ) const;

};

ImageScalingOptions::ScalingScheme scalingSchemeFromString(std::string const& name);
std::string scalingSchemeToString(ImageScalingOptions::ScalingScheme scheme);

}}} // namespace lsst::afw::fits

#endif // ifndef LSST_AFW_fitsCompression_h_INCLUDED