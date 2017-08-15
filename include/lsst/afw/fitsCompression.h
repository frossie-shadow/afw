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

/// FITS BITPIX header value by C++ type
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

/// Abstract base class for an array of pixel values
///
/// For writing pixels out with the C-based cfitsio API, it's helpful to
/// have a type-agnostic holder of the pixel values (because although we can
/// template on the array type, the output type can be configured at run-time,
/// which means we can't use ndarray::Array as the carrier).
/// This is essentially a C++-ish void* array.
class PixelArrayBase {
  public:
    virtual ~PixelArrayBase() {}

    /// Return a void* array of the pixels
    virtual void * getData() const = 0;

    /// Return the number of pixels
    std::size_t getNum() const { return _num; }

  protected:
    PixelArrayBase(std::size_t num) : _num(num) {}

  private:
    std::size_t _num;  // Number of pixel values
};

/// Typed array of pixel values
template <typename T>
class PixelArray : public PixelArrayBase {
  public:

    /// Construct from an ndarray::Array
    ///
    /// Populating the array involves a copy. This is necessary when changing
    /// types. We don't currently include the magic necessary to avoid the copy
    /// when the type is unchanged.
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

/// Create a PixelArray suitable for an image with the nominated BITPIX
///
/// @param[in] bitpix  Bits per pixel (0,8,16,32,64,-32,-64) for output.
/// @param[in] array  Array with pixel values to convert.
/// @return A new PixelArray.
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
///
/// Tile compression is a feature provided by cfitsio, where contiguous parts
/// of the image ("tiles", e.g., rows, multiple rows, blocks or the entire image)
/// are compressed separately. The aim of this struct is to hold the parameters
/// used to configure the compression.
///
/// Floating-point images can be losslessly compressed (quantizeLevel=0.0) using
/// (only) the GZIP or GZIP_SHUFFLE compression schemes, but the compression
/// factor so achieved is modest (e.g., ~ 10% compression). Better compression
/// factors can be achieved if floating-point images are first quantised into
/// integer images. This can be done by cfitsio (through the quantizeLevel parameter)
/// or through use of the ImageScalingOptions.
///
/// The Compression is specified by:
/// * the compression scheme (what algorithm is used to do the actual compression)
/// * the tile size
/// * the quantization level (for quantization applied by cfitsio; for floating-point images)
///
/// Due to bugs, cfitsio may require setting the quantizeLevel to a value other than
/// zero when compressing integer data, but in this case it should have no effect.
struct ImageCompressionOptions {
    /// Compression algorithms
    ///
    /// cfitsio's compression schemes are #defines; these have a namespace.
    ///
    /// We deliberately don't support HCOMPRESS: it doesn't appear to be useful to us (e.g., lossy)
    /// and it requires extra configuration.
    enum CompressionScheme {
        NONE,              ///< No compression
        GZIP,              ///< Standard GZIP compression
        GZIP_SHUFFLE,      ///< GZIP compression with shuffle (most-significant byte first)
        RICE,              ///< RICE compression
        PLIO,              ///< PLIO compression
    };
    typedef ndarray::Array<long, 1> Tiles;

    CompressionScheme scheme;  ///< Compresion algorithm to use
    Tiles tiles;           ///< Tile size; a dimension with 0 means infinite (e.g., to specify one row: 0,1)
    float quantizeLevel;   ///< quantization level: 0.0 = none requires use of GZIP or GZIP_SHUFFLE

    /// Custom compression
    explicit ImageCompressionOptions(
        CompressionScheme scheme_,
        Tiles tiles_,
        int quantizeLevel_=0.0
    ) : scheme(scheme_), tiles(tiles_), quantizeLevel(quantizeLevel_) {}

    /// Compression by rows or entire image
    ///
    /// @param[in] scheme_  Compression algorithm to use
    /// @param[in] rows  Number of rows per tile (0 = entire image)
    /// @param[in] quantizeLevel_  cfitsio quantization level
    explicit ImageCompressionOptions(
        CompressionScheme scheme_,
        int rows=1,
        int quantizeLevel_=0.0
    );

    /// Default compression for a particular style of image
    template <typename T>
    explicit ImageCompressionOptions(image::Image<T> const& image) :
        ImageCompressionOptions(image.getBBox().getArea() > 0 ? GZIP_SHUFFLE : NONE) {}
    template <typename T>
    explicit ImageCompressionOptions(image::Mask<T> const& mask) :
        ImageCompressionOptions(mask.getBBox().getArea() > 0 ? GZIP_SHUFFLE : NONE) {}

    // Disable compression for int64: cfitsio won't compress them
    explicit ImageCompressionOptions(image::Image<std::int64_t> const& image) :
        ImageCompressionOptions(NONE) {}
    explicit ImageCompressionOptions(image::Mask<std::int64_t> const& mask) :
        ImageCompressionOptions(NONE) {}
    explicit ImageCompressionOptions(image::Image<std::uint64_t> const& image) :
        ImageCompressionOptions(NONE) {}
    explicit ImageCompressionOptions(image::Mask<std::uint64_t> const& mask) :
        ImageCompressionOptions(NONE) {}

};

/// Interpret compression scheme expressed in string
ImageCompressionOptions::CompressionScheme compressionSchemeFromString(std::string const& name);

/// Provide string version of compression scheme
std::string compressionSchemeToString(ImageCompressionOptions::CompressionScheme scheme);

/// Convert compression scheme from cfitsio to ImageCompressionOptions::CompressionScheme
ImageCompressionOptions::CompressionScheme compressionSchemeFromCfitsio(int cfitsio);

/// Convert ImageCompressionOptions::CompressionScheme to cfitsio
int compressionSchemeToCfitsio(ImageCompressionOptions::CompressionScheme scheme);


/// Scale to apply to image
///
/// Images are scaled to the type implied by the provided BITPIX
/// using the provided scale and zero-point:
///
///    value in memory = BZERO + BSCALE * value in FITS
///
/// In addition to scaling, a random field of values distributed [0,1) may be
/// added before quantisation ("fuzz"); this preserves the expectation value of
/// the floating-point image, while increasing the variance by 1/12.
struct ImageScale {
    int bitpix;  ///< Bits per pixel; negative means floating-point: 8,16,32,64,-32,-64
    double bscale;  ///< Scale to apply when reading from FITS
    double bzero;  ///< Zero-point to apply when reading from FITS

    /// Constructor
    ImageScale(int bitpix_, double bscale_, double bzero_) :
      bitpix(bitpix_), bscale(bscale_), bzero(bzero_) {}

    /// Convert to an array of pixel values to write to FITS
    ///
    /// @param[in] image  Image to scale
    /// @param[in] fuzz  Add random values before quantising?
    /// @param[in] seed  Seed for random number generator
    /// @return Array of pixel values, appropriately scaled.
    template <typename T>
    std::shared_ptr<detail::PixelArrayBase> toFits(
        ndarray::Array<T const, 2, 2> const& image,
        bool fuzz,
        unsigned long seed=1
    ) const;

    /// Convert to an array
    ///
    /// Use of this method is generally not necessary, since cfitsio automatically
    /// applies the scaling on read.
    template <typename T>
    ndarray::Array<T, 2, 2> fromFits(ndarray::Array<T, 2, 2> const& image) const;
};


/// Options for scaling image pixels
///
/// Scaling (quantisation) of floating-point images is important in order to achieve
/// respectable compression factors. Unfortunately, scaling a floating-point image
/// means losing some information, in two ways: values are quantised so values in
/// between the quanta will be rounded up or down; and values outside the range
/// of supported values in the integer image will be clipped. This implementation
/// is based on the successful implementation used by Pan-STARRS.
///
/// cfitsio provides a facility for scaling image pixels on write (see
/// ImageCompressionOptions.quantizeLevel), but that facility is limited
/// (it provides a single scaling options, and doesn't use our masks when
/// collecting statistics). This scaling facility provides multiple scaling
/// options, and could be extended to support more in the future (e.g., logarithmic
/// or asinh scaling).
///
/// Scaling is specified by:
/// * the scaling scheme: the algorithm used to determining the scales (and therefore
///   the dynamic range).
/// * the bits per pixel: 8,16,32,64,-32,-64 (negative means floating-point but you
///   probably don't want to use those here); larger values give more dynamic range,
///   produce larger images (which effect may be negated by compression).
/// * fuzz: whether to add a random field of values distributed [0,1) before
///   quantisation; this preserves the expectation value of the floating-point image,
///   while increasing the variance by 1/12.
/// * seed: seed for the random number generator used by the fuzz.
/// * maskPlanes: a list of mask planes to be ignored when doing statistics.
/// * quantizeLevel: for the STDEV_* algorithms, specifies the ratio of the quantum
///   size to the image standard deviation.
/// * quantizePad: for the STDEV_POSITIVE and STDEV_NEGATIVE algorithms, specifies
///   how many standard deviations to allow on the short side.
/// * bscale, bzero: for the MANUAL algorithm, specifies the BSCALE and BZERO to use.
///
/// Scaling schemes are:
/// * NONE: no scaling or quantisation at all. The image goes out the way it came in.
/// * RANGE: scale based on the dynamic range in the image. This preserves dynamic
///   range at the cost of sensitivity to low-level fluctuations.
/// * STDEV_POSITIVE: the scale is set to sample the standard deviation of the image
///   (bscale = stdev/quantizeLevel) and the dynamic range extends principally
///   positive (allowing only quantizePad standard deviations negative).
/// * STDEV_NEGATIVE: similar to STDEV_POSITIVE, but the dynamic range extends
///   principally negative.
/// * STDEV_BOTH: the scale is set similar to STDEV_POSITIVE, but the dynamic range
///   is shared between the positive and negative sides.
/// * MANUAL: the scale is set manually. We do what we're told, no more, no less.
///
/// Perhaps this one class could/should have been polymorphic, with different
/// subclasses for different algorithms? But I went with a C-like approach keying
/// off the enum, probably because I was influenced by Pan-STARRS' C code.
class ImageScalingOptions {
  public:
    enum ScalingScheme {
        NONE,  ///< No scaling
        RANGE,  ///< Scale to preserve dynamic range
        STDEV_POSITIVE,  ///< Scale based on the standard deviation. dynamic range positive
        STDEV_NEGATIVE,  ///< Scale based on the standard deviation, dynamic range negative
        STDEV_BOTH,  ///< Scale based on the standard deviation, dynamic range positive+negative
        MANUAL,  ///< Scale set manually
    };
    ScalingScheme scheme;  ///< Scaling algorithm to use
    int bitpix;  ///< Bits per pixel (0, 8,16,32,64,-32,-64)
    bool fuzz;  ///< Fuzz the values when quantising floating-point values?
    unsigned long seed;  ///< Seed for random number generator when fuzzing
    std::vector<std::string> maskPlanes;  ///< Mask planes to ignore when doing statistics
    float quantizeLevel;  ///< Divisor of the standard deviation for STDEV_* scaling
    float quantizePad;  ///< Number of stdev to allow on the low side (for STDEV_POSITIVE/NEGATIVE)
    double bscale, bzero;  ///< Manually specified BSCALE and BZERO (for MANUAL scaling)

    /// Default Ctor
    ///
    /// Scaling is disabled by default.
    explicit ImageScalingOptions()
      : scheme(NONE), bitpix(0), fuzz(false), seed(1), quantizeLevel(4.0), quantizePad(5.0),
        bscale(std::numeric_limits<double>::quiet_NaN()), bzero(std::numeric_limits<double>::quiet_NaN()) {}

    /// General purpose Ctor
    ///
    /// @param[in] scheme_  Scaling algorithm to use
    /// @param[in] bitpix_  Bits per pixel (8,16,32,64,-32,-64)
    /// @param[in] maskPlanes_  Mask planes to ignore when doing statistics
    /// @param[in] seed_  Seed for random number generator when fuzzing
    /// @param[in] quantizeLevel_  Divisor of the standard deviation for STDEV_* scaling
    /// @param[in] quantizePad_  Number of stdev to allow on the low side (for STDEV_POSITIVE/NEGATIVE)
    /// @param[in] fuzz_  Fuzz the values when quantising floating-point values?
    /// @param[in] bscale_  Manually specified BSCALE (for MANUAL scaling)
    /// @param[in] bzero_  Manually specified BZERO (for MANUAL scaling)
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

    /// Manual scaling Ctor
    ///
    /// @param[in] bitpix_  Bits per pixel (8,16,32,64,-32,-64)
    /// @param[in] bscale_  Manually specified BSCALE
    /// @param[in] bzero_  Manually specified BZERO
    ImageScalingOptions(int bitpix_, double bscale_=1.0, double bzero_=0.0)
      : scheme(MANUAL), bitpix(bitpix_), fuzz(false), seed(1), maskPlanes({}),
        quantizeLevel(4.0), quantizePad(5.0), bscale(bscale_), bzero(bzero_) {}

    ///{
    /// Determine the scaling for a particular image
    ///
    /// @param[in] image  Image for which to determine scaling
    /// @param[in] mask  Mask for image (used to measuring statistics)
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
    ///}

  private:
    /// Convert image,mask to arrays
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

    /// Determine scale using the RANGE algorithm
    ///
    /// @param[in] image  Image for which to determine scaling
    /// @param[in] mask  Mask for image (used to measuring statistics)
    /// @param[in] isUnsigned  Is the output using an unsigned integer?
    template <typename T, int N>
    ImageScale determineFromRange(
        ndarray::Array<T const, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask,
        bool isUnsigned=false
    ) const;

    /// Determine scale using the STDEV algorithm
    ///
    /// @param[in] image  Image for which to determine scaling
    /// @param[in] mask  Mask for image (used to measuring statistics)
    template <typename T, int N>
    ImageScale determineFromStdev(
        ndarray::Array<T const, N, N> const& image,
        ndarray::Array<bool, N, N> const& mask
    ) const;

};

/// Interpret scaling scheme expressed in string
ImageScalingOptions::ScalingScheme scalingSchemeFromString(std::string const& name);

/// Provide string version of compression scheme
std::string scalingSchemeToString(ImageScalingOptions::ScalingScheme scheme);

}}} // namespace lsst::afw::fits

#endif // ifndef LSST_AFW_fitsCompression_h_INCLUDED
