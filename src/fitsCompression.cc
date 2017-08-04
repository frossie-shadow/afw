// -*- lsst-c++ -*-

#include "fitsio.h"
extern "C" {
#include "fitsio2.h"
}

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/Random.h"

#include "lsst/afw/fitsCompression.h"

namespace lsst {
namespace afw {
namespace fits {

ImageCompressionOptions::CompressionScheme compressionSchemeFromString(std::string const& name)
{
    if (name == "NONE") return ImageCompressionOptions::NONE;
    if (name == "GZIP") return ImageCompressionOptions::GZIP;
    if (name == "GZIP_SHUFFLE") return ImageCompressionOptions::GZIP_SHUFFLE;
    if (name == "RICE") return ImageCompressionOptions::RICE;
    if (name == "HCOMPRESS") throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                                               "HCOMPRESS is unsupported");
    if (name == "PLIO") return ImageCompressionOptions::PLIO;
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Unrecognised compression scheme: " + name);
}

std::string compressionSchemeToString(ImageCompressionOptions::CompressionScheme scheme)
{
    switch (scheme) {
      case ImageCompressionOptions::NONE: return "NONE";
      case ImageCompressionOptions::GZIP: return "GZIP";
      case ImageCompressionOptions::GZIP_SHUFFLE: return "GZIP_SHUFFLE";
      case ImageCompressionOptions::RICE: return "RICE";
      case ImageCompressionOptions::PLIO: return "PLIO";
      default:
        std::ostringstream os;
        os << "Unrecognized compression scheme: " << scheme;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}

ImageCompressionOptions::CompressionScheme compressionSchemeFromCfitsio(int cfitsio)
{
    switch (cfitsio) {
      case 0: return ImageCompressionOptions::NONE;
      case RICE_1: return ImageCompressionOptions::RICE;
      case GZIP_1: return ImageCompressionOptions::GZIP;
      case GZIP_2: return ImageCompressionOptions::GZIP_SHUFFLE;
      case PLIO_1: return ImageCompressionOptions::PLIO;
      case HCOMPRESS_1: throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                                          "Unsupported compression scheme: HCOMPRESS_1");
      default:
        std::ostringstream os;
        os << "Unrecognized cfitsio compression: " << cfitsio;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}

int compressionSchemeToCfitsio(ImageCompressionOptions::CompressionScheme scheme)
{
    switch (scheme) {
      case ImageCompressionOptions::NONE: return 0;
      case ImageCompressionOptions::GZIP: return GZIP_1;
      case ImageCompressionOptions::GZIP_SHUFFLE: return GZIP_2;
      case ImageCompressionOptions::RICE: return RICE_1;
      case ImageCompressionOptions::PLIO: return PLIO_1;
      default:
        std::ostringstream os;
        os << "Unrecognized compression scheme: " << scheme;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}


ImageCompressionOptions::ImageCompressionOptions(
    ImageCompressionOptions::CompressionScheme scheme_,
    int quantizeLevel_
) : scheme(scheme_), quantizeLevel(quantizeLevel_) {
    tiles = ndarray::allocate(MAX_COMPRESS_DIM);
    tiles[0] = 1;
    for (int ii = 1; ii < MAX_COMPRESS_DIM; ++ii) tiles[ii] = 1;
}


ImageScalingOptions::ScalingScheme scalingSchemeFromString(std::string const& name)
{
    if (name == "NONE") return ImageScalingOptions::NONE;
    if (name == "RANGE") return ImageScalingOptions::RANGE;
    if (name == "STDEV_POSITIVE") return ImageScalingOptions::STDEV_POSITIVE;
    if (name == "STDEV_NEGATIVE") return ImageScalingOptions::STDEV_NEGATIVE;
    if (name == "STDEV_BOTH") return ImageScalingOptions::STDEV_BOTH;
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Unrecognized scaling scheme: " + name);
}


std::string scalingSchemeToString(ImageScalingOptions::ScalingScheme scheme)
{
    switch (scheme) {
      case ImageScalingOptions::NONE: return "NONE";
      case ImageScalingOptions::RANGE: return "RANGE";
      case ImageScalingOptions::STDEV_POSITIVE: return "STDEV_POSITIVE";
      case ImageScalingOptions::STDEV_NEGATIVE: return "STDEV_NEGATIVE";
      case ImageScalingOptions::STDEV_BOTH: return "STDEV_BOTH";
      default:
        std::ostringstream os;
        os << "Unrecognized scaling scheme: " << scheme;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}


namespace {

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
#if 0
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
#endif


template <typename T, int N>
std::pair<T, T> calculateMedianStdev(
    ndarray::Array<T const, N, N> const& image,
    ndarray::Array<bool const, N, N> const& mask
) {
    std::size_t num = 0;
    auto const& flatMask = ndarray::flatten<1>(mask);
    for (auto mm = flatMask.begin(); mm != flatMask.end(); ++mm) {
        if (!*mm) ++num;
    }
    ndarray::Array<T, 1, 1> array = ndarray::allocate(num);
    auto const& flatImage = ndarray::flatten<1>(image);
    auto mm = ndarray::flatten<1>(mask).begin();
    auto aa = array.begin();
    for (auto ii = flatImage.begin(); ii != flatImage.end(); ++ii, ++mm, ++aa) {
        if (*mm) continue;
        *aa = *ii;
    }

    // Quartiles; from https://stackoverflow.com/a/11965377/834250
    auto const q1 = num/4;
    auto const q2 = num/2;
    auto const q3 = q1 + q2;
    std::nth_element(array.begin(), array.begin() + q1, array.end());
    std::nth_element(array.begin() + q1 + 1, array.begin() + q2, array.end());
    std::nth_element(array.begin() + q2 + 1, array.begin() + q3, array.end());

    T const median = num % 2 ? array[num/2] : 0.5*(array[num/2] + array[num/2 - 1]);
    // No, we're not doing any interpolation for the lower and upper quartiles.
    // We're estimating the noise, so it doesn't need to be super precise.
    T const lq = array[q1];
    T const uq = array[q3];
    return std::make_pair(median, 0.741*(uq - lq));
}

} // anonymous namespace

template <typename T, int N>
ImageScale ImageScalingOptions::determineFromRange(
    ndarray::Array<T const, N, N> const& image,
    ndarray::Array<bool const, N, N> const& mask
) {
    double const range = std::pow(2.0, bitpix); // Range of values for target BITPIX
    T min = std::numeric_limits<T>::max(), max = std::numeric_limits<T>::min();
    auto mm = ndarray::flatten<1>(mask).begin();
    auto const& flatImage = ndarray::flatten<1>(image);
    for (auto ii = flatImage.begin(); ii != flatImage.end(); ++ii, ++mm) {
        if (*mm) continue;
        if (!std::isfinite(*ii)) continue;
        if (*ii > max) max = *ii;
        if (*ii < min) min = *ii;
    }
    if (min == max) return ImageScale(bitpix, 1.0, min);
    double bscale = (max - min)/range;
    return ImageScale(bitpix, bscale, min + 0.5*range*bscale);
}

template <typename T, int N>
ImageScale ImageScalingOptions::determineFromStdev(
    ndarray::Array<T const, N, N> const& image,
    ndarray::Array<bool const, N, N> const& mask
) {
    auto stats = calculateMedianStdev(image, mask);
    auto const median = stats.first, stdev = stats.second;

    double imageVal;                    // Value on image
    long diskVal;                       // Corresponding quantized value
    switch (scheme) {
      case ImageScalingOptions::STDEV_POSITIVE:
        // Put (mean - N sigma) at the lowest possible value: predominantly positive images
        imageVal = median - quantizePad * stdev;
        diskVal = - (1L << (bitpix - 1));
        break;
      case ImageScalingOptions::STDEV_NEGATIVE:
        // Put (mean + N sigma) at the highest possible value: predominantly negative images
        imageVal = median + quantizePad * stdev;
        diskVal = (1L << (bitpix - 1)) - 1;
        break;
      case ImageScalingOptions::STDEV_BOTH:
        // Put mean right in the middle: images with an equal abundance of positive and negative values
        imageVal = median;
        diskVal = 0;
        break;
      default:
        std::abort(); // Programming error: should never get here
    }

    double const bscale = stdev/quantizeLevel;
    return ImageScale(bitpix, bscale, imageVal - bscale*diskVal);
}


template <typename T, int N>
ImageScale ImageScalingOptions::determineImpl(
    ndarray::Array<T const, N, N> const& image,
    ndarray::Array<bool const, N, N> const& mask
) {
    switch (scheme) {
      case NONE: return ImageScale(bitpix, 1.0, 0.0);
      case RANGE: return determineFromRange(image, mask);
      case ImageScalingOptions::STDEV_POSITIVE:
      case ImageScalingOptions::STDEV_NEGATIVE:
      case ImageScalingOptions::STDEV_BOTH:
        return determineFromStdev(image, mask);
      default:
        std::abort();  // should never get here
    }
}

template <typename DiskT, typename MemT>
ndarray::Array<DiskT, 2, 2> ImageScale::toDisk(
    ndarray::Array<MemT const, 2, 2> const& image,
    bool fuzz,
    unsigned long seed
) {
    // Having the user choose the right template is simpler than boost::mpl::for_each
    if (Bitpix<DiskT>::value != bitpix) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Output type doesn't match bitpix");
    }
    ndarray::Array<DiskT, 2, 2> disk(image.getShape());

    if (bscale == 1.0 && bzero == 0.0 && !fuzz) {
        disk.deep() = image;
        return disk;
    }

    math::Random rng(math::Random::MT19937, seed);

    // Note: BITPIX=8 treated differently, since it uses unsigned values; the rest use signed */
    double const min = bitpix == 8 ? 0 : -std::pow(2.0, bitpix - 1);
    double const max = bitpix == 8 ? 255 : (std::pow(2.0, bitpix - 1) - 1.0);

    double const scale = 1.0/bscale;
    auto outIter = ndarray::flatten<1>(disk).begin();
    auto const& flatImage = ndarray::flatten<1>(image);
    for (auto inIter = flatImage.begin(); inIter != flatImage.end(); ++inIter, ++outIter) {
        double value = (*inIter - bzero)*scale;
        if (!std::isfinite(value)) {
            // This choice of "max" for non-finite and overflow pixels is mainly cosmetic --- it has to be
            // something, and "min" would produce holes in the cores of bright stars.
            *outIter = max;
            continue;
        }
        if (fuzz && std::floor(value) != value) {
            // Add random factor [0.0,1.0): adds a variance of 1/12,
            // but preserves the expectation value given the floor()
            value += rng.uniform();
        }
        // Check for underflow and overflow; set either to max
        *outIter = (value < min || value > max ? max : std::floor(value));
    }

    return disk;
}


template <typename MemT, typename DiskT>
ndarray::Array<MemT, 2, 2> ImageScale::fromDisk(ndarray::Array<DiskT const, 2, 2> const& image) {
    ndarray::Array<MemT, 2, 2> memory = ndarray::allocate(image.getShape());
    memory.deep() = bscale*image + bzero;
    return memory;
}


// Explicit instantiation
#define INSTANTIATE2(MEMTYPE, DISKTYPE) \
template ndarray::Array<DISKTYPE, 2, 2> ImageScale::toDisk<DISKTYPE, MEMTYPE>( \
    ndarray::Array<MEMTYPE const, 2, 2> const&, bool, unsigned long); \
template ndarray::Array<MEMTYPE, 2, 2> ImageScale::fromDisk<MEMTYPE, DISKTYPE>( \
    ndarray::Array<DISKTYPE const, 2, 2> const&); \
template ndarray::Array<DISKTYPE, 2, 2> ImageScalingOptions::apply<DISKTYPE, MEMTYPE>( \
    image::Image<MEMTYPE> const&, std::shared_ptr<image::Mask<image::MaskPixel>>);

#define INSTANTIATE1(MEMTYPE) \
template ImageScale ImageScalingOptions::determine<MEMTYPE>( \
    image::Image<MEMTYPE> const&, std::shared_ptr<image::Mask<image::MaskPixel>>); \
template ImageScale ImageScalingOptions::determineImpl<MEMTYPE, 2>( \
    ndarray::Array<MEMTYPE const, 2, 2> const& image, ndarray::Array<bool const, 2, 2> const& mask); \
INSTANTIATE2(MEMTYPE, std::int16_t); \
INSTANTIATE2(MEMTYPE, std::int32_t); \
INSTANTIATE2(MEMTYPE, std::int64_t);

INSTANTIATE1(std::uint16_t);
INSTANTIATE1(std::int32_t);
INSTANTIATE1(std::uint64_t);
INSTANTIATE1(float);
INSTANTIATE1(double);

}}} // namespace lsst::afw::fits