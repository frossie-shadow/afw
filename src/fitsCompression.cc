// -*- lsst-c++ -*-

#include "fitsio.h"
extern "C" {
#include "fitsio2.h"
}

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/Statistics.h"
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
math::Statistics calculateStats(
    image::Image<T> const& image,
    int flags,
    std::shared_ptr<image::Mask<image::MaskPixel>> mask=std::shared_ptr<image::Mask<image::MaskPixel>>(),
    std::vector<std::string> const& maskPlanes={}
) {
    if (mask) {
        math::StatisticsControl ctrl;
        ctrl.setAndMask(mask->getPlaneBitMask(maskPlanes));
        return math::makeStatistics(image, flags, ctrl);
    }
    return math::makeStatistics(image, *mask, flags);
}

} // anonymous namespace

template <typename T>
ImageScale ImageScalingOptions::determineFromRange(
    image::Image<T> const& image,
    std::shared_ptr<image::Mask<image::MaskPixel>> mask
) {
    double const range = std::pow(2.0, bitpix); // Range of values for target BITPIX
    auto const& stats = calculateStats(image, math::MIN | math::MAX, mask, maskPlanes);
    auto const min = stats.getValue(math::MIN), max = stats.getValue(math::MAX);
    if (!std::isfinite(min) || !std::isfinite(max)) {
        /// XXX warning or exception?
        return ImageScale(bitpix, 1.0, 0.0);
    }
    if (min == max) return ImageScale(bitpix, 1.0, min);
    double bscale = (max - min)/range;
    return ImageScale(bitpix, bscale, min + 0.5*range*bscale);
}

template <typename T>
ImageScale ImageScalingOptions::determineFromStdev(
    image::Image<T> const& image,
    std::shared_ptr<image::Mask<image::MaskPixel>> mask
) {
    auto const& stats = calculateStats(image, math::MEANCLIP | math::STDEVCLIP, mask, maskPlanes);
    auto const mean = stats.getValue(math::MEANCLIP), stdev = stats.getValue(math::STDEVCLIP);
    if (!std::isfinite(mean) || !std::isfinite(stdev)) {
        /// XXX warning or exception?
        return ImageScale(bitpix, 1.0, 0.0);
    }

    double imageVal;                    // Value on image
    long diskVal;                       // Corresponding quantized value
    switch (scheme) {
      case ImageScalingOptions::STDEV_POSITIVE:
        // Put (mean - N sigma) at the lowest possible value: predominantly positive images
        imageVal = mean - quantizePad * stdev;
        diskVal = - (1L << (bitpix - 1));
        break;
      case ImageScalingOptions::STDEV_NEGATIVE:
        // Put (mean + N sigma) at the highest possible value: predominantly negative images
        imageVal = mean + quantizePad * stdev;
        diskVal = (1L << (bitpix - 1)) - 1;
        break;
      case ImageScalingOptions::STDEV_BOTH:
        // Put mean right in the middle: images with an equal abundance of positive and negative values
        imageVal = mean;
        diskVal = 0;
        break;
      default:
        std::abort(); // Programming error: should never get here
    }

    double const bscale = stdev/quantizeLevel;
    return ImageScale(bitpix, bscale, imageVal - bscale*diskVal);
}



template <typename T>
ImageScale ImageScalingOptions::determine(
    image::Image<T> const& image,
    std::shared_ptr<image::Mask<image::MaskPixel>> mask
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
image::Image<DiskT> ImageScale::toDisk(image::Image<MemT> const& image, bool fuzz, unsigned long seed) {
    // Having the user choose the right template is simpler than boost::mpl::for_each
    if (FitsBitPix<DiskT>::CONSTANT != bitpix) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "Output type doesn't match bitpix");
    }

    if (bscale == 1.0 && bzero == 0.0 && !fuzz) {
        return image::Image<DiskT>(image, true);
    }

    image::Image<DiskT> disk(image.getBBox());
    math::Random rng(math::Random::MT19937, seed);

    // Note: BITPIX=8 treated differently, since it uses unsigned values; the rest use signed */
    double const min = bitpix == 8 ? 0 : -std::pow(2.0, bitpix - 1);
    double const max = bitpix == 8 ? 255 : (std::pow(2.0, bitpix - 1) - 1.0);

    double const scale = 1.0/bscale;
    for (int y = 0; y < image.getHeight(); ++y) {
        auto outIter = disk.row_begin(y);
        for (auto inIter = image.row_begin(y), inEnd = image.row_end(y); inIter != inEnd;
             ++inIter, ++outIter) {
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
    }

    return disk;
}


template <typename MemT, typename DiskT>
image::Image<MemT> ImageScale::fromDisk(image::Image<DiskT> const& image) {
    image::Image<MemT> memory(image.getBBox());
    memory.getArray().deep() = bscale*image.getArray() + bzero;
    return memory;
}


// Explicit instantiation
#define INSTANTIATE2(MEMTYPE, DISKTYPE) \
template image::Image<DISKTYPE> ImageScale::toDisk<DISKTYPE, MEMTYPE>( \
    image::Image<MEMTYPE> const&, bool, unsigned long); \
template image::Image<MEMTYPE> ImageScale::fromDisk<MEMTYPE, DISKTYPE>(image::Image<DISKTYPE> const&); \
template image::Image<DISKTYPE> ImageScalingOptions::apply<DISKTYPE, MEMTYPE>( \
    image::Image<MEMTYPE> const&, std::shared_ptr<image::Mask<image::MaskPixel>>);

#define INSTANTIATE1(MEMTYPE) \
template ImageScale ImageScalingOptions::determine<MEMTYPE>( \
    image::Image<MEMTYPE> const&, std::shared_ptr<image::Mask<image::MaskPixel>>); \
INSTANTIATE2(MEMTYPE, std::uint8_t); \
INSTANTIATE2(MEMTYPE, std::int16_t); \
INSTANTIATE2(MEMTYPE, std::int32_t); \
INSTANTIATE2(MEMTYPE, std::int64_t);

INSTANTIATE1(std::uint8_t);
INSTANTIATE1(std::int16_t);
INSTANTIATE1(std::uint16_t);
INSTANTIATE1(std::int32_t);
INSTANTIATE1(std::uint32_t);
INSTANTIATE1(std::uint64_t);
INSTANTIATE1(float);
INSTANTIATE1(double);


}}} // namespace lsst::afw::fits