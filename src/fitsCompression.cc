// -*- lsst-c++ -*-

#include "fitsio.h"
extern "C" {
#include "fitsio2.h"
}

#include "lsst/pex/exceptions.h"

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

}}} // namespace lsst::afw::fits