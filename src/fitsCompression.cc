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
    if (name == "NONE") return ImageCompressionOptions::COMPRESS_NONE;
    if (name == "GZIP") return ImageCompressionOptions::COMPRESS_GZIP;
    if (name == "GZIP_SORTED") return ImageCompressionOptions::COMPRESS_GZIP_SORTED;
    if (name == "RICE") return ImageCompressionOptions::COMPRESS_RICE;
    if (name == "HCOMPRESS") throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                                               "HCOMPRESS is unsupported");
    if (name == "PLIO") return ImageCompressionOptions::COMPRESS_PLIO;
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Unrecognised compression scheme: " + name);
}

std::string compressionSchemeToString(ImageCompressionOptions::CompressionScheme scheme)
{
    switch (scheme) {
      case ImageCompressionOptions::COMPRESS_NONE: return "NONE";
      case ImageCompressionOptions::COMPRESS_GZIP: return "GZIP";
      case ImageCompressionOptions::COMPRESS_GZIP_SORTED: return "GZIP_SORTED";
      case ImageCompressionOptions::COMPRESS_RICE: return "RICE";
      case ImageCompressionOptions::COMPRESS_PLIO: return "PLIO";
      default:
        std::ostringstream os;
        os << "Unrecognized compression scheme: " << scheme;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}

ImageCompressionOptions::CompressionScheme compressionSchemeFromCfitsio(int cfitsio)
{
    switch (cfitsio) {
      case 0: return ImageCompressionOptions::COMPRESS_NONE;
      case RICE_1: return ImageCompressionOptions::COMPRESS_RICE;
      case GZIP_1: return ImageCompressionOptions::COMPRESS_GZIP;
      case GZIP_2: return ImageCompressionOptions::COMPRESS_GZIP_SORTED;
      case PLIO_1: return ImageCompressionOptions::COMPRESS_PLIO;
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
      case ImageCompressionOptions::COMPRESS_NONE: return 0;
      case ImageCompressionOptions::COMPRESS_GZIP: return GZIP_1;
      case ImageCompressionOptions::COMPRESS_GZIP_SORTED: return GZIP_2;
      case ImageCompressionOptions::COMPRESS_RICE: return RICE_1;
      case ImageCompressionOptions::COMPRESS_PLIO: return PLIO_1;
      default:
        std::ostringstream os;
        os << "Unrecognized compression scheme: " << scheme;
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, os.str());
    }
}

}}} // namespace lsst::afw::fits