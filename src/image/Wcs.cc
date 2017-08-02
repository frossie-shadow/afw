// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <limits>

#include "boost/format.hpp"

#include "wcslib/wcs.h"
#include "wcslib/wcsfix.h"
#include "wcslib/wcshdr.h"

#include "lsst/daf/base.h"
#include "lsst/daf/base/Citizen.h"
#include "lsst/afw/formatters/Utils.h"
#include "lsst/afw/formatters/WcsFormatter.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/SpherePoint.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/aggregates.h"

namespace except = lsst::pex::exceptions;

using namespace std;

typedef lsst::daf::base::PropertySet PropertySet;
typedef lsst::daf::base::PropertyList PropertyList;
typedef lsst::afw::geom::Point2D GeomPoint;
typedef std::shared_ptr<lsst::afw::geom::SpherePoint> CoordPtr;

// The amount of space allocated to strings in wcslib
const int STRLEN = 72;

namespace lsst {
namespace afw {
namespace image {

// Set internal params for wcslib
void Wcs::_setWcslibParams() {
    _wcsfixCtrl =        // ctrl for wcsfix
            2;           // Translate "H" to "h"
    _wcshdrCtrl =        // ctrl for wcspih
            2;           // Report each rejected keyrecord and the reason why it was rejected
    _relax =             // relax parameter for wcspih;
            WCSHDR_all;  // Accept all extensions recognized by the parser
}

const int lsstToFitsPixels = +1;
const int fitsToLsstPixels = -1;

//
// Constructors
//

Wcs::Wcs()
        : daf::base::Citizen(typeid(this)),
          _wcsInfo(NULL),
          _nWcsInfo(0),
          _relax(0),
          _wcsfixCtrl(0),
          _wcshdrCtrl(0),
          _nReject(0),
          _coordSystem(coord::UNKNOWN)  // set by _initWcs
{
    _setWcslibParams();
    _initWcs();
}

Wcs::Wcs(std::shared_ptr<daf::base::PropertySet const> const& fitsMetadata)
        : daf::base::Citizen(typeid(this)),
          _wcsInfo(NULL),
          _nWcsInfo(0),
          _relax(0),
          _wcsfixCtrl(0),
          _wcshdrCtrl(0),
          _nReject(0),
          _coordSystem(coord::UNKNOWN)  // set by _initWcs
{
    _setWcslibParams();

    initWcsLibFromFits(fitsMetadata);
    _initWcs();
}

void Wcs::_initWcs() {
    // first four characters of CTYPE1 (name of first axis)
    std::string ctype1 = std::string(_wcsInfo->ctype[0]).substr(0, 4);

    if (_wcsInfo) {
        if (ctype1[0] == 'G') {
            _coordSystem = coord::GALACTIC;
            _skyAxesSwapped = (ctype1[2] == 'A');  // GLAT instead of GLON
        } else if (ctype1[0] == 'E') {
            _coordSystem = coord::ECLIPTIC;
            _skyAxesSwapped = (ctype1[2] == 'A');  // ELAT instead of ELON
        } else {
            _coordSystem = coord::makeCoordEnum(_wcsInfo->radesys);
            _skyAxesSwapped = (ctype1[0] == 'D');  // DEC instead of RA
        }

        // tell WCSlib that values have been updated
        _wcsInfo->flag = 0;
        // and then tell it to do its internal magic.
        int status = wcsset(_wcsInfo);
        if (status != 0) {
            throw LSST_EXCEPT(except::RuntimeError,
                              (boost::format("Failed to setup wcs structure with wcsset. Status %d: %s") %
                               status % wcs_errmsg[status])
                                      .str());
        }
    }
}

Wcs::Wcs(GeomPoint const& crval, GeomPoint const& crpix, Eigen::Matrix2d const& CD, std::string const& ctype1,
         std::string const& ctype2, double equinox, std::string const& raDecSys, std::string const& cunits1,
         std::string const& cunits2)
        : daf::base::Citizen(typeid(this)),
          _wcsInfo(NULL),
          _nWcsInfo(0),
          _relax(0),
          _wcsfixCtrl(0),
          _wcshdrCtrl(0),
          _nReject(0),
          _coordSystem(coord::UNKNOWN)  // set by _initWcs
{
    _setWcslibParams();
    initWcsLib(crval, crpix, CD, ctype1, ctype2, equinox, raDecSys, cunits1, cunits2);
    _initWcs();
}

void Wcs::initWcsLibFromFits(std::shared_ptr<daf::base::PropertySet const> const& header) {
    /*
     * Access control for the input header
     *
     * We want to hack up the input, and in order to do so we need to do a deep copy on it.
     * We only want to do that copy once, and would like to avoid doing it altogether.
     */
    class HeaderAccess {
    public:
        // Return a readable version of the metadata
        std::shared_ptr<daf::base::PropertySet const> const& toRead() { return _constHeader; }
        // Return a writable version of the metadata
        std::shared_ptr<daf::base::PropertySet> const& toWrite() {
            if (!_hackHeader) {
                _hackHeader = _constHeader->deepCopy();
                _constHeader = _hackHeader;
            }
            return _hackHeader;
        }

        // Ctor
        HeaderAccess(std::shared_ptr<daf::base::PropertySet const> const& header)
                : _constHeader(header), _hackHeader() {}

    private:
        std::shared_ptr<daf::base::PropertySet const> _constHeader;
        std::shared_ptr<daf::base::PropertySet> _hackHeader;
    };

    HeaderAccess access(header);

    // Some headers (e.g. SDSS ones from FNAL) have EQUINOX as a string.  Fix this,
    // as wcslib 4.4.4 refuses to handle it.  Furthermore, some headers (e.g. Skymapper)
    // use a string but prepend 'J':  "J2000.0" -- if we don't handle this, we deduce
    // an equinox of 0.0
    {
        std::string const& key = "EQUINOX";
        if (access.toRead()->exists(key) && access.toRead()->typeOf(key) == typeid(std::string)) {
            auto equinoxStr = access.toRead()->getAsString(key).c_str();
            double equinox = ::atof(equinoxStr[0] == 'J' ? equinoxStr + 1 : equinoxStr);
            access.toWrite()->set(key, equinox);
        }
    }

    //Check header isn't empty
    int nCards = lsst::afw::formatters::countFitsHeaderCards(*(access.toRead()));
    if (nCards <= 0) {
        string msg = "Could not parse FITS WCS: no header cards found";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    // While the standard does not insist on CRVAL and CRPIX being present, it
    // is almost certain their absence indicates a problem.
    // Check for CRPIX
    if (!access.toRead()->exists("CRPIX1") && !access.toRead()->exists("CRPIX1a")) {
        string msg = "Neither CRPIX1 not CRPIX1a found";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    if (!access.toRead()->exists("CRPIX2") && !access.toRead()->exists("CRPIX2a")) {
        string msg = "Neither CRPIX2 not CRPIX2a found";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    // And the same for CRVAL
    if (!access.toRead()->exists("CRVAL1") && !access.toRead()->exists("CRVAL1a")) {
        string msg = "Neither CRVAL1 not CRVAL1a found";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    if (!access.toRead()->exists("CRVAL2") && !access.toRead()->exists("CRVAL2a")) {
        string msg = "Neither CRVAL2 not CRVAL2a found";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }
    /*
     * According to Greisen and Calabretta (A&A 395, 1061–1075 (2002)) it's illegal to mix PCi_j and CDi_j
     * headers; unfortunately Subaru puts both in its headers.  It actually uses PC001002 instead of PC1_2
     * (dating to a proposed FITS standard from 1996) and at least sometimes fails to include CDELT[12],
     * so the CD and PC matrices are inconsistent
     *
     * If we detect any part of a CD matrix, delete all PC matrices
     */
    if (access.toRead()->exists("CD1_1") || access.toRead()->exists("CD1_2") ||
        access.toRead()->exists("CD2_1") || access.toRead()->exists("CD2_2")) {
        for (int i = 1; i <= 2; ++i) {
            for (int j = 1; j <= 2; ++j) {
                std::string key = (boost::format("PC%i_%i") % j % i).str();
                if (access.toRead()->exists(key)) {
                    double const val = access.toRead()->getAsDouble(key);
                    access.toWrite()->remove(key);
                    access.toWrite()->add("X_" + key, val);
                }

                key = (boost::format("PC%03d%03d") % j % i).str();
                if (access.toRead()->exists(key)) {
                    double const val = access.toRead()->getAsDouble(key);
                    access.toWrite()->remove(key);
                    access.toWrite()->add("X_" + key, val);
                }
            }
        }
    }

    //Pass the header into wcslib's formatter to extract & setup the Wcs. First need
    //to convert to a C style string, so the compile doesn't complain about constness
    std::string metadataStr = formatters::formatFitsProperties(*(access.toRead()));
    // We own the data, and wcslib is slack about constness, so no qualms with casting away const
    char* hdrString = const_cast<char*>(metadataStr.c_str());
    // printf("wcspih string:\n%s\n", hdrString);

    nCards = formatters::countFitsHeaderCards(*(access.toRead()));  // we may have dropped some
    int pihStatus = wcspih(hdrString, nCards, _relax, _wcshdrCtrl, &_nReject, &_nWcsInfo, &_wcsInfo);

    if (pihStatus != 0) {
        throw LSST_EXCEPT(except::RuntimeError,
                          (boost::format("Could not parse FITS WCS: wcspih status = %d (%s)") % pihStatus %
                           wcs_errmsg[pihStatus])
                                  .str());
    }

    // Run wcsfix on _wcsInfo to try and fix any problems it knows about.
    const int* naxes = NULL;  // should be {NAXIS1, NAXIS2, ...} to check cylindrical projections
    int stats[NWCSFIX];       // status returns from wcsfix
    int fixStatus = wcsfix(_wcsfixCtrl, naxes, _wcsInfo, stats);
    if (fixStatus != 0) {
        std::stringstream errStream;
        errStream << "Could not parse FITS WCS: wcsfix failed " << std::endl;
        for (int ii = 0; ii < NWCSFIX; ++ii) {
            if (stats[ii] >= 0) {
                errStream << "\t" << ii << ": " << stats[ii] << " " << wcsfix_errmsg[stats[ii]] << std::endl;
            } else {
                errStream << "\t" << ii << ": " << stats[ii] << std::endl;
            }
        }
    }

    // The Wcs standard requires a default value for RADESYS if the keyword
    // doesn't exist in header, but wcslib doesn't set it. So we do so here. This code
    // conforms to Calabretta & Greisen 2002 \S 3.1
    if (!(access.toRead()->exists("RADESYS") || access.toRead()->exists("RADESYSa"))) {
        // If RADECSYS exists, use that (counter to Calabretta & Greisen 2002 \S 3.1, but commonly used).
        // If equinox exist and < 1984, use FK4. If >= 1984, use FK5
        if (access.toRead()->exists("RADECSYS")) {
            strncpy(_wcsInfo->radesys, access.toRead()->getAsString("RADECSYS").c_str(), STRLEN);
        } else if (access.toRead()->exists("EQUINOX") || access.toRead()->exists("EQUINOXa")) {
            std::string const EQUINOX = access.toRead()->exists("EQUINOX") ? "EQUINOX" : "EQUINOXa";
            double const equinox = access.toRead()->getAsDouble(EQUINOX);
            if (equinox < 1984) {
                strncpy(_wcsInfo->radesys, "FK4", STRLEN);
            } else {
                strncpy(_wcsInfo->radesys, "FK5", STRLEN);
            }
        } else {
            // If Equinox doesn't exist, default to ICRS
            strncpy(_wcsInfo->radesys, "ICRS", STRLEN);
        }
    }
    // strip trailing whitespace
    {
        for (int i = strlen(_wcsInfo->radesys) - 1; i >= 0; i--) {
            if (isspace(_wcsInfo->radesys[i])) {
                _wcsInfo->radesys[i] = '\0';
            }
        }
    }
    //
    // If there are no CDi_j cards in the header, set CDi_j from PCi_j
    // CDi_j == CDELTi*PCi_j
    //
    if ((_wcsInfo->altlin & 2) == 0) {  // no CDi_j cards were found in the header
        double const* cdelt = _wcsInfo->cdelt;
        double const* pc = _wcsInfo->pc;
        double* cd = _wcsInfo->cd;

        cd[0] = cdelt[0] * pc[0];  // 1_1
        cd[1] = cdelt[0] * pc[1];  // 1_2
        cd[2] = cdelt[1] * pc[2];  // 2_1
        cd[3] = cdelt[1] * pc[3];  // 2_2
    }
}

void Wcs::initWcsLib(GeomPoint const& crval, GeomPoint const& crpix, Eigen::Matrix2d const& CD,
                     std::string const& ctype1, std::string const& ctype2, double equinox,
                     std::string const& raDecSys, std::string const& cunits1, std::string const& cunits2) {
    // Check CD is a valid size
    if ((CD.rows() != 2) || (CD.cols() != 2)) {
        string msg = "CD is not a 2x2 matrix";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    // Check that cunits are legitimate values
    bool isValid = (cunits1 == "deg");
    isValid |= (cunits1 == "arcmin");
    isValid |= (cunits1 == "arcsec");
    isValid |= (cunits1 == "mas");

    if (!isValid) {
        string msg = "CUNITS1 must be one of {deg|arcmin|arcsec|mas}";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    isValid = (cunits2 == "deg");
    isValid |= (cunits2 == "arcmin");
    isValid |= (cunits2 == "arcsec");
    isValid |= (cunits2 == "mas");

    if (!isValid) {
        string msg = "CUNITS2 must be one of {deg|arcmin|arcsec|mas}";
        throw LSST_EXCEPT(except::InvalidParameterError, msg);
    }

    // Initialise the wcs struct
    _wcsInfo = static_cast<struct wcsprm*>(malloc(sizeof(struct wcsprm)));
    if (_wcsInfo == NULL) {
        throw LSST_EXCEPT(except::MemoryError, "Cannot allocate WCS info");
    }

    _wcsInfo->flag = -1;
    int status = wcsini(true, 2, _wcsInfo);  // 2 indicates a naxis==2, a two dimensional image
    if (status != 0) {
        throw LSST_EXCEPT(except::MemoryError,
                          (boost::format("Failed to allocate memory with wcsini. Status %d: %s") % status %
                           wcs_errmsg[status])
                                  .str());
    }

    // Set crval, crpix and CD. Internally to the class, we use fits units for consistency with
    // wcslib.
    _wcsInfo->crval[0] = crval.getX();
    _wcsInfo->crval[1] = crval.getY();
    _wcsInfo->crpix[0] = crpix.getX() + lsstToFitsPixels;
    _wcsInfo->crpix[1] = crpix.getY() + lsstToFitsPixels;

    // Set the CD matrix
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            _wcsInfo->cd[(2 * i) + j] = CD(i, j);
        }
    }

    // Specify that we have a CD matrix, but no PC or CROTA
    _wcsInfo->altlin = 2;
    _wcsInfo->flag = 0;  // values have been updated

    // This is a work around for what I think is a bug in wcslib. ->types is neither
    // initialised or set to NULL by default, so if I try to delete a Wcs object,
    // wcslib then attempts to free non-existent space, and the code can crash.
    _wcsInfo->types = NULL;

    // Set the coordinate system
    strncpy(_wcsInfo->ctype[0], ctype1.c_str(), STRLEN);
    strncpy(_wcsInfo->ctype[1], ctype2.c_str(), STRLEN);
    strncpy(_wcsInfo->radesys, raDecSys.c_str(), STRLEN);
    _wcsInfo->equinox = equinox;

    // Set the units
    strncpy(_wcsInfo->cunit[0], cunits1.c_str(), STRLEN);
    strncpy(_wcsInfo->cunit[1], cunits2.c_str(), STRLEN);

    _nWcsInfo = 1;  // Specify that we have only one coordinate representation

    // Tell wcslib that we are need to set up internal values
    status = wcsset(_wcsInfo);
    if (status != 0) {
        throw LSST_EXCEPT(except::RuntimeError,
                          (boost::format("Failed to setup wcs structure with wcsset. Status %d: %s") %
                           status % wcs_errmsg[status])
                                  .str());
    }
}

Wcs::Wcs(Wcs const& rhs)
        : daf::base::Citizen(typeid(this)),
          _wcsInfo(NULL),
          _nWcsInfo(rhs._nWcsInfo),
          _relax(rhs._relax),
          _wcsfixCtrl(rhs._wcsfixCtrl),
          _wcshdrCtrl(rhs._wcshdrCtrl),
          _nReject(rhs._nReject),
          _coordSystem(coord::UNKNOWN)  // set by _initWcs
{
    if (rhs._nWcsInfo > 0) {
        _wcsInfo = static_cast<struct wcsprm*>(calloc(rhs._nWcsInfo, sizeof(struct wcsprm)));
        if (_wcsInfo == NULL) {
            throw LSST_EXCEPT(pex::exceptions::MemoryError, "Cannot allocate WCS info");
        }

        int alloc = 1;  // Unconditionally allocate memory when calling
        for (int i = 0; i != rhs._nWcsInfo; ++i) {
            _wcsInfo[i].flag = -1;
            int status = wcscopy(alloc, &rhs._wcsInfo[i], &_wcsInfo[i]);
            if (status != 0) {
                wcsvfree(&i, &_wcsInfo);
                throw LSST_EXCEPT(pex::exceptions::MemoryError,
                                  (boost::format("Could not copy WCS: wcscopy status = %d : %s") % status %
                                   wcs_errmsg[status])
                                          .str());
            }
        }
    }
    _initWcs();
}

bool Wcs::operator==(Wcs const& other) const {
    if (&other == this) return true;
    // We do a bidirectional test with a virtual member function in case one of us is a derived
    // class with members we don't know about here.
    // This is not the most efficient possible implementation, but I think it's the easiest one
    // with which to ensure correctness, and I think that's more important in this case.
    return this->_isSubset(other) && other._isSubset(*this);
}

// convenience functions and a macro for implementing isSubset
namespace {

inline bool compareArrays(double const* a, double const* b, int n) {
    for (int i = 0; i < n; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

template <typename T>
inline bool compareStringArrays(T a, T b, int n) {
    for (int i = 0; i < n; ++i)
        if (strcmp(a[i], b[i]) != 0) return false;
    return true;
}

#define CHECK_NULLS(a, b)                 \
    do {                                  \
        if ((a) == NULL) {                \
            if ((b) == NULL) return true; \
            return false;                 \
        }                                 \
        if ((b) == NULL) return false;    \
    } while (false)

}  // anonymous

bool Wcs::_isSubset(Wcs const& rhs) const {
    CHECK_NULLS(_wcsInfo, rhs._wcsInfo);
    CHECK_NULLS(_wcsInfo->crval, rhs._wcsInfo->crval);
    CHECK_NULLS(_wcsInfo->cd, rhs._wcsInfo->cd);
    CHECK_NULLS(_wcsInfo->crpix, rhs._wcsInfo->crpix);
    CHECK_NULLS(_wcsInfo->cunit, rhs._wcsInfo->cunit);
    CHECK_NULLS(_wcsInfo->ctype, rhs._wcsInfo->ctype);
    return _nWcsInfo == rhs._nWcsInfo && _coordSystem == rhs._coordSystem &&
           _wcsInfo->naxis == rhs._wcsInfo->naxis && _wcsInfo->equinox == rhs._wcsInfo->equinox &&
           _wcsInfo->altlin == rhs._wcsInfo->altlin &&
           compareArrays(_wcsInfo->crval, rhs._wcsInfo->crval, 2) &&
           compareArrays(_wcsInfo->crpix, rhs._wcsInfo->crpix, 2) &&
           compareArrays(_wcsInfo->cd, rhs._wcsInfo->cd, 4) &&
           compareStringArrays(_wcsInfo->cunit, rhs._wcsInfo->cunit, 2) &&
           compareStringArrays(_wcsInfo->ctype, rhs._wcsInfo->ctype, 2) &&
           skyToPixel(_wcsInfo->crval[0] * geom::degrees, _wcsInfo->crval[1] * geom::degrees) ==
                   rhs.skyToPixel(_wcsInfo->crval[0] * geom::degrees, _wcsInfo->crval[1] * geom::degrees) &&
           *pixelToSky(_wcsInfo->crpix[0], _wcsInfo->crpix[1]) ==
                   *rhs.pixelToSky(_wcsInfo->crpix[0], _wcsInfo->crpix[1]);
}

Wcs::~Wcs() {
    if (_wcsInfo != NULL) {
        wcsvfree(&_nWcsInfo, &_wcsInfo);
    }
}

std::shared_ptr<Wcs> Wcs::clone(void) const { return std::shared_ptr<Wcs>(new Wcs(*this)); }

//
// Accessors
//

CoordPtr Wcs::getSkyOrigin() const {
    assert(_wcsInfo);
    return std::make_shared<geom::SpherePoint>(_wcsInfo->crval[0] * geom::degrees,
                                               _wcsInfo->crval[1] * geom::degrees);
}

GeomPoint Wcs::getPixelOrigin() const {
    assert(_wcsInfo);
    // Convert from fits units back to lsst units
    double p1 = _wcsInfo->crpix[0] + fitsToLsstPixels;
    double p2 = _wcsInfo->crpix[1] + fitsToLsstPixels;
    return geom::Point2D(p1, p2);
}

Eigen::Matrix2d Wcs::getCDMatrix() const {
    assert(_wcsInfo);
    int const naxis = _wcsInfo->naxis;

    // If naxis != 2, I'm not sure if any of what follows is correct
    assert(naxis == 2);

    Eigen::Matrix2d C;

    for (int i = 0; i < naxis; ++i) {
        for (int j = 0; j < naxis; ++j) {
            C(i, j) = _wcsInfo->cd[(i * naxis) + j];
        }
    }

    return C;
}

void Wcs::flipImage(int flipLR, int flipTB, geom::Extent2I dimensions) const {
    assert(_wcsInfo);

    int const naxis = _wcsInfo->naxis;

    // If naxis != 2, I'm not sure if any of what follows is correct
    assert(naxis == 2);
    // Origin is at (1,1).  Adjust to avoid off by one.
    if (flipLR) {
        _wcsInfo->cd[0] = -_wcsInfo->cd[0];
        _wcsInfo->cd[2] = -_wcsInfo->cd[2];
        _wcsInfo->crpix[0] = -_wcsInfo->crpix[0] + dimensions.getX() + 1;
    }
    if (flipTB) {
        _wcsInfo->cd[1] = -_wcsInfo->cd[1];
        _wcsInfo->cd[3] = -_wcsInfo->cd[3];
        _wcsInfo->crpix[1] = -_wcsInfo->crpix[1] + dimensions.getY() + 1;
    }

    // tells libwcs to invalidate cached data, since transformation has been modified
    _wcsInfo->flag = 0;
}

void Wcs::rotateImageBy90(int nQuarter, geom::Extent2I dimensions) const {
    assert(_wcsInfo);

    while (nQuarter < 0) {
        nQuarter += 4;
    }

    int const naxis = _wcsInfo->naxis;

    // If naxis != 2, I'm not sure if any of what follows is correct
    assert(naxis == 2);
    double a = _wcsInfo->cd[0];
    double b = _wcsInfo->cd[1];
    double c = _wcsInfo->cd[2];
    double d = _wcsInfo->cd[3];
    double crpx = _wcsInfo->crpix[0];
    double crpy = _wcsInfo->crpix[1];
    // Origin is at (1,1).  Adjust to avoid off by one.
    // E.g. CRPIX one pixel off the UR of the grid should go to
    //(0,0) for nQuarter=2
    switch (nQuarter % 4) {
        case 0:
            break;
        case 1:
            _wcsInfo->cd[0] = -b;
            _wcsInfo->cd[1] = a;
            _wcsInfo->cd[2] = -d;
            _wcsInfo->cd[3] = c;
            _wcsInfo->crpix[0] = -crpy + dimensions.getY() + 1;
            _wcsInfo->crpix[1] = crpx;
            break;
        case 2:
            _wcsInfo->cd[0] = -a;
            _wcsInfo->cd[1] = -b;
            _wcsInfo->cd[2] = -c;
            _wcsInfo->cd[3] = -d;
            _wcsInfo->crpix[0] = -crpx + dimensions.getX() + 1;
            _wcsInfo->crpix[1] = -crpy + dimensions.getY() + 1;
            break;
        case 3:
            _wcsInfo->cd[0] = b;
            _wcsInfo->cd[1] = -a;
            _wcsInfo->cd[2] = d;
            _wcsInfo->cd[3] = -c;
            _wcsInfo->crpix[0] = crpy;
            _wcsInfo->crpix[1] = -crpx + dimensions.getX() + 1;
            break;
    }

    // tells libwcs to invalidate cached data, since transformation has been modified
    _wcsInfo->flag = 0;
}

std::shared_ptr<PropertyList> Wcs::getFitsMetadata() const {
    return formatters::WcsFormatter::generatePropertySet(*this);
}

bool Wcs::isFlipped() const {
    // We calculate the determinant of the CD (i.e the rotation and scaling)
    // matrix. If this determinant is positive, then the image can be rotated
    // to a position where increasing the right ascension and declination
    // increases the horizontal and vertical pixel position. In this case the
    // image is flipped.
    assert(_wcsInfo);
    double det = (_wcsInfo->cd[0] * _wcsInfo->cd[3]) - (_wcsInfo->cd[1] * _wcsInfo->cd[2]);

    if (det == 0) {
        throw(LSST_EXCEPT(except::RuntimeError, "Wcs CD matrix is singular"));
    }

    return (det > 0);
}

static double square(double x) { return x * x; }

double Wcs::pixArea(GeomPoint pix00) const {
    //
    // Figure out the (0, 0), (0, 1), and (1, 0) ra/dec coordinates of the corners of a square drawn in pixel
    // It'd be better to centre the square at sky00, but that would involve another conversion between sky and
    // pixel coordinates so I didn't bother
    //
    const double side = 1;  // length of the square's sides in pixels

    // Work in 3-space to avoid RA wrapping and pole issues.
    geom::Point3D v0 = pixelToSky(pix00)->getVector();

    // Step by "side" in x and y pixel directions...
    GeomPoint px(pix00);
    GeomPoint py(pix00);
    px.shift(geom::Extent2D(side, 0));
    py.shift(geom::Extent2D(0, side));
    // Push the points through the WCS, and find difference in 3-space.
    geom::Extent3D dx = pixelToSky(px)->getVector() - v0;
    geom::Extent3D dy = pixelToSky(py)->getVector() - v0;

    // Compute |cross product| = area of parallelogram with sides dx,dy
    // FIXME -- this is slightly incorrect -- it's making the small-angle
    // approximation, taking the distance *through* the unit sphere
    // rather than over its surface.
    // This is in units of ~radians^2
    double area = sqrt(square(dx[1] * dy[2] - dx[2] * dy[1]) + square(dx[2] * dy[0] - dx[0] * dy[2]) +
                       square(dx[0] * dy[1] - dx[1] * dy[0]));

    return area / square(side) * square(180. / geom::PI);
}

geom::Angle Wcs::pixelScale() const { return sqrt(pixArea(getPixelOrigin())) * geom::degrees; }

GeomPoint Wcs::skyToPixelImpl(geom::Angle sky1, geom::Angle sky2) const {
    assert(_wcsInfo);

    double skyTmp[2];
    double imgcrd[2];
    double phi, theta;
    double pixTmp[2];
    /*
     printf("_skyCoordsReversed: %c\n", (_skyCoordsReversed ? 'T' : 'F'));
     printf("wcsinfo.lat: %i,  lng: %i\n", _wcsInfo->lat, _wcsInfo->lng);
     */
    // WCSLib is smart enough to notice and handle crazy SDSS CTYPE1 = DEC--TAN,
    // by recording the indices of the long and lat coordinates.
    skyTmp[_wcsInfo->lng] = sky1.asDegrees();
    skyTmp[_wcsInfo->lat] = sky2.asDegrees();

    int stat[1];
    int status = 0;
    status = wcss2p(_wcsInfo, 1, 2, skyTmp, &phi, &theta, imgcrd, pixTmp, stat);
    if (status == 9) {
        throw LSST_EXCEPT(except::DomainError,
                          (boost::format("sky coordinates %s, %s degrees is not valid for this WCS") %
                           sky1.asDegrees() % sky2.asDegrees())
                                  .str());
    }
    if (status > 0) {
        throw LSST_EXCEPT(except::RuntimeError,
                          (boost::format("Error: wcslib returned a status code of %d at sky %s, %s deg: %s") %
                           status % sky1.asDegrees() % sky2.asDegrees() % wcs_errmsg[status])
                                  .str());
    }

    // wcslib assumes 1-indexed coords
    return geom::Point2D(pixTmp[0] + PixelZeroPos + fitsToLsstPixels,
                         pixTmp[1] + PixelZeroPos + fitsToLsstPixels);
}

GeomPoint Wcs::skyToPixel(geom::SpherePoint const& coord) const {
    return skyToPixelImpl(coord.getLongitude(), coord.getLatitude());
}

GeomPoint Wcs::skyToPixel(geom::Angle sky1, geom::Angle sky2) const { return skyToPixelImpl(sky1, sky2); }

GeomPoint Wcs::skyToIntermediateWorldCoord(geom::SpherePoint const& coord) const {
    assert(_wcsInfo);

    double skyTmp[2];
    double imgcrd[2];
    double phi, theta;
    double pixTmp[2];

    skyTmp[_wcsInfo->lng] = coord.getLongitude().asDegrees();
    skyTmp[_wcsInfo->lat] = coord.getLatitude().asDegrees();

    // Estimate pixel coordinates
    int stat[1];
    int status = 0;
    imgcrd[0] = imgcrd[1] = -1e6;
    /*
     printf("  skyTmp[] = (%.3f, %.3f)\n", skyTmp[0], skyTmp[1]);
     printf("  _wcsInfo->lng,lat = %i, %i\n", _wcsInfo->lng, _wcsInfo->lat);
     */
    status = wcss2p(_wcsInfo, 1, 2, skyTmp, &phi, &theta, imgcrd, pixTmp, stat);
    if (status > 0) {
        throw LSST_EXCEPT(except::RuntimeError,
                          (boost::format("Error: wcslib returned a status code of %d at sky %s, %s deg: %s") %
                           status % skyTmp[0] % skyTmp[1] % wcs_errmsg[status])
                                  .str());
    }
    /*
     printf("->iwc (%.3f, %.3f)\n", imgcrd[0], imgcrd[1]);
     printf("-> pix (%.2f, %.2f)\n", pixTmp[0], pixTmp[1]);
     std::shared_ptr<geom::SpherePoint> crval = getSkyOrigin();
     printf("(crval is (%.3f, %.3f))\n", crval->getLongitude().asDegrees(), crval->getLatitude().asDegrees());
     */
    return GeomPoint(imgcrd[0], imgcrd[1]);
}

double Wcs::getEquinox() const {
    if (_wcsInfo == nullptr) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return _wcsInfo->equinox;
}

bool Wcs::isSameSkySystem(Wcs const& wcs) const {
    if (_isIcrs() && wcs._isIcrs()) {
        return true;
    }
    return (getCoordSystem() == wcs.getCoordSystem()) && (getEquinox() == wcs.getEquinox());
}

void Wcs::pixelToSkyImpl(double pixel1, double pixel2, geom::Angle skyTmp[2]) const {
    assert(_wcsInfo);

    // wcslib assumes 1-indexed coordinates
    double pixTmp[2] = {pixel1 - PixelZeroPos + lsstToFitsPixels, pixel2 - PixelZeroPos + lsstToFitsPixels};
    double imgcrd[2];
    double phi, theta;

    double sky[2];
    int status = 0;
    status = wcsp2s(_wcsInfo, 1, 2, pixTmp, imgcrd, &phi, &theta, sky, &status);
    if (status > 0) {
        throw LSST_EXCEPT(except::RuntimeError,
                          (boost::format("Error: wcslib returned a status code of %d at pixel %s, %s: %s") %
                           status % pixel1 % pixel2 % wcs_errmsg[status])
                                  .str());
    }
    // FIXME -- _wcsInfo.lat, _wcsInfo.lng ?
    skyTmp[0] = sky[0] * geom::degrees;
    skyTmp[1] = sky[1] * geom::degrees;
}

CoordPtr Wcs::pixelToSky(GeomPoint const& pixel) const { return pixelToSky(pixel.getX(), pixel.getY()); }

CoordPtr Wcs::pixelToSky(double pixel1, double pixel2) const {
    assert(_wcsInfo);

    geom::Angle skyTmp[2];
    pixelToSkyImpl(pixel1, pixel2, skyTmp);
    return std::make_shared<geom::SpherePoint>(skyTmp[0], skyTmp[1]);
}

void Wcs::pixelToSky(double pixel1, double pixel2, geom::Angle& sky1, geom::Angle& sky2) const {
    geom::Angle skyTmp[2];
    // HACK -- we shouldn't need to initialize these -- pixelToSkyImpl() sets them unless an
    // exception is thrown -- but be safe.
    skyTmp[0] = 0. * geom::radians;
    skyTmp[1] = 0. * geom::radians;
    pixelToSkyImpl(pixel1, pixel2, skyTmp);
    sky1 = skyTmp[0];
    sky2 = skyTmp[1];
}

geom::AffineTransform Wcs::linearizePixelToSky(geom::SpherePoint const& coord, geom::AngleUnit skyUnit) const {
    return linearizePixelToSkyInternal(skyToPixel(coord), coord, skyUnit);
}
geom::AffineTransform Wcs::linearizePixelToSky(GeomPoint const& pix, geom::AngleUnit skyUnit) const {
    return linearizePixelToSkyInternal(pix, *pixelToSky(pix), skyUnit);
}

geom::AffineTransform Wcs::linearizePixelToSkyInternal(GeomPoint const& pix00, geom::SpherePoint const& coord,
                                                       geom::AngleUnit skyUnit) const {
    //
    // Figure out the (0, 0), (0, 1), and (1, 0) ra/dec coordinates of the corners of a square drawn in pixel
    // It'd be better to centre the square at sky00, but that would involve another conversion between sky and
    // pixel coordinates so I didn't bother
    //
    const double side = 10;  // length of the square's sides in pixels
    GeomPoint const sky00 = coord.getPosition(skyUnit);
    typedef std::pair<geom::Angle, geom::Angle> AngleAngle;
    AngleAngle const dsky10 =
            coord::getTangentPlaneOffset(*pixelToSky(pix00 + geom::Extent2D(side, 0)), coord);
    AngleAngle const dsky01 =
            coord::getTangentPlaneOffset(*pixelToSky(pix00 + geom::Extent2D(0, side)), coord);

    Eigen::Matrix2d m;
    m(0, 0) = dsky10.first.asAngularUnits(skyUnit) / side;
    m(0, 1) = dsky01.first.asAngularUnits(skyUnit) / side;
    m(1, 0) = dsky10.second.asAngularUnits(skyUnit) / side;
    m(1, 1) = dsky01.second.asAngularUnits(skyUnit) / side;

    Eigen::Vector2d sky00v;
    sky00v << sky00.getX(), sky00.getY();
    Eigen::Vector2d pix00v;
    pix00v << pix00.getX(), pix00.getY();
    // return geom::AffineTransform(m, geom::Extent2D(sky00v - m * pix00v));
    return geom::AffineTransform(m, (sky00v - m * pix00v));
}

geom::AffineTransform Wcs::linearizeSkyToPixel(geom::SpherePoint const& coord, geom::AngleUnit skyUnit) const {
    return linearizeSkyToPixelInternal(skyToPixel(coord), coord, skyUnit);
}

geom::AffineTransform Wcs::linearizeSkyToPixel(GeomPoint const& pix, geom::AngleUnit skyUnit) const {
    return linearizeSkyToPixelInternal(pix, *pixelToSky(pix), skyUnit);
}

geom::AffineTransform Wcs::linearizeSkyToPixelInternal(GeomPoint const& pix00, geom::SpherePoint const& coord,
                                                       geom::AngleUnit skyUnit) const {
    geom::AffineTransform inverse = linearizePixelToSkyInternal(pix00, coord, skyUnit);
    return inverse.invert();
}

geom::LinearTransform Wcs::getLinearTransform() const { return geom::LinearTransform(getCDMatrix()); }

// -------- table-based persistence -------------------------------------------------------------------------

class WcsFactory : public table::io::PersistableFactory {
public:
    explicit WcsFactory(std::string const& name) : table::io::PersistableFactory(name) {}

    virtual std::shared_ptr<table::io::Persistable> read(InputArchive const& archive,
                                                         CatalogVector const& catalogs) const;
};

namespace {

// Read-only singleton struct containing the schema and keys that a simple Wcs is mapped
// to in record persistence.
struct WcsPersistenceHelper {
    table::Schema schema;
    table::PointKey<double> crval;
    table::PointKey<double> crpix;
    table::Key<table::Array<double> > cd;
    table::Key<std::string> ctype1;
    table::Key<std::string> ctype2;
    table::Key<double> equinox;
    table::Key<std::string> radesys;
    table::Key<std::string> cunit1;
    table::Key<std::string> cunit2;

    static WcsPersistenceHelper const& get() {
        static WcsPersistenceHelper instance;
        return instance;
    };

    // No copying
    WcsPersistenceHelper(const WcsPersistenceHelper&) = delete;
    WcsPersistenceHelper& operator=(const WcsPersistenceHelper&) = delete;

    // No moving
    WcsPersistenceHelper(WcsPersistenceHelper&&) = delete;
    WcsPersistenceHelper& operator=(WcsPersistenceHelper&&) = delete;

private:
    WcsPersistenceHelper()
            : schema(),
              crval(table::PointKey<double>::addFields(schema, "crval", "celestial reference point", "deg")),
              crpix(table::PointKey<double>::addFields(schema, "crpix", "pixel reference point", "pixel")),
              cd(schema.addField<table::Array<double> >(
                      "cd", "linear transform matrix, ordered (1_1, 2_1, 1_2, 2_2)", 4)),
              ctype1(schema.addField<std::string>("ctype1", "coordinate type", 72)),
              ctype2(schema.addField<std::string>("ctype2", "coordinate type", 72)),
              equinox(schema.addField<double>("equinox", "equinox of coordinates")),
              radesys(schema.addField<std::string>("radesys", "coordinate system for equinox", 72)),
              cunit1(schema.addField<std::string>("cunit1", "coordinate units", 72)),
              cunit2(schema.addField<std::string>("cunit2", "coordinate units", 72)) {
        schema.getCitizen().markPersistent();
    }
};

std::string getWcsPersistenceName() { return "Wcs"; }

WcsFactory registration(getWcsPersistenceName());

}  // anonymous

std::string Wcs::getPersistenceName() const { return getWcsPersistenceName(); }

std::string Wcs::getPythonModule() const { return "lsst.afw.image"; }

void Wcs::write(OutputArchiveHandle& handle) const {
    WcsPersistenceHelper const& keys = WcsPersistenceHelper::get();
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    std::shared_ptr<afw::table::BaseRecord> record = catalog.addNew();
    record->set(keys.crval.getX(), _wcsInfo[0].crval[0]);
    record->set(keys.crval.getY(), _wcsInfo[0].crval[1]);
    record->set(keys.crpix, getPixelOrigin());
    Eigen::Matrix2d cdIn = getCDMatrix();
    Eigen::Map<Eigen::Matrix2d> cdOut((*record)[keys.cd].getData());
    cdOut = cdIn;
    record->set(keys.ctype1, std::string(_wcsInfo[0].ctype[0]));
    record->set(keys.ctype2, std::string(_wcsInfo[0].ctype[1]));
    record->set(keys.equinox, _wcsInfo[0].equinox);
    record->set(keys.radesys, std::string(_wcsInfo[0].radesys));
    record->set(keys.cunit1, std::string(_wcsInfo[0].cunit[0]));
    record->set(keys.cunit2, std::string(_wcsInfo[0].cunit[1]));
    handle.saveCatalog(catalog);
}

bool Wcs::_mayBePersistable() const {
    if (_wcsInfo[0].naxis != 2) return false;
    if (std::strcmp(_wcsInfo[0].cunit[0], "deg") != 0) return false;
    if (std::strcmp(_wcsInfo[0].cunit[1], "deg") != 0) return false;

    return true;
}

bool Wcs::isPersistable() const {
    if (!_mayBePersistable()) {
        return false;
    }
    // The current table persistence only works for TAN and TAN-SIP projections
    return false;
}

Wcs::Wcs(afw::table::BaseRecord const& record)
        : daf::base::Citizen(typeid(this)),
          _wcsInfo(NULL),
          _nWcsInfo(0),
          _relax(0),
          _wcsfixCtrl(0),
          _wcshdrCtrl(0),
          _nReject(0),
          _coordSystem(static_cast<afw::coord::CoordSystem>(-1))  // set by _initWcs
{
    WcsPersistenceHelper const& keys = WcsPersistenceHelper::get();
    if (!record.getSchema().contains(keys.schema)) {
        throw LSST_EXCEPT(afw::table::io::MalformedArchiveError, "Incorrect schema for Wcs persistence");
    }
    _setWcslibParams();
    Eigen::Matrix2d cd = Eigen::Map<Eigen::Matrix2d const>(record[keys.cd].getData());
    initWcsLib(record.get(keys.crval), record.get(keys.crpix), cd, record.get(keys.ctype1),
               record.get(keys.ctype2), record.get(keys.equinox), record.get(keys.radesys),
               record.get(keys.cunit1), record.get(keys.cunit2));
    _initWcs();
}

std::shared_ptr<table::io::Persistable> WcsFactory::read(InputArchive const& inputs,
                                                         CatalogVector const& catalogs) const {
    WcsPersistenceHelper const& keys = WcsPersistenceHelper::get();
    LSST_ARCHIVE_ASSERT(catalogs.size() >= 1u);
    LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
    LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys.schema);
    std::shared_ptr<Wcs> result(new Wcs(catalogs.front().front()));
    return result;
}

// ----------------------------------------------------------------------------------------------------------

//
// Mutators
//

void Wcs::shiftReferencePixel(double dx, double dy) {
    assert(_wcsInfo);
    _wcsInfo->crpix[0] += dx;
    _wcsInfo->crpix[1] += dy;

    // tells libwcs to invalidate cached data, since transformation has been modified
    _wcsInfo->flag = 0;
}

/*
 * Now WCSA, pixel coordinates, but allowing for X0 and Y0
 */
namespace detail {

/**
 * @internal Define a trivial WCS that maps the lower left corner (LLC) pixel of an image to a given value
 */
std::shared_ptr<daf::base::PropertyList> createTrivialWcsAsPropertySet(
        std::string const& wcsName,  ///< @internal Name of desired WCS
        int const x0,                ///< @internal Column coordinate of LLC pixel
        int const y0                 ///< @internal Row coordinate of LLC pixel
        ) {
    std::shared_ptr<daf::base::PropertyList> wcsMetaData(new daf::base::PropertyList);

    wcsMetaData->set("CRVAL1" + wcsName, x0, "Column pixel of Reference Pixel");
    wcsMetaData->set("CRVAL2" + wcsName, y0, "Row pixel of Reference Pixel");
    wcsMetaData->set("CRPIX1" + wcsName, 1, "Column Pixel Coordinate of Reference");
    wcsMetaData->set("CRPIX2" + wcsName, 1, "Row Pixel Coordinate of Reference");
    wcsMetaData->set("CTYPE1" + wcsName, "LINEAR", "Type of projection");
    wcsMetaData->set("CTYPE2" + wcsName, "LINEAR", "Type of projection");
    wcsMetaData->set("CUNIT1" + wcsName, "PIXEL", "Column unit");
    wcsMetaData->set("CUNIT2" + wcsName, "PIXEL", "Row unit");

    return wcsMetaData;
}
/**
 * @internal Return a Point2I(x0, y0) given a PropertySet containing a suitable WCS (e.g. "A")
 *
 * The WCS must have CRPIX[12] == 1 and CRVAL[12] must be present.  If this is true, the WCS
 * cards are removed from the metadata
 *
 * @param wcsName the WCS to search (E.g. "A")
 * @param metadata the metadata, maybe containing the WCS
 */
geom::Point2I getImageXY0FromMetadata(std::string const& wcsName, daf::base::PropertySet* metadata) {
    int x0 = 0;  // Our value of X0
    int y0 = 0;  // Our value of Y0

    try {
        //
        // Only use WCS if CRPIX[12] == 1 and CRVAL[12] is present
        //
        if (metadata->getAsDouble("CRPIX1" + wcsName) == 1 &&
            metadata->getAsDouble("CRPIX2" + wcsName) == 1) {
            x0 = metadata->getAsInt("CRVAL1" + wcsName);
            y0 = metadata->getAsInt("CRVAL2" + wcsName);
            //
            // OK, we've got it.  Remove it from the header
            //
            metadata->remove("CRVAL1" + wcsName);
            metadata->remove("CRVAL2" + wcsName);
            metadata->remove("CRPIX1" + wcsName);
            metadata->remove("CRPIX2" + wcsName);
            metadata->remove("CTYPE1" + wcsName);
            metadata->remove("CTYPE1" + wcsName);
            metadata->remove("CUNIT1" + wcsName);
            metadata->remove("CUNIT2" + wcsName);
        }
    } catch (pex::exceptions::NotFoundError&) {
        ;  // OK, not present
    }

    return geom::Point2I(x0, y0);
}

/**
 * @internal Strip keywords from the input metadata that are related to the generated Wcs
 *
 * It isn't entirely obvious that this is enough --- e.g. if the input metadata has deprecated
 * WCS keywords such as CDELT[12] they won't be stripped.  Well, actually we catch CDELT[12], LTV[12], and
 * PC00[12]00[12]
 * but there may be others
 *
 * @param metadata Metadata to be stripped
 * @param wcs A Wcs with (implied) keywords
 */
int stripWcsKeywords(std::shared_ptr<daf::base::PropertySet> const& metadata,
                     std::shared_ptr<Wcs const> const& wcs) {
    std::shared_ptr<daf::base::PropertySet> wcsMetadata = wcs->getFitsMetadata();
    std::vector<std::string> paramNames = wcsMetadata->paramNames();
    paramNames.push_back("CDELT1");
    paramNames.push_back("CDELT2");
    paramNames.push_back("LTV1");
    paramNames.push_back("LTV2");
    paramNames.push_back("PC001001");
    paramNames.push_back("PC001002");
    paramNames.push_back("PC002001");
    paramNames.push_back("PC002002");
    for (std::vector<std::string>::const_iterator ptr = paramNames.begin(); ptr != paramNames.end(); ++ptr) {
        metadata->remove(*ptr);
    }

    return 0;  // would be ncard if remove returned a status
}
}

// -------------------------------------------------------------------------------------------------
//
// XYTransformFromWcsPair

XYTransformFromWcsPair::XYTransformFromWcsPair(std::shared_ptr<Wcs const> dst, std::shared_ptr<Wcs const> src)
        : XYTransform(), _dst(dst), _src(src), _isSameSkySystem(dst->isSameSkySystem(*src)) {}

std::shared_ptr<geom::XYTransform> XYTransformFromWcsPair::clone() const {
    return std::make_shared<XYTransformFromWcsPair>(_dst->clone(), _src->clone());
}

geom::Point2D XYTransformFromWcsPair::forwardTransform(Point2D const& pixel) const {
    if (_isSameSkySystem) {
        // high performance branch; no coordinate conversion required
        afw::geom::Angle sky0, sky1;
        _src->pixelToSky(pixel[0], pixel[1], sky0, sky1);
        return _dst->skyToPixel(sky0, sky1);
    } else {
        std::shared_ptr<afw::geom::SpherePoint const> coordPtr{_src->pixelToSky(pixel)};
        return _dst->skyToPixel(*coordPtr);
    }
}

geom::Point2D XYTransformFromWcsPair::reverseTransform(Point2D const& pixel) const {
    if (_isSameSkySystem) {
        // high performance branch; no coordinate conversion required
        afw::geom::Angle sky0, sky1;
        _dst->pixelToSky(pixel[0], pixel[1], sky0, sky1);
        return _src->skyToPixel(sky0, sky1);
    } else {
        std::shared_ptr<afw::geom::SpherePoint const> coordPtr{_dst->pixelToSky(pixel)};
        return _src->skyToPixel(*coordPtr);
    }
}

std::shared_ptr<geom::XYTransform> XYTransformFromWcsPair::invert() const {
    // just swap src, dst
    return std::make_shared<XYTransformFromWcsPair>(_src, _dst);
}
}
}
}  // end image
