// -*- LSST-C++ -*-
/** \file
 * \brief Implementation of PsfFormatter class
 *
 * \version $Revision: 2151 $
 * \date $Date$
 *
 * Contact: Kian-Tat Lim (ktl@slac.stanford.edu)
 *
 * \ingroup afw
 */
#include <stdexcept>
#include <string>
#include <vector>

#include "boost/serialization/nvp.hpp"

#include "lsst/afw/detection/PsfFormatter.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/daf/persistence/FormatterImpl.h"
#include "lsst/daf/persistence/LogicalLocation.h"
#include "lsst/daf/persistence/BoostStorage.h"
#include "lsst/daf/persistence/XmlStorage.h"
#include "lsst/pex/exceptions.h"
#include "lsst/log/Log.h"
#include "lsst/pex/policy/Policy.h"

BOOST_CLASS_EXPORT(lsst::afw::detection::Psf)

namespace {
LOG_LOGGER _log = LOG_GET("afw.detection.PsfFormatter");
}

namespace afwDetect = lsst::afw::detection;
namespace afwMath = lsst::afw::math;
namespace dafBase = lsst::daf::base;
namespace dafPersist = lsst::daf::persistence;
namespace pexPolicy = lsst::pex::policy;

using boost::serialization::make_nvp;

/** Register this Formatter subclass through a static instance of
 * FormatterRegistration.
 */
dafPersist::FormatterRegistration
afwDetect::PsfFormatter::registration("Psf", typeid(afwDetect::Psf), createInstance);

/** Constructor.
 * \param[in] policy Policy for configuring this Formatter
 */
afwDetect::PsfFormatter::PsfFormatter(
    pexPolicy::Policy::Ptr policy) :
    dafPersist::Formatter(typeid(this)), _policy(policy) {}

/** Minimal destructor.
 */
afwDetect::PsfFormatter::~PsfFormatter(void) {}

void afwDetect::PsfFormatter::write(
    dafBase::Persistable const* persistable,
    dafPersist::Storage::Ptr storage,
    dafBase::PropertySet::Ptr) {
    LOGL_DEBUG(_log, "PsfFormatter write start");
    afwDetect::Psf const* ps = dynamic_cast<afwDetect::Psf const*>(persistable);
    if (ps == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Persisting non-Psf");
    }
    if (typeid(*storage) == typeid(dafPersist::BoostStorage)) {
        LOGL_DEBUG(_log, "PsfFormatter write BoostStorage");
        dafPersist::BoostStorage* boost =
            dynamic_cast<dafPersist::BoostStorage*>(storage.get());
        boost->getOArchive() & ps;
        LOGL_DEBUG(_log, "PsfFormatter write end");
        return;
    }
    else if (typeid(*storage) == typeid(dafPersist::XmlStorage)) {
        LOGL_DEBUG(_log, "PsfFormatter write XmlStorage");
        dafPersist::XmlStorage* xml =
            dynamic_cast<dafPersist::XmlStorage*>(storage.get());
        xml->getOArchive() & make_nvp("psf", ps);
        LOGL_DEBUG(_log, "PsfFormatter write end");
        return;
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Unrecognized Storage for Psf");
}

dafBase::Persistable* afwDetect::PsfFormatter::read(
    dafPersist::Storage::Ptr storage, dafBase::PropertySet::Ptr) {
    LOGL_DEBUG(_log, "PsfFormatter read start");
    afwDetect::Psf* ps;
    if (typeid(*storage) == typeid(dafPersist::BoostStorage)) {
        LOGL_DEBUG(_log, "PsfFormatter read BoostStorage");
        dafPersist::BoostStorage* boost =
            dynamic_cast<dafPersist::BoostStorage*>(storage.get());
        boost->getIArchive() & ps;
        LOGL_DEBUG(_log, "PsfFormatter read end");
        return ps;
    }
    else if (typeid(*storage) == typeid(dafPersist::XmlStorage)) {
        LOGL_DEBUG(_log, "PsfFormatter read XmlStorage");
        dafPersist::XmlStorage* xml =
            dynamic_cast<dafPersist::XmlStorage*>(storage.get());
        xml->getIArchive() & make_nvp("psf", ps);
        LOGL_DEBUG(_log, "PsfFormatter read end");
        return ps;
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Unrecognized Storage for Psf");
}

void afwDetect::PsfFormatter::update(dafBase::Persistable* ,
                                   dafPersist::Storage::Ptr,
                                   dafBase::PropertySet::Ptr) {
    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Unexpected call to update for Psf");
}

/** Serialize a Psf to a Boost archive.  Handles text or XML
 * archives, input or output.
 */
template <class Archive>
void afwDetect::PsfFormatter::delegateSerialize(
        Archive& ar,                    ///< Boost archive
        unsigned int const,             ///< Version of the Psf class
        dafBase::Persistable* persistable ///< persistable Pointer to the Psf as a Persistable
                                               ) {
    LOGL_DEBUG(_log, "PsfFormatter delegateSerialize start");
    afwDetect::Psf* ps = dynamic_cast<afwDetect::Psf*>(persistable);
    if (ps == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Serializing non-Psf");
    }
#if 0                                   // not present in baseclass
    ar & make_nvp("width", ps->_width) & make_nvp("height", ps->_height);
    ar & make_nvp("k", ps->_kernel);
#endif

    LOGL_DEBUG(_log, "PsfFormatter delegateSerialize end");
}

/** Factory method for PsfFormatter.
 * \param[in] policy Policy for configuring the PsfFormatter
 * \return Shared pointer to a new instance
 */
dafPersist::Formatter::Ptr afwDetect::PsfFormatter::createInstance(
    pexPolicy::Policy::Ptr policy) {
    return dafPersist::Formatter::Ptr(new afwDetect::PsfFormatter(policy));
}
