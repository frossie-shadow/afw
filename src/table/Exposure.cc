// -*- lsst-c++ -*-
#include <memory>
#include <typeinfo>
#include <string>

#include "lsst/daf/base/PropertySet.h"
#include "lsst/daf/base/PropertyList.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/FitsWriter.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/detail/Access.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/ApCorrMap.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/image/VisitInfo.h"

namespace lsst {
namespace afw {
namespace table {

//-----------------------------------------------------------------------------------------------------------
//----- Private ExposureTable/Record classes ---------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// These private derived classes are what you actually get when you do ExposureTable::make; like the
// private classes in BaseTable.cc, it's more convenient to have an extra set of trivial derived
// classes than to do a lot of friending.

namespace {

int const EXPOSURE_TABLE_CURRENT_VERSION = 2;                   // current version of ExposureTable
std::string const EXPOSURE_TABLE_VERSION_KEY = "EXPTABLE_VER";  // FITS header key for ExposureTable version

// Field names used to store the archive IDs of different components (used in
// multiple places, so we define them here instead of as multiple string
// literals).
std::string const WCS_FIELD_NAME = "wcs";
std::string const PSF_FIELD_NAME = "psf";
std::string const CALIB_FIELD_NAME = "calib";
std::string const VISIT_INFO_FIELD_NAME = "visitInfo";
std::string const AP_CORR_MAP_FIELD_NAME = "apCorrMap";
std::string const VALID_POLYGON_FIELD_NAME = "validPolygon";

int getTableVersion(daf::base::PropertySet &metadata) {
    return metadata.exists(EXPOSURE_TABLE_VERSION_KEY) ? metadata.get<int>(EXPOSURE_TABLE_VERSION_KEY) : 1;
}

/**
 * @internal Helper class for for persisting ExposureRecord
 *
 * Contains keys for columns beyond BaseRecord, a schema mapper and and helper functions
 */
struct PersistenceHelper {
    Schema schema;
    Key<int> wcs;
    Key<int> psf;
    Key<int> calib;
    Key<int> apCorrMap;
    Key<int> validPolygon;
    Key<int> visitInfo;

    // Create a SchemaMapper that maps an ExposureRecord to a BaseRecord with IDs for Wcs, Psf, etc.
    SchemaMapper makeWriteMapper(Schema const &inputSchema) const {
        std::vector<Schema> inSchemas;
        inSchemas.push_back(PersistenceHelper().schema);
        inSchemas.push_back(inputSchema);
        // don't need front; it's an identity mapper
        SchemaMapper result = SchemaMapper::join(inSchemas).back();
        result.editOutputSchema().setAliasMap(inputSchema.getAliasMap());
        return result;
    }

    // Create a SchemaMapper that maps a BaseRecord to an ExposureRecord with IDs for WCS, Psf, etc.
    SchemaMapper makeReadMapper(Schema const &inputSchema) const {
        SchemaMapper result = SchemaMapper::removeMinimalSchema(inputSchema, schema);
        result.editOutputSchema().setAliasMap(inputSchema.getAliasMap());
        return result;
    }

    // Write psf, wcs, etc. from an ExposureRecord to an archive
    template <typename OutputArchiveIsh>
    void writeRecord(ExposureRecord const &input, BaseRecord &output, SchemaMapper const &mapper,
                     OutputArchiveIsh &archive, bool permissive) const {
        output.assign(input, mapper);
        output.set(psf, archive.put(input.getPsf(), permissive));
        output.set(wcs, archive.put(input.getWcs(), permissive));
        output.set(calib, archive.put(input.getCalib(), permissive));
        output.set(apCorrMap, archive.put(input.getApCorrMap(), permissive));
        output.set(validPolygon, archive.put(input.getValidPolygon(), permissive));
        output.set(visitInfo, archive.put(input.getVisitInfo(), permissive));
    }

    // Read psf, wcs, etc. from an archive to an ExposureRecord
    void readRecord(BaseRecord const &input, ExposureRecord &output, SchemaMapper const &mapper,
                    io::InputArchive const &archive) const {
        output.assign(input, mapper);
        if (psf.isValid()) {
            output.setPsf(archive.get<detection::Psf>(input.get(psf)));
        }
        if (wcs.isValid()) {
            output.setWcs(archive.get<image::Wcs>(input.get(wcs)));
        }
        if (calib.isValid()) {
            output.setCalib(archive.get<image::Calib>(input.get(calib)));
        }
        if (apCorrMap.isValid()) {
            output.setApCorrMap(archive.get<image::ApCorrMap>(input.get(apCorrMap)));
        }
        if (validPolygon.isValid()) {
            output.setValidPolygon(archive.get<geom::polygon::Polygon>(input.get(validPolygon)));
        }
        if (visitInfo.isValid()) {
            output.setVisitInfo(archive.get<image::VisitInfo>(input.get(visitInfo)));
        }
    }

    // No copying
    PersistenceHelper(const PersistenceHelper &) = delete;
    PersistenceHelper &operator=(const PersistenceHelper &) = delete;

    // No moving
    PersistenceHelper(PersistenceHelper &&) = delete;
    PersistenceHelper &operator=(PersistenceHelper &&) = delete;

    // Construct a PersistenceHelper using the most modern schema.
    PersistenceHelper()
            : schema(),
              wcs(schema.addField<int>(WCS_FIELD_NAME, "archive ID for Wcs object")),
              psf(schema.addField<int>(PSF_FIELD_NAME, "archive ID for Psf object")),
              calib(schema.addField<int>(CALIB_FIELD_NAME, "archive ID for Calib object")),
              apCorrMap(schema.addField<int>(AP_CORR_MAP_FIELD_NAME, "archive ID for ApCorrMap object")),
              validPolygon(schema.addField<int>(VALID_POLYGON_FIELD_NAME, "archive ID for Polygon object")),
              visitInfo(schema.addField<int>(VISIT_INFO_FIELD_NAME, "archive ID for VisitInfo object")) {}

    // Add a field to this->schema, saving its key in 'key', if and only if 'name' is a field in 'oldSchema'
    void addIfPresent(Schema const &oldSchema, Key<int> &key, std::string const &name) {
        try {
            auto item = oldSchema.find<int>(name);
            key = schema.addField(item.field);
        } catch (pex::exceptions::NotFoundError &) {
        }
    }

    // Construct a PersistenceHelper from a possibly old on-disk schema
    PersistenceHelper(Schema const &oldSchema) {
        addIfPresent(oldSchema, wcs, WCS_FIELD_NAME);
        addIfPresent(oldSchema, psf, PSF_FIELD_NAME);
        addIfPresent(oldSchema, calib, CALIB_FIELD_NAME);
        addIfPresent(oldSchema, apCorrMap, AP_CORR_MAP_FIELD_NAME);
        addIfPresent(oldSchema, validPolygon, VALID_POLYGON_FIELD_NAME);
        addIfPresent(oldSchema, visitInfo, VISIT_INFO_FIELD_NAME);
        assert(oldSchema.contains(schema));
    }
};

}  // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ExposureFitsWriter ---------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// A custom FitsWriter for Exposure - this sets the AFW_TYPE key to EXPOSURE, which should ensure
// we use ExposureFitsReader to read it, and sets EXPOSURE_TABLE_VERSION_KEY to the current version:
// EXPOSURE_TABLE_CURRENT_VERSION

namespace {

class ExposureFitsWriter : public io::FitsWriter {
public:
    ExposureFitsWriter(Fits *fits, std::shared_ptr<io::OutputArchive> archive, int flags)
            : io::FitsWriter(fits, flags), _doWriteArchive(false), _archive(archive), _helper() {
        if (!_archive) {
            _doWriteArchive = true;
            _archive.reset(new io::OutputArchive());
        }
    }

protected:
    virtual void _writeTable(std::shared_ptr<BaseTable const> const &table, std::size_t nRows);

    virtual void _writeRecord(BaseRecord const &r);

    virtual void _finish() {
        if (_doWriteArchive) _archive->writeFits(*_fits);
    }

    bool _doWriteArchive;
    std::shared_ptr<io::OutputArchive> _archive;
    std::shared_ptr<BaseRecord> _record;
    PersistenceHelper _helper;
    SchemaMapper _mapper;
};

void ExposureFitsWriter::_writeTable(std::shared_ptr<BaseTable const> const &t, std::size_t nRows) {
    std::shared_ptr<ExposureTable const> inTable = std::dynamic_pointer_cast<ExposureTable const>(t);
    if (!inTable) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                          "Cannot use a ExposureFitsWriter on a non-Exposure table.");
    }
    _mapper = _helper.makeWriteMapper(inTable->getSchema());
    std::shared_ptr<BaseTable> outTable = BaseTable::make(_mapper.getOutputSchema());
    io::FitsWriter::_writeTable(outTable, nRows);
    _fits->writeKey("AFW_TYPE", "EXPOSURE", "Tells lsst::afw to load this as an Exposure table.");
    _fits->writeKey(EXPOSURE_TABLE_VERSION_KEY, EXPOSURE_TABLE_CURRENT_VERSION, "Exposure table version");
    _record = outTable->makeRecord();
}

void ExposureFitsWriter::_writeRecord(BaseRecord const &r) {
    ExposureRecord const &record = static_cast<ExposureRecord const &>(r);
    _helper.writeRecord(record, *_record, _mapper, *_archive, false);
    io::FitsWriter::_writeRecord(*_record);
}

}  // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ExposureFitsReader ---------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

// FitsColumnReader that reads a Persistable subclass T (Wcs, Psf, or Calib here) by using an int
// column to retrieve the object from an InputArchive and attach it to an ExposureRecord via
// the Setter member function pointer.
template <typename T, void (ExposureRecord::*Setter)(std::shared_ptr<T const>)>
class PersistableObjectColumnReader : public io::FitsColumnReader {
public:
    static void setup(std::string const &name, io::FitsSchemaInputMapper &mapper) {
        auto item = mapper.find(name);
        if (item) {
            if (mapper.hasArchive()) {
                std::unique_ptr<io::FitsColumnReader> reader(new PersistableObjectColumnReader(item->column));
                mapper.customize(std::move(reader));
            }
            mapper.erase(item);
        }
    }

    PersistableObjectColumnReader(int column) : _column(column) {}

    virtual void readCell(BaseRecord &record, std::size_t row, fits::Fits &fits,
                          std::shared_ptr<io::InputArchive> const &archive) const {
        int id = 0;
        fits.readTableScalar<int>(row, _column, id);
        std::shared_ptr<T> value = archive->get<T>(id);
        (static_cast<ExposureRecord &>(record).*(Setter))(value);
    }

private:
    bool _noHeavy;
    int _column;
};

namespace {

class ExposureFitsReader : public io::FitsReader {
public:
    ExposureFitsReader() : afw::table::io::FitsReader("EXPOSURE") {}

    virtual std::shared_ptr<BaseTable> makeTable(io::FitsSchemaInputMapper &mapper,
                                                 std::shared_ptr<daf::base::PropertyList> metadata,
                                                 int ioFlags, bool stripMetadata) const {
        // We rely on the table version stored in the metadata when loading an ExposureCatalog
        // persisted on its own.  This is not as flexible in terms of backwards compatibility
        // as the code that loads ExposureCatalogs persisted as part of something else, but
        // we happen to know there are no ExposureCatalogs sitting on disk with with versions
        // older than what this routine supports.
        auto tableVersion = getTableVersion(*metadata);
        PersistableObjectColumnReader<detection::Psf, &ExposureRecord::setPsf>::setup("psf", mapper);
        PersistableObjectColumnReader<image::Wcs, &ExposureRecord::setWcs>::setup("wcs", mapper);
        PersistableObjectColumnReader<image::Calib, &ExposureRecord::setCalib>::setup("calib", mapper);
        PersistableObjectColumnReader<image::ApCorrMap, &ExposureRecord::setApCorrMap>::setup("apCorrMap",
                                                                                              mapper);
        PersistableObjectColumnReader<geom::polygon::Polygon, &ExposureRecord::setValidPolygon>::setup(
                "validPolygon", mapper);
        if (tableVersion > 1) {
            PersistableObjectColumnReader<image::VisitInfo, &ExposureRecord::setVisitInfo>::setup("visitInfo",
                                                                                                  mapper);
        }
        std::shared_ptr<ExposureTable> table = ExposureTable::make(mapper.finalize());
        table->setMetadata(metadata);
        return table;
    }

    virtual bool usesArchive(int ioFlags) const { return true; }
};

static ExposureFitsReader const exposureFitsReader;

}  // anonymous

//-----------------------------------------------------------------------------------------------------------
//----- ExposureTable/Record member function implementations -----------------------------------------------
//-----------------------------------------------------------------------------------------------------------

geom::Box2I ExposureRecord::getBBox() const {
    return geom::Box2I(get(ExposureTable::getBBoxMinKey()), get(ExposureTable::getBBoxMaxKey()));
}

void ExposureRecord::setBBox(geom::Box2I const &bbox) {
    set(ExposureTable::getBBoxMinKey(), bbox.getMin());
    set(ExposureTable::getBBoxMaxKey(), bbox.getMax());
}

bool ExposureRecord::contains(geom::SpherePoint const &coord, bool includeValidPolygon) const {
    if (!getWcs()) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "ExposureRecord does not have a Wcs; cannot call contains()");
    }

    // If there is no valid polygon set to false
    if (includeValidPolygon && !getValidPolygon()) {
        includeValidPolygon = false;
    }

    try {
        geom::Point2D point = getWcs()->skyToPixel(coord);
        if (includeValidPolygon)
            return (geom::Box2D(getBBox()).contains(point) && getValidPolygon()->contains(point));
        else
            return geom::Box2D(getBBox()).contains(point);
    } catch (pex::exceptions::DomainError &) {
        // Wcs can throw if the given coordinate is outside the region
        // where the Wcs is valid.
        return false;
    }
}

bool ExposureRecord::contains(geom::Point2D const &point, image::Wcs const &wcs,
                              bool includeValidPolygon) const {
    return contains(*wcs.pixelToSky(point), includeValidPolygon);
}

ExposureRecord::ExposureRecord(std::shared_ptr<ExposureTable> const &table) : BaseRecord(table) {}

void ExposureRecord::_assign(BaseRecord const &other) {
    try {
        ExposureRecord const &s = dynamic_cast<ExposureRecord const &>(other);
        _psf = s._psf;
        _wcs = s._wcs;
        _calib = s._calib;
        _apCorrMap = s._apCorrMap;
        _validPolygon = s._validPolygon;
        _visitInfo = s._visitInfo;
    } catch (std::bad_cast &) {
    }
}

std::shared_ptr<ExposureTable> ExposureTable::make(Schema const &schema) {
    if (!checkSchema(schema)) {
        throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterError,
                "Schema for Exposure must contain at least the keys defined by makeMinimalSchema().");
    }
    return std::shared_ptr<ExposureTable>(new ExposureTable(schema));
}

ExposureTable::ExposureTable(Schema const &schema) : BaseTable(schema) {}

ExposureTable::ExposureTable(ExposureTable const &other) : BaseTable(other) {}

ExposureTable::MinimalSchema::MinimalSchema() {
    id = schema.addField<RecordId>("id", "unique ID");
    bboxMin = PointKey<int>::addFields(schema, "bbox_min", "bbox minimum point", "pixel");
    bboxMax = PointKey<int>::addFields(schema, "bbox_max", "bbox maximum point", "pixel");
    schema.getCitizen().markPersistent();
}

ExposureTable::MinimalSchema &ExposureTable::getMinimalSchema() {
    static MinimalSchema it;
    return it;
}

std::shared_ptr<io::FitsWriter> ExposureTable::makeFitsWriter(fits::Fits *fitsfile, int flags) const {
    return std::make_shared<ExposureFitsWriter>(fitsfile, std::shared_ptr<io::OutputArchive>(), flags);
}

std::shared_ptr<io::FitsWriter> ExposureTable::makeFitsWriter(fits::Fits *fitsfile,
                                                              std::shared_ptr<io::OutputArchive> archive,
                                                              int flags) const {
    return std::make_shared<ExposureFitsWriter>(fitsfile, archive, flags);
}

std::shared_ptr<BaseTable> ExposureTable::_clone() const {
    return std::shared_ptr<ExposureTable>(new ExposureTable(*this));
}

std::shared_ptr<BaseRecord> ExposureTable::_makeRecord() {
    return std::shared_ptr<ExposureRecord>(new ExposureRecord(getSelf<ExposureTable>()));
}

//-----------------------------------------------------------------------------------------------------------
//----- ExposureCatalogT member function implementations ----------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

template <typename RecordT>
void ExposureCatalogT<RecordT>::writeToArchive(io::OutputArchiveHandle &handle, bool permissive) const {
    PersistenceHelper helper{};
    SchemaMapper mapper = helper.makeWriteMapper(this->getSchema());
    BaseCatalog outputCat = handle.makeCatalog(mapper.getOutputSchema());
    outputCat.reserve(this->size());
    for (const_iterator i = this->begin(); i != this->end(); ++i) {
        helper.writeRecord(*i, *outputCat.addNew(), mapper, handle, permissive);
    }
    handle.saveCatalog(outputCat);
}

template <typename RecordT>
ExposureCatalogT<RecordT> ExposureCatalogT<RecordT>::readFromArchive(io::InputArchive const &archive,
                                                                     BaseCatalog const &catalog) {
    // Helper constructor will infer which components are available
    // (effectively the version, but more flexible).
    PersistenceHelper helper{catalog.getSchema()};
    SchemaMapper mapper = helper.makeReadMapper(catalog.getSchema());
    ExposureCatalogT<ExposureRecord> result(mapper.getOutputSchema());
    result.reserve(catalog.size());
    for (BaseCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
        helper.readRecord(*i, *result.addNew(), mapper, archive);
    }
    return result;
}

template <typename RecordT>
ExposureCatalogT<RecordT> ExposureCatalogT<RecordT>::subsetContaining(geom::SpherePoint const &coord,
                                                                      bool includeValidPolygon) const {
    ExposureCatalogT result(this->getTable());
    for (const_iterator i = this->begin(); i != this->end(); ++i) {
        if (i->contains(coord, includeValidPolygon)) {
            result.push_back(i);
        }
    }
    return result;
}

template <typename RecordT>
ExposureCatalogT<RecordT> ExposureCatalogT<RecordT>::subsetContaining(geom::Point2D const &point,
                                                                      image::Wcs const &wcs,
                                                                      bool includeValidPolygon) const {
    ExposureCatalogT result(this->getTable());
    for (const_iterator i = this->begin(); i != this->end(); ++i) {
        if (i->contains(point, wcs, includeValidPolygon)) {
            result.push_back(i);
        }
    }
    return result;
}

//-----------------------------------------------------------------------------------------------------------
//----- Explicit instantiation ------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------

template class CatalogT<ExposureRecord>;
template class CatalogT<ExposureRecord const>;

template class SortedCatalogT<ExposureRecord>;
template class SortedCatalogT<ExposureRecord const>;

template class ExposureCatalogT<ExposureRecord>;
template class ExposureCatalogT<ExposureRecord const>;
}
}
}  // namespace lsst::afw::table
