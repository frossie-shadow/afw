// -*- lsst-c++ -*-
#ifndef LSST_AFW_TABLE_slots_h_INCLUDED
#define LSST_AFW_TABLE_slots_h_INCLUDED

#include "lsst/afw/table/aggregates.h"

namespace lsst {

namespace daf { namespace base {
class PropertySet;
}} // namespace daf::base

namespace afw {

namespace fits {
class Fits;
}; // namespace fits

namespace table {

class SlotDefinition {
public:

    explicit SlotDefinition(std::string const & name) : _name(name) {}

    std::string getName() const { return _name; }

    std::string getAlias(int version) const {
        return (version > 0 ? "slot_" : "slot.") + _name;
    }

protected:
    std::string _name;
};

class FluxSlotDefinition : public SlotDefinition {
public:

    typedef double MeasValue;
    typedef double ErrValue;
    typedef Key<double> MeasKey;
    typedef Key<double> ErrKey;

    explicit FluxSlotDefinition(std::string const & name) : SlotDefinition(name) {}

    bool isValid() const { return _measKey.isValid(); }

    MeasKey getMeasKey() const { return _measKey; }

    ErrKey getErrKey() const { return _errKey; }

    Key<Flag> getFlagKey() const { return _flagKey; }

    void setKeys(std::string const & alias, Schema const & schema);

private:
    MeasKey _measKey;
    ErrKey _errKey;
    Key<Flag> _flagKey;
};

class CentroidSlotDefinition : public SlotDefinition {
public:

    typedef geom::Point2D MeasValue;
    typedef Eigen::Matrix<float,2,2> ErrValue;
    typedef Point2DKey MeasKey;
    typedef CovarianceMatrixKey<float,2> ErrKey;

    explicit CentroidSlotDefinition(std::string const & name) : SlotDefinition(name) {}

    bool isValid() const { return _measKey.isValid(); }

    MeasKey getMeasKey() const { return _measKey; }

    ErrKey getErrKey() const { return _errKey; }

    Key<Flag> getFlagKey() const { return _flagKey; }

    void setKeys(std::string const & alias, Schema const & schema);

private:
    MeasKey _measKey;
    ErrKey _errKey;
    Key<Flag> _flagKey;
};

class ShapeSlotDefinition : public SlotDefinition {
public:

    typedef geom::ellipses::Quadrupole MeasValue;
    typedef Eigen::Matrix<float,3,3> ErrValue;
    typedef QuadrupoleKey MeasKey;
    typedef CovarianceMatrixKey<float,3> ErrKey;

    explicit ShapeSlotDefinition(std::string const & name) : SlotDefinition(name) {}

    bool isValid() const { return _measKey.isValid(); }

    MeasKey getMeasKey() const { return _measKey; }

    ErrKey getErrKey() const { return _errKey; }

    Key<Flag> getFlagKey() const { return _flagKey; }

    void setKeys(std::string const & alias, Schema const & schema);

private:
    MeasKey _measKey;
    ErrKey _errKey;
    Key<Flag> _flagKey;
};


struct SlotSuite {
    FluxSlotDefinition defPsfFlux;
    FluxSlotDefinition defApFlux;
    FluxSlotDefinition defInstFlux;
    FluxSlotDefinition defModelFlux;
    CentroidSlotDefinition defCentroid;
    ShapeSlotDefinition defShape;

    void handleAliasChange(std::string const & alias, Schema const & schema);

    explicit SlotSuite(Schema const & schema);
};


//@{
/**
 *  Utilities for version 0 tables and measurement algorithms that fill them.
 *
 *  These are deprecated and should not be used by new code; they will be removed when
 *  the new measurement framework in meas_base is complete and the old on in meas_algorithms
 *  is retired.
 */

#ifndef SWIG

template <typename MeasTagT, typename ErrTagT>
struct Measurement {
    typedef MeasTagT MeasTag;  ///< the tag (template parameter) type used for the measurement
    typedef ErrTagT ErrTag;    ///< the tag (template parameter) type used for the uncertainty
    typedef typename Field<MeasTag>::Value MeasValue; ///< the value type used for the measurement
    typedef typename Field<ErrTag>::Value ErrValue;   ///< the value type used for the uncertainty
    typedef Key<MeasTag> MeasKey;  ///< the Key type for the actual measurement
    typedef Key<ErrTag> ErrKey;    ///< the Key type for the error on the measurement
};

/// A collection of types useful for flux measurement algorithms.
struct Flux : public Measurement<double, double> {}; //pgee temporary

/// A collection of types useful for centroid measurement algorithms.
struct Centroid : public Measurement< Point<double>, Covariance< Point<float> > > {};

/// A collection of types useful for shape measurement algorithms.
struct Shape : public Measurement< Moments<double>, Covariance< Moments<float> > > {};

/**
 *  @brief A three-element tuple of measurement, uncertainty, and flag keys.
 *
 *  Most measurement should have more than one flag key to indicate different kinds of failures.
 *  This flag key should usually be set to be a logical OR of all of them, so it is set whenever
 *  a measurement cannot be fully trusted.
 */
template <typename MeasurementT>
struct KeyTuple {
    typename MeasurementT::MeasKey meas; ///< Key used for the measured value.
    typename MeasurementT::ErrKey err;   ///< Key used for the uncertainty.
    Key<Flag> flag;                      ///< Failure bit; set if the measurement did not fully succeed.

    /// Default-constructor; all keys will be invalid.
    KeyTuple() {}

    /// Main constructor.
    KeyTuple(
        typename MeasurementT::MeasKey const & meas_,
        typename MeasurementT::ErrKey const & err_,
        Key<Flag> const & flag_
    ) : meas(meas_), err(err_), flag(flag_) {}

};

/// Convenience function to setup fields for centroid measurement algorithms.
KeyTuple<Centroid> addCentroidFields(Schema & schema, std::string const & name, std::string const & doc);

/// Convenience function to setup fields for shape measurement algorithms.
KeyTuple<Shape> addShapeFields(Schema & schema, std::string const & name, std::string const & doc);

/// Convenience function to setup fields for flux measurement algorithms.
KeyTuple<Flux> addFluxFields(Schema & schema, std::string const & name, std::string const & doc);

#endif // !SWIG

//@

}}} // lsst::afw::table

#endif // !LSST_AFW_TABLE_slots_h_INCLUDED