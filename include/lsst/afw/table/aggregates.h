// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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
#ifndef AFW_TABLE_aggregates_h_INCLUDED
#define AFW_TABLE_aggregates_h_INCLUDED

#include "lsst/afw/table/FunctorKey.h"
#include "lsst/afw/table/Schema.h"

namespace lsst {
namespace afw {

namespace geom {

template <typename T, int N>
class Point;

namespace ellipses {

class Quadrupole;

}  // namespace ellipses
}  // namespace geom

namespace table {

/**
 *  A FunctorKey used to get or set a geom::Point from an (x,y) pair of int or double Keys.
 */
template <typename T>
class PointKey : public FunctorKey<lsst::afw::geom::Point<T, 2> > {
public:
    /**
     *  Add a pair of _x, _y fields to a Schema, and return a PointKey that points to them.
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_x", "_y", will be appended to this
     *                         to form the full field names.
     *  @param[in]     doc     String used as the documentation for the fields.
     *  @param[in]     unit    String used as the unit for all fields.
     */
    static PointKey addFields(Schema& schema, std::string const& name, std::string const& doc,
                              std::string const& unit);

    /// Default constructor; instance will not be usable unless subsequently assigned to.
    PointKey() : _x(), _y() {}

    /// Construct from a pair of Keys
    PointKey(Key<T> const& x, Key<T> const& y) : _x(x), _y(y) {}

    /**
     *  Construct from a subschema, assuming x and y subfields
     *
     *  If a schema has "a_x" and "a_y" fields, this constructor allows you to construct
     *  a PointKey via:
     *
     *      PointKey<T> k(schema["a"]);
     */
    PointKey(SubSchema const& s) : _x(s["x"]), _y(s["y"]) {}

    /// Get a Point from the given record
    virtual geom::Point<T, 2> get(BaseRecord const& record) const;

    /// Set a Point in the given record
    virtual void set(BaseRecord& record, geom::Point<T, 2> const& value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying x and y Keys
    bool operator==(PointKey<T> const& other) const { return _x == other._x && _y == other._y; }
    bool operator!=(PointKey<T> const& other) const { return !(*this == other); }
    //@}

    /// Return True if both the x and y Keys are valid.
    bool isValid() const { return _x.isValid() && _y.isValid(); }

    /// Return the underlying x Key
    Key<T> getX() const { return _x; }

    /// Return the underlying y Key
    Key<T> getY() const { return _y; }

private:
    Key<T> _x;
    Key<T> _y;
};

typedef PointKey<int> Point2IKey;
typedef PointKey<double> Point2DKey;

/**
 *  A FunctorKey used to get or set celestial coordinates from a pair of Angle keys.
 *
 *  Coords are always assumed to be stored and returned in the ICRS system.
 */
class CoordKey : public FunctorKey<geom::SpherePoint> {
public:
    /**
     *  Add a pair of _ra, _dec fields to a Schema, and return a CoordKey that points to them.
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_ra", "_dec", will be appended
     *                         to this to form the full field names.
     *  @param[in]     doc     String used as the documentation for the fields.
     */
    static CoordKey addFields(afw::table::Schema& schema, std::string const& name, std::string const& doc);

    /// Default constructor; instance will not be usable unless subsequently assigned to.
    CoordKey() : _ra(), _dec() {}

    /// Construct from a pair of Keys
    CoordKey(Key<geom::Angle> const& ra, Key<geom::Angle> const& dec) : _ra(ra), _dec(dec) {}

    /**
     *  Construct from a subschema, assuming ra and dec subfields.
     *
     *  If a schema has "a_ra" and "a_dec" fields, this constructor allows you to
     *  construct a CoordKey via:
     *
     *      CoordKey k(schema["a"]);
     */
    CoordKey(SubSchema const& s) : _ra(s["ra"]), _dec(s["dec"]) {}

    /// Get an IcrsCoord from the given record
    virtual geom::SpherePoint get(BaseRecord const& record) const;

    /// Set an IcrsCoord in the given record
    virtual void set(BaseRecord& record, geom::SpherePoint const& value) const;

    //@{
    /// Compare CoordKeys for equality using the constituent `ra` and `dec` Keys
    bool operator==(CoordKey const& other) const { return _ra == other._ra && _dec == other._dec; }
    bool operator!=(CoordKey const& other) const { return !(*this == other); }
    //@}

    bool isValid() const { return _ra.isValid() && _dec.isValid(); }

    //@{
    /// Return a constituent Key
    Key<geom::Angle> getRa() const { return _ra; }
    Key<geom::Angle> getDec() const { return _dec; }
    //@}

private:
    Key<geom::Angle> _ra;
    Key<geom::Angle> _dec;
};

//@{
/// Compare CoordKeys for equality using the constituent Keys
bool operator==(CoordKey const& lhs, CoordKey const& rhs);
bool operator!=(CoordKey const& lhs, CoordKey const& rhs);
//@}

/// Enum used to set units for geometric FunctorKeys
enum class CoordinateType { PIXEL, CELESTIAL };

/**
 *  A FunctorKey used to get or set a geom::ellipses::Quadrupole from a tuple of constituent Keys.
 */
class QuadrupoleKey : public FunctorKey<lsst::afw::geom::ellipses::Quadrupole> {
public:
    /**
     *
     *  Add a set of quadrupole subfields to a schema and return a QuadrupoleKey that points to them.
     *
     *  @param[in,out] schema     Schema to add fields to.
     *  @param[in]     name       Name prefix for all fields; ("_xx", "_yy", "_xy") will be appended to
     *                            this to form the full field names. In celestial coordinates, we use "x"
     *                            as a synonym for "RA" and "y" for "dec".
     *  @param[in]     doc        String used as the documentation for the fields.
     *  @param[in]     coordType  Type of coordinates in use (PIXEL or CELESTIAL).
     */
    static QuadrupoleKey addFields(Schema& schema, std::string const& name, std::string const& doc,
                                   CoordinateType coordType = CoordinateType::PIXEL);

    /// Default constructor; instance will not be usable unless subsequently assigned to.
    QuadrupoleKey() : _ixx(), _iyy(), _ixy() {}

    /// Construct from individual Keys
    QuadrupoleKey(Key<double> const& ixx, Key<double> const& iyy, Key<double> const& ixy)
            : _ixx(ixx), _iyy(iyy), _ixy(ixy) {}

    /**
     *  Construct from a subschema with appropriate subfields
     *
     *  If the schema has "a_xx", "a_yy" and "a_xy" fields this constructor enables you to
     *  construct a QuadrupoleKey via:
     *
     *      QuadrupoleKey k(schema["a"], coordType);
     */
    QuadrupoleKey(SubSchema const& s) : _ixx(s["xx"]), _iyy(s["yy"]), _ixy(s["xy"]) {}

    /// Get a Quadrupole from the given record
    virtual geom::ellipses::Quadrupole get(BaseRecord const& record) const;

    /// Set a Quadrupole in the given record
    virtual void set(BaseRecord& record, geom::ellipses::Quadrupole const& value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Ixx, Iyy, Ixy Keys
    bool operator==(QuadrupoleKey const& other) const {
        return _ixx == other._ixx && _iyy == other._iyy && _ixy == other._ixy;
    }
    bool operator!=(QuadrupoleKey const& other) const { return !(*this == other); }
    //@}

    /// Return True if all the constituent Keys are valid.
    bool isValid() const { return _ixx.isValid() && _iyy.isValid() && _ixy.isValid(); }

    //@{
    /// Return a constituent Key
    Key<double> getIxx() const { return _ixx; }
    Key<double> getIyy() const { return _iyy; }
    Key<double> getIxy() const { return _ixy; }
    //@}

private:
    Key<double> _ixx;
    Key<double> _iyy;
    Key<double> _ixy;
};

/**
 *  A FunctorKey used to get or set a geom::ellipses::Ellipse from an (xx,yy,xy,x,y) tuple of Keys.
 */
class EllipseKey : public FunctorKey<lsst::afw::geom::ellipses::Ellipse> {
public:
    /**
     *  Add a set of _xx, _yy, _xy, _x, _y fields to a Schema, and return an EllipseKey that points to them.
     *
     *  @param[in,out] schema  Schema to add fields to.
     *  @param[in]     name    Name prefix for all fields; "_xx", "_yy", "_xy", "_x" ,"_y", will be
     *                         appended to this to form the full field names.
     *  @param[in]     doc     String used as the documentation for the fields.
     *  @param[in]     unit    String used as the unit for x and y fields; "<unit>^2" will be used for
     *                         xx, yy, and xy fields.
     */
    static EllipseKey addFields(Schema& schema, std::string const& name, std::string const& doc,
                                std::string const& unit);

    /// Default constructor; instance will not be usable unless subsequently assigned to.
    EllipseKey() : _qKey(), _pKey() {}

    /// Construct from individual Keys
    EllipseKey(QuadrupoleKey const& qKey, PointKey<double> const& pKey) : _qKey(qKey), _pKey(pKey) {}

    /**
     *  Construct from a subschema, assuming (xx, yy, xy, x, y) subfields
     *
     *  If a schema has "a_xx", "a_yy", "a_xy", "a_x", and "a_y" fields, this constructor allows you to
     *  construct an EllipseKey via:
     *
     *      EllipseKey k(schema["a"]);
     */
    EllipseKey(SubSchema const& s) : _qKey(s), _pKey(s) {}

    /// Get an Ellipse from the given record
    virtual geom::ellipses::Ellipse get(BaseRecord const& record) const;

    /// Set an Ellipse in the given record
    virtual void set(BaseRecord& record, geom::ellipses::Ellipse const& value) const;

    //@{
    /// Compare the FunctorKey for equality with another, using the underlying Ixx, Iyy, Ixy Keys
    bool operator==(EllipseKey const& other) const { return _qKey == other._qKey && _pKey == other._pKey; }
    bool operator!=(EllipseKey const& other) const { return !(*this == other); }
    //@}

    /// Return True if all the constituent Keys are valid.
    bool isValid() const { return _qKey.isValid() && _pKey.isValid(); }

    //@{
    /// Return constituent FunctorKeys
    QuadrupoleKey getCore() const { return _qKey; }
    PointKey<double> getCenter() const { return _pKey; }
    //@}

private:
    QuadrupoleKey _qKey;
    PointKey<double> _pKey;
};

template <typename T, int N>
class CovarianceMatrixKey : public FunctorKey<Eigen::Matrix<T, N, N> > {
public:
    typedef std::vector<Key<T> > SigmaKeyArray;
    typedef std::vector<Key<T> > CovarianceKeyArray;
    typedef std::vector<std::string> NameArray;

    /**
     *  Add covariance matrix fields to a Schema, and return a CovarianceMatrixKey to manage them.
     *
     *  @param[out] schema    Schema to add fields to.
     *  @param[in]  prefix    String used to form the first part of all field names.  Suffixes of
     *                        the form '_xSigma' and '_x_y_Cov' will be added to form the full
     *                        field names (using names={'x', 'y'} as an example).
     *  @param[in]  unit      Unit for for sigma values (covariance matrix elements will be
     *                        unit^2).
     *  @param[in]  names     Vector of strings containing the names of the quantities the
     *                        covariance matrix represents the uncertainty of.
     *  @param[in]  diagonalOnly   If true, only create fields for the Sigma values.
     */
    static CovarianceMatrixKey addFields(Schema& schema, std::string const& prefix, NameArray const& names,
                                         std::string const& unit, bool diagonalOnly = false);

    /**
     *  Add covariance matrix fields to a Schema, and return a CovarianceMatrixKey to manage them.
     *
     *  @param[out] schema    Schema to add fields to.
     *  @param[in]  prefix    String used to form the first part of all field names.  Suffixes of
     *                        the form '_xSigma' and '_x_y_Cov' will be added to form the full
     *                        field names (using names={'x', 'y'} as an example).
     *  @param[in]  units     Vector of units for for sigma values (covariance matrix elements will
     *                        have "{units[i]} {units[j]}" or "{units[i]}^2", depending on whether
     *                        units[i] == units[j].
     *  @param[in]  names     Vector of strings containing the names of the quantities the
     *                        covariance matrix represents the uncertainty of.
     *  @param[in]  diagonalOnly   If true, only create fields for the Sigma values.
     */
    static CovarianceMatrixKey addFields(Schema& schema, std::string const& prefix, NameArray const& names,
                                         NameArray const& units, bool diagonalOnly = false);

    /// Construct an invalid instance; must assign before subsequent use.
    CovarianceMatrixKey();

    /**
     *  Construct a from arrays of per-element Keys
     *
     *  The sigma array Keys should point to the square root of the diagonal of the
     *  covariance matrix.  The cov array Keys should point to the off-diagonal elements
     *  of the lower-triangle, packed first in rows, then in columns (or equivalently,
     *  in the upper-triangle, packed first in columns, then in rows).  For a 4x4 matrix,
     *  the order is is:
     *
     *      sigma[0]^2   cov[0]       cov[1]       cov[3]
     *      cov[0]       sigma[1]^2   cov[2]       cov[4]
     *      cov[1]       cov[2]       sigma[2]^2   cov[5]
     *      cov[3]       cov[4]       cov[5]       sigma[3]^2
     *
     *  The cov array may also be empty, to indicate that no off-diagonal elements are
     *  stored, and should be set to zero.  If not empty, the size of the cov matrix
     *  must be exactly n*(n-1)/2, where n is the size of the sigma matrix.
     */
    explicit CovarianceMatrixKey(SigmaKeyArray const& sigma,
                                 CovarianceKeyArray const& cov = CovarianceKeyArray());

    /**
     *  Construct from a subschema and an array of names for each parameter of the matrix.
     *
     *  The field names should match the following convention:
     *   - diagonal elements should have names like "p1Sigma", where "p1" is the name of the parameter,
     *     and should contain the square root of the variance in that parameter.
     *   - off-diagonal elements hould have names like "p1_p2_Cov", where "p1" and "p2" are names of
     *     parameters.
     *  For example, for the covariance matrix of a position, we'd look for "xSigma", "ySigma", and
     *  "x_y_Cov".
     */
    CovarianceMatrixKey(SubSchema const& s, NameArray const& names);

    /// Get a covariance matrix from the given record
    virtual Eigen::Matrix<T, N, N> get(BaseRecord const& record) const;

    /// Set a covariance matrix in the given record (uses only the lower triangle of the given matrix)
    virtual void set(BaseRecord& record, Eigen::Matrix<T, N, N> const& value) const;

    /// Return the element in row i and column j
    T getElement(BaseRecord const& record, int i, int j) const;

    /// Set the element in row i and column j
    void setElement(BaseRecord& record, int i, int j, T value) const;

    /**
     *  Return True if all the constituent sigma Keys are valid
     *
     *  Note that if the only one or more off-diagonal keys are invalid, we assume that means those terms
     *  are zero, not that the whole FunctorKey is invalid.
     */
    bool isValid() const;

    //@{
    /// Compare the FunctorKey for equality with another, using its constituent Keys
    bool operator==(CovarianceMatrixKey const& other) const;
    bool operator!=(CovarianceMatrixKey const& other) const { return !(*this == other); }
    //@}

private:
    SigmaKeyArray _sigma;
    CovarianceKeyArray _cov;
};
}
}
}  // namespace lsst::afw::table

#endif  // !AFW_TABLE_aggregates_h_INCLUDED
