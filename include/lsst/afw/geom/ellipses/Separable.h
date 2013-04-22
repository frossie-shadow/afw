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

#ifndef LSST_AFW_GEOM_ELLIPSES_Separable_h_INCLUDED
#define LSST_AFW_GEOM_ELLIPSES_Separable_h_INCLUDED

/**
 *  \file
 *  @brief Definitions and inlines for Separable.
 *
 *  \note Do not include directly; use the main ellipse header file.
 */

#include "lsst/afw/geom/ellipses/EllipseCore.h"
#include "lsst/afw/geom/ellipses/Convolution.h"
#include "lsst/afw/geom/ellipses/Transformer.h"
#include "lsst/afw/geom/ellipses/GridTransform.h"

namespace lsst { namespace afw { namespace geom { namespace ellipses {

/**
 *  @brief An ellipse core with a complex ellipticity and radius parameterization.
 *
 *  
 */
template <typename Ellipticity_, typename Radius_>
class Separable : public EllipseCore {
public:

    enum ParameterEnum { E1=0, E2=1, RADIUS=2 }; ///< Definitions for elements of a core vector.

    typedef Ellipticity_ Ellipticity;
    typedef Radius_ Radius;

    double const getE1() const { return _ellipticity.getE1(); }
    void setE1(double e1) { _ellipticity.setE1(e1); }

    double const getE2() const { return _ellipticity.getE2(); }
    void setE2(double e2) { _ellipticity.setE2(e2); }

    Radius const & getRadius() const { return _radius; }
    Radius & getRadius() { return _radius; }
    void setRadius(double radius) { _radius = radius; }
    void setRadius(Radius const & radius) { _radius = radius; }

    Ellipticity const & getEllipticity() const { return _ellipticity; }
    Ellipticity & getEllipticity() { return _ellipticity; }

    /// @brief Deep copy the ellipse core.
    PTR(Separable) clone() const { return boost::static_pointer_cast<Separable>(_clone()); }

    /// Return a string that identifies this parametrization.
    virtual std::string getName() const;

    /**
     *  @brief Put the parameters into a "standard form", and throw InvalidEllipseParameters
     *         if they cannot be normalized.
     */
    virtual void normalize();

    /// @brief Standard assignment.
    Separable & operator=(Separable const & other);

    /// @brief Converting assignment.
    Separable & operator=(EllipseCore const & other) { EllipseCore::operator=(other); return *this; }

    /// @brief Construct a circle with the given radius.
    explicit Separable(double radius=1.0);

    /// @brief Construct from parameter values.
    Separable(double e1, double e2, double radius=Radius(), bool normalize=true);

    /// @brief Construct from parameter values.
    explicit Separable(std::complex<double> const & complex, 
                       double radius=Radius(), bool normalize=true);

    /// @brief Construct from parameter values.
    explicit Separable(Ellipticity const & ellipticity, double radius=Radius(), bool normalize=true);

    /// @brief Construct from a parameter vector.
    explicit Separable(EllipseCore::ParameterVector const & vector, bool normalize=false);

    /// @brief Copy constructor.
    Separable(Separable const & other) : _ellipticity(other._ellipticity), _radius(other._radius) {}

    /// @brief Converting copy constructor.
    Separable(EllipseCore const & other) { *this = other; }

#ifndef SWIG
    /// @brief Converting copy constructor.
    Separable(EllipseCore::Transformer const & transformer) {
        transformer.apply(*this);
    }

    /// @brief Converting copy constructor.
    Separable(EllipseCore::Convolution const & convolution) {
        convolution.apply(*this);
    }
#endif
protected:

    virtual PTR(EllipseCore) _clone() const { return boost::make_shared<Separable>(*this); }

    virtual void readParameters(double const * iter);
    virtual void writeParameters(double * iter) const;

    virtual void _assignToQuadrupole(double & ixx, double & iyy, double & ixy) const;
    virtual void _assignFromQuadrupole(double ixx, double iyy, double ixy);

    virtual void _assignToAxes(double & a, double & b, double & theta) const;
    virtual void _assignFromAxes(double a, double b, double theta);

    virtual Jacobian _dAssignToQuadrupole(double & ixx, double & iyy, double & ixy) const;
    virtual Jacobian _dAssignFromQuadrupole(double ixx, double iyy, double ixy);

    virtual Jacobian _dAssignToAxes(double & a, double & b, double & theta) const;
    virtual Jacobian _dAssignFromAxes(double a, double b, double theta);

private:

    static EllipseCore::Registrar<Separable> registrar;

    Ellipticity _ellipticity;
    Radius _radius;
};

}}}} // namespace lsst::afw::geom::ellipses

#endif // !LSST_AFW_GEOM_ELLIPSES_Separable_h_INCLUDED
