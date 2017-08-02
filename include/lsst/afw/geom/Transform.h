// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2017 LSST Corporation.
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

#ifndef LSST_AFW_GEOM_TRANSFORM_H
#define LSST_AFW_GEOM_TRANSFORM_H

#include <memory>
#include <vector>

#include "astshim.h"
#include "Eigen/Core"
#include "ndarray.h"

#include "lsst/afw/geom/Endpoint.h"

namespace lsst {
namespace afw {
namespace geom {

/**
Transform LSST spatial data, such as Point2D and SpherePoint, using an AST transform.

This class contains two Endpoints, to specify the "from" and "to" LSST data type,
and an ast::FrameSet or ast::Mapping to specify the transformation.
In the case of a FrameSet the transformation is from the `BASE` frame to the `CURRENT` frame.
The endpoints convert the data between the LSST Form (e.g. Point2D) and the form used by astshim.

Depending on the ast::FrameSet or ast::Mapping used to define it, a Transform may
provide either a forward transform, an inverse transform, or both. In particular, the
@ref getInverse "inverse" of a forward-only transform is an inverse-only transform. The
@ref hasForward and @ref hasInverse methods can be used to check which transforms are available.

Unless otherwise stated, all constructors and methods may throw `std::runtime_error` to indicate
internal errors within AST.

@note You gain some safety by constructing a Transform from an ast::FrameSet,
since the base and current frames in the FrameSet can be checked against by the appropriate endpoint.

@note "In place" versions of `applyForward` and `applyInverse` are not available
because data must be copied when converting from LSST data types to the type used by astshim,
so it didn't seem worth the bother.
*/
template <class FromEndpoint, class ToEndpoint>
class Transform {
public:
    using FromArray = typename FromEndpoint::Array;
    using FromPoint = typename FromEndpoint::Point;
    using ToArray = typename ToEndpoint::Array;
    using ToPoint = typename ToEndpoint::Point;

    Transform(Transform const &) = default;
    Transform(Transform &&) = default;
    Transform &operator=(Transform const &) = delete;
    Transform &operator=(Transform &&) = delete;

    /**
    Construct a Transform from a deep copy of an ast::Mapping

    The internal FrameSet consists of a frame constructed by each endpoint
    connected by the mapping.

    @param[in] mapping  ast::Mapping describing the desired transformation
    @param[in] simplify  Simplify the mapping? This combines component mappings
        and removes redundant components where possible.
    */
    explicit Transform(ast::Mapping const &mapping, bool simplify = true);

    /**
    Construct a Transform from a deep copy of a FrameSet.

    The result transforms from the "base" frame to the "current" frame.
    The "from" endpoint is used to normalize the "base" frame
    and the "to" endpoint is used to normalize the "current" frame.

    This is pickier than the constructor that takes an ast::Mapping in that:
    - SphereEndpoint must be associated with an ast::SkyFrame and the SkyFrame axes
      are swapped if necessary to the standard order: longitude, latitude.
    - Point2Endpoint must be associated with an ast::Frame (not a subclass),
      because Frame is the only kind of Frame that is sure to be Cartesian.

    @param[in] frameSet  ast::FrameSet describing the desired transformation in the usual way:
                         from "base" frame to "current" frame
    @param[in] simplify  Simplify the frame set? This simplifies each mapping
                         in the frame set by combining component mappings and removing
                         redundant components where possible. However it
                         does not remove any frames.
    */
    explicit Transform(ast::FrameSet const &frameSet, bool simplify = true);

    virtual ~Transform(){};

    /**
     * Test if this method has a forward transform.
     *
     * @exceptsafe Provides basic exception safety.
     */
    bool hasForward() const { return _frameSet->hasForward(); }

    /**
     * Test if this method has an inverse transform.
     *
     * @exceptsafe Provides basic exception safety.
     */
    bool hasInverse() const { return _frameSet->hasInverse(); }

    /**
    Get the "from" endpoint
    */
    FromEndpoint getFromEndpoint() const { return _fromEndpoint; }

    /**
    Get the contained frameset
    */
    std::shared_ptr<const ast::FrameSet> getFrameSet() const { return _frameSet; }

    /**
    Get the "to" endpoint
    */
    ToEndpoint getToEndpoint() const { return _toEndpoint; }

    /**
    Transform one point in the forward direction ("from" to "to")
    */
    ToPoint applyForward(FromPoint const &point) const;

    /**
    Transform an array of points in the forward direction ("from" to "to")

    The first dimension of the array must match the number of input axes, and the data order is
    values for the first axis, then values for the next axis, and so on, e.g. for 2 axes:
        x0, x1, x2, ..., y0, y1, y2...
    */
    ToArray applyForward(FromArray const &array) const;

    /**
    Transform one point in the inverse direction ("to" to "from")
    */
    FromPoint applyInverse(ToPoint const &point) const;

    /**
    Transform an array of points in the inverse direction ("to" to "from")

    The first dimension of the array must match the number of output axes, and the data order is
    values for the first axis, then values for the next axis, and so on, e.g. for 2 axes:
        x0, x1, x2, ..., y0, y1, y2...
    */
    FromArray applyInverse(ToArray const &array) const;

    /**
     * The inverse of this Transform.
     *
     * @returns a Transform whose `applyForward` is equivalent to this Transform's
     *          `applyInverse`, and vice versa.
     *
     * @exceptsafe Provides basic exception safety.
     */
    Transform<ToEndpoint, FromEndpoint> getInverse() const;

    /**
     * The Jacobian matrix of this Transform.
     *
     * The matrix is defined only if this object has a forward transform.
     *
     * @param x the position at which the Jacobian shall be evaluated
     * @returns a matrix `J` with `getToEndpoint().getNAxes()` rows and
     *          `getFromEndpoint().getNAxes()` columns. `J(i,j)` shall be the
     *          rate of change of the `i`th output coordinate with respect to
     *          the `j`th input coordinate, or `NaN` if the derivative cannot
     *          be calculated.
     *
     * @exceptsafe Provides basic exception safety.
     *
     * @note The derivatives may be estimated by sampling and interpolating
     *       this Transform in the neighborhood of `x`. If the implementation
     *       requires interpolation, computation of the Jacobian may require
     *       hundreds of evaluations of @ref applyForward.
     */
    Eigen::MatrixXd getJacobian(FromPoint const &x) const;

    /**
     * Concatenate two Transforms.
     *
     * @tparam NextToEndpoint the "to" Endpoint of `next`
     * @param next the Transform to apply after this one
     * @returns a Transform that first applies this transform to its input, and then
     *          `next` to the result. Its inverse shall first apply the
     *          inverse of `next`, and then the inverse of this transform.
     *
     * @throws pex::exceptions::InvalidParameterError Thrown if
     *         `getToEndpoint()` and `next.getFromEndpoint()` do not
     *         have the same number of axes.
     * @exceptsafe Provides basic exception safety.
     *
     * More than two Transforms can be combined in series. For example:
     *
     *     auto pixelsToSky = pixelsToFP.then(fpToPupil).then(pupilToSky);
     */
    template <class NextToEndpoint>
    Transform<FromEndpoint, NextToEndpoint> then(Transform<ToEndpoint, NextToEndpoint> const &next) const;

    /**
     * Return a short version of the class name with no punctuation
     *
     * Used as the Python class name and for persistence as a string
     *
     * Returns "Transform" + fromEndpoint.getClassPrefix() + "To" + toEndpoint.getClassPrefix(),
     * for example "TransformPoint2ToSpherePoint" or "TransformPoint2ToGeneric".
     */
    static std::string getShortClassName();

    /**
     * Deserialize a Transform of this type from an input stream
     *
     * @param[in] is  input stream from which to deserialize this Transform
     */
    static Transform<FromEndpoint, ToEndpoint> readStream(std::istream & is);

    /// Deserialize a Transform of this type from a string, using the same format as readStream
    static Transform<FromEndpoint, ToEndpoint> readString(std::string & str);

    /**
     * Serialize this Transform to an output stream
     *
     * Version 1 format is as follows:
     * - The version number (an integer)
     * - A space
     * - The short class name, as obtained from getShortClassName
     * - A space
     * - The contained ast::FrameSet written using FrameSet.show(os, false)
     *
     * @param[out] os  outpu stream to which to serialize this Transform
     */
    virtual void writeStream(std::ostream & os) const;

    /// Serialize this Transform to a string, using the same format as writeStream
    virtual std::string writeString() const;

protected:
    /**
    Construct a Transform from a shared pointer to a FrameSet

    @note The FrameSet may be modified by normalizing the base and current frame.
    */
    explicit Transform(std::shared_ptr<ast::FrameSet> &&frameSet);

private:
    FromEndpoint _fromEndpoint;
    std::shared_ptr<const ast::FrameSet> _frameSet;
    ToEndpoint _toEndpoint;
};

/**
Print a Transform to an ostream

The format is "Transform<_fromEndpoint_, _toEndpoint_>"
where _fromEndpoint_ and _toEndpoint_ are the appropriate endpoint printed to the ostream;
for example "Transform<GenericEndpoint(4), Point2Endpoint()>"
*/
template <class FromEndpoint, class ToEndpoint>
std::ostream &operator<<(std::ostream &os, Transform<FromEndpoint, ToEndpoint> const &transform);

}  // namespace geom
}  // namespace afw
}  // namespace lsst

#endif
