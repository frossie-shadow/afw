/*
 * LSST Data Management System
 * Copyright 2014 LSST Corporation.
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

#include <exception>
#include <memory>
#include <sstream>
#include <type_traits>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/cameraGeom/TransformMap.h"

namespace lsst {
namespace afw {
namespace cameraGeom {

namespace {

/**
 * Represent a set of camera Transforms as a FrameSet.
 *
 * @tparam Map Any type satisfying the STL map API and mapping CameraSys to
 *             shared_ptr<Transform>.
 *
 * @param root Common coordinate system serving as root of FrameSet's tree.
 * @param transforms A map whose keys are coordinate systems, and whose
 *                   values point to Transforms that convert from `root` to the
 *                   corresponding key. All Transforms must be invertible.
 * @return a newly allocated FrameSet corresponding to the method arguments
 *
 * @throws lsst::pex::exceptions::InvalidParameterError Thrown if
 *         `transforms` contains `root`, or if any Transform is not
 *         invertible.
 */
template <class Map>
std::unique_ptr<ast::FrameSet> makeTransforms(CameraSys const &root, Map const &transforms) {
    ast::Frame rootFrame(2, "Ident=" + root.getSysName());
    auto result = std::unique_ptr<ast::FrameSet>(new ast::FrameSet(rootFrame));

    for (auto keyValue : transforms) {
        CameraSys const key = keyValue.first;
        std::shared_ptr<TransformMap::Transform> const value = keyValue.second;

        if (key == root) {
            std::ostringstream buffer;
            buffer << "Cannot specify a Transform from " << root << " to itself.";
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, buffer.str());
        }
        if (!value->hasForward()) {
            std::ostringstream buffer;
            buffer << *value << " from " << key << " has no forward transform.";
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, buffer.str());
        }
        if (!value->hasInverse()) {
            std::ostringstream buffer;
            buffer << *value << " from " << key << " has no inverse transform.";
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, buffer.str());
        }

        auto toFrame = value->getFrameSet()->getFrame(ast::FrameSet::CURRENT);
        toFrame->setIdent(key.getSysName());
        result->addFrame(ast::FrameSet::BASE, *(value->getFrameSet()->getMapping()), *toFrame);
    }
    return result;
}

/**
 * Identify the frames in to a FrameSet constructed by makeTransforms.
 *
 * @tparam Map Any type satisfying the STL map API and mapping CameraSys to
 *             shared_ptr<Transform>.
 *
 * @param reference Coordinate system to which `transforms` converts.
 * @param transforms A map whose keys are coordinate systems, and whose
 *                   values point to Transforms that convert from `reference`
 *                   to the corresponding key. All Transforms must be
 *                   invertible.
 * @return a map from `reference` and each key in `transforms` to the
 *         corresponding frame number in `makeTransforms(reference, transforms)`.
 *
 * @warning Does not perform input validation.
 */
template <class Map>
std::unordered_map<CameraSys, int> makeTranslator(CameraSys const &reference, Map const &transforms) {
    std::unordered_map<CameraSys, int> result({std::make_pair(reference, 1)});
    int nFrames = 1;

    for (auto keyValue : transforms) {
        CameraSys const key = keyValue.first;
        result.emplace(key, ++nFrames);
    }
    return result;
}

}  // namespace

geom::Point2Endpoint TransformMap::_pointConverter;

TransformMap::TransformMap(CameraSys const &reference,
                           std::map<CameraSys, std::shared_ptr<Transform>> const &transforms)
        : _transforms(makeTransforms(reference, transforms)),
          _frameIds(makeTranslator(reference, transforms)) {}

TransformMap::TransformMap(CameraSys const &reference,
                           std::unordered_map<CameraSys, std::shared_ptr<Transform>> const &transforms)
        : _transforms(makeTransforms(reference, transforms)),
          _frameIds(makeTranslator(reference, transforms)) {}

// TransformMap is immutable, so we can just copy the shared_ptr
TransformMap::TransformMap(TransformMap const &other) = default;

// Cannot do any move optimizations without breaking immutability
TransformMap::TransformMap(TransformMap const &&other) : TransformMap(other) {}

geom::Point2D TransformMap::transform(geom::Point2D const &point, CameraSys const &fromCoordSys,
                                      CameraSys const &toCoordSys) const {
    auto mapping = _getMapping(fromCoordSys, toCoordSys);
    return _pointConverter.pointFromData(mapping->applyForward(_pointConverter.dataFromPoint(point)));
}

std::vector<geom::Point2D> TransformMap::transform(std::vector<geom::Point2D> const &pointList,
                                                   CameraSys const &fromCoordSys,
                                                   CameraSys const &toCoordSys) const {
    auto mapping = _getMapping(fromCoordSys, toCoordSys);
    return _pointConverter.arrayFromData(mapping->applyForward(_pointConverter.dataFromArray(pointList)));
}

bool TransformMap::contains(CameraSys const &system) const noexcept { return _frameIds.count(system) > 0; }

std::shared_ptr<TransformMap::Transform> TransformMap::getTransform(CameraSys const &fromCoordSys,
                                                                    CameraSys const &toCoordSys) const {
    return std::make_shared<Transform>(*_getMapping(fromCoordSys, toCoordSys));
}

int TransformMap::_getFrame(CameraSys const &system) const {
    try {
        return _frameIds.at(system);
    } catch (std::out_of_range const &e) {
        std::ostringstream buffer;
        buffer << "Unsupported coordinate system: " << system;
        std::throw_with_nested(LSST_EXCEPT(pex::exceptions::InvalidParameterError, buffer.str()));
    }
}

std::shared_ptr<ast::Mapping const> TransformMap::_getMapping(CameraSys const &fromCoordSys,
                                                              CameraSys const &toCoordSys) const {
    return _transforms->getMapping(_getFrame(fromCoordSys), _getFrame(toCoordSys));
}

size_t TransformMap::size() const noexcept { return _frameIds.size(); }
}  // namespace cameraGeom
}  // namespace afw
}  // namespace lsst
