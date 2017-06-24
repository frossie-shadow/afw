// -*- LSST-C++ -*-

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

#ifndef LSST_AFW_IMAGE_DETAILS_WCSUTILS_H
#define LSST_AFW_IMAGE_DETAILS_WCSUTILS_H

#include <memory>
#include <string>

#include "lsst/daf/base.h"

namespace lsst {
namespace afw {
namespace image {
namespace detail {

std::shared_ptr<lsst::daf::base::PropertyList> createTrivialWcsAsPropertySet(std::string const& wcsName,
                                                                             int const x0 = 0,
                                                                             int const y0 = 0);

geom::Point2I getImageXY0FromMetadata(std::string const& wcsName, lsst::daf::base::PropertySet* metadata);

}  // namespace detail
}  // namespace image
}  // namespace afw
}  // namespace lsst

#endif  // LSST_AFW_IMAGE_DETAILS_WCSUTILS_H
