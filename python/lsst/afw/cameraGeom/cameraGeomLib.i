// -*- lsst-++ -*-

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
 
#define CAMERA_GEOM_LIB_I

%define cameraGeomLib_DOCSTRING
"
Python bindings for classes describing the the geometry of a mosaic camera
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.afw.cameraGeom", docstring=cameraGeomLib_DOCSTRING) cameraGeomLib

%{
#include "lsst/afw/geom/ellipses.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/afw/utils.i" 
#if defined(IMPORT_IMAGE_I)
%import  "lsst/afw/image/imageLib.i"
%import  "lsst/afw/geom/ellipses/ellipsesLib.i"
#endif

%lsst_exceptions();

%include "lsst/afw/cameraGeom/cameraGeomPtrs.i"

%template(AmpSet) std::vector<boost::shared_ptr<lsst::afw::cameraGeom::Amp> >;

%include "lsst/afw/cameraGeom/Id.i"
%include "lsst/afw/cameraGeom/FpPoint.i"
%include "lsst/afw/cameraGeom/Orientation.i"
%include "lsst/afw/cameraGeom/Detector.i"
%include "lsst/afw/cameraGeom/Amp.i"
%include "lsst/afw/cameraGeom/DetectorMosaic.i"
%include "lsst/afw/cameraGeom/Ccd.h"
%include "lsst/afw/cameraGeom/Raft.h"
%include "lsst/afw/cameraGeom/Camera.h"

%include "lsst/afw/cameraGeom/Distortion.i"

%inline %{
    lsst::afw::cameraGeom::Ccd::Ptr
    cast_Ccd(lsst::afw::cameraGeom::Detector::Ptr detector) {
        return boost::shared_dynamic_cast<lsst::afw::cameraGeom::Ccd>(detector);
    }

    lsst::afw::cameraGeom::Raft::Ptr
    cast_Raft(lsst::afw::cameraGeom::Detector::Ptr detector) {
        return boost::shared_dynamic_cast<lsst::afw::cameraGeom::Raft>(detector);
    }

%}

%definePythonIterator(lsst::afw::cameraGeom::Ccd);
