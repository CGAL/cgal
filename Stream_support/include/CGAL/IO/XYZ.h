// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_XYZ_H
#define CGAL_IO_XYZ_H


#include <fstream>
#include <iostream>
#include <vector>
#include <type_traits>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {

/**
   \ingroup PkgStreamSupportIoFuncsXYZ

   \brief reads points (positions + normals, if available), using the \ref IOStreamXYZ.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<OutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param is input stream.
   \param output output iterator over points.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa \ref IOStreamXYZ
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(std::istream& is,
              OutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values());




/**
   \ingroup PkgStreamSupportIoFuncsXYZ

   \brief reads points (positions + normals, if available), using the \ref IOStreamXYZ.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It must be a model of `DefaultConstructible` and defaults to `value_type_traits<OutputIterator>::%type`.
   It can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.
   \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

   \param fname input file name.
   \param output output iterator over points.
   \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{normal_map}
       \cgalParamDescription{a property map associating normals to the elements of the point range}
       \cgalParamType{a model of `WritablePropertyMap` with value type `geom_traits::Vector_3`}
       \cgalParamDefault{If this parameter is omitted, normals in the input stream are ignored.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \returns `true` if reading was successful, `false` otherwise.

   \sa \ref IOStreamXYZ
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
           typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_XYZ(const std::string& fname,
              OutputIterator output,
              const CGAL_NP_CLASS& np = parameters::default_values());
} // namespace IO
} // namespace CGAL


#include <CGAL/IO/XYZ/read_xyz_points.h>
#include <CGAL/IO/XYZ/write_xyz_points.h>

#endif // CGAL_IO_XYZ_H
