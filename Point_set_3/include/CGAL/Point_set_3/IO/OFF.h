// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_IO_OFF_H
#define CGAL_POINT_SET_IO_OFF_H

#include <CGAL/license/Point_set_3.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>

#include <fstream>
#include <string>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

template <typename Point, typename Vector>
class Point_set_3;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {

/*!
  \ingroup PkgPointSet3IOOFF

  \brief reads the content of an intput stream in the \ref IOStreamOFF into a point set.

  If normals are present in the input (NOFF), a normal map will be created and filled.

  \param is the input stream
  \param point_set the point set

  \return `true` if the reading was successful, `false` otherwise.
 */
template <typename Point, typename Vector>
bool read_OFF(std::istream& is,
              CGAL::Point_set_3<Point, Vector>& point_set)
{
  point_set.add_normal_map();

  bool out = CGAL::IO::read_OFF(is, point_set.index_back_inserter(),
                                CGAL::parameters::point_map(point_set.point_push_map())
                                                 .normal_map(point_set.normal_push_map()));

  bool has_normals = false;
  for(typename CGAL::Point_set_3<Point, Vector>::const_iterator it=point_set.begin(); it!=point_set.end(); ++it)
  {
    if(point_set.normal(*it) != CGAL::NULL_VECTOR)
    {
      has_normals = true;
      break;
    }
  }

  if(!has_normals)
    point_set.remove_normal_map();

  return out;
}

/*!
  \ingroup PkgPointSet3IOOFF

  \brief reads the content of an input file in the \ref IOStreamOFF into a point set.

  If normals are present in the input (NOFF), a normal map will be created and filled.

  \param fname the path to the input file
  \param point_set the point set

  \return `true` if the reading was successful, `false` otherwise.
*/
template <typename Point, typename Vector>
bool read_OFF(const std::string& fname, CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ifstream is(fname);
  return read_OFF(is, point_set);
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3,
              \link PkgPointSet3IO `CGAL::IO::read_OFF()` \endlink  should be used instead.
 */
template <typename Point, typename Vector>
CGAL_DEPRECATED bool read_off_point_set(std::istream& is,  ///< input stream.
                                        CGAL::Point_set_3<Point, Vector>& point_set)  ///< point set.
{
  return IO::read_OFF(is, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

/*!
  \ingroup PkgPointSet3IOOFF

  \brief writes the content of a point set into an output stream in the \ref IOStreamOFF.

  If it exists, the normal map associated to the point set is inserted in the stream.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param point_set the point set
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`the precision of the stream `os``}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the writing was successful, `false` otherwise.
 */
template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const CGAL::Point_set_3<Point, Vector>& point_set,
               const CGAL_BGL_NP_CLASS& np)
{
  if(point_set.has_normal_map())
    return Point_set_processing_3::internal::write_OFF_PSP(os, point_set,
                                                           np.point_map(point_set.point_map())
                                                             .normal_map(point_set.normal_map()));

  return Point_set_processing_3::internal::write_OFF_PSP(os, point_set,
                                                         np.point_map(point_set.point_map()));
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector>
bool write_OFF(std::ostream& os, const CGAL::Point_set_3<Point, Vector>& point_set)
{
  return write_OFF(os, point_set, parameters::all_default());
}

/// \endcond

/*!
  \ingroup PkgPointSet3IOOFF

  \brief writes the content of a point set into an output file in the \ref IOStreamOFF.

  If it exists, the normal map associated to the point set is inserted in the file.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the path to the output file
  \param point_set the point set
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the writing was successful, `false` otherwise.
*/
template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& fname, const CGAL::Point_set_3<Point, Vector>& point_set, const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname);
  return write_OFF(os, point_set, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector>
bool write_OFF(const std::string& fname, const CGAL::Point_set_3<Point, Vector>& point_set)
{
  std::ofstream os(fname);
  return write_OFF(os, point_set, parameters::all_default());
}

/// \endcond

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgPointSet3IODeprecated

  \deprecated This function is deprecated since \cgal 5.3,
              \link PkgPointSet3IO `CGAL::IO::write_OFF()` \endlink  should be used instead.
 */
template <typename Point, typename Vector>
CGAL_DEPRECATED bool write_off_point_set(std::ostream& os, ///< output stream.
                                         const CGAL::Point_set_3<Point, Vector>& point_set) ///< point set
{
  return IO::write_OFF(os, point_set);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_POINT_SET_IO_OFF_H
