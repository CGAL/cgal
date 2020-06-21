// Copyright (c) 2016 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_IO
#define CGAL_POINT_SET_3_IO

#include <CGAL/license/Point_set_3.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO/LAS.h>
#include <CGAL/Point_set_3/IO/OFF.h>
#include <CGAL/Point_set_3/IO/PLY.h>
#include <CGAL/Point_set_3/IO/XYZ.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <fstream>
#include <string>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
  \ingroup PkgPointSet3IO

  \brief Reads the point set from an input stream that can be either:

  - \link IOStreamXYZ XYZ \endlink
  - \link IOStreamOFF OFF \endlink
  - \link IOStreamPLY PLY \endlink
  - \link IOStreamLAS LAS \endlink

  The format is detected from the stream. If the stream contains
  normal vectors, the normal map is added to the point set. For PLY
  input, all point properties found in the header are added.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param is the input stream
  \param ps the point set

  \return `is`

  \relates Point_set_3
*/
template <typename Point, typename Vector>
std::istream& operator>>(std::istream& is,
                         CGAL::Point_set_3<Point, Vector>& ps)
{
  // Check format identifier on first line
  std::string line;
  if(!getline(is, line))
    return is;

  is.seekg(0);
  if(line.find("OFF") == 0 || line.find("NOFF") == 0)
    CGAL::read_OFF(is, ps);
  else if(line.find("ply") == 0)
    CGAL::read_PLY(is, ps);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(line.find("LASF") == 0)
    CGAL::read_LAS(is, ps);
#endif // LAS
  else
    CGAL::read_XYZ(is, ps);

  return is;
}

/*!
  \ingroup PkgPointSet3IO

  \brief Reads the point set from an input file that can be either:

  - \link IOStreamXYZ XYZ \endlink
  - \link IOStreamOFF OFF \endlink
  - \link IOStreamPLY PLY \endlink
  - \link IOStreamLAS LAS \endlink

  The format is detected from the filename extension. If the file contains
  normal vectors, the normal map is added to the point set. For PLY
  input, all point properties found in the header are added.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param filename the path to the input file
  \param ps the point set

  \return `true` if the reading was successful, `false` otherwise.
 */
template <typename Point,
          typename Vector>
bool read_point_set(const std::string& filename,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  const std::string ext = IO::internal::get_file_extension(filename);

  if(ext == "xyz")
    return read_XYZ(filename, ps);
  else if(ext == "off")
    return read_OFF(filename, ps);
  else if(ext =="ply")
    return read_PLY(filename, ps);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return read_LAS(filename, ps);
#endif

  return false;
}

template <typename Point,
          typename Vector>
bool read_point_set(const char* filename,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  return read_point_set(std::string(filename), ps);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
  \ingroup PkgPointSet3IO

  \brief Inserts the point set in an output stream in ASCII PLY format.
         All properties are inserted in their instantiation order.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`

  \param os the output stream
  \param ps the point set

  \return `os`

  \relates Point_set_3
*/
template <typename Point, typename Vector>
std::ostream& operator<<(std::ostream& os,
                         const CGAL::Point_set_3<Point, Vector>& ps)
{
  write_ply_point_set(os, ps);
  return os;
}

/*!
  \ingroup PkgPointSet3IO

  \brief Inserts the point set in an output file that can be either:

  - \link IOStreamXYZ XYZ \endlink
  - \link IOStreamOFF OFF \endlink
  - \link IOStreamPLY PLY \endlink
  - \link IOStreamLAS LAS \endlink

  The format is detected from the filename extension.

  \tparam Point a `CGAL::Point_3`
  \tparam Vector a `CGAL::Vector_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param filename the path to the output file
  \param ps the point set
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
bool write_point_set(const std::string& filename,
                     CGAL::Point_set_3<Point, Vector>& ps,
                     const CGAL_BGL_NP_CLASS& np)
{
  const std::string ext = IO::internal::get_file_extension(filename);

  if(ext == "xyz")
    return write_XYZ(filename, ps, np);
  else if(ext == "off")
    return write_OFF(filename, ps, np);
  else if(ext == "ply")
    return write_PLY(filename, ps, np);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return write_LAS(filename, ps, np);
#endif

  return false;
}

template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_point_set(const char* filename,
                     CGAL::Point_set_3<Point, Vector>& ps,
                     const CGAL_BGL_NP_CLASS& np)
{
  return write_point_set(std::string(filename), ps, np);
}

template <typename Point, typename Vector>
bool write_point_set(const char* filename, CGAL::Point_set_3<Point, Vector>& ps)
{
  return write_point_set(filename, ps, parameters::all_default());
}

template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_point_set(const std::string& filename, CGAL::Point_set_3<Point, Vector>& ps, const CGAL_BGL_NP_CLASS& np)
{
  return write_point_set(filename.c_str(), ps, np);
}

template <typename Point, typename Vector>
bool write_point_set(const std::string& filename, CGAL::Point_set_3<Point, Vector>& ps)
{
  return write_point_set(filename.c_str(), ps, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_POINT_SET_3_IO
