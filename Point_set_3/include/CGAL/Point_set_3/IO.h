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

#include <CGAL/IO/helpers.h>
#include <CGAL/Point_set_3/IO/LAS.h>
#include <CGAL/Point_set_3/IO/OFF.h>
#include <CGAL/Point_set_3/IO/PLY.h>
#include <CGAL/Point_set_3/IO/XYZ.h>

#include <fstream>
#include <string>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {

template <typename Point, typename Vector>
class Point_set_3;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

/*!
  \ingroup PkgPointSet3IO

  \brief reads the point set from an input stream.

  Supported file formats are the following:
  - \ref IOStreamOFF (`.off`)
  - \ref IOStreamPLY (`.ply`)
  - \ref IOStreamLAS (`.las`)
  - \ref IOStreamXYZ (`.xyz`)

  The format is detected from the stream. If the stream contains
  normal vectors, the normal map is added to the point set. For PLY
  input, all point properties found in the header are added.

  \attention When reading a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ifstream`.

  \param is input stream
  \param ps point set

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
    CGAL::IO::read_OFF(is, ps);
  else if(line.find("ply") == 0)
    CGAL::IO::read_PLY(is, ps);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(line.find("LASF") == 0)
    CGAL::IO::read_LAS(is, ps);
#endif // LAS
  else
    CGAL::IO::read_XYZ(is, ps);

  return is;
}

namespace IO {

/*!
  \ingroup PkgPointSet3IO

  \brief reads the point set from an input file.

  Supported file formats are the following:
  - \ref IOStreamOFF (`.off`)
  - \ref IOStreamPLY (`.ply`)
  - \ref IOStreamLAS (`.las`)
  - \ref IOStreamXYZ (`.xyz`)

  The format is detected from the filename extension (letter case is not important).
  If the file contains normal vectors, the normal map is added to the point set.
  For PLY input, all point properties found in the header are added.

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname name of the input file
  \param ps the point set
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be read in binary (`true`) or in ASCII (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{This parameter is only relevant for `PLY` writing: the `OFF` and `XYZ` formats
                       are always ASCII, and the `LAS` format is always binary.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the reading was successful, `false` otherwise.
 */
template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_point_set(const std::string& fname,
                    CGAL::Point_set_3<Point, Vector>& ps,
                    const CGAL_BGL_NP_CLASS& np)
{
  const std::string ext = internal::get_file_extension(fname);

  if(ext == "xyz" || ext == "pwn")
    return read_XYZ(fname, ps);
  else if(ext == "off")
    return read_OFF(fname, ps);
  else if(ext =="ply")
    return read_PLY(fname, ps, np);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return read_LAS(fname, ps);
#endif

  return false;
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector>
bool read_point_set(const std::string& fname, CGAL::Point_set_3<Point, Vector>& ps)
{
  return read_point_set(fname, ps, parameters::all_default());
}
/// \endcond

} // namespace IO


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
  \ingroup PkgPointSet3IO

  \brief writes the point set in an output stream in the \ref IOStreamPLY.

  All properties are inserted in their instantiation order.

  \attention When writing a binary file, the flag `std::ios::binary` flag must be set during the creation of the `ofstream`.

  \param os the output stream
  \param ps the point set

  \return `os`
*/
template <typename Point, typename Vector>
std::ostream& operator<<(std::ostream& os,
                         const CGAL::Point_set_3<Point, Vector>& ps)
{
  IO::write_PLY(os, ps);
  return os;
}

namespace IO {

/*!
  \ingroup PkgPointSet3IO

  \brief writes the point set in an output file.

  Supported file formats are the following:
  - \ref IOStreamOFF (`.off`)
  - \ref IOStreamPLY (`.ply`)
  - \ref IOStreamLAS (`.las`)
  - \ref IOStreamXYZ (`.xyz`)

  The format is detected from the filename extension (letter case is not important).

  \tparam Point the point type of the `Point_set_3`
  \tparam Vector the vector type of the `Point_set_3`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname name of the output file
  \param ps the point set
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{use_binary_mode}
      \cgalParamDescription{indicates whether data should be written in binary (`true`) or in ASCII (`false`)}
      \cgalParamType{Boolean}
      \cgalParamDefault{`true`}
      \cgalParamExtra{This parameter is only relevant for `PLY` writing: the `OFF` and `XYZ` formats
                      are always ASCII, and the `LAS` format is always binary.}
    \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the writing was successful, `false` otherwise.
*/
template <typename Point, typename Vector, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_point_set(const std::string& fname,
                     CGAL::Point_set_3<Point, Vector>& ps,
                     const CGAL_BGL_NP_CLASS& np)
{
  const std::string ext = internal::get_file_extension(fname);

  if(ext == "xyz")
    return write_XYZ(fname, ps, np);
  else if(ext == "off")
    return write_OFF(fname, ps, np);
  else if(ext == "ply")
    return write_PLY(fname, ps, np);
#ifdef CGAL_LINKED_WITH_LASLIB
  else if(ext == "las")
    return write_LAS(fname, ps);
#endif

  return false;
}

/// \cond SKIP_IN_MANUAL

template <typename Point, typename Vector>
bool write_point_set(const std::string& fname, CGAL::Point_set_3<Point, Vector>& ps)
{
  return write_point_set(fname, ps, parameters::all_default());
}

/// \endcond

} // namespace IO

} // namespace CGAL

#endif // CGAL_POINT_SET_3_IO
