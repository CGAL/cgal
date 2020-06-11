// Copyright (c) 2020  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno

#ifndef CGAL_POINT_SET_IO_H
#define CGAL_POINT_SET_IO_H

#include <string>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/polygon_soup_io.h>


namespace CGAL {

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point,
          typename Vector>
bool read_point_set(const std::string& fname,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  const std::string ext = IO::internal::get_file_extension(fname);

  if (ext=="xyz") {
    return read_XYZ(fname, ps);
  }

  if (ext == "off") {
    return read_OFF(fname, ps);
  }

  if (ext =="ply") {
    return read_PLY(fname, ps);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (ext == "las") {
    return read_LAS(fname, ps);
  }
#endif
  return false;
}

template <typename Point,
          typename Vector>
bool read_point_set(const char* fname,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  return read_point_set(std::string(fname),
                        ps);
}

/*!
  \ingroup PkgPointSet3IO
 */
template <typename Point,
          typename Vector>
bool write_point_set(const std::string& fname,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  const std::string ext = IO::internal::get_file_extension(fname);
  if (ext == "xyz") {
    return write_XYZ(fname, ps);
  }

  if (ext == "off") {
    return write_OFF(fname, ps);
  }

  if (ext == "ply") {
    return write_PLY(fname, ps);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (ext == "las") {
    return write_LAS(fname, ps);
  }
#endif
  return false;
}


}

#endif // CGAL_POINT_SET_IO_H
