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


namespace CGAL {

template <typename Point,
          typename Vector>
bool read_point_set(const std::string& fname,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  if (fname.find(".xyz") != std::string::npos) {
    return read_XYZ(fname, ps);
  }

  if (fname.find(".off") != std::string::npos) {
    return read_OFF(fname, ps);
  }

  if (fname.find(".ply") != std::string::npos) {
    return read_PLY(fname, ps);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (fname.find(".las") != std::string::npos) {
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

template <typename Point,
          typename Vector>
bool write_point_set(const std::string& fname,
                    CGAL::Point_set_3<Point, Vector>& ps)
{
  if (fname.find(".xyz") != std::string::npos) {
    return write_XYZ(fname, ps);
  }

  if (fname.find(".off") != std::string::npos) {
    return write_OFF(fname, ps);
  }

  if (fname.find(".ply") != std::string::npos) {
    return write_PLY(fname, ps);
  }

#ifdef CGAL_LINKED_WITH_LASLIB
  if (fname.find(".las") != std::string::npos) {
    return write_LAS(fname, ps);
  }
#endif
  return false;
}
}

#endif // CGAL_POINT_SET_IO_H
