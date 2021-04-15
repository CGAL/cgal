// Copyright (c) 2012,2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_ERROR_CODE_H
#define CGAL_MESH_ERROR_CODE_H

#include <CGAL/license/Mesh_3.h>

#include <string>
#include <sstream>

namespace CGAL {

enum Mesh_error_code {
  CGAL_MESH_3_NO_ERROR = 0,
  CGAL_MESH_3_MAXIMAL_NUMBER_OF_VERTICES_REACHED,
  CGAL_MESH_3_STOPPED
};

inline
std::string mesh_error_string(const Mesh_error_code& error_code) {
  switch(error_code) {
  case CGAL_MESH_3_NO_ERROR:
    return "no error";
  case CGAL_MESH_3_MAXIMAL_NUMBER_OF_VERTICES_REACHED:
    return "the maximal number of vertices has been reached";
  case CGAL_MESH_3_STOPPED:
    return "the meshing process was stopped";
  default:
    std::stringstream str("");
    str << "unknown error (error_code="
        << error_code
        << ")";
    return str.str();
  }
}

}

#endif // CGAL_MESH_ERROR_CODE_H
