// Copyright (c) 2012,2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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
