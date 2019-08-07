// Copyright (c) 2019 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_GOCAD_H
#define CGAL_IO_GOCAD_H

#include <iostream>
#include <CGAL/IO/GOCAD/GOCAD_internals.h>
#include <string>

namespace CGAL{
template <typename FaceGraph>
bool
read_gocad(FaceGraph& facegraph, std::istream& in, std::string& name, std::string& color)
{
  return GOCAD_internal::read_gocad(facegraph, in, name, color);
}

template <typename FaceGraph>
bool
write_gocad(FaceGraph& facegraph, std::ostream& os, const std::string& name)
{
  return GOCAD_internal::write_gocad(facegraph, os, name);
}

}//end CGAL

#endif //CGAL_IO_GOCAD_H
