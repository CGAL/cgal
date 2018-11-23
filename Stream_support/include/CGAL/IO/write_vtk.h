// Copyright (c) 2018
// GeometryFactory( France) All rights reserved. 
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
//
// Author(s)     : Stephane Tayeb

#ifndef CGAL_WRITE_VTK_IO_H
#define CGAL_WRITE_VTK_IO_H
#include <fstream>
#include <vector>
template <class FT> 
void
write_vector(std::ostream& os, 
             const std::vector<FT>& vect) 
{
  const char* buffer = reinterpret_cast<const char*>(&(vect[0]));
  std::size_t size = vect.size()*sizeof(FT);
  
  os.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t)); // number of bytes encoded
  os.write(buffer, vect.size()*sizeof(FT));                     // encoded data
}
#endif 
