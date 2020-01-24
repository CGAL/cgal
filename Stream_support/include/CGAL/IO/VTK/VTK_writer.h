// Copyright (c) 2018
// GeometryFactory (France) All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb

#ifndef CGAL_IO_VTK_WRITER_H
#define CGAL_IO_VTK_WRITER_H

#include <fstream>
#include <iostream>
#include <vector>

namespace CGAL {
namespace IO {
namespace internal {

template <class FT>
void write_vector(std::ostream& os,
                  const std::vector<FT>& vect)
{
  const char* buffer = reinterpret_cast<const char*>(&(vect[0]));
  std::size_t size = vect.size()*sizeof(FT);

  os.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t)); // number of bytes encoded
  os.write(buffer, vect.size()*sizeof(FT)); // encoded data
}

} // namespace internal
} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_VTK_WRITER_H
