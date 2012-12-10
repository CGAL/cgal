// Copyright (c) 2011 GeometryFactory, Sophia Antipolis (France)
//  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author: Laurent Rineau

#ifndef CGAL_GL_CHECK_ERROR_H
#define CGAL_GL_CHECK_ERROR_H

#include <iostream>
#include <CGAL/glu.h>

namespace CGAL {

inline bool check_gl_error(const char* filename, long line)
{
  GLenum error = glGetError();
  if(error != GL_NO_ERROR) {
    std::cerr << "GL errors! file " << filename << ", line:" << line << "\n";
    do {
      std::cerr << gluErrorString(error) << std::endl;
      error = glGetError();
    }
    while(error != GL_NO_ERROR);
    std::cerr << "end of errors\n";
    return true;
  }
  return false;
}

} // end namespace CGAL

#endif // CGAL_GL_CHECK_ERROR_H
