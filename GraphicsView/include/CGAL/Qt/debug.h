// Copyright (c) 2008  GeometryFactory Sarl (France).
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
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_DEBUG_H
#define CGAL_QT_DEBUG_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/auto_link/Qt.h>
#include <CGAL/export/Qt.h>

#include <QString>

namespace CGAL {
namespace Qt {

/**
 *  Must be used like that:
 *     CGAL::Qt:traverse_resources(":/cgal"); // view CGAL resources
 *  or
 *     CGAL::Qt:traverse_resources(":"); // view all resources
 *  and displays the resources tree on std::cerr.
 */
CGAL_QT_EXPORT void traverse_resources(const QString& name,
                                        const QString& dirname = QString(),
                                        int indent = 0);

/**
 * Call this in the end of an OpenGL implementation to check if it returns errors. 
 */
template <typename QtOpenGLFunctions>
void opengl_check_errors(QtOpenGLFunctions* gl,
                         unsigned int line) {
  GLenum error = gl->glGetError();
  while (error != GL_NO_ERROR)
  {
    if(error == GL_INVALID_ENUM)
      std::cerr << "An unacceptable value is specified for an enumerated argument." << "@" << line << std::endl;
    if(error == GL_INVALID_VALUE)
      std::cerr << "A numeric argument is out of range." << "@" << line << std::endl;
    if(error == GL_INVALID_OPERATION)
      std::cerr << "The specified operation is not allowed in the current state." << "@" << line << std::endl;
    if(error == GL_INVALID_FRAMEBUFFER_OPERATION)
      std::cerr << "The framebuffer object is not complete." << "@" << line << std::endl;
    if(error == GL_OUT_OF_MEMORY)
      std::cerr << "There is not enough memory left to execute the command." << "@" << line << std::endl;
#ifdef GL_STACK_UNDERFLOW
    if(error == GL_STACK_UNDERFLOW)
      std::cerr << "An attempt has been made to perform an operation that would cause an internal stack to underflow." << "@" << line << std::endl;
#endif
#ifdef GL_STACK_OVERFLOW
    if(error == GL_STACK_OVERFLOW)
      std::cerr << "An attempt has been made to perform an operation that would cause an internal stack to overflow." << "@" << line << std::endl;
#endif
    error = gl->glGetError();
  }
}

} // namespace Qt
} // namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/debug_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_QT_DEBUG_H
