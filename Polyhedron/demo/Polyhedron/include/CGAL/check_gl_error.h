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
#include <QOpenGLFunctions>
#include <iostream>
#include <string>
/*

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

inline bool gl_check_compilation(GLuint shader)
{
    GLint result;
    glGetShaderiv(shader,GL_COMPILE_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Vertex compilation OK"<<std::endl;
        return true;
    } else {
        int maxLength;
        int length;
        glGetShaderiv(shader,GL_INFO_LOG_LENGTH,&maxLength);
        std::string log;
        glGetShaderInfoLog(shader,maxLength,&length,&log[0]);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
        return false;
    }
}
inline bool gl_check_link(GLuint *program)
{
    GLint result;
    glGetProgramiv(*program,GL_LINK_STATUS,&result);
    if(result == GL_TRUE){
        std::cout<<"Link OK"<<std::endl;
        return true;
    } else {
        int maxLength;
        int length;
        glGetProgramiv(*program,GL_INFO_LOG_LENGTH,&maxLength);
        std::string log;
        glGetProgramInfoLog(*program,maxLength,&length,&log[0]);
        std::cout<<"link error : Length = "<<length<<", log ="<<log<<std::endl;
        return false;
    }
}

} // end namespace CGAL
*/
#endif // CGAL_GL_CHECK_ERROR_H
