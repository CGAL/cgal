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
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
   
#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/debug.h>
#include <QDir>
#include <iostream>
#include <CGAL/gl.h>
#include <qopenglfunctions.h>
namespace CGAL {
namespace Qt {


CGAL_INLINE_FUNCTION
void traverse_resources(const QString& name, const QString& dirname, int indent)
{
  std::cerr << qPrintable(QString(indent, ' '))
            << qPrintable(name);
  QString fullname = 
    dirname.isEmpty() ?
    name :
    dirname + "/" + name;
  QDir dir(fullname);
  if(dir.exists()) {
    std::cerr << "/\n";
    Q_FOREACH(QString path, dir.entryList())
    {
      traverse_resources(path, fullname, indent + 2);
    }
  }
  else {
    std::cerr << "\n";
  }
}

CGAL_INLINE_FUNCTION
void opengl_check_errors(unsigned int line)
{
 GLenum error = ::glGetError();
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
   if(error == GL_STACK_UNDERFLOW)
     std::cerr << "An attempt has been made to perform an operation that would cause an internal stack to underflow." << "@" << line << std::endl;
   if(error == GL_STACK_OVERFLOW)
     std::cerr << "An attempt has been made to perform an operation that would cause an internal stack to overflow." << "@" << line << std::endl;
   error = ::glGetError();
 }
}
} // namesapce Qt
} // namespace CGAL
