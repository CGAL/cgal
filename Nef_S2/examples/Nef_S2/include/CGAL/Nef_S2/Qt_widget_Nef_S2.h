// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_QT_WIDGET_NEF_S2_H
#define CGAL_QT_WIDGET_NEF_S2_H

#include <CGAL/Nef_S2/Sphere_geometry_OGL.h>
#include "CGAL/Nef_S2/Qt_widget_OpenGL.h"
#include <boost/thread/mutex.hpp>

namespace CGAL {

class Qt_widget_Nef_S2 : public Qt_widget_OpenGL {
    Q_OBJECT

 public:
        template <typename Nef_polyhedron>
            void addObject(const typename Nef_polyhedron::Const_decorator& N, const std::string & label=std::string()) {
                CGAL::OGL::OGL_base_object* O = new CGAL::OGL::Unit_sphere(CGAL::OGL::NefS2_to_UnitSphere<Nef_polyhedron>::convert(N));
                mutex.lock();
                object_.push_back(O);
                label_.push_back(label);
                mutex.unlock();
            }

  Qt_widget_Nef_S2();

  virtual ~Qt_widget_Nef_S2(); 

 protected:
  boost::mutex mutex;
  
};

} //namespace CGAL
#endif // CGAL_QT_WIDGET_NEF_S2_H
