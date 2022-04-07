// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_OGL_BASE_OBJECT_H
#define CGAL_OGL_BASE_OBJECT_H

#include <CGAL/Simple_cartesian.h>

namespace CGAL {

namespace OGL {

  class OGL_base_object {
  public:

    typedef CGAL::Simple_cartesian<double>       Double_kernel;
    typedef Double_kernel::Point_3               Double_point;
    typedef Double_kernel::Vector_3              Double_vector;
    typedef Double_kernel::Segment_3             Double_segment;
    typedef Double_kernel::Aff_transformation_3  Affine_3;

    virtual void draw() const  = 0;
    virtual void init() = 0;
    virtual void toggle(int) = 0;
    virtual void set_style(int) = 0;
    virtual ~OGL_base_object() {}
  };
}

} //namespace CGAL
#endif // CGAL_OGL_BASE_OBJECT_H
