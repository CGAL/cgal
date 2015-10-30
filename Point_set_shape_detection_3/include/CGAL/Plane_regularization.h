// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s)     : 
//

/**
* \ingroup PkgPointSetShapeDetection3
* \file CGAL/Plane_regularization.h
*
*/


#ifndef CGAL_PLANE_REGULARIZATION_H
#define CGAL_PLANE_REGULARIZATION_H

template <typename GeomTraits>
class Plane_regularization
{
public:

  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Plane_3 Plane;

private:

public:

  Plane_regularization ()
  {

  }

  virtual ~Plane_regularization ()
  {

  }

  
  

};



#endif // CGAL_PLANE_REGULARIZATION_H
