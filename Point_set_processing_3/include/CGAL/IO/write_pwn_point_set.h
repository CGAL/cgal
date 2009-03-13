// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_WRITE_PWN_POINT_SET_H
#define CGAL_WRITE_PWN_POINT_SET_H

#include <CGAL/point_set_processing_assertions.h>

#include <stdio.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE


/// Save points to .pwn file (position + normal, ASCII).
/// @return true on success.
template <typename InputIterator> ///< InputIterator value_type must be
                                  ///< a model of the PointWithNormal_3 concept.
bool write_pwn_point_set(const char* pFilename,  ///< output file
                         InputIterator first,    ///< first input point
                         InputIterator beyond)   ///< past-the-end input point
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_with_normal;

  typedef typename Point_with_normal::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;

  CGAL_precondition(pFilename != NULL);
  CGAL_point_set_processing_precondition(first != beyond);

  FILE *pFile = fopen(pFilename,"wt");
  if(pFile == NULL)
    return false;

  for(InputIterator it = first; it != beyond; it++)
  {
    const Point& p = *it;
    const Vector& n = it->normal();
    fprintf(pFile,"%lf %lf %lf %lf %lf %lf\n",
      p.x(),p.y(),p.z(),
      n.x(),n.y(),n.z());
  }

  fclose(pFile);
  return true;
}


CGAL_END_NAMESPACE

#endif // CGAL_WRITE_PWN_POINT_SET_H
