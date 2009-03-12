// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_SURFACE_RECONSTRUCTION_WRITE_OFF_POINT_CLOUD_H
#define CGAL_SURFACE_RECONSTRUCTION_WRITE_OFF_POINT_CLOUD_H

#include <CGAL/surface_reconstruction_assertions.h>

#include <stdio.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE


/// Save points (positions + normals) to a .off file (ASCII).
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be a model of the PointWithNormal_3 concept.
///
/// @return true on success.
template <typename InputIterator>
bool surface_reconstruction_write_off_point_cloud(const char* pFilename,
                                                  InputIterator first,    ///< first input point
                                                  InputIterator beyond)   ///< past-the-end input point
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_with_normal;

  typedef typename Point_with_normal::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;

  CGAL_precondition(pFilename != NULL);
  CGAL_surface_reconstruction_precondition(first != beyond);

  FILE *pFile = fopen(pFilename,"wt");
  if(pFile == NULL)
    return false;

  // Write header
  const int num_input_vertices = std::distance(first, beyond);
  fprintf(pFile,"NOFF\n");
  fprintf(pFile,"%d 0 0\n", num_input_vertices);

  // Write positions + normals
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

/// Save points (positions + optionally normals) to a .off file (ASCII).
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be a model of PointWithNormal_3 if
/// write_normals is true, else a model of Kernel::Point_3.
///
/// @return true on success.
template <typename InputIterator>
bool surface_reconstruction_write_off_point_cloud(const char* pFilename,
                                                  InputIterator first,    ///< first input point
                                                  InputIterator beyond,   ///< past-the-end input point
                                                  bool write_normals)
{
  if(write_normals)
  {
    return surface_reconstruction_write_off_point_cloud(pFilename, first, beyond);
  }
  else
  {
    // model of Kernel::Point_3
    typedef typename std::iterator_traits<InputIterator>::value_type Point;

    CGAL_precondition(pFilename != NULL);
    CGAL_surface_reconstruction_precondition(first != beyond);

    FILE *pFile = fopen(pFilename,"wt");
    if(pFile == NULL)
      return false;

    // Write header
    const int num_input_vertices = std::distance(first, beyond);
    fprintf(pFile,"OFF\n");
    fprintf(pFile,"%d 0 0\n", num_input_vertices);

    // Write positions
    for(InputIterator it = first; it != beyond; it++)
    {
      const Point& p = *it;
      fprintf(pFile,"%lf %lf %lf\n",
                    p.x(),p.y(),p.z());
    }

    fclose(pFile);
    return true;
  }
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_WRITE_OFF_POINT_CLOUD_H
