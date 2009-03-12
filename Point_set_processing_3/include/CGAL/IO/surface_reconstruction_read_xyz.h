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

#ifndef CGAL_SURFACE_RECONSTRUCTION_READ_XYZ_H
#define CGAL_SURFACE_RECONSTRUCTION_READ_XYZ_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <deque>
#include <stdio.h>

CGAL_BEGIN_NAMESPACE


/// Read points (positions + normals, if available) from a .xyz file (ASCII).
///
/// @commentheading Template Parameters:
/// @param OutputIterator value_type must be a model of the PointWithNormal_3 concept.
///
/// @return true on success.
template <typename OutputIterator>
bool surface_reconstruction_read_xyz(const char* pFilename,
                                     OutputIterator output)
{
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename value_type_traits<OutputIterator>::type Point_with_normal;

  typedef typename Point_with_normal::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;

  CGAL_precondition(pFilename != NULL);

  FILE *pFile = fopen(pFilename,"rt");
  if(pFile == NULL)
  {
    std::cerr << "Error: cannot open " << pFilename;
    return false;
  }

  // scan points
  int lineNumber = 0;
  char pLine[4096];
  while(fgets(pLine,4096,pFile))
  {
    lineNumber++;

    // Read position + normal...
    double x,y,z;
    double nx,ny,nz;
    long pointsCount;
    if(sscanf(pLine,"%lg %lg %lg %lg %lg %lg",&x,&y,&z,&nx,&ny,&nz) == 6)
    {
      Point point(x,y,z);
      Vector normal(nx,ny,nz);
      *output = Point_with_normal(point,normal);
      output++;
    }
    // ...or read only position...
    else if(sscanf(pLine,"%lg %lg %lg",&x,&y,&z) == 3)
    {
      Point point(x,y,z);
      *output = point;
      output++;
    }
    // ...or skip number of points (optional)
    else if (lineNumber == 1 && sscanf(pLine,"%ld",&pointsCount) == 1)
    {
      continue;
    }
    else
    {
      std::cerr << "Error line " << lineNumber << " of " << pFilename;
      return false;
    }
  }

  fclose(pFile);
  return true;
}

/// Read points (positions only) from a .xyz file (ASCII).
/// Normals are ignored.
///
/// @commentheading Template Parameters:
/// @param OutputIterator value_type must be Point_3.
///
/// @return true on success.
template <typename OutputIterator>
bool surface_reconstruction_read_xyz(const char* pFilename,
                                     OutputIterator output,
                                     bool read_normals /* must be false */)
{
  CGAL_surface_reconstruction_assertion(read_normals == false);

  typedef typename value_type_traits<OutputIterator>::type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  typedef typename CGAL::Point_with_normal_3<Kernel> Point_with_normal;

  // Read file in temporary Point_with_normal_3 container
  std::deque<Point_with_normal> pwns;
  if (surface_reconstruction_read_xyz(pFilename, std::back_inserter(pwns)))
  {
    // copy to Point_3 container, removing normals
    std::copy(pwns.begin(), pwns.end(), output);

    return true;
  }
  else
    return false;
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_READ_XYZ_H
