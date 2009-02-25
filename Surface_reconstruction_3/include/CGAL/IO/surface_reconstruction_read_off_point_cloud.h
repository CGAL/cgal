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

#ifndef CGAL_SURFACE_RECONSTRUCTION_READ_OFF_POINT_CLOUD_H
#define CGAL_SURFACE_RECONSTRUCTION_READ_OFF_POINT_CLOUD_H

#include <CGAL/value_type_traits.h>
#include <CGAL/surface_reconstruction_assertions.h>

#include <stdio.h>

CGAL_BEGIN_NAMESPACE


/// Read points (positions + normals, if available) from a .off file (ASCII).
///
/// @commentheading Template Parameters:
/// @param OutputIterator value_type must be a model of the PointWithNormal_3 concept.
///
/// @return true on success.
template <typename OutputIterator>
bool surface_reconstruction_read_off_point_cloud(const char* pFilename, OutputIterator output)
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
  long pointsCount = 0, facesCount = 0, edgesCount = 0; // number of items in file
  int pointsRead = 0; // current number of points read
  int lineNumber = 0; // current line number
  char pLine[4096]; // current line buffer
  while(fgets(pLine,4096,pFile))
  {
    lineNumber++;

    // Read file signature on first line
    if (lineNumber == 1)
    {
      char signature[4096];
      if ( (sscanf(pLine,"%s",signature) != 1)
        || (strcmp(signature, "OFF") != 0 && strcmp(signature, "NOFF") != 0) )
      {
        // if unsupported file format
        std::cerr << "Incorrect file format line " << lineNumber << " of " << pFilename;
        return false;
      }
    }

    // Read number of points on 2nd line
    else if (lineNumber == 2)
    {
      if (sscanf(pLine,"%ld %ld %ld",&pointsCount,&facesCount,&edgesCount) != 3)
      {
        std::cerr << "Error line " << lineNumber << " of " << pFilename;
        return false;
      }
    }

    // Read 3D points on next lines
    else if (pointsRead < pointsCount)
    {
      // Read position + normal...
      double x,y,z;
      double nx,ny,nz;
      if(sscanf(pLine,"%lg %lg %lg %lg %lg %lg",&x,&y,&z,&nx,&ny,&nz) == 6)
      {
        Point point(x,y,z);
        Vector normal(nx,ny,nz);
        *output = Point_with_normal(point,normal);
        output++;
        pointsRead++;
      }
      // ...or read only position...
      else if(sscanf(pLine,"%lg %lg %lg",&x,&y,&z) == 3)
      {
        Point point(x,y,z);
        *output = point;
        output++;
        pointsRead++;
      }
      // ...or skip comment line
    }
    // Skip remaining lines
  }

  fclose(pFile);
  return true;
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_READ_OFF_POINT_CLOUD_H
