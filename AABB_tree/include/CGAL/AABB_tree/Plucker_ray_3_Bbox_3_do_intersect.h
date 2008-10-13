// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : 
#ifndef CGAL_PLUCKER_RAY_3_BBOX_3_DO_INTERSECT_H
#define CGAL_PLUCKER_RAY_3_BBOX_3_DO_INTERSECT_H

CGAL_BEGIN_NAMESPACE

template <class K>
bool do_intersect_type_0(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;
  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero < xmin || ftzero < ymin || ftzero < zmin)
    return false;

  //if it lies on correct side of the silhouette(6 edges) then it intersects
  //MM*
  if((diry*xmax) > (dirx*ymin))//DH
    return false;
  if((diry*xmin) < (dirx*ymax))//BF
    return false;

  //*MM
  if((diry*zmax) > (dirz*ymin))//HE
    return false;
  if((diry*zmin) < (dirz*ymax))//CB
    return false;

  //M*M
  if((dirx*zmax) > (dirz*xmin))//EF
    return false;
  if((dirx*zmin) < (dirz*xmax))//DC
    return false;

  return true;
}


template <class K>
bool do_intersect_type_1(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero < xmin || ftzero < ymin || ftzero > zmax)
    return false;

  //MM*
  if((diry*xmax) > (dirx*ymin))//DH
    return false;
  if((diry*xmin) < (dirx*ymax))//BF
    return false;

  //*MP
  if((diry*zmax) > (dirz*ymax))//GF
    return false;
  if((diry*zmin) < (dirz*ymin))//DA
    return false;

  //M*P
  if((dirx*zmax) > (dirz*xmax))//HG
    return false;
  if((dirx*zmin) < (dirz*xmin))//AB
    return false;

  return true;
}

template <class K>
bool do_intersect_type_2(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero < xmin || ftzero > ymax || ftzero < zmin)
    return false;

  //MP*
  if((diry*xmax) < (dirx*ymax))//CG
    return false;
  if((diry*xmin) > (dirx*ymin))//AE
    return false;

  //*PM
  if((diry*zmax) < (dirz*ymax))//GF
    return false;
  if((diry*zmin) > (dirz*ymin))//DA
    return false;

  //M*M
  if((dirx*zmax) > (dirz*xmin))//EF
    return false;
  if((dirx*zmin) < (dirz*xmax))//DC
    return false;

  return true;
}

template <class K>
bool do_intersect_type_3(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero < xmin || ftzero > ymax || ftzero > zmax)
    return false;

  //MP*
  if((diry*xmax) < (dirx*ymax))//CG
    return false;
  if((diry*xmin) > (dirx*ymin))//AE
    return false;

  //*PP
  if((diry*zmax) < (dirz*ymin))//HE
    return false;
  if((diry*zmin) > (dirz*ymax))//CB
    return false;

  //M*P
  if((dirx*zmax) > (dirz*xmax))//HG
    return false;
  if((dirx*zmin) < (dirz*xmin))//AB
    return false;

  return true;
}

template <class K>
bool do_intersect_type_4(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero > xmax || ftzero < ymin || ftzero < zmin)
    return false;

  //PM*
  if((diry*xmax) > (dirx*ymax))//CG
    return false;
  if((diry*xmin) < (dirx*ymin))//AE
    return false;

  //*MM
  if((diry*zmax) > (dirz*ymin))//HE
    return false;
  if((diry*zmin) < (dirz*ymax))//CB
    return false;

  //P*M
  if((dirx*zmax) < (dirz*xmax))//HG
    return false;
  if((dirx*zmin) > (dirz*xmin))//AB
    return false;

  return true;
}

template <class K>
bool do_intersect_type_5(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero > xmax || ftzero < ymin || ftzero > zmax)
    return false;

  //PM*
  if((diry*xmax) > (dirx*ymax))//CG
    return false;
  if((diry*xmin) < (dirx*ymin))//AE
    return false;

  //*MP
  if((diry*zmax) > (dirz*ymax))//GF
    return false;
  if((diry*zmin) < (dirz*ymin))//DA
    return false;

  //P*P;
  if((dirx*zmax) < (dirz*xmin))//EF
    return false;
  if((dirx*zmin) > (dirz*xmax))//DC
    return false;

  return true;
}

template <class K>
bool do_intersect_type_6(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero > xmax || ftzero > ymax || ftzero < zmin)
    return false;

  //PP*
  if((diry*xmax) < (dirx*ymin))//DH
    return false;
  if((diry*xmin) > (dirx*ymax))//BF
    return false;

  //*PM
  if((diry*zmax) < (dirz*ymax))//GF
    return false;
  if((diry*zmin) > (dirz*ymin))//DA
    return false;

  //P*M
  if((dirx*zmax) < (dirz*xmax))//HG
    return false;
  if((dirx*zmin) > (dirz*xmin))//AB
    return false;

  return true;
}

template <class K>
bool do_intersect_type_7(const typename K::Ray_3& ray, 
		                     const CGAL::Bbox_3& bbox)
{
	typedef typename K::FT FT;
	typedef typename K::Point_3 Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Ray_3 Ray;

  const Point source = ray.source();

  const FT xmin = bbox.xmin()-source.x();
  const FT xmax = bbox.xmax()-source.x();
  const FT ymin = bbox.ymin()-source.y();
  const FT ymax = bbox.ymax()-source.y();
  const FT zmin = bbox.zmin()-source.z();
  const FT zmax = bbox.zmax()-source.z();

  const Vector direction1 = ray.to_vector();

  const FT dirx = direction1.x();
  const FT diry = direction1.y();
  const FT dirz = direction1.z();

  const FT ftzero = (FT)0.0;

  if(ftzero > xmax || ftzero > ymax || ftzero > zmax)
    return false;

  //PP*
  if((diry*xmax) < (dirx*ymin))//DH
    return false;
  if((diry*xmin) > (dirx*ymax))//BF
    return false;

  //*PP
  if((diry*zmax) < (dirz*ymin))//HE
    return false;
  if((diry*zmin) > (dirz*ymax))//CB
    return false;

  //P*P
  if((dirx*zmax) < (dirz*xmin))//EF
    return false;
  if((dirx*zmin) > (dirz*xmax))//DC
    return false;

  return true;
}


CGAL_END_NAMESPACE

#endif  // CGAL_PLUCKER_RAY_3_BBOX_3_DO_INTERSECT_H
