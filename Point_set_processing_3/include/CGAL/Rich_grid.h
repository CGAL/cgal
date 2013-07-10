// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
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
// Author(s) : Shihao Wu, Cl¨¦ment Jamin 

#ifndef CGAL_RICH_GRID_H
#define CGAL_RICH_GRID_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>
#include <CGAL/value_type_traits.h>

#include <iterator>
#include <algorithm>
#include <cmath>
#include <ctime>


//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Rich Grid section
// ----------------------------------------------------------------------------
//namespace rich_grid_internel{

namespace rich_grid_internel{
 
/// The Rich_point class represents a 3D point with inedxes of neighbor points;
/// - a position,
/// - a index.
/// - self point set neighbors.
/// - other point set neighbors.
///
/// @heading Parameters:
/// @param Kernel       Geometric traits class.
template <typename Kernel>
class Rich_point
{
public:
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::FT FT;

public:
  Rich_point(const Point& p = CGAL::ORIGIN,
             const int& i = 0,
             const Vector& v = CGAL::NULL_VECTOR
             ):pt(p),index(i),normal(v){} 

  Point pt;
  Vector normal;
  unsigned int index;
  std::vector<unsigned int> neighbors;
  std::vector<unsigned int> original_neighbors;//is not necessary
};

template <typename Kernel>
class Rich_box
{
public:
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;

  Rich_box(){init();}

  void init()
  {
    FT max_double = (FT)(std::numeric_limits<double>::max)();
    min_x = min_y = min_z = max_double;
    max_x = max_y = max_z = -max_double;
  }

  void add_point(const Point& p)
  {
    if (p.x() < min_x) min_x = p.x();
    if (p.y() < min_y) min_y = p.y();
    if (p.z() < min_z) min_z = p.z();
    if (p.x() > max_x) max_x = p.x();
    if (p.y() > max_y) max_y = p.y();
    if (p.z() > max_z) max_z = p.z();
  }

  Point get_min(){return Point(min_x, min_y, min_z);}
  Point get_max(){return Point(max_x, max_y, max_z);}

private:
  FT min_x, min_y, min_z, max_x, max_y, max_z;
};


template <typename Kernel>
class Rich_grid 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  typedef std::vector<Rich_point<Kernel>*> Point_Pointer_vector;
  typedef typename Point_Pointer_vector::iterator iterator;

public:

  Rich_grid() {}

  void init(std::vector<Rich_point<Kernel> > &vert,
            Rich_box<Kernel>& box,
            const FT _radius); 

  // Travel for the point set itself 
  void travel_itself(void (*self)(iterator starta, iterator enda, FT radius),
    void (*other)(iterator starta, iterator enda, 
    iterator startb, iterator endb, FT radius));

  // Travel other self between two point set(original and samples) 
  void travel_others(Rich_grid &points, 
    void (*travel_others)(iterator starta, iterator enda, 
    iterator startb, iterator endb, FT radius));


  void static __cdecl find_original_neighbors(iterator starta, iterator enda, 
    iterator startb, iterator endb, FT radius);
  void static  __cdecl find_self_neighbors(iterator start, 
                                           iterator end, FT radius);
  void static  __cdecl find_other_neighbors(iterator starta, iterator enda, 
    iterator startb, iterator endb, FT radius);

private:
  
  std::vector<Rich_point<Kernel>*> rich_points;  
  std::vector<int> index;   
  int xside, yside, zside;
  FT radius;

  int cell(int x, int y, int z) { return x + xside*(y + yside*z); }
  bool isEmpty(int cell) { return index[cell+1] == index[cell]; }
  iterator startV(int origin) { return rich_points.begin() + index[origin]; }  
  iterator endV(int origin) { return rich_points.begin() + index[origin+1]; }
};


template <typename Kernel>
class XSort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.x() < b->pt.x();
  }
};

template <typename Kernel>
class YSort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.y() < b->pt.y();
  }
};

template <typename Kernel>
class ZSort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.z() < b->pt.z();
  }
};


// divide spoints into some grids
// and each grid has their points index in the index vector of sample.
template <typename Kernel>
void Rich_grid<Kernel>::init(std::vector<Rich_point<Kernel> > &vert, 
                      Rich_box<Kernel>& box, const typename Kernel::FT _radius) 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;

  radius = _radius;

  Point min = box.get_min();
  Point max = box.get_max();

  rich_points.resize(vert.size());
  for(int i = 0; i < rich_points.size(); i++)
  {
    Point& pt = vert[i].pt;  
    rich_points[i] = &vert[i];
  }

  xside = (int)ceil((max.x() - min.x())/radius);
  yside = (int)ceil((max.y() - min.y())/radius);
  zside = (int)ceil((max.z() - min.z())/radius);

  xside = (xside > 0) ? xside : 1;
  yside = (yside > 0) ? yside : 1;
  zside = (zside > 0) ? zside : 1;

  index.resize(xside*yside*zside + 1, -1);  

  std::sort(rich_points.begin(), rich_points.end(), ZSort<Kernel>()); 

  unsigned int startz = 0;
  for(unsigned int z = 0; z < zside; z++) 
  {
    int endz = startz;
    FT maxz = min.z() + (z+1)*radius;
    while(endz < rich_points.size() && rich_points[endz]->pt.z()< maxz)
      ++endz; 

    sort(rich_points.begin()+startz, 
         rich_points.begin()+endz, YSort<Kernel>());

    int starty = startz;
    for(int y = 0; y < yside; y++) 
    {
      int endy = starty;        
      FT maxy = min.y() + (y+1)*radius;
      while(endy < endz && rich_points[endy]->pt.y() < maxy)
        ++endy;

      sort(rich_points.begin()+starty, 
           rich_points.begin()+endy, XSort<Kernel>());

      int startx = starty;
      for(int x = 0; x < xside; x++) 
      {
        int endx = startx;
        index[x + xside*y + xside*yside*z] = endx;          
        FT maxx = min.x() + (x+1)*radius;
        while(endx < endy && rich_points[endx]->pt.x() < maxx)
          ++endx;

        startx = endx;
      }
      starty = endy;
    }
    startz = endz;
  }

  //compute the last grid's range
  index[xside*yside*zside] = startz;
}

/// define how to travel in the same gird 
template <typename Kernel>
void Rich_grid<Kernel>::travel_itself(
  void (*self)(iterator starta, iterator enda, 
               const typename Kernel::FT radius),
  void (*other)(iterator starta, iterator enda, 
  iterator startb, iterator endb, FT radius)
) 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;

  static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
    0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

  static int diagonals[14*2] = { 0, 0, //remove this to avoid self intesextion
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7,
    2, 3, 1, 3, 1, 2,                       
    1, 4, 2, 5, 3, 6 };

  for(int z = 0; z < zside; z++) {
    for(int y = 0; y < yside; y++) {
      for(int x = 0; x < xside; x++) {
        int origin = cell(x, y, z);
        self(startV(origin), endV(origin), radius);   
        // compute between other girds
        for(int d = 2; d < 28; d += 2) { // skipping self
          int *cs = corner + 3*diagonals[d];
          int *ce = corner + 3*diagonals[d+1];
          if((x + cs[0] < xside) && (y + cs[1] < yside) && (z + cs[2] < zside) 
          && (x + ce[0] < xside) && (y + ce[1] < yside) && (z + ce[2] < zside)) 
          {
            origin = cell(x+cs[0], y+cs[1], z+cs[2]);
            int dest = cell(x+ce[0], y+ce[1], z+ce[2]);
            other(startV(origin), endV(origin), 
              startV(dest),   endV(dest), radius);        
          }
        }   
      }
    }
  }
}

/// define how to travel in other gird 
template <typename Kernel>
void Rich_grid<Kernel>::travel_others(
  Rich_grid &points, 
  void (*travel_others)(iterator starta, iterator enda, 
                        iterator startb, iterator endb, 
                        const typename Kernel::FT radius)
) 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;

  static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
    0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

  static int diagonals[14*2] = { 0, 0, //remove this to avoid self intesextion
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7,
    2, 3, 1, 3, 1, 2,                       
    1, 4, 2, 5, 3, 6 };

  for(int z = 0; z < zside; z++) {
    for(int y = 0; y < yside; y++) {
      for(int x = 0; x < xside; x++) {     
        int origin = cell(x, y, z);  


        if(!isEmpty(origin) && !points.isEmpty(origin)) 
          travel_others(startV(origin), endV(origin), 
          points.startV(origin), points.endV(origin), radius);  


        for(int d = 2; d < 28; d += 2) { //skipping self
          int *cs = corner + 3*diagonals[d];
          int *ce = corner + 3*diagonals[d+1];
          if((x+cs[0] < xside) && (y+cs[1] < yside) && (z+cs[2] < zside) &&
            (x+ce[0] < xside) && (y+ce[1] < yside) && (z+ce[2] < zside)) {

              origin   = cell(x+cs[0], y+cs[1], z+cs[2]);

              int dest = cell(x+ce[0], y+ce[1], z+ce[2]);

              if(!isEmpty(origin) && !points.isEmpty(dest))        // Locally 
                travel_others(startV(origin), endV(origin), 
                points.startV(dest),   points.endV(dest), radius); 

              if(!isEmpty(dest) && !points.isEmpty(origin))  
                travel_others(startV(dest), endV(dest), 
                points.startV(origin),   points.endV(origin), radius);        
          }
        }      
      }
    }
  }
}

/// grid travel function to find the neighbors in the original point set
template <typename Kernel>
void Rich_grid<Kernel>::find_original_neighbors(
  iterator starta, 
  iterator enda, 
  iterator startb, 
  iterator endb,
  FT radius
)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  FT radius2 = radius*radius;

  for(Rich_grid::iterator dest = starta; dest != enda; dest++) 
  {
    Rich_point<Kernel> &v = *(*dest);

    Point &p = v.pt;
    for(Rich_grid::iterator origin = startb; origin != endb; origin++)
    {
      Rich_point<Kernel> &t = *(*origin);
      Point &q = t.pt;

      FT dist2 = CGAL::squared_distance(p, q);

      if(dist2 < radius2) 
      {                          
        v.original_neighbors.push_back((*origin)->index);
      }
    }
  }
}

/// grid travel function to find the neighbors in the same point set
template <typename Kernel>
void Rich_grid<Kernel>::find_self_neighbors(
  iterator start, iterator end, FT radius)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  FT radius2 = radius*radius;
  for(iterator dest = start; dest != end; dest++)
  {
    Rich_point<Kernel> &v = *(*dest);
    Point &p = v.pt;

    for(iterator origin = dest+1; origin != end; origin++)
    {
      Rich_point<Kernel> &t = *(*origin);
      Point &q = t.pt;

      FT dist2 = CGAL::squared_distance(p, q);

      if(dist2 < radius2) 
      {   
        v.neighbors.push_back((*origin)->index);
        t.neighbors.push_back((*dest)->index);
      }
    }
  }
}

/// grid travel function to find the neighbors in the same point set
template <typename Kernel>
void Rich_grid<Kernel>::find_other_neighbors(
  iterator starta, iterator enda, 
  iterator startb, iterator endb, FT radius)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  FT radius2 = radius*radius;
  for(iterator dest = starta; dest != enda; dest++)
  {
    Rich_point<Kernel> &v = *(*dest);
    Point &p = v.pt;

    for(iterator origin = startb; origin != endb; origin++)
    {
      Rich_point<Kernel> &t = *(*origin);
      Point &q = t.pt;

      FT dist2 = CGAL::squared_distance(p, q);

      if(dist2 < radius2) 
      {   
        v.neighbors.push_back((*origin)->index);
        t.neighbors.push_back((*dest)->index);
      }
    }
  }
}




/// Compute ball neighbors for each point in the same point set.
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return 
template <typename Kernel>
void compute_ball_neighbors_one_self(
  std::vector<Rich_point<Kernel> >& points,
  Rich_box<Kernel>& box,
  const typename Kernel::FT radius)
{
  typedef typename Kernel::FT FT;
  CGAL_point_set_processing_precondition(radius > 0);

  for (unsigned int i = 0; i < points.size(); i++)
  {
    points[i].neighbors.clear();
  }

  Rich_grid<Kernel> points_grid;
  points_grid.init(points, box, radius);
  points_grid.travel_itself(Rich_grid<Kernel>::find_self_neighbors, 
                            Rich_grid<Kernel>::find_other_neighbors);
}

/// Compute ball neighbors for each (sample)points in the other point set
/// 
/// \pre `radius > 0`
///
/// @tparam Kernel Geometric traits class.
///
/// @return 
template <typename Kernel>
void compute_ball_neighbors_one_to_another(
  std::vector<Rich_point<Kernel> >& samples, ///< sample point set
  std::vector<Rich_point<Kernel> >& original,///< original point set
  Rich_box<Kernel>& box, ///< bounding box
  const typename Kernel::FT radius ///< neighbor radius
)
{
  typedef typename Kernel::FT FT;
  if (radius < FT(0.0))
  {
    return;
  }

  for (unsigned int i = 0; i < samples.size(); i++)
  {
    samples[i].original_neighbors.clear();
  }

  Rich_grid<Kernel> samples_grid;
  samples_grid.init(samples, box, radius);

  // here can be optimized by initial the original grid just one time.
  Rich_grid<Kernel> original_grid;
  original_grid.init(original, box, radius);

  samples_grid.travel_others(original_grid, 
                             Rich_grid<Kernel>::find_original_neighbors);
}


} //namespace internal

} //namespace CGAL

#endif // CGAL_RICH_GRID_H
