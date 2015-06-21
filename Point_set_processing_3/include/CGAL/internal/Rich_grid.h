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
// Author(s) : Shihao Wu, Clement Jamin, Pierre Alliez  

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

namespace CGAL {

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Rich Grid section
// ----------------------------------------------------------------------------
//namespace rich_grid_internal{

namespace rich_grid_internal{
 
/// The Rich_point class represents a 3D point with inedxes of neighbor points;
/// - a position,
/// - an index.
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
             int i = 0,
             const Vector& n = CGAL::NULL_VECTOR
             ):pt(p), index(i), normal(n){} 
  
  ~Rich_point()
  {
    neighbors.clear();
    original_neighbors.clear();
  }

public:
  Point pt;
  unsigned int index;
  Vector normal;
  std::vector<unsigned int> neighbors;
  std::vector<unsigned int> original_neighbors;//it's not necessary
};

template <typename Kernel>
class Rich_grid 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  typedef std::vector<Rich_point<Kernel>*> Point_ptr_vector;
  typedef typename Point_ptr_vector::iterator iterator;

public:

  Rich_grid() {}

  void init(std::vector<Rich_point<Kernel> > &vert,
            CGAL::Bbox_3 bbox,
            const FT _radius); 

  // Travel for the point set itself 
  void travel_itself(void (*self)(iterator starta, iterator enda, FT radius),
                    void (*other)(iterator starta, iterator enda, 
                    iterator startb, iterator endb, FT radius));

  // Travel other self between two point set(original and samples) 
  void travel_others(Rich_grid &points, 
                     void (*travel_others)(iterator starta, 
                                           iterator enda, 
                                           iterator startb, 
                                           iterator endb, 
                                           FT radius));


  // functions for neighborhood searching
  static void find_original_neighbors(iterator starta, 
                                              iterator enda, 
                                              iterator startb,
                                              iterator endb, 
                                              FT radius);

  static void find_self_neighbors(iterator start, 
                                           iterator end, 
                                           FT radius);
 
  static void find_other_neighbors(iterator starta, 
                                            iterator enda, 
                                            iterator startb, 
                                            iterator endb, 
                                            FT radius);

private:
  
  std::vector<Rich_point<Kernel>*> rich_points;  
  std::vector<int> indices;   
  int x_side, y_side, z_side;
  FT radius;

  int cell(int x, int y, int z) { return x + x_side * (y + y_side * z); }
  bool is_empty(int cell) { return indices[cell+1] == indices[cell]; }

  iterator get_start_iter(int origin) 
  { 
    return rich_points.begin() + indices[origin]; 
  }  

  iterator get_end_iter(int origin) 
  { 
    return rich_points.begin() + indices[origin+1]; 
  }
};


template <typename Kernel>
class X_Sort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.x() < b->pt.x();
  }
};

template <typename Kernel>
class Y_Sort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.y() < b->pt.y();
  }
};

template <typename Kernel>
class Z_Sort {
public:
  bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
    return a->pt.z() < b->pt.z();
  }
};


// divide spoints into some grids
// and each grid has their points index in the index vector of sample.
template <typename Kernel>
void Rich_grid<Kernel>::init(std::vector<Rich_point<Kernel> > &vert, 
                             CGAL::Bbox_3 bbox,
                             const typename Kernel::FT _radius) 
{
  typedef typename Kernel::FT FT;

  rich_points.resize(vert.size());
  for(size_t i = 0; i < rich_points.size(); ++i)
  {
    rich_points[i] = &vert[i];
  }

  radius = _radius;

  x_side = (unsigned int)ceil((bbox.xmax() - bbox.xmin()) / radius);
  y_side = (unsigned int)ceil((bbox.ymax() - bbox.ymin()) / radius);
  z_side = (unsigned int)ceil((bbox.zmax() - bbox.zmin()) / radius);

  x_side = (x_side > 0) ? x_side : 1;
  y_side = (y_side > 0) ? y_side : 1;
  z_side = (z_side > 0) ? z_side : 1;

  indices.resize(x_side * y_side * z_side + 1, -1);  

  std::sort(rich_points.begin(), rich_points.end(), Z_Sort<Kernel>()); 

  unsigned int start_z = 0;
  for(int z = 0; z < z_side; z++) 
  {
    unsigned int end_z = start_z;
    FT max_z = bbox.zmin() + FT(z+1)*radius;
    while(end_z < rich_points.size() && rich_points[end_z]->pt.z() < max_z)
      ++end_z; 

    sort(rich_points.begin() + start_z, 
         rich_points.begin() + end_z, Y_Sort<Kernel>());

    unsigned int start_y = start_z;
    for(int y = 0; y < y_side; y++) 
    {
      unsigned int end_y = start_y;        
      FT max_y = bbox.ymin() + FT(y+1) * radius;
      while(end_y < end_z && rich_points[end_y]->pt.y() < max_y)
        ++end_y;

      sort(rich_points.begin() + start_y, 
           rich_points.begin() + end_y, X_Sort<Kernel>());

      unsigned int start_x = start_y;
      for(int x = 0; x < x_side; x++) 
      {
        unsigned int end_x = start_x;
        indices[x + x_side * y + x_side * y_side * z] = end_x;          
        FT max_x = bbox.xmin() + FT(x+1) * radius;
        while(end_x < end_y && rich_points[end_x]->pt.x() < max_x)
          ++end_x;

        start_x = end_x;
      }
      start_y = end_y;
    }
    start_z = end_z;
  }

  //compute the last grid's range
  indices[x_side * y_side * z_side] = start_z;
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
  static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
    0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

  static int diagonals[14*2] = { 0, 0, //remove this to avoid self intesextion
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7,
    2, 3, 1, 3, 1, 2,                       
    1, 4, 2, 5, 3, 6 };

  for(int z = 0; z < z_side; z++) {
    for(int y = 0; y < y_side; y++) {
      for(int x = 0; x < x_side; x++) {
        int origin = cell(x, y, z);
        self(get_start_iter(origin), get_end_iter(origin), radius);   
        // compute between other girds
        for(int d = 2; d < 28; d += 2) { // skipping self
          int *cs = corner + 3*diagonals[d];
          int *ce = corner + 3*diagonals[d+1];
          if((x + cs[0] < x_side) && (y + cs[1] < y_side) && (z + cs[2] < z_side) 
          && (x + ce[0] < x_side) && (y + ce[1] < y_side) && (z + ce[2] < z_side)) 
          {
            origin = cell(x+cs[0], y+cs[1], z+cs[2]);
            int dest = cell(x+ce[0], y+ce[1], z+ce[2]);
            other(get_start_iter(origin), get_end_iter(origin), 
                  get_start_iter(dest),   get_end_iter(dest), radius);        
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
  static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
    0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

  static int diagonals[14*2] = { 0, 0, //remove this to avoid self intesextion
    0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7,
    2, 3, 1, 3, 1, 2,                       
    1, 4, 2, 5, 3, 6 };

  for(int z = 0; z < z_side; z++) {
    for(int y = 0; y < y_side; y++) {
      for(int x = 0; x < x_side; x++) {     
        int origin = cell(x, y, z);  


        if(!is_empty(origin) && !points.is_empty(origin)) 
          travel_others(get_start_iter(origin), get_end_iter(origin), 
          points.get_start_iter(origin), points.get_end_iter(origin), radius);  


        for(int d = 2; d < 28; d += 2) { //skipping self
          int *cs = corner + 3*diagonals[d];
          int *ce = corner + 3*diagonals[d+1];
          if((x+cs[0] < x_side) && (y+cs[1] < y_side) && (z+cs[2] < z_side) &&
            (x+ce[0] < x_side) && (y+ce[1] < y_side) && (z+ce[2] < z_side)) {

              origin   = cell(x+cs[0], y+cs[1], z+cs[2]);

              int dest = cell(x+ce[0], y+ce[1], z+ce[2]);

              if(!is_empty(origin) && !points.is_empty(dest))        
                  travel_others(get_start_iter(origin), 
                                get_end_iter(origin), 
                                points.get_start_iter(dest),   
                                points.get_end_iter(dest), 
                                radius); 

              if(!is_empty(dest) && !points.is_empty(origin))  
                  travel_others(get_start_iter(dest), 
                                get_end_iter(dest), 
                                points.get_start_iter(origin),   
                                points.get_end_iter(origin), 
                                radius);        
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

  iterator dest;
  for(dest = starta; dest != enda; dest++) 
  {
    Rich_point<Kernel> &v = *(*dest);

    Point &p = v.pt;

    iterator origin;
    for(origin = startb; origin != endb; origin++)
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
  std::vector<Rich_point<Kernel> >& points, ///< sample point set
  CGAL::Bbox_3 bbox, ///< bounding box
  const typename Kernel::FT radius)
{
  CGAL_point_set_processing_precondition(radius > 0);

  for (unsigned int i = 0; i < points.size(); ++i)
  {
    points[i].neighbors.clear();
  }

  Rich_grid<Kernel> points_grid;
  points_grid.init(points, bbox, radius);
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
  CGAL::Bbox_3 bbox, ///< bounding box
  const typename Kernel::FT radius ///< neighbor radius
)
{
  typedef typename Kernel::FT FT;
  if (radius < FT(0.0))
  {
    return;
  }

  for (unsigned int i = 0; i < samples.size(); ++i)
  {
    samples[i].original_neighbors.clear();
  }

  Rich_grid<Kernel> samples_grid;
  samples_grid.init(samples, bbox, radius);

  // here can be initial the original grid just one time ?
  Rich_grid<Kernel> original_grid;
  original_grid.init(original, bbox, radius);

  samples_grid.travel_others(original_grid, 
                             Rich_grid<Kernel>::find_original_neighbors);
}


} //namespace rich_grid_internal

/// \endcond

} //namespace CGAL

#endif // CGAL_RICH_GRID_H
