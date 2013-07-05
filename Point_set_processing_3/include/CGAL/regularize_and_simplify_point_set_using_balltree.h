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
// Author(s) : Shihao Wu, Clement Jamin 

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_BALL_TREE_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_BALL_TREE_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>
#include <ctime>


//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>

/// \cond SKIP_IN_MANUAL
// ----------------------------------------------------------------------------
// Testing code section
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Rich Grid section
// ----------------------------------------------------------------------------
//namespace rich_grid_internel
//{
//}

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace regularize_and_simplify_internal{

  class Timer
  {
  public:

    void start(const std::string& str)
    {
      std::cout << std::endl;
      starttime = clock();
      //std::cout << "@@@@@ Time Count Strat For: " << str << std::endl;
      _str = str;
    }

    void end()
    {
      stoptime = clock();
      timeused = stoptime - starttime;
      std::cout << /*endl <<*/ "@@@@ finish	" << _str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << std::endl;
      std::cout << std::endl;
    }

  private:
    int starttime, mid_start, mid_end, stoptime, timeused;
    std::string _str;
  };


  // ----------------------------------------------------------------------------
  // Ball Tree section
  // ----------------------------------------------------------------------------
  template <typename Kernel>
  class Rich_point
  {
  public:
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;

  public:
    Rich_point(){}
    Rich_point(const Point& p, const int& i):pt(p),index(i){} 

    Point pt;
    unsigned int index;
    std::vector<unsigned int> neighbors;
    std::vector<unsigned int> original_neighbors;//need more memory, but make the code a littel easier to read.
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
      FT max_double = std::numeric_limits<FT>::max();
      min_x = min_y = min_z = max_double;
      max_x = max_y = max_z = -max_double;
    }

    void add_point(const Point& p)
    {
      if (p[0] < min_x) min_x = p[0];
      if (p[1] < min_y) min_y = p[1];
      if (p[2] < min_z) min_z = p[2];
      if (p[0] > max_x) max_x = p[0];
      if (p[1] > max_y) max_y = p[1];
      if (p[2] > max_z) max_z = p[2];
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

    void init(std::vector<Rich_point<Kernel>> &vert, Rich_box<Kernel>& box, const FT _radius); 

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
    void static  __cdecl find_self_neighbors(iterator start, iterator end, FT radius);
    void static  __cdecl find_other_neighbors(iterator starta, iterator enda, 
      iterator startb, iterator endb, FT radius);

  private:
  
    std::vector<Rich_point<Kernel>*> rich_points;  
    std::vector<int> index;    // the start index of each grid in the sample points which is order by Zsort
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
      return a->pt[0] < b->pt[0];
    }
  };

  template <typename Kernel>
  class YSort {
  public:
    bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
      return a->pt[1] < b->pt[1];
    }
  };

  template <typename Kernel>
  class ZSort {
  public:
    bool operator()(const Rich_point<Kernel> *a, const Rich_point<Kernel> *b) {
      return a->pt[2] < b->pt[2];
    }
  };


  // divide spoints into some grids
  // and each grid has their points index in the index vector of sample.
  template <typename Kernel>
  void Rich_grid<Kernel>::init(std::vector<Rich_point<Kernel>> &vert, Rich_box<Kernel>& box, const typename Kernel::FT _radius) 
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

    std::sort(rich_points.begin(), rich_points.end(), ZSort<Kernel>()); //this would be very slow 

    unsigned int startz = 0;
    for(unsigned int z = 0; z < zside; z++) 
    {
      int endz = startz;
      FT maxz = min.z() + (z+1)*radius;
      while(endz < rich_points.size() && rich_points[endz]->pt.z()< maxz)
        ++endz; 

      sort(rich_points.begin()+startz, rich_points.begin()+endz, YSort<Kernel>());

      int starty = startz;
      for(int y = 0; y < yside; y++) 
      {
        int endy = starty;        
        FT maxy = min.y() + (y+1)*radius;
        while(endy < endz && rich_points[endy]->pt.y() < maxy)
          ++endy;

        sort(rich_points.begin()+starty, rich_points.begin()+endy, XSort<Kernel>());

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
    index[xside*yside*zside] = startz;  // in order to compute the last grid's range
  }

  /// define how to travel in the same gird 
  template <typename Kernel>
  void Rich_grid<Kernel>::travel_itself(void (*self)(iterator starta, iterator enda, const typename Kernel::FT radius),
    void (*other)(iterator starta, iterator enda, 
    iterator startb, iterator endb, FT radius)) 
  {
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;

    static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
      0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

    static int diagonals[14*2] = { 0, 0, //remove this line to avoid self intesextion
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
            if((x + cs[0] < xside) && (y + cs[1] < yside) && (z + cs[2] < zside) &&
              (x + ce[0] < xside) && (y + ce[1] < yside) && (z + ce[2] < zside)) 
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
  void Rich_grid<Kernel>::travel_others(Rich_grid &points, 
    void (*travel_others)(iterator starta, iterator enda, 
    iterator startb, iterator endb, const typename Kernel::FT radius)) 
  {
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;

    static int corner[8*3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
      0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

    static int diagonals[14*2] = { 0, 0, //remove this line to avoid self intesextion
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

                if(!isEmpty(origin) && !points.isEmpty(dest))           // Locally 
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
    iterator endb, FT radius)
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
  void Rich_grid<Kernel>::find_self_neighbors(iterator start, iterator end, FT radius)
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
  void Rich_grid<Kernel>::find_other_neighbors(iterator starta, iterator enda, 
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




  /// Compute ball neighbors for each point in the same point set(sample or original points) 
  /// 
  /// \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return 
  template <typename Kernel>
  void compute_ball_neighbors_one_self(
    std::vector<Rich_point<Kernel>>& points,
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
    points_grid.travel_itself(Rich_grid<Kernel>::find_self_neighbors, Rich_grid<Kernel>::find_other_neighbors);
  }

  /// Compute ball neighbors for each point(sample points) in the other point set(original points) 
  /// 
  /// \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return 
  template <typename Kernel>
  void compute_ball_neighbors_one_to_another(
    std::vector<Rich_point<Kernel>>& samples, ///< sample point set
    std::vector<Rich_point<Kernel>>& original,///< original point set
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

    // can be speed up by initial the original grid just one time(as paramter of this function)
    Rich_grid<Kernel> original_grid;
    original_grid.init(original, box, radius);

    samples_grid.travel_others(original_grid, Rich_grid<Kernel>::find_original_neighbors);
  }

  // ----------------------------------------------------------------------------
  // WLOP algorithm section
  // ----------------------------------------------------------------------------
  /// Compute average term for each sample points
  /// According to their ball neighborhood original points
  /// 
  /// \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return average term vector
  template <typename Kernel>
    typename Kernel::Vector_3
    compute_average_term(
    const typename Kernel::Point_3& query, ///< 3D point to project
    const std::vector<Rich_point<Kernel>> neighbor_original_points, ///< neighbor original points
    const typename Kernel::FT radius, ///<accept neighborhood radius
    const std::vector<typename Kernel::FT>& density_weight_set ///<if user need density
    )
  {
    CGAL_point_set_processing_precondition(radius > 0);
    CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);
    bool is_density_weight_set_empty = density_weight_set.empty();

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    FT radius2 = radius * radius;

    //Compute average term
    FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
    FT iradius16 = -(FT)4.0/radius2;
    Vector average = CGAL::NULL_VECTOR; 
    for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
    {
      const Point& np = neighbor_original_points[i].pt;
      unsigned int idx_of_original = neighbor_original_points[i].index;

      FT dist2 = CGAL::squared_distance(query, np);
      weight = exp(dist2 * iradius16);

      if(!is_density_weight_set_empty)
      {
        weight *= density_weight_set[idx_of_original];
      }

      average_weight_sum += weight;
      average = average + (np - CGAL::ORIGIN) * weight;
    }

    // output
    return average/average_weight_sum;
  }

  /// Compute repulsion term for each sample points
  /// According to their Ball neighborhood sample points
  /// 
  /// \pre \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return repulsion term vector
  template <typename Kernel>
    typename Kernel::Vector_3
    compute_repulsion_term(
    const typename Kernel::Point_3& query, ///< 3D point to project
    const std::vector<Rich_point<Kernel>> neighbor_sample_points, ///< neighbor sample points
    const typename Kernel::FT radius, ///<accept neighborhood radius
    const std::vector<typename Kernel::FT>& density_weight_set ///< if user need density
    )
  {
    CGAL_point_set_processing_precondition(radius > 0);
    CGAL_point_set_processing_precondition(neighbor_sample_points.size() >= 1);

    bool is_density_weight_set_empty = density_weight_set.empty();

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    FT radius2 = radius * radius;

    //Compute average term
    FT weight = (FT)0.0, repulsion_weight_sum = (FT)0.0;
    FT iradius16 = -(FT)4.0/radius2;

    Vector repulsion = CGAL::NULL_VECTOR; 
    for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
    {
      const Point& np = neighbor_sample_points[i].pt;
      unsigned int idx_of_sample = neighbor_sample_points[i].index;

      Vector diff = query - np;

      FT dist2 = CGAL::squared_distance(query, np);
      FT dist = std::sqrt(dist2);

      weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0)/dist, 2);
      if(!is_density_weight_set_empty)
      {
        weight *= density_weight_set[idx_of_sample];
      }

      repulsion_weight_sum += weight;
      repulsion = repulsion + diff * weight;
    }

    // output
    return repulsion/repulsion_weight_sum;
  }




  /// Compute density weight for each original points,
  /// according to their ball neighborhood original points
  /// 
  /// \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return density weight
  template <typename Kernel>
    typename Kernel::FT
    compute_density_weight_for_original_point(
    const typename Kernel::Point_3& query, ///< 3D point to project
    const std::vector<typename Kernel::Point_3> neighbor_original_points, ///< neighbor original points
    const typename Kernel::FT radius
    )
  {
    CGAL_point_set_processing_precondition(radius > 0);

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    //Compute density weight
    FT radius2 = radius * radius;
    FT density_weight = (FT)1.0;
    FT iradius16 = -(FT)4.0/radius2;

    for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
    {
      const Point& np = neighbor_original_points[i];
      FT dist2 = CGAL::squared_distance(query, np);
      density_weight += std::exp(dist2 * iradius16);
    }

    // output
    return FT(1.0) / density_weight;
  }


  /// Compute density weight for sample point,
  /// according to their ball neighborhood sample points
  /// 
  /// \pre `radius > 0`
  ///
  /// @tparam Kernel Geometric traits class.
  ///
  /// @return density weight
  template <typename Kernel>
    typename Kernel::FT
    compute_density_weight_for_sample_point(
    const typename Kernel::Point_3& query, ///< 3D point to project
    const std::vector<typename Kernel::Point_3> neighbor_sample_points, ///< neighbor sample points
    const typename Kernel::FT radius
    )
  {
    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    //Compute density weight
    FT radius2 = radius * radius;
    FT density_weight = (FT)1.0;
    FT iradius16 = -(FT)4.0/radius2;

    for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
    {
      const Point& np = neighbor_sample_points[i];
      FT dist2 = CGAL::squared_distance(query, np);
      density_weight += std::exp(dist2 * iradius16);
    }

    // output
    //return std::sqrt(density_weight); 
    return density_weight;
  }
}



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------
namespace CGAL {

  /// \ingroup PkgPointSetProcessing
  /// WLOP Algorithm: The simplification algorithm can produces a set of 
  ///	denoised, outlier-free and evenly distributed particles over the original dense point cloud,
  ///	so as to improve the reliability of current normal orientation algorithm. 
  ///	The core of the algorithm is a Weighted Locally Optimal projection operator with a density uniformization term. 
  /// More deatail please see: http://web.siat.ac.cn/~huihuang/WLOP/WLOP_page.html
  ///
  /// @tparam ForwardIterator iterator over input points.
  /// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
  ///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
  /// @tparam Kernel Geometric traits class.
  ///        It can be omitted and deduced automatically from PointPMap value_type.
  ///
  /// @return iterator of the first point to downsampled points.

  // This variant requires all parameters.
  template <typename ForwardIterator,
    typename PointPMap,
    typename Kernel
  >
  ForwardIterator
  regularize_and_simplify_point_set_using_balltree(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  const typename Kernel::FT neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const bool need_compute_density, ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  const Kernel& /*kernel*/) ///< geometric traits.
  {
    CGAL_point_set_processing_precondition(neighbor_radius > 0);
    regularize_and_simplify_internal::Timer time;

    // basic geometric types
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;
    typedef typename regularize_and_simplify_internal::Rich_point<Kernel> Rich_point;

    // precondition: at least one element in the container.
    // to fix: should have at least three distinct points
    // but this is costly to check
    CGAL_point_set_processing_precondition(first != beyond);
    CGAL_point_set_processing_precondition(retain_percentage >= 0 && retain_percentage <= 100);

    // Random shuffle
    std::random_shuffle (first, beyond);

    // Computes original(input) and sample points size 
    std::size_t nb_points_original = std::distance(first, beyond);
    std::size_t nb_points_sample = (std::size_t)(FT(nb_points_original) * (retain_percentage/100.0));
    std::size_t first_index_to_sample = nb_points_original - nb_points_sample;

    // The first point iter of original and sample points
    ForwardIterator it;// point iterator
    ForwardIterator first_original_point = first;
    ForwardIterator first_sample_point = first;
    std::advance(first_sample_point, first_index_to_sample);

    //Copy sample points
    std::vector<Point> sample_points(nb_points_sample);
    unsigned int i; // sample point index
    for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
      sample_points[i] = get(point_pmap, it);

    //Copy original points(Maybe not the best choice)
    std::vector<Point> original_points(nb_points_original);
    for(it = first_original_point, i = 0; it != beyond; ++it, i++)
      original_points[i] = get(point_pmap, it);


    // Initilization
    std::vector<Rich_point> original_rich_point_set(nb_points_original);
    std::vector<Rich_point> sample_rich_point_set(nb_points_sample);

    regularize_and_simplify_internal::Rich_box<Kernel> box;
    for (it = first_original_point, i = 0; it != beyond ; ++it, i++)
    {  
      Point& p0 = get(point_pmap,it);
      Rich_point rp(p0, i);
      original_rich_point_set[i] = rp;
      box.add_point(rp.pt);
    }

    // Compute original density weight for original points if user needed
    std::vector<FT> original_density_weight_set;
    if (need_compute_density)
    {
      time.start("Buile Ball Tree For Original");
      regularize_and_simplify_internal::compute_ball_neighbors_one_self(original_rich_point_set, box, neighbor_radius);
      time.end();

      time.start("Compute Density For Original");
      for (it = first_original_point, i = 0; it != beyond ; ++it, i++)
      {
        std::vector<Point> original_neighbors;
        std::vector<unsigned int>& neighors_indexes = original_rich_point_set[i].neighbors;
        for (unsigned int j = 0; j < neighors_indexes.size(); j++)
        {
          original_neighbors.push_back(original_points[neighors_indexes[j]]);
        }

        FT density = regularize_and_simplify_internal::compute_density_weight_for_original_point<Kernel>(get(point_pmap,it),
                     original_neighbors, neighbor_radius);

        original_density_weight_set.push_back(density);
        original_rich_point_set[i].neighbors.clear();
      }
      time.end();
    }


    for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
    {
      // Build Ball Tree For Sample Neighbor
      time.start("Build Ball Tree For Sample Neighbor");
      for (i=0 ; i < sample_points.size(); i++)
      {
        Point& p0 = sample_points[i];
        Rich_point rp(p0, i);
        sample_rich_point_set[i] = rp;
      }
      regularize_and_simplify_internal::compute_ball_neighbors_one_self(sample_rich_point_set, box, neighbor_radius);
      time.end();

      // Compute sample density weight for sample points if user needed
      std::vector<FT> sample_density_weight_set;
      time.start("Compute Density For Sample");
      if (need_compute_density)
      {
        for (i=0 ; i < sample_points.size(); i++)
        {
          std::vector<Point> sample_neighbors;
          std::vector<unsigned int>& neighors_indexes = sample_rich_point_set[i].neighbors;
          for (unsigned int j = 0; j < neighors_indexes.size(); j++)
          {
            sample_neighbors.push_back(sample_points[neighors_indexes[j]]);
          }

          FT density = regularize_and_simplify_internal::compute_density_weight_for_sample_point<Kernel>
                       (sample_points[i], sample_neighbors, neighbor_radius);
          sample_density_weight_set.push_back(density);
        }
      }
      time.end();

      // Build Ball Tree For Sample-Original Neighbor
      time.start("Build Ball Tree For Sample-Original Neighbor");
      regularize_and_simplify_internal::compute_ball_neighbors_one_to_another(sample_rich_point_set,
                                                                original_rich_point_set, box, neighbor_radius);
      time.end();

      // Compute average term and repulsion term for each sample points separately,
      // then update each sample points
      std::vector<Vector> average_set(nb_points_sample);
      std::vector<Vector> repulsion_set(nb_points_sample);
      time.start("Compute Average Term");
      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        std::vector<Rich_point> rich_original_neighbors;
        std::vector<unsigned int>& neighors_indexes = sample_rich_point_set[i].original_neighbors;

        if (neighors_indexes.empty())
        {
          average_set[i] = p - CGAL::ORIGIN;
          continue;
        }

        for (unsigned int j = 0; j < neighors_indexes.size(); j++)
        {
          unsigned int idx_of_original = neighors_indexes[j];
          Rich_point rp(original_points[idx_of_original], idx_of_original);
          rich_original_neighbors.push_back(rp);
        }

        average_set[i] = regularize_and_simplify_internal::compute_average_term<Kernel>(p, rich_original_neighbors, neighbor_radius, original_density_weight_set); 
      }
      time.end();

      time.start("Compute Repulsion Term");
      for (i = 0; i < sample_points.size(); i++)
      {
        std::vector<Rich_point> rich_sample_neighbors;
        std::vector<unsigned int>& neighors_indexes = sample_rich_point_set[i].neighbors;
        
        if (neighors_indexes.empty())
        {
          repulsion_set[i] = CGAL::NULL_VECTOR;
          continue;
        }

        for (unsigned int j = 0; j < neighors_indexes.size(); j++)
        {
          unsigned int idx_of_sample = neighors_indexes[j];
          Rich_point rp(sample_points[idx_of_sample], idx_of_sample);
          rich_sample_neighbors.push_back(rp);
        }

        Point& p = sample_points[i];
        repulsion_set[i] = regularize_and_simplify_internal::compute_repulsion_term<Kernel>(p, rich_sample_neighbors, neighbor_radius, sample_density_weight_set);
      }
      time.end();

      for (i = 0; i < sample_points.size(); i++)
      {
        Point& p = sample_points[i];
        p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
      }

      std::cout << "iterate:	" << iter_n + 1 <<  "	"<< std::endl;
    }

    //Copy back modified sample points to original points for output
    for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
    {
      Point& original_p = get(point_pmap, it);
      const Point& sample_p = sample_points[i];
      original_p = sample_p;
    }

    return first_sample_point;
  }

  /// @cond SKIP_IN_MANUAL
  // This variant deduces the kernel from the iterator type.
  template <typename ForwardIterator,
    typename PointPMap
  >
  ForwardIterator
  regularize_and_simplify_point_set_using_balltree(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density  ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  ) 
  {
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    return regularize_and_simplify_point_set_using_balltree(
      first,beyond,
      point_pmap,
      retain_percentage,
      neighbor_radius,
      iter_number,
      need_compute_density,
      Kernel());
  }
  /// @endcond

  /// @cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Dereference_property_map.
  template <typename ForwardIterator 
  >
  ForwardIterator
  regularize_and_simplify_point_set_using_balltree(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double retain_percentage, ///< percentage of points to retain
  double neighbor_radius, ///< size of neighbors.
  const unsigned int iter_number, ///< number of iterations.
  const bool need_compute_density ///< if needed to compute density to generate more rugularized result, 
  ///  especially when the density of input is uneven.
  ) 
  {
    return regularize_and_simplify_point_set_using_balltree(
      first,beyond,
      make_dereference_property_map(first),
      retain_percentage, neighbor_radius, iter_number, need_compute_density);
  }
  /// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
