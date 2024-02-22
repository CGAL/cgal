// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  NN

#ifndef CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H
#define CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H

// #include <memory>

// #include <CGAL/license/Ball_merge_surface_reconstruction.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <boost/iterator/zip_iterator.hpp>

#include <queue>

namespace CGAL {

template <typename Traits, typename ConcurrencyTag>
class Ball_merge_surface_reconstruction {

  typedef CGAL::Triangulation_vertex_base_with_info_3<int, Traits> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<unsigned int, Traits, CGAL::Delaunay_triangulation_cell_base_3<Traits>> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb, ConcurrencyTag> Tds;
  typedef CGAL::Delaunay_triangulation_3<Traits, Tds> Delaunay;
  typedef typename Delaunay::Point Point;
  typedef typename Delaunay::Cell_handle Cell_handle;
  static constexpr int NOT_VISITED=-1;
  int group = 1, gcount = 0, maxg = 0, maingroup, max1 = 0;
  double par, bbdiaglen;
  Bbox_3 bbox;
  std::vector<int> cell_groups;

public:
  int secondgroup;
  int option = 0,tlen=0;
  Delaunay dt3;

private:

  double distance(const Point& p1, const Point& p2) const
  { return sqrt(squared_distance(p1,p2)); }

  //Function that computes the Intersection Ratio (IR) given two cell handles
  // AF This could become a filtered predicate if robustness might be an issue
  // AF Find a better name and an enum for the return type
  bool _function(Cell_handle fh1, Cell_handle fh2) const
  {
    Point p1 = dt3.dual(fh1); //p1 is the circumcenter of the first cell
    Point p2 = dt3.dual(fh2); //p2 is the circumcenter of the second cell
    double d = distance(p1, p2); //distance between the circumcenters
    double a = distance(fh1->vertex(0)->point(), p1);//circumradius of circumsphere of the first cell
    double b = distance(fh2->vertex(0)->point(), p2);//circumradius of circumsphere of the second cell
    if ((a + b - d) > a * par || (a + b - d) > b * par)//Check whether the IR is less than a user given parameter
      return true;
    return false;
  }

  //Function to recursively group all the cells that are mergeable from the cell fh with user specified parameter
  void regroup(const Delaunay& dt3, Cell_handle fh)
  {
    std::queue<Cell_handle> q; //A queue to do the grouping
    q.push(fh); //The cell itself is mergeable
    while (!q.empty()){ //A traversal
      fh = q.front();
      q.pop();
      cell_groups[fh->info()] = group; //A label
      ++gcount;
      for (int i = 0; i < 4; ++i){ //For each neighbor
        if ((!dt3.is_infinite(fh->neighbor(i)))&&(cell_groups[fh->neighbor(i)->info()] == 0 && _function(fh, fh->neighbor(i)) == 1)){//If it is unlabeled and can be merged - using the function that computes IR          fh->neighbor(i)->info() = group;//If mergeable, then change the label
          cell_groups[fh->neighbor(i)->info()] = group;//If mergeable, then change the label
          q.push(fh->neighbor(i));//Push it into the traversal list
        }
      }
    }
  }

  //Function to check whether a triangle made by three points has an edge greater than bbdiag/(user defined paramter - 200 by default) - used for filtering long triangles in local ballmerge
  int bblen(const Point& a, const Point& b, const Point& c) const
  {
    if (distance(a, b) * tlen > bbdiaglen || distance(c, b) * tlen > bbdiaglen || distance(a, c) * tlen > bbdiaglen)
      return 0;
    return 1;
  }

public:

  void operator()(const std::vector<Point>& points, double par_, double tlen_)
  {
    dt3.clear();
    par = par_;
    tlen = tlen_;
    bbox = bbox_3(points.begin(), points.end());
 
    bbdiaglen = distance(Point(bbox.xmin(), bbox.ymin(), bbox.zmin()), Point(bbox.xmax(), bbox.ymax(), bbox.zmax()));

    std::vector<int> vids(points.size());
    std::iota(vids.begin(), vids.end(), 0);

    if constexpr (std::is_same_v<ConcurrencyTag, Parallel_tag>) {
      typename Delaunay::Lock_data_structure locking_ds(bbox, 50);
      dt3.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), vids.begin())),
                 boost::make_zip_iterator(boost::make_tuple(points.end(), vids.end())),
                 &locking_ds);
    } else {
      dt3.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), vids.begin())),
                 boost::make_zip_iterator(boost::make_tuple(points.end(), vids.end())));//Sequential Delaunay computation
    }
    unsigned int cell_id=0;
    for (Cell_handle ch : dt3.all_cell_handles())
    {
      ch->info()=cell_id++;
    }

    cell_groups.assign(dt3.number_of_cells(), 0);
    if (option == 1){//If the user opted for global algorithm
      for (Cell_handle cell : dt3.finite_cell_handles())
      {
        if (cell_groups[cell->info()] == 0){//If the cell label is not altered
          cell_groups[cell->info()] = group;//Assign a label
          regroup(dt3, cell);//Group all the mergeable cells
          if (maxg < gcount){
            //Remember the largest group - based on the number of tetrahedra in the group
            max1 = maxg;
            secondgroup = maingroup;//To remember the second largest group
            maxg = gcount;
            maingroup = group;
          }
          else if (max1 < gcount && gcount != maxg){
            max1 = gcount;
            secondgroup = group;//To remember the second largest group
          }
          gcount = 0;
          ++group;//Update the label for next cell
        }
      }
    }

  }

// void flood_from_infinity(std::vector<bool> &outside, int group) const
// {
//   std::vector<Cell_handle> queue;
//     queue.push_back(dt3.infinite_cell());
//     while(!queue.empty())
//     {
//       Cell_handle cell = queue.back();
//       queue.pop_back();
//       if (outside[cell->info()]) continue;
//       outside[cell->info()]=true;
//       for(int i=0;i<4;++i)
//       {
//         if(cell_groups[cell->neighbor(i)->info()] == group) continue;
//         if (!outside[cell->neighbor(i)->info()])
//           queue.push_back(cell->neighbor(i));
//       }
//     }
// }

void set_triangle_indices_hull1(std::vector<std::array<int,3>>& meshFaceIndices) const
{
  std::vector<bool> visited(dt3.number_of_cells(),false);
  for (Cell_handle cell : dt3.finite_cell_handles()){
    for (int i = 0; i < 4; ++i){
      if (option == 1){
        //If global, write the cell details of the largest group to the PLY file
        if (cell_groups[cell->info()] == maingroup && (cell_groups[cell->neighbor(i)->info()] != maingroup || dt3.is_infinite(cell->neighbor(i)))){
          //Write the triangles between cells if they have have different labels and one of them is labeled as the same as the largest group
          std::array<int, 3> indices;
          for (int j = 0; j < 3; ++j){
            indices[j] = cell->vertex((i + 1 + j) % 4)->info();
          }
          if (i%2==1)
              std::swap(indices[0], indices[1]);
          meshFaceIndices.push_back(indices);
        }
      }
      else if (!visited[cell->neighbor(i)->info()])
      {
        if (bblen(cell->vertex((i + 1) % 4)->point(), cell->vertex((i + 2) % 4)->point(), cell->vertex((i + 3) % 4)->point())){
          //If the triangle crosses our bbdiagonal based criteria
          if (dt3.is_infinite(cell->neighbor(i))||!_function(cell, cell->neighbor(i)) == 1){
            //If the cells cannot be merged, then write the triangle between these two cells to the PLY file
            visited[cell->info()]=true;
            std::array<int,3> indices;
            for (int j = 0; j < 3; ++j){
              indices[j] = cell->vertex((i + 1 + j) % 4)->info();
            }
            meshFaceIndices.push_back(indices);
          }
        }
      }
    }
  }
}

  void set_triangle_indices_hull2(std::vector<std::array<int,3>>& meshFaceIndices) const
  {
  for (Cell_handle cell : dt3.finite_cell_handles())
    {
      for (int i = 0; i < 4; ++i)
        if (cell_groups[cell->info()] == secondgroup && (dt3.is_infinite(cell->neighbor(i)) || cell_groups[cell->neighbor(i)->info()] != secondgroup))
        {
          //Write the triangles between cells if they have have different labels and one of them is labeled as the same as the second largest group
          std::array<int,3> indices;
          for (int j = 0; j < 3; j++)
            indices[j] = cell->vertex((i + 1 + j) % 4)->info();
          if (i%2==1)
            std::swap(indices[0], indices[1]);
          meshFaceIndices.push_back(indices);
        }
    }
  }
};

/// \ingroup PkgBallMergeRef
template <class Concurrency_tag, class Traits>
void ball_merge_surface_reconstruction_local(const std::vector<CGAL::Point_3<Traits>>& points,
                                             std::vector<std::array<int, 3> >& out_triangles,
                                             double parameter, double tlen=200.)
{
  CGAL::Ball_merge_surface_reconstruction<Traits, Concurrency_tag> bmsr;
  bmsr.option=0;
  bmsr(points, parameter, tlen);
  bmsr.set_triangle_indices_hull1(out_triangles);
}

/// \ingroup PkgBallMergeRef
template <class Concurrency_tag, class Traits>
void ball_merge_surface_reconstruction_global(const std::vector<CGAL::Point_3<Traits>>& points,
                                              std::vector<std::array<int, 3> >& out_triangles1,
                                              std::vector<std::array<int, 3> >& out_triangles2,
                                              double parameter)
{
  CGAL::Ball_merge_surface_reconstruction<Traits, Concurrency_tag> bmsr;
  bmsr.option=1;
  bmsr(points, parameter, 0);
  bmsr.set_triangle_indices_hull1(out_triangles1);
  bmsr.set_triangle_indices_hull2(out_triangles2);
}


} // namespace CGAL



#endif // CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H
