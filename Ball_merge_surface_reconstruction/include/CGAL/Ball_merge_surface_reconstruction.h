// Copyright (c) 2023 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: ( GPL-3.0-or-later OR LicenseRef-Commercial ) AND MIT
//
//
// Author(s)     :  Amal Dev Parakkat, Andreas Fabri, Sébastien Loriot
//
// This file incorporates work covered by the following copyright and permission notice:
//
//     MIT License
//
//     Copyright (c) 2023 Parakkat, Amal Dev and Ohrhallinger, Stefan and Eisemann, Elmar and Memari, Pooran
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy
//     of this software and associated documentation files (the "Software"), to deal
//     in the Software without restriction, including without limitation the rights
//     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//     copies of the Software, and to permit persons to whom the Software is
//     furnished to do so, subject to the following conditions:
//
//     The above copyright notice and this permission notice shall be included in all
//     copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//     SOFTWARE.
//
//
// @inproceedings{parakkat2024ballmerge,
//   title={BallMerge: High-quality Fast Surface Reconstruction via Voronoi Balls},
//   author={Parakkat, Amal Dev and Ohrhallinger, Stefan and Eisemann, Elmar and Memari, Pooran},
//   booktitle={Computer Graphics Forum},
//   year={2024}
// }
//


#ifndef CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H
#define CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H

#include <CGAL/license/Ball_merge_surface_reconstruction.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Container_helper.h>
#include <boost/iterator/zip_iterator.hpp>

#include <queue>

namespace CGAL {

namespace internal{

template <typename Traits, typename ConcurrencyTag>
class Ball_merge_surface_reconstruction {

  typedef CGAL::Triangulation_vertex_base_with_info_3<std::size_t, Traits> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<unsigned int, Traits, CGAL::Delaunay_triangulation_cell_base_3<Traits>> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb, ConcurrencyTag> Tds;
  typedef CGAL::Delaunay_triangulation_3<Traits, Tds> Delaunay;
  typedef typename Delaunay::Point Point;
  typedef typename Delaunay::Cell_handle Cell_handle;
  int group = 1, gcount = 0, maxg = 0, maingroup, max1 = 0;
  double par, bbdiaglen;
  Bbox_3 bbox;
  std::vector<int> cell_groups;

public:
  int secondgroup;
  int option = 0,eta=0;
  Delaunay dt3;

private:

  double distance(const Point& p1, const Point& p2) const
  { return std::sqrt(typename Traits::Compute_squared_distance_3()(p1,p2)); }

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
    if (distance(a, b) * eta > bbdiaglen || distance(c, b) * eta > bbdiaglen || distance(a, c) * eta > bbdiaglen)
      return 0;
    return 1;
  }

public:

  template <class PointRange>
  void operator()(const PointRange& points, double par_, double eta_)
  {
    dt3.clear();
    par = par_;
    eta = eta_;
    bbox = bbox_3(points.begin(), points.end());

    bbdiaglen = distance(Point(bbox.xmin(), bbox.ymin(), bbox.zmin()), Point(bbox.xmax(), bbox.ymax(), bbox.zmax()));

    std::vector<std::size_t> vids(points.size());
    std::iota(vids.begin(), vids.end(), 0);

    if constexpr (std::is_same<ConcurrencyTag, Parallel_tag>::value) {
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

  template <class TripleIndexRange>
  void set_triangle_indices_hull1(TripleIndexRange& meshFaceIndices) const
  {
    using TripleIndex = typename std::iterator_traits<typename TripleIndexRange::iterator>::value_type;
    std::vector<bool> visited(dt3.number_of_cells(),false);
    for (Cell_handle cell : dt3.finite_cell_handles()){
      for (int i = 0; i < 4; ++i){
        if (option == 1){
          //If global, write the cell details of the largest group to the PLY file
          if (cell_groups[cell->info()] == maingroup && (cell_groups[cell->neighbor(i)->info()] != maingroup || dt3.is_infinite(cell->neighbor(i)))){
            //Write the triangles between cells if they have have different labels and one of them is labeled as the same as the largest group
            TripleIndex indices;
            CGAL::internal::resize(indices, 3);
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
            if (dt3.is_infinite(cell->neighbor(i)) || !_function(cell, cell->neighbor(i))){
              //If the cells cannot be merged, then write the triangle between these two cells to the PLY file
              visited[cell->info()]=true;
              TripleIndex indices;
              CGAL::internal::resize(indices, 3);
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


  template <class TripleIndexRange>
  void set_triangle_indices_hull2(TripleIndexRange& meshFaceIndices) const
  {
    using TripleIndex = typename std::iterator_traits<typename TripleIndexRange::iterator>::value_type;
    for (Cell_handle cell : dt3.finite_cell_handles())
    {
      for (int i = 0; i < 4; ++i)
        if (cell_groups[cell->info()] == secondgroup && (dt3.is_infinite(cell->neighbor(i)) || cell_groups[cell->neighbor(i)->info()] != secondgroup))
        {
          //Write the triangles between cells if they have have different labels and one of them is labeled as the same as the second largest group
          TripleIndex indices;
          CGAL::internal::resize(indices, 3);
          for (int j = 0; j < 3; j++)
            indices[j] = cell->vertex((i + 1 + j) % 4)->info();
          if (i%2==1)
            std::swap(indices[0], indices[1]);
          meshFaceIndices.emplace_back(indices);
        }
    }
  }

};

} // end of internal namespace


/*!
 *
 * \ingroup PkgBallMergeRef
 *
 * creates a triangle soup approximating the surface.
 * After creation, resulting triangular faces will be stored in `out_triangles`, allowing us to use it independently for further processing.
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * \tparam Traits a model of `DelaunayTriangulationTraits_3`, with default using the value type of `PointRange` plugged in `Kernel_traits`
 * \tparam PointRange a model of `RandomAccessContainer`
 * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
 *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
 *                          being constructible from `std::size_t`.
 *
 * \param points is the input points
 * \param out_triangles is the output parameter storing triangles approximating the surface
 * \param delta is the value of \f$ \delta \f$
 * \param eta is the optional parameter \f$ \eta \f$
 *
 */
template <class Concurrency_tag = Sequential_tag, class Traits = Default, class PointRange, class TripleIndexRange>
void ball_merge_surface_reconstruction_local(const PointRange& points,
                                             TripleIndexRange& out_triangles,
                                             double delta, double eta=200.)
{
  using Point_3 = std::remove_const_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Traits_ = typename Default::Get<Traits,typename Kernel_traits<Point_3>::type>::type;

  CGAL::internal::Ball_merge_surface_reconstruction<Traits_, Concurrency_tag> bmsr;
  bmsr.option=0;
  bmsr(points, delta, eta);
  bmsr.set_triangle_indices_hull1(out_triangles);
}

/*! \ingroup PkgBallMergeRef
 *
 * creates a watertight mesh approximating the surface.
 * After creation, resulting meshes will be stored in `out_triangles1` and `out_triangles2`, allowing us to use it independently for further processing.
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * \tparam Traits a model of `DelaunayTriangulationTraits_3`, with default using the value type of `PointRange` plugged in `Kernel_traits`
 * @tparam PointRange a model of the concepts `RandomAccessContainer`
 * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
 *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
 *                          being constructible from `std::size_t`.
 *
 * \param points is the input points representing a single connected component of a watertight surface
 * \param out_triangles1 is the output parameter storing the first resulting mesh
 * \param out_triangles2 is the output parameter storing the second resulting mesh
 * \param delta is the value of \f$ \delta \f$
 */
template <class Concurrency_tag= Sequential_tag, class Traits = Default, class PointRange, class TripleIndexRange>
void ball_merge_surface_reconstruction_global(const PointRange& points,
                                              TripleIndexRange& out_triangles1,
                                              TripleIndexRange& out_triangles2,
                                              double delta)
{
  using Point_3 = std::remove_const_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Traits_ = typename Default::Get<Traits,typename Kernel_traits<Point_3>::type>::type;
  CGAL::internal::Ball_merge_surface_reconstruction<Traits_, Concurrency_tag> bmsr;
  bmsr.option=1;
  bmsr(points, delta, 0);
  bmsr.set_triangle_indices_hull1(out_triangles1);
  bmsr.set_triangle_indices_hull2(out_triangles2);
}


} // namespace CGAL



#endif // CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H
