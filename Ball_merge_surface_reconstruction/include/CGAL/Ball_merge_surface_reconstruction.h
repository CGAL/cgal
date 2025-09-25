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
// The original code from the publication mentioned below has been used as the base for the CGAL package,
// the original version is currently not publicly available.
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
#include <CGAL/Named_function_parameters.h>

#include <queue>

namespace CGAL {

/*!
 *
 * \ingroup PkgBallMergeRef
 *
 * Class that can be used for running the ball merge surface reconstruction algorithms. Once input points are passed,
 * it is possible to run a reconstruction algorithm with different parameters without rebuilding the internal Delaunay triangulation.
 *
 * \tparam Traits a model of `DelaunayTriangulationTraits_3`, with default using the value type of `PointRange` plugged in `Kernel_traits`
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 */
template <typename Traits, typename ConcurrencyTag>
class Ball_merge_surface_reconstruction
{
  enum Rec_option { LOCAL=0, GLOBAL=1};

#ifndef DOXYGEN_RUNNING
  typedef CGAL::Triangulation_vertex_base_with_info_3<std::size_t, Traits> Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<unsigned int, Traits, CGAL::Delaunay_triangulation_cell_base_3<Traits>> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb, ConcurrencyTag> Tds;
  typedef CGAL::Delaunay_triangulation_3<Traits, Tds> Delaunay;
  typedef typename Delaunay::Point Point;
  typedef typename Delaunay::Cell_handle Cell_handle;
  int m_group = 1, m_gcount = 0, m_maingroup=0, m_secondgroup=0;
  double m_delta, m_eta, m_bbdiaglen;

  std::vector<int> m_cell_groups;
  Rec_option m_option = GLOBAL;

  Delaunay m_dt3;

  double distance(const Point& p1, const Point& p2) const
  { return std::sqrt(typename Traits::Compute_squared_distance_3()(p1,p2)); }

  //Function that computes the Intersection Ratio (IR) given two cell handles
  // AF This could become a filtered predicate if robustness might be an issue
  // AF Find a better name and an enum for the return type
  bool _function(Cell_handle fh1, Cell_handle fh2) const
  {
    Point p1 = m_dt3.dual(fh1); //p1 is the circumcenter of the first cell
    Point p2 = m_dt3.dual(fh2); //p2 is the circumcenter of the second cell
    double d = distance(p1, p2); //distance between the circumcenters
    double a = distance(fh1->vertex(0)->point(), p1);//circumradius of circumsphere of the first cell
    double b = distance(fh2->vertex(0)->point(), p2);//circumradius of circumsphere of the second cell
    if ((a + b - d) > a * m_delta || (a + b - d) > b * m_delta)//Check whether the IR is less than a user given parameter
      return true;
    return false;
  }

  //Function to recursively group all the cells that are mergeable from the cell fh with user specified parameter
  void regroup(Cell_handle fh)
  {
    std::queue<Cell_handle> q; //A queue to do the grouping
    q.push(fh); //The cell itself is mergeable
    while (!q.empty()){ //A traversal
      fh = q.front();
      q.pop();
      m_cell_groups[fh->info()] = m_group; //A label
      ++m_gcount;
      for (int i = 0; i < 4; ++i){ //For each neighbor
        if ((!m_dt3.is_infinite(fh->neighbor(i)))&&(m_cell_groups[fh->neighbor(i)->info()] == 0 && _function(fh, fh->neighbor(i)) == 1)){
          //If it is unlabeled and can be merged - using the function that computes IR
          m_cell_groups[fh->neighbor(i)->info()] = m_group;//If mergeable, then change the label
          q.push(fh->neighbor(i));//Push it into the traversal list
        }
      }
    }
  }

  //Function to check whether a triangle made by three points has an edge greater than bbdiag/(user defined paramter - 200 by default) - used for filtering long triangles in local ballmerge
  int bblen(const Point& a, const Point& b, const Point& c) const
  {
    if (distance(a, b) * m_eta > m_bbdiaglen || distance(c, b) * m_eta > m_bbdiaglen || distance(a, c) * m_eta > m_bbdiaglen)
      return 0;
    return 1;
  }

public:


  void set_parameters(double delta_, double eta_, Rec_option option)
  {
    m_delta = delta_;
    m_eta = eta_;
    m_option = option;
  }

  void run_reconstruction()
  {
    m_group = 1;
    m_gcount = 0;
    m_maingroup=-1;
    m_secondgroup=-1;

    int maxg = 0;
    int max1 = 0;

    m_cell_groups.clear();
    m_cell_groups.assign(m_dt3.number_of_cells(), 0);
    if (m_option == GLOBAL){//If the user opted for global algorithm
      for (Cell_handle cell : m_dt3.finite_cell_handles())
      {
        if (m_cell_groups[cell->info()] == 0){//If the cell label is not altered
          m_cell_groups[cell->info()] = m_group;//Assign a label
          regroup(cell);//Group all the mergeable cells
          if (maxg < m_gcount){
            //Remember the largest group - based on the number of tetrahedra in the group
            max1 = maxg;
            m_secondgroup = m_maingroup;//To remember the second largest group
            maxg = m_gcount;
            m_maingroup = m_group;
          }
          else if (max1 < m_gcount && m_gcount != maxg){
            max1 = m_gcount;
            m_secondgroup = m_group;//To remember the second largest group
          }
          m_gcount = 0;
          ++m_group;//Update the label for next cell
        }
      }
    }


  }

  template <class TripleIndexRange>
  void set_triangle_indices_hull1(TripleIndexRange& meshFaceIndices) const
  {
    using TripleIndex = typename std::iterator_traits<typename TripleIndexRange::iterator>::value_type;

    if (m_option == GLOBAL)
    {
      for (Cell_handle cell : m_dt3.finite_cell_handles()){
        for (int i = 0; i < 4; ++i){
          //If global, write the cell details of the largest group to the PLY file
          if (m_cell_groups[cell->info()] == m_maingroup && (m_cell_groups[cell->neighbor(i)->info()] != m_maingroup || m_dt3.is_infinite(cell->neighbor(i)))){
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
      }
    }
    else
    {
      std::vector<bool> visited(m_dt3.number_of_cells(),false);
      for (Cell_handle cell : m_dt3.finite_cell_handles()){
        for (int i = 0; i < 4; ++i){
          if (!visited[cell->neighbor(i)->info()])
          {
            if (bblen(cell->vertex((i + 1) % 4)->point(), cell->vertex((i + 2) % 4)->point(), cell->vertex((i + 3) % 4)->point())){
              //If the triangle crosses our bbdiagonal based criteria
              if (m_dt3.is_infinite(cell->neighbor(i)) || !_function(cell, cell->neighbor(i))){
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
  }


  template <class TripleIndexRange>
  void set_triangle_indices_hull2(TripleIndexRange& meshFaceIndices) const
  {
    CGAL_assertion(m_option==GLOBAL);

    using TripleIndex = typename std::iterator_traits<typename TripleIndexRange::iterator>::value_type;
    for (Cell_handle cell : m_dt3.finite_cell_handles())
    {
      for (int i = 0; i < 4; ++i)
        if (m_cell_groups[cell->info()] == m_secondgroup && (m_dt3.is_infinite(cell->neighbor(i)) || m_cell_groups[cell->neighbor(i)->info()] != m_secondgroup))
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
#endif
public:
  /*!
   * default constructor
   */
  Ball_merge_surface_reconstruction()
  {}

  /*!
   * constructs the class and calls \link build_triangulation `build_triangulation(points)` \endlink.
   *
   * \tparam PointRange a model of `RandomAccessContainer`, with `Traits::Point_3` as value type
   * \param points is the input points
   */
  template <class PointRange>
  Ball_merge_surface_reconstruction(const PointRange& points);


  /*!
   * sets the input points for the triangulation, and builds the internal triangulation.
   * If called several times, only the last point range will be considered for the reconstructions.
   *
   * \tparam PointRange a model of `RandomAccessContainer`, with `Traits::Point_3` as value type
   * \param points is the input points
   */
  template <class PointRange>
  void build_triangulation(const PointRange& points)
  {
    m_dt3.clear();

    Bbox_3 bbox = bbox_3(points.begin(), points.end());

    m_bbdiaglen = distance(Point(bbox.xmin(), bbox.ymin(), bbox.zmin()), Point(bbox.xmax(), bbox.ymax(), bbox.zmax()));

    std::vector<std::size_t> vids(points.size());
    std::iota(vids.begin(), vids.end(), 0);

    if constexpr (std::is_same<ConcurrencyTag, Parallel_tag>::value) {
      typename Delaunay::Lock_data_structure locking_ds(bbox, 50);
      m_dt3.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), vids.begin())),
                 boost::make_zip_iterator(boost::make_tuple(points.end(), vids.end())),
                 &locking_ds);
    } else {
      m_dt3.insert(boost::make_zip_iterator(boost::make_tuple(points.begin(), vids.begin())),
                 boost::make_zip_iterator(boost::make_tuple(points.end(), vids.end())));//Sequential Delaunay computation
    }
    unsigned int cell_id=0;
    for (Cell_handle ch : m_dt3.all_cell_handles())
    {
      ch->info()=cell_id++;
    }
  }

#ifndef DOXYGEN_RUNNING
  /*!
   *
   * creates a triangle soup approximating the surface with sample points passed to `build_triangulation()`,
   * and puts the resulting triangle faces in `out_triangles`.
   *
   * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
   *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
   *                          being constructible from `std::size_t`.
   *
   * \param out_triangles is the output parameter storing triangles approximating the surface
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \cgalNamedParamsBegin
   *
   *   \cgalParamNBegin{delta}
   *     \cgalParamDescription{the value of \f$ \delta \f$}
   *     \cgalParamType{double}
   *     \cgalParamDefault{1.7}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{eta}
   *     \cgalParamDescription{the value of \f$ \eta \f$}
   *     \cgalParamType{double}
   *     \cgalParamDefault{200.}
   *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
   *   \cgalParamNEnd
   *
   * \cgalNamedParamsEnd
   *
   */
  template <class NamedParameters = parameters::Default_named_parameters,
            class TripleIndexRange>
  void local_reconstruction(TripleIndexRange& out_triangles,
                            const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    if (m_dt3.dimension()!=3) return;

    const double delta = choose_parameter(get_parameter(np, internal_np::delta), 1.7);
    const double eta = choose_parameter(get_parameter(np, internal_np::eta), 200.);

    set_parameters(delta, eta, LOCAL);
    set_triangle_indices_hull1(out_triangles);
  }
#endif

  /*!
   *
   * creates two watertight meshes approximating the surface with sample points passed to `build_triangulation()`,
   * and puts the resulting triangle faces in `out_triangles1` and `out_triangles2`. Output triangle faces are triple of indices
   * refering to the position of the input points passed to `build_triangulation()`.
   * \note As this function creates two shells (outer and inner in arbitrary order) the input point set must only be sampled on a single connected component.
   *
   * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
   *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
   *                          being constructible from `std::size_t`.
   *
   * \param out_triangles1 is the output parameter storing the first resulting mesh
   * \param out_triangles2 is the output parameter storing the second resulting mesh
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \return the delta value used for the reconstruction
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{delta}
   *     \cgalParamDescription{the value of \f$ \delta \f$}
   *     \cgalParamType{double}
   *     \cgalParamDefault{Value iteratively estimated, as described in \ref secBMSRPar}
   *   \cgalParamNEnd
   *
   * \cgalNamedParamsEnd
   *
   * \pre `build_triangulation()` should have been called before calling this function.
   */
  template <class NamedParameters = parameters::Default_named_parameters,
            class TripleIndexRange>
  double reconstruction(TripleIndexRange& out_triangles1,
                        TripleIndexRange& out_triangles2,
                        const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    if (m_dt3.dimension()!=3) return 0;

    if constexpr (!is_default_parameter<NamedParameters, internal_np::delta_t>::value)
    {
      const double delta = get_parameter(np, internal_np::delta);

      set_parameters(delta, 0, GLOBAL);
      run_reconstruction();
      set_triangle_indices_hull1(out_triangles1);
      set_triangle_indices_hull2(out_triangles2);

      return delta;
    }
    else
    {
      double delta = 2.;
      std::size_t nb_pt_prev=0;
      bool first=true;

      do
      {
        set_parameters(delta, 0, GLOBAL);
        run_reconstruction();

        std::vector<std::array<std::size_t, 3> > out_triangles;
        set_triangle_indices_hull1(out_triangles);

        std::set<std::size_t> indices_used;
        for (auto t : out_triangles)
        {
          indices_used.insert(t[0]);
          indices_used.insert(t[1]);
          indices_used.insert(t[2]);
        }

        std::size_t nb_pt=indices_used.size();

        if (!first && nb_pt_prev > nb_pt &&  (nb_pt_prev-nb_pt)  < 0.1*m_dt3.number_of_vertices())
        {
          delta+=0.01;
          break;
        }

        delta-=0.01;

        nb_pt_prev=nb_pt;
        first = false;
      }
      while(delta>1);

      set_triangle_indices_hull1(out_triangles1);
      set_triangle_indices_hull2(out_triangles2);

      return delta;
    }
  }

};

#ifndef DOXYGEN_RUNNING
/*!
 *
 * \ingroup PkgBallMergeRef
 *
 * creates a triangle soup approximating the surface, and puts the resulting triangle faces in `out_triangles`.
 * Output triangle faces are triple of indices refering to the position of the input points in `points`.
 *
 * \tparam Traits a model of `DelaunayTriangulationTraits_3`, with default using the value type of `PointRange` plugged in `Kernel_traits`
 * \tparam PointRange a model of `RandomAccessContainer`, with `Traits::Point_3` as value type
 * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
 *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
 *                          being constructible from `std::size_t`.
 *
 * \param points is the input points
 * \param out_triangles is the output parameter storing triangles approximating the surface
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{delta}
 *     \cgalParamDescription{the value of \f$ \delta \f$}
 *     \cgalParamType{double}
 *     \cgalParamDefault{1.7}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{eta}
 *     \cgalParamDescription{the value of \f$ \eta \f$}
 *     \cgalParamType{double}
 *     \cgalParamDefault{200.}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `DelaunayTriangulationTraits_3`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the value type of `PointRange`, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class NamedParameters = parameters::Default_named_parameters,
          class PointRange,
          class TripleIndexRange>
void ball_merge_surface_reconstruction_local(const PointRange& points,
                                             TripleIndexRange& out_triangles,
                                             const NamedParameters& np = parameters::default_values())
{
  using Point_3 = std::remove_const_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Traits_ = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t, NamedParameters,
                                                               typename Kernel_traits<Point_3>::type>::type;
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters, Sequential_tag>::type;

  Ball_merge_surface_reconstruction<Traits_, Concurrency_tag> bmsr;
  bmsr.build_triangulation(points);
  bmsr.local_reconstruction(out_triangles, np);
}
#endif

/*! \ingroup PkgBallMergeRef
 *
 * creates two watertight meshes approximating the surface, and puts the resulting triangle faces in `out_triangles1` and `out_triangles2`.
 * As this function creates two shells (outer and inner) the input point set shall only be sampled on a single connected component.
 * Output triangle faces are triple of indices refering to the position of the input points in `points`.
 *
 * \tparam PointRange a model of the concepts `RandomAccessContainer`
 * \tparam TripleIndexRange a model of `BackInsertionSequence` with `value_type`
 *                          being a model of `RandomAccessContainer` and `BackInsertionSequence` with `value_type`
 *                          being constructible from `std::size_t`.
 *
 * \param points is the input points representing a single connected component of a watertight surface
 * \param out_triangles1 is the output parameter storing the first resulting mesh
 * \param out_triangles2 is the output parameter storing the second resulting mesh
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \return the delta value used for the reconstruction
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{delta}
 *     \cgalParamDescription{the value of \f$ \delta \f$}
 *     \cgalParamType{double}
 *     \cgalParamDefault{Value iteratively estimated, as described in \ref secBMSRPar}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `DelaunayTriangulationTraits_3`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the value type of `PointRange`, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class NamedParameters = parameters::Default_named_parameters,
          class PointRange,
          class TripleIndexRange>
double ball_merge_surface_reconstruction(const PointRange& points,
                                         TripleIndexRange& out_triangles1,
                                         TripleIndexRange& out_triangles2,
                                         const NamedParameters& np = parameters::default_values())
{
  using Point_3 = std::remove_const_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Traits_ = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t, NamedParameters,
                                                               typename Kernel_traits<Point_3>::type>::type;
  using Concurrency_tag = typename internal_np::Lookup_named_param_def<internal_np::concurrency_tag_t, NamedParameters, Sequential_tag>::type;

  Ball_merge_surface_reconstruction<Traits_, Concurrency_tag> bmsr;
  bmsr.build_triangulation(points);
  double delta = bmsr.reconstruction(out_triangles1, out_triangles2, np);

  return delta;
}


} // namespace CGAL



#endif // CGAL_BALL_MERGE_SURFACE_RECONSTRUCTION_H
