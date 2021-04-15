// Copyright (c) 1999-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 Clement Jamin

#ifndef CGAL_REGULAR_TRIANGULATION_3_H
#define CGAL_REGULAR_TRIANGULATION_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <CGAL/point_generators_3.h>
# include <tbb/parallel_for.h>
# include <thread>
# include <tbb/enumerable_thread_specific.h>
# include <tbb/concurrent_vector.h>
#endif

#include <CGAL/Triangulation_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/internal/Has_nested_type_Bare_point.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/result_of.h>

#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/internal/info_check.h>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>
#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#ifdef CGAL_TRIANGULATION_3_PROFILING
# include <CGAL/Mesh_3/Profiling_tools.h>
#endif

#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
#include <CGAL/point_generators_3.h>
#endif

#include <boost/bind.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/utility/result_of.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <thread>
#include <utility>
#include <vector>

namespace CGAL {

/************************************************
   *
   * Regular_triangulation_3 class
   *
   ************************************************/

template < class Gt, class Tds_ = Default, class Lock_data_structure_ = Default >
class Regular_triangulation_3
  : public Triangulation_3<
             Gt,
             typename Default::Get<Tds_, Triangulation_data_structure_3 <
                                           Regular_triangulation_vertex_base_3<Gt>,
                                           Regular_triangulation_cell_base_3<Gt> > >::type,
             Lock_data_structure_>
{
private:
  typedef typename Default::Get<Tds_, Triangulation_data_structure_3 <
                                        Regular_triangulation_vertex_base_3<Gt>,
                                        Regular_triangulation_cell_base_3<Gt> >
                               >::type                            Tds;

  typedef Regular_triangulation_3<Gt, Tds_, Lock_data_structure_> Self;

public:
  typedef Triangulation_3<Gt, Tds, Lock_data_structure_>          Tr_Base;

  typedef Gt                                    Geom_traits;
  typedef Tds                                   Triangulation_data_structure;

  typedef Geom_traits                           Traits;
  typedef typename Tr_Base::Concurrency_tag     Concurrency_tag;
  typedef typename Tr_Base::Lock_data_structure Lock_data_structure;

  typedef typename Tr_Base::Vertex_handle       Vertex_handle;
  typedef typename Tr_Base::Cell_handle         Cell_handle;
  typedef typename Tr_Base::Vertex              Vertex;
  typedef typename Tr_Base::Cell                Cell;
  typedef typename Tr_Base::Facet               Facet;
  typedef typename Tr_Base::Edge                Edge;

  typedef typename Tr_Base::size_type           size_type;
  typedef typename Tr_Base::Locate_type         Locate_type;
  typedef typename Tr_Base::Cell_iterator       Cell_iterator;
  typedef typename Tr_Base::Facet_iterator      Facet_iterator;
  typedef typename Tr_Base::Edge_iterator       Edge_iterator;
  typedef typename Tr_Base::Facet_circulator    Facet_circulator;

  typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Tr_Base::All_cells_iterator       All_cells_iterator;

  // Traits are not supposed to define Bare_point, but leaving below
  // for backward compatibility
  typedef typename boost::mpl::eval_if_c<
    internal::Has_nested_type_Bare_point<Gt>::value,
    typename internal::Bare_point_type<Gt>,
    boost::mpl::identity<typename Gt::Point_3>
  >::type                                          Bare_point;
  typedef typename Gt::Weighted_point_3            Weighted_point;

  typedef typename Gt::Segment_3                   Segment;
  typedef typename Gt::Triangle_3                  Triangle;
  typedef typename Gt::Tetrahedron_3               Tetrahedron;

  // types for dual:
  typedef typename Gt::Line_3        Line;
  typedef typename Gt::Ray_3         Ray;
  typedef typename Gt::Plane_3       Plane;
  typedef typename Gt::Object_3      Object;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_true                   Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_false                  Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Tr_Base::geom_traits;
#endif
  using Tr_Base::adjacent_vertices;
  using Tr_Base::adjacent_vertices_threadsafe;
  using Tr_Base::cw;
  using Tr_Base::ccw;
  using Tr_Base::construct_point;
  using Tr_Base::coplanar_orientation;
  using Tr_Base::dimension;
  using Tr_Base::find_conflicts;
  using Tr_Base::finite_facets_begin;
  using Tr_Base::finite_facets_end;
  using Tr_Base::finite_vertices_begin;
  using Tr_Base::finite_vertices_end;
  using Tr_Base::finite_cells_begin;
  using Tr_Base::finite_cells_end;
  using Tr_Base::finite_edges_begin;
  using Tr_Base::finite_edges_end;
  using Tr_Base::incident_facets;
  using Tr_Base::insert_in_conflict;
  using Tr_Base::infinite_vertex;
  using Tr_Base::is_infinite;
  using Tr_Base::is_valid;
  using Tr_Base::is_valid_finite;
  using Tr_Base::locate;
  using Tr_Base::mirror_vertex;
  using Tr_Base::mirror_index;
  using Tr_Base::next_around_edge;
  using Tr_Base::number_of_vertices;
  using Tr_Base::orientation;
  using Tr_Base::point;
  using Tr_Base::side_of_segment;
  using Tr_Base::side_of_edge;
  using Tr_Base::tds;
  using Tr_Base::vertex_triple_index;

  Regular_triangulation_3(const Gt& gt = Gt(), Lock_data_structure *lock_ds = nullptr)
    : Tr_Base(gt, lock_ds), hidden_point_visitor(this)
  { }

  Regular_triangulation_3(Lock_data_structure *lock_ds, const Gt& gt = Gt())
    : Tr_Base(lock_ds, gt), hidden_point_visitor(this)
  { }

  Regular_triangulation_3(const Regular_triangulation_3& rt)
    : Tr_Base(rt), hidden_point_visitor(this)
  {
    CGAL_triangulation_postcondition(is_valid());
  }

  Regular_triangulation_3(Regular_triangulation_3&& rt)
    noexcept(noexcept(Tr_Base(std::move(rt))))
    : Tr_Base(std::move(rt)), hidden_point_visitor(this)
  {
    CGAL_triangulation_postcondition(is_valid());
  }

  ~Regular_triangulation_3() = default;

  void swap(Regular_triangulation_3& tr)
  {
    // The 'vertices' and 'hidden_points' members of
    // 'hidden_point_visitor' should be empty as they are only filled
    // (and cleared) during the insertion of a point.  Hidden points
    // are not stored there, but rather in cells. Thus, the only thing
    // that must be set is the triangulation pointer, and it is
    // already correctly set. There is nothing to do about
    // 'hidden_point_visitor'.
    Tr_Base::swap(tr);
  }

  Regular_triangulation_3& operator=(const Regular_triangulation_3& tr)
  {
    Regular_triangulation_3 copy(tr);
    copy.swap(*this);
    return *this;
  }

  Regular_triangulation_3& operator=(Regular_triangulation_3&& tr)
    noexcept(noexcept(Regular_triangulation_3(std::move(tr))))
  {
    Tr_Base::operator=(std::move(tr));
    return *this;
  }

  //insertion
  template < typename InputIterator >
  Regular_triangulation_3(InputIterator first, InputIterator last,
                          const Gt& gt = Gt(), Lock_data_structure *lock_ds = nullptr)
    : Tr_Base(gt, lock_ds), hidden_point_visitor(this)
  {
    insert(first, last);
  }

  template < typename InputIterator >
  Regular_triangulation_3(InputIterator first, InputIterator last,
                          Lock_data_structure *lock_ds, const Gt& gt = Gt())
    : Tr_Base(gt, lock_ds), hidden_point_visitor(this)
  {
    insert(first, last);
  }

private:
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
  std::vector<Vertex_handle>
  add_temporary_points_on_far_sphere(const size_t num_points)
  {
    std::vector<Vertex_handle> far_sphere_vertices;

    const size_t MIN_NUM_POINTS_FOR_FAR_SPHERE_POINTS = 1000000;
    if(num_points >= MIN_NUM_POINTS_FOR_FAR_SPHERE_POINTS)
    {
      // Add temporary vertices on a "far sphere" to reduce contention on
      // the infinite vertex

      // Get bbox
      const Bbox_3& bbox = *this->get_bbox();
      // Compute radius for far sphere
      const double& xdelta = bbox.xmax() - bbox.xmin();
      const double& ydelta = bbox.ymax() - bbox.ymin();
      const double& zdelta = bbox.zmax() - bbox.zmin();
      const double radius = 1.3 * 0.5 * std::sqrt(xdelta*xdelta +
                                                  ydelta*ydelta +
                                                  zdelta*zdelta);

      // WARNING - TODO @fixme this code has to be fixed because Vector_3 is not
      // required by the traits concept
      const typename Gt::Vector_3 center(bbox.xmin() + 0.5*xdelta,
                                         bbox.ymin() + 0.5*ydelta,
                                         bbox.zmin() + 0.5*zdelta);
      Random_points_on_sphere_3<Bare_point> random_point(radius);
      const int NUM_PSEUDO_INFINITE_VERTICES = static_cast<int>(
                                                 std::thread::hardware_concurrency() * 3.5);
      typename Gt::Construct_weighted_point_3 cwp =
          geom_traits().construct_weighted_point_3_object();

      std::vector<Weighted_point> points_on_far_sphere;
      for(int i = 0 ; i < NUM_PSEUDO_INFINITE_VERTICES ; ++i, ++random_point)
        points_on_far_sphere.push_back(cwp(*random_point + center));

      // Spatial sorting can only be applied to bare points, so we need an adaptor
      typedef typename Geom_traits::Construct_point_3 Construct_point_3;
      typedef typename boost::result_of<const Construct_point_3(const Weighted_point&)>::type Ret;
      typedef boost::function_property_map<Construct_point_3, Weighted_point, Ret> fpmap;
      typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, fpmap> Search_traits_3;

      spatial_sort(points_on_far_sphere.begin(), points_on_far_sphere.end(),
                   Search_traits_3(
                     boost::make_function_property_map<Weighted_point, Ret, Construct_point_3>(
                       geom_traits().construct_point_3_object()), geom_traits()));

      typename std::vector<Weighted_point>::const_iterator it_p =
          points_on_far_sphere.begin();
      typename std::vector<Weighted_point>::const_iterator it_p_end =
          points_on_far_sphere.end();

      for(; it_p != it_p_end ; ++it_p)
      {
        Locate_type lt;
        Cell_handle c, hint;
        int li, lj;

        c = locate(*it_p, lt, li, lj, hint);
        Vertex_handle v = insert(*it_p, lt, c, li, lj);
        hint = (v == Vertex_handle() ? c : v->cell());

        far_sphere_vertices.push_back(v);
      }
    }

    return far_sphere_vertices;
  }

  void remove_temporary_points_on_far_sphere(
      const std::vector<Vertex_handle>& far_sphere_vertices)
  {
    if(!far_sphere_vertices.empty())
    {
      // Remove the temporary vertices on far sphere
      remove(far_sphere_vertices.begin(), far_sphere_vertices.end());
    }
  }
#endif // CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE

public:
#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        typename boost::enable_if<
                          boost::is_convertible<
                          typename std::iterator_traits<InputIterator>::value_type,
                          Weighted_point> >::type* = nullptr)
#else
  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  {
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
    static Profile_branch_counter_3 bcounter(
          "early withdrawals / late withdrawals / successes [Regular_tri_3::insert]");
#endif

#ifdef CGAL_TRIANGULATION_3_PROFILING
    WallClockTimer t;
#endif

    size_type n = number_of_vertices();
    std::vector<Weighted_point> points(first, last);

    // Spatial sorting can only be applied to bare points, so we need an adaptor
    // @todo Unary_function_to_property_map makes a copy (get() returns a value_type) but
    // we could hope to get a const& to the bare point. Unfortunately, the lazy
    // kernel creates temporaries and prevent it.
    typedef typename Geom_traits::Construct_point_3 Construct_point_3;
    typedef typename boost::result_of<const Construct_point_3(const Weighted_point&)>::type Ret;
    typedef boost::function_property_map<Construct_point_3, Weighted_point, Ret> fpmap;
    typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, fpmap> Search_traits_3;

    spatial_sort(points.begin(), points.end(),
                 Search_traits_3(
                   boost::make_function_property_map<Weighted_point, Ret, Construct_point_3>(
                     geom_traits().construct_point_3_object()), geom_traits()));

    // Parallel
#ifdef CGAL_LINKED_WITH_TBB
    if(this->is_parallel())
    {
      size_t num_points = points.size();
      Cell_handle hint;

#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
      std::vector<Vertex_handle> far_sphere_vertices =
          add_temporary_points_on_far_sphere(num_points);
#endif

      size_t i = 0;
      // Insert "num_points_seq" points sequentially
      // (or more if dim < 3 after that)
      size_t num_points_seq = (std::min)(num_points, (size_t)100);
      while (i < num_points_seq || (dimension() < 3 && i < num_points))
      {
        Locate_type lt;
        Cell_handle c;
        int li, lj;

        c = locate(points[i], lt, li, lj, hint);
        Vertex_handle v = insert (points[i], lt, c, li, lj);

        hint = (v == Vertex_handle() ? c : v->cell());
        ++i;
      }

      tbb::enumerable_thread_specific<Vertex_handle> tls_hint(hint->vertex(0));
      tbb::parallel_for(tbb::blocked_range<size_t>(i, num_points),
                        Insert_point<Self>(*this, points, tls_hint));

#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
      remove_temporary_points_on_far_sphere(far_sphere_vertices);
#endif
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      Cell_handle hint;
      for(typename std::vector<Weighted_point>::const_iterator p = points.begin(),
           end = points.end(); p != end; ++p)
      {
        Locate_type lt;
        Cell_handle c;
        int li, lj;
        c = locate(*p, lt, li, lj, hint);

        Vertex_handle v = insert(*p, lt, c, li, lj);

        hint = v == Vertex_handle() ? c : v->cell();
      }
    }
#ifdef CGAL_TRIANGULATION_3_PROFILING
    std::cerr << "Points inserted in " << t.elapsed() << " seconds." << std::endl;
#endif
    return number_of_vertices() - n;
  }

#ifndef CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
private:

  //top stands for tuple-or-pair
  template <class Info>
  const Weighted_point& top_get_first(const std::pair<Weighted_point,Info>& pair) const { return pair.first; }

  template <class Info>
  const Info& top_get_second(const std::pair<Weighted_point,Info>& pair) const { return pair.second; }

  template <class Info>
  const Weighted_point& top_get_first(const boost::tuple<Weighted_point,Info>& tuple) const { return boost::get<0>(tuple); }

  template <class Info>
  const Info& top_get_second(const boost::tuple<Weighted_point,Info>& tuple) const { return boost::get<1>(tuple); }

  // Functor to go from an index of a container of Weighted_point to
  // the corresponding Bare_point
  template<class Construct_bare_point, class Container>
  struct Index_to_Bare_point
  {
      typename boost::result_of<const Construct_bare_point(const Weighted_point&)>::type
      operator()(const std::size_t& i) const
    {
      return cp(c[i]);
    }

    Index_to_Bare_point(const Container& c, const Construct_bare_point& cp)
      : c(c), cp(cp) { }

    const Container& c;
    const Construct_bare_point cp;
  };

  template <class Tuple_or_pair,class InputIterator>
  std::ptrdiff_t insert_with_info(InputIterator first,InputIterator last)
  {
    size_type n = number_of_vertices();
    std::vector<std::size_t> indices;
    std::vector<Weighted_point> points;
    std::vector<typename Triangulation_data_structure::Vertex::Info> infos;
    std::size_t index=0;
    for(InputIterator it=first;it!=last;++it)
    {
      Tuple_or_pair pair = *it;
      points.push_back(top_get_first(pair));
      infos.push_back(top_get_second(pair));
      indices.push_back(index++);
    }

    // We need to sort the points and their info at the same time through
    // the `indices` vector AND spatial sort can only handle Gt::Point_3.
    typedef typename Geom_traits::Construct_point_3 Construct_point_3;
    typedef Index_to_Bare_point<Construct_point_3,
        std::vector<Weighted_point> > Access_bare_point;
    typedef typename boost::result_of<const Construct_point_3(const Weighted_point&)>::type Ret;
    typedef boost::function_property_map<Access_bare_point, std::size_t, Ret> fpmap;
    typedef CGAL::Spatial_sort_traits_adapter_3<Gt, fpmap> Search_traits_3;

    Access_bare_point accessor(points, geom_traits().construct_point_3_object());
    spatial_sort(indices.begin(), indices.end(),
                 Search_traits_3(
                   boost::make_function_property_map<
                     std::size_t, Ret, Access_bare_point>(accessor),
                   geom_traits()));

#ifdef CGAL_LINKED_WITH_TBB
    if(this->is_parallel())
    {
      size_t num_points = points.size();
      Cell_handle hint;

#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
      std::vector<Vertex_handle> far_sphere_vertices =
          add_temporary_points_on_far_sphere(num_points);
#endif

      size_t i = 0;
      // Insert "num_points_seq" points sequentially
      // (or more if dim < 3 after that)
      size_t num_points_seq = (std::min)(num_points, (size_t)100);
      while (i < num_points_seq || (dimension() < 3 && i < num_points))
      {
        Locate_type lt;
        Cell_handle c;
        int li, lj;
        c = locate(points[indices[i]], lt, li, lj, hint);

        Vertex_handle v = insert(points[indices[i]], lt, c, li, lj);
        if(v != Vertex_handle())
        {
          v->info() = infos[indices[i]];
          hint = v->cell();
        }
        else
          hint = c;

        ++i;
      }

      tbb::enumerable_thread_specific<Vertex_handle> tls_hint(hint->vertex(0));
      tbb::parallel_for(tbb::blocked_range<size_t>(i, num_points),
                        Insert_point_with_info<Self>(*this, points, infos, indices, tls_hint));

#ifdef CGAL_CONCURRENT_TRIANGULATION_3_ADD_TEMPORARY_POINTS_ON_FAR_SPHERE
      remove_temporary_points_on_far_sphere(far_sphere_vertices);
#endif
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      Cell_handle hint;
      for(typename std::vector<std::size_t>::const_iterator
           it = indices.begin(), end = indices.end();
           it != end; ++it)
      {
        Locate_type lt;
        Cell_handle c;
        int li, lj;
        c = locate(points[*it], lt, li, lj, hint);

        Vertex_handle v = insert(points[*it], lt, c, li, lj);
        if(v!=Vertex_handle())
        {
          v->info()=infos[*it];
          hint=v->cell();
        }
        else
        {
          hint = c;
        }
      }
    }

    return number_of_vertices() - n;
  }

public:

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first,
                        InputIterator last,
                        typename boost::enable_if<
                        boost::is_convertible<
                        typename std::iterator_traits<InputIterator>::value_type,
                        std::pair<Weighted_point,typename internal::Info_check<typename Triangulation_data_structure::Vertex>::type>
                        >
                        >::type* = nullptr)
  {
    return insert_with_info<
             std::pair<Weighted_point,
                       typename internal::Info_check<
                         typename Triangulation_data_structure::Vertex>::type>
           >(first,last);
  }

  template <class  InputIterator_1,class InputIterator_2>
  std::ptrdiff_t
  insert(boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > first,
         boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > last,
         typename boost::enable_if<
           boost::mpl::and_<
           typename boost::is_convertible< typename std::iterator_traits<InputIterator_1>::value_type, Weighted_point >,
           typename boost::is_convertible< typename std::iterator_traits<InputIterator_2>::value_type, typename internal::Info_check<typename Triangulation_data_structure::Vertex>::type >
         > >::type* =nullptr)
  {
    return insert_with_info<
             boost::tuple<Weighted_point,
                          typename internal::Info_check<
                            typename Triangulation_data_structure::Vertex>::type>
            >(first,last);
  }
#endif //CGAL_TRIANGULATION_3_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO


  Vertex_handle insert(const Weighted_point& p, Vertex_handle hint,
                       bool *could_lock_zone = nullptr)
  {
    return insert(p,
                  hint == Vertex_handle() ? this->infinite_cell() : hint->cell(),
                  could_lock_zone);
  }

  Vertex_handle insert(const Weighted_point& p,
                       Cell_handle start = Cell_handle(), bool *could_lock_zone = nullptr);

  Vertex_handle insert(const Weighted_point& p, Locate_type lt,
                       Cell_handle c, int li, int, bool *could_lock_zone = nullptr);

  template <class CellIt>
  Vertex_handle insert_in_hole(const Weighted_point& p,
                               CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i);

  template <class CellIt>
  Vertex_handle insert_in_hole(const Weighted_point& p,
                               CellIt cell_begin, CellIt cell_end,
                               Cell_handle begin, int i, Vertex_handle newv);

  template <class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(const Weighted_point& p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 OutputIteratorInternalFacets ifit,
                 bool *could_lock_zone = nullptr,
                 const Facet *this_facet_must_be_in_the_cz = nullptr,
                 bool *the_facet_is_in_its_cz = nullptr) const
  {
    CGAL_triangulation_precondition(dimension() >= 2);

    std::vector<Cell_handle> cells;
    cells.reserve(32);
    std::vector<Facet> facets;
    facets.reserve(64);

    if(dimension() == 2)
    {
      Conflict_tester_2 tester(p, this);
      if(! tester(c))
        return make_triple(bfit, cit, ifit);

      ifit = Tr_Base::find_conflicts(c, tester,
                                     make_triple(std::back_inserter(facets),
                                                 std::back_inserter(cells),
                                                 ifit),
                                     could_lock_zone,
                                     this_facet_must_be_in_the_cz,
                                     the_facet_is_in_its_cz).third;
    }
    else
    {
      Conflict_tester_3 tester(p, this);
      if(! tester(c))
        return make_triple(bfit, cit, ifit);

      ifit = Tr_Base::find_conflicts(c, tester,
                                     make_triple(std::back_inserter(facets),
                                                 std::back_inserter(cells),
                                                 ifit),
                                     could_lock_zone,
                                     this_facet_must_be_in_the_cz,
                                     the_facet_is_in_its_cz).third;
    }

    // Reset the conflict flag on the boundary.
    for(typename std::vector<Facet>::iterator fit = facets.begin();
                                              fit != facets.end(); ++fit)
    {
      fit->first->neighbor(fit->second)->tds_data().clear();
      *bfit++ = *fit;
    }

    // Reset the conflict flag in the conflict cells.
    for(typename std::vector<Cell_handle>::iterator ccit = cells.begin();
                                                    ccit != cells.end(); ++ccit)
    {
      (*ccit)->tds_data().clear();
      *cit++ = *ccit;
    }
    return make_triple(bfit, cit, ifit);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Weighted_point& p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
                 bool *could_lock_zone = nullptr) const
  {
    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           Emptyset_iterator> t = find_conflicts(p, c, bfit, cit,
                                                 Emptyset_iterator(),
                                                 could_lock_zone);
    return std::make_pair(t.first, t.second);
  }

  // Returns the vertices on the interior of the conflict hole.
  template <class OutputIterator>
  OutputIterator vertices_inside_conflict_zone(const Weighted_point&p, Cell_handle c,
                                               OutputIterator res) const
  {
    CGAL_triangulation_precondition(dimension() >= 2);

    // Get the facets on the boundary of the hole, and the cells of the hole
    std::vector<Cell_handle> cells;
    std::vector<Facet> facets;
    find_conflicts(p, c, std::back_inserter(facets),
                   std::back_inserter(cells), Emptyset_iterator());

    // Put all vertices on the hole in 'vertices'
    const int d = dimension();
    std::set<Vertex_handle> vertices;
    for(typename std::vector<Cell_handle>::const_iterator it = cells.begin(),
        end = cells.end(); it != end; ++it)
    {
      for(int i = 0; i <= d; ++i)
        vertices.insert((*it)->vertex(i));
    }
    // Then extract the vertices of the boundary and remove them from
    // 'vertices'
    if(dimension() == 3)
    {
      for(typename std::vector<Facet>::const_iterator i = facets.begin();
           i != facets.end(); ++i)
      {
        vertices.erase(i->first->vertex((i->second+1)&3));
        vertices.erase(i->first->vertex((i->second+2)&3));
        vertices.erase(i->first->vertex((i->second+3)&3));
      }
    }
    else
    {
      for(typename std::vector<Facet>::const_iterator i = facets.begin();
           i != facets.end(); ++i)
      {
        vertices.erase(i->first->vertex(cw(i->second)));
        vertices.erase(i->first->vertex(ccw(i->second)));
      }
    }

    return std::copy(vertices.begin(), vertices.end(), res);
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  // Returns the vertices on the boundary of the conflict hole.
  template <class OutputIterator>
  OutputIterator vertices_in_conflict(const Weighted_point&p, Cell_handle c,
                                      OutputIterator res) const
  {
    return vertices_on_conflict_zone_boundary(p, c, res);
  }
#endif // CGAL_NO_DEPRECATED_CODE

  // Returns the vertices on the boundary of the conflict hole.
  template <class OutputIterator>
  OutputIterator vertices_on_conflict_zone_boundary(const Weighted_point&p,
                                                    Cell_handle c,
                                                    OutputIterator res) const
  {
    CGAL_triangulation_precondition(dimension() >= 2);

    // Get the facets on the boundary of the hole.
    std::vector<Facet> facets;
    find_conflicts(p, c, std::back_inserter(facets),
                   Emptyset_iterator(), Emptyset_iterator());

    // Then extract uniquely the vertices.
    std::set<Vertex_handle> vertices;
    if(dimension() == 3)
    {
      for(typename std::vector<Facet>::const_iterator i = facets.begin();
                                                      i != facets.end(); ++i)
      {
        vertices.insert(i->first->vertex((i->second+1)&3));
        vertices.insert(i->first->vertex((i->second+2)&3));
        vertices.insert(i->first->vertex((i->second+3)&3));
      }
    } else {
      for(typename std::vector<Facet>::const_iterator i = facets.begin();
                                                      i != facets.end(); ++i)
      {
        vertices.insert(i->first->vertex(cw(i->second)));
        vertices.insert(i->first->vertex(ccw(i->second)));
      }
    }

    return std::copy(vertices.begin(), vertices.end(), res);
  }

    // In parallel operations, we need to be able to check the health of the 'hint' vertex handle,
    // which might be invalided by other threads. One way to do that is the 'is_vertex()' function
    // of the TDS, but it runs in O(sqrt(n)) complexity. When we are using our TDS, we can use
    // a lower level function from the compact container, which runs in constant time.
    BOOST_MPL_HAS_XXX_TRAIT_DEF(Is_CGAL_TDS_3)

    template <typename TDS_,
              bool has_TDS_tag = has_Is_CGAL_TDS_3<TDS_>::value>
    struct Is_CGAL_TDS_3 : public CGAL::Tag_false
    { };

    template <typename TDS_>
    struct Is_CGAL_TDS_3<TDS_, true> : public CGAL::Boolean_tag<TDS_::Is_CGAL_TDS_3::value>
    { };

    template <typename TDS_,
              bool is_CGAL_TDS_3 = Is_CGAL_TDS_3<TDS_>::value>
    struct Vertex_validity_checker
    {
      bool operator()(const typename TDS_::Vertex_handle vh_, const TDS_& tds_) {
        return tds_.is_vertex(vh_);
      }
    };

    template <typename TDS_>
    struct Vertex_validity_checker<TDS_, true /* is_CGAL_TDS_3 */>
    {
      bool operator()(const typename TDS_::Vertex_handle vh_, const TDS_& tds_) {
        return tds_.vertices().is_used(vh_);
      }
    };

  void remove(Vertex_handle v);
  // Concurrency-safe
  // See Triangulation_3::remove for more information
  bool remove(Vertex_handle v, bool *could_lock_zone);

  template < typename InputIterator >
  size_type remove(InputIterator first, InputIterator beyond)
  {
    CGAL_triangulation_precondition(!this->does_repeat_in_range(first, beyond));
    size_type n = number_of_vertices();

#ifdef CGAL_TRIANGULATION_3_PROFILING
    WallClockTimer t;
#endif

    // Parallel
#ifdef CGAL_LINKED_WITH_TBB
    if(this->is_parallel())
    {
      // TODO: avoid that by asking for ramdom-access iterators?
      std::vector<Vertex_handle> vertices(first, beyond);
      tbb::concurrent_vector<Vertex_handle> vertices_to_remove_sequentially;

      tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size()),
                        Remove_point<Self>(*this, vertices, vertices_to_remove_sequentially));

      // Do the rest sequentially
      for(typename tbb::concurrent_vector<Vertex_handle>::const_iterator
            it = vertices_to_remove_sequentially.begin(),
            it_end = vertices_to_remove_sequentially.end()
            ; it != it_end
            ; ++it)
      {
        remove(*it);
      }
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      while(first != beyond)
      {
        remove(*first);
        ++first;
      }
    }

#ifdef CGAL_TRIANGULATION_3_PROFILING
    std::cerr << "Points removed in " << t.elapsed() << " seconds." << std::endl;
#endif
    return n - number_of_vertices();
  }

  template <class OutputItCells>
  void remove_and_give_new_cells(Vertex_handle v, OutputItCells cit)
  {
    Self tmp;
    Vertex_remover<Self> remover(tmp);
    Tr_Base::remove_and_give_new_cells(v, remover, cit);

    CGAL_triangulation_expensive_postcondition(is_valid());
  }

  // Displacement works only for regular triangulation
  // without hidden points at any time
  Vertex_handle move_if_no_collision(Vertex_handle v, const Weighted_point& p);
  Vertex_handle move(Vertex_handle v, const Weighted_point& p);

  // REMOVE CLUSTER - works only when Regular has no hidden point at all
  // "regular as Delaunay"
  template < typename InputIterator >
  size_type remove_cluster(InputIterator first, InputIterator beyond)
  {
    Self tmp;
    Vertex_remover<Self> remover(tmp);
    return Tr_Base::remove(first, beyond, remover);
  }

protected:

  Oriented_side side_of_oriented_power_sphere(const Weighted_point& p0,
                                              const Weighted_point& p1,
                                              const Weighted_point& p2,
                                              const Weighted_point& p3,
                                              const Weighted_point& p,
                                              bool perturb = false) const;

  Oriented_side side_of_oriented_power_circle(const Weighted_point& p0,
                                              const Weighted_point& p1,
                                              const Weighted_point& p2,
                                              const Weighted_point& p,
                                              bool perturb = false) const;

  Bounded_side side_of_bounded_power_circle(const Weighted_point& p0,
                                            const Weighted_point& p1,
                                            const Weighted_point& p2,
                                            const Weighted_point& p,
                                            bool perturb = false) const;

  Bounded_side side_of_bounded_power_segment(const Weighted_point& p0,
                                             const Weighted_point& p1,
                                             const Weighted_point& p,
                                             bool perturb = false) const;

public:
  // Queries
  Bounded_side side_of_power_sphere(Cell_handle c, const Weighted_point& p,
                                    bool perturb = false) const;

  Bounded_side side_of_power_circle(const Facet& f, const Weighted_point& p,
                                    bool /* perturb */ = false) const
  {
    return side_of_power_circle(f.first, f.second, p);
  }

  Bounded_side side_of_power_circle(Cell_handle c, int i, const Weighted_point& p,
                                    bool perturb = false) const;

  Bounded_side side_of_power_segment(Cell_handle c, const Weighted_point& p,
                                     bool perturb = false) const;

  // Undocumented, needed for Mesh_3 (because of Periodic_3_mesh_3)
  bool greater_or_equal_power_distance(const Bare_point& p,
                                       const Weighted_point& q,
                                       const Weighted_point& r) const;

  Vertex_handle nearest_power_vertex_in_cell(const Bare_point& p,
                                             Cell_handle c)  const;

  Vertex_handle nearest_power_vertex(const Bare_point& p, Cell_handle c = Cell_handle()) const;

  bool is_Gabriel(Cell_handle c, int i) const;
  bool is_Gabriel(Cell_handle c, int i, int j) const;
  bool is_Gabriel(const Facet& f)const ;
  bool is_Gabriel(const Edge& e) const;
  bool is_Gabriel(Vertex_handle v) const;

  // Dual functions
  Bare_point dual(Cell_handle c) const;
  Object dual(Cell_handle c, int i) const;
  Object dual(const Facet& facet) const;

  void dual_segment(Cell_handle c, int i, Bare_point& p, Bare_point&q) const;
  void dual_segment(const Facet& facet, Bare_point& p, Bare_point&q) const;
  void dual_segment_exact(const Facet& facet, Bare_point& p, Bare_point&q) const;
  void dual_ray(Cell_handle c, int i, Ray& ray) const;
  void dual_ray(const Facet& facet, Ray& ray) const;
  void dual_ray_exact(const Facet& facet, Ray& ray) const;

  template < class Stream>
  Stream& draw_dual(Stream& os) const;

  bool is_valid(bool verbose = false, int level = 0) const;

protected:
  bool less_power_distance(const Bare_point& p,
                           const Weighted_point& q,
                           const Weighted_point& r)  const
  {
    return geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
  }

  Bare_point construct_weighted_circumcenter(const Weighted_point& p,
                                             const Weighted_point& q,
                                             const Weighted_point& r,
                                             const Weighted_point& s) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
  }

  Bare_point construct_weighted_circumcenter(const Weighted_point& p,
                                             const Weighted_point& q,
                                             const Weighted_point& r) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
  }

  Segment construct_segment(const Bare_point& p, const Bare_point& q) const
  {
    return geom_traits().construct_segment_3_object()(p, q);
  }

  Line construct_perpendicular_line(const Plane& pl, const Bare_point& p) const
  {
    return geom_traits().construct_perpendicular_line_3_object()(pl, p);
  }

  Plane construct_plane(const Bare_point& p, const Bare_point& q, const Bare_point& r) const
  {
    return geom_traits().construct_plane_3_object()(p, q, r);
  }

  Ray construct_ray(const Bare_point& p, const Line& l) const
  {
    return geom_traits().construct_ray_3_object()(p, l);
  }

  Object construct_object(const Bare_point& p) const
  {
    return geom_traits().construct_object_3_object()(p);
  }

  Object construct_object(const Segment& s) const
  {
    return geom_traits().construct_object_3_object()(s);
  }

  Object construct_object(const Ray& r) const
  {
    return geom_traits().construct_object_3_object()(r);
  }

  Vertex_handle nearest_power_vertex(const Bare_point& p,
                                     Vertex_handle v,
                                     Vertex_handle w) const
  {
    // In case of equality, v is returned.
    CGAL_triangulation_precondition(v != w);
    if(is_infinite(v))
      return w;

    if(is_infinite(w))
      return v;

    return less_power_distance(p, w->point(), v->point()) ? w : v;
  }

  Oriented_side power_test(const Weighted_point& p, const Weighted_point& q) const
  {
    CGAL_triangulation_precondition(this->equal(p, q));
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(p, q);
  }

  Oriented_side power_test(const Weighted_point& p, const Weighted_point& q,
                           const Weighted_point& r) const
  {
    CGAL_triangulation_precondition(this->collinear(p, q, r));
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(p, q, r);
  }

  Oriented_side power_test(const Weighted_point& p, const Weighted_point& q,
                           const Weighted_point& r, const Weighted_point& s) const
  {
    CGAL_triangulation_precondition(this->coplanar(p, q, r, s));
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(p, q, r, s);
  }

  Oriented_side power_test(const Weighted_point& p, const Weighted_point& q,
                           const Weighted_point& r, const Weighted_point& s,
                           const Weighted_point& t) const
  {
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(p, q, r, s, t);
  }

  bool in_conflict_3(const Weighted_point& p, const Cell_handle c) const
  {
    return side_of_power_sphere(c, p, true) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_2(const Weighted_point& p, const Cell_handle c, int i) const
  {
    return side_of_power_circle(c, i, p, true) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_1(const Weighted_point& p, const Cell_handle c) const
  {
    return side_of_power_segment(c, p, true) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_0(const Weighted_point& p, const Cell_handle c) const
  {
    return power_test(c->vertex(0)->point(), p) == ON_POSITIVE_SIDE;
  }

  bool in_conflict(const Weighted_point& p, const Cell_handle c) const
  {
    switch(dimension())
    {
      case 0:
        return in_conflict_0(p, c);
      case 1:
        return in_conflict_1(p, c);
      case 2:
        return in_conflict_2(p, c, 3);
      case 3:
        return in_conflict_3(p, c);
    }
    return true;
  }

  class Conflict_tester_3
  {
    const Weighted_point& p;
    const Self *t;

  public:
    Conflict_tester_3(const Weighted_point& pt, const Self *tr)
      : p(pt), t(tr)
    {
    }

    bool operator()(const Cell_handle c) const
    {
      return t->in_conflict_3(p, c);
    }

    bool test_initial_cell(const Cell_handle c) const
    {
      return operator()(c);
    }
    Oriented_side compare_weight(const Weighted_point& wp1,
                                 const Weighted_point& wp2) const
    {
      return t->power_test(wp1, wp2);
    }
  };

  class Conflict_tester_2
  {
    const Weighted_point& p;
    const Self *t;

  public:
    Conflict_tester_2(const Weighted_point& pt, const Self *tr)
      : p(pt), t(tr)
    {
    }

    bool operator()(const Cell_handle c) const
    {
      return t->in_conflict_2(p, c, 3);
    }

    bool test_initial_cell(const Cell_handle c) const
    {
      return operator()(c);
    }

    Oriented_side compare_weight(const Weighted_point& wp1,
                                 const Weighted_point& wp2) const
    {
      return t->power_test(wp1, wp2);
    }
  };

  class Conflict_tester_1
  {
    const Weighted_point& p;
    const Self *t;

  public:
    Conflict_tester_1(const Weighted_point& pt, const Self *tr)
      : p(pt), t(tr)
    {
    }

    bool operator()(const Cell_handle c) const
    {
      return t->in_conflict_1(p, c);
    }

    bool test_initial_cell(const Cell_handle c) const
    {
      return operator()(c);
    }

    Oriented_side compare_weight(const Weighted_point& wp1,
                                 const Weighted_point& wp2) const
    {
      return t->power_test(wp1, wp2);
    }
  };

  class Conflict_tester_0
  {
    const Weighted_point& p;
    const Self *t;

  public:
    Conflict_tester_0(const Weighted_point& pt, const Self *tr)
      : p(pt), t(tr)
    {
    }

    bool operator()(const Cell_handle c) const
    {
      return t->in_conflict_0(p, c);
    }

    bool test_initial_cell(const Cell_handle c) const
    {
      return operator()(c);
    }

    int compare_weight(const Weighted_point& wp1,
                       const Weighted_point& wp2) const
    {
      return t->power_test(wp1, wp2);
    }
  };

  // Sequential version
  // "dummy" is here to allow the specialization (see below)
  // See http://groups.google.com/group/comp.lang.c++.moderated/browse_thread/thread/285ab1eec49e1cb6
  template<typename Concurrency_tag_, typename dummy = void>
  class Hidden_point_visitor
  {
    Self *t;
    mutable std::vector<Vertex_handle> vertices;
    mutable std::vector<Weighted_point> hidden_points;

  public:
    Hidden_point_visitor(Self *tr) : t(tr) {}

    template <class InputIterator>
    void process_cells_in_conflict(InputIterator start, InputIterator end) const
    {
      int dim = t->dimension();
      while(start != end)
      {
        std::copy((*start)->hidden_points_begin(),
                  (*start)->hidden_points_end(),
                  std::back_inserter(hidden_points));

        for(int i=0; i<=dim; i++)
        {
          Vertex_handle v = (*start)->vertex(i);
          if(v->cell() != Cell_handle())
          {
            vertices.push_back(v);
            v->set_cell(Cell_handle());
          }
        }
        start ++;
      }
    }

    void reinsert_vertices(Vertex_handle v)
    {
      Cell_handle hc = v->cell();
      for(typename std::vector<Vertex_handle>::iterator
           vi = vertices.begin(); vi != vertices.end(); ++vi)
      {
        if((*vi)->cell() != Cell_handle())
          continue;

        hc = t->locate((*vi)->point(), hc);
        hide_point(hc, (*vi)->point());
        t->tds().delete_vertex(*vi);
      }

      vertices.clear();
      for(typename std::vector<Weighted_point>::iterator
           hp = hidden_points.begin(); hp != hidden_points.end(); ++hp)
      {
        hc = t->locate(*hp, hc);
        hide_point (hc, *hp);
      }
      hidden_points.clear();
    }

    Vertex_handle replace_vertex(Cell_handle c, int index, const Weighted_point& p)
    {
      Vertex_handle v = c->vertex(index);
      hide_point(c, v->point());
      v->set_point(p);
      return v;
    }

    void hide_point(Cell_handle c, const Weighted_point& p)
    {
      c->hide_point(p);
    }
  };

#ifdef CGAL_LINKED_WITH_TBB
  // Parallel version specialization
  template<typename dummy>
  class Hidden_point_visitor<Parallel_tag, dummy>
  {
    typedef Hidden_point_visitor<Parallel_tag> HPV;

    Self *t;
    mutable tbb::enumerable_thread_specific<std::vector<Vertex_handle> >  vertices;
    mutable tbb::enumerable_thread_specific<std::vector<Weighted_point> > hidden_points;

  public:
    Hidden_point_visitor(Self *tr) : t(tr) {}

    template <class InputIterator>
    void process_cells_in_conflict(InputIterator start, InputIterator end) const
    {
      int dim = t->dimension();
      while(start != end)
      {
        std::copy((*start)->hidden_points_begin(),
                  (*start)->hidden_points_end(),
                  std::back_inserter(hidden_points.local()));

        for(int i=0; i<=dim; i++)
        {
          Vertex_handle v = (*start)->vertex(i);
          if(v->cell() != Cell_handle())
          {
            vertices.local().push_back(v);
            v->set_cell(Cell_handle());
          }
        }
        start ++;
      }
    }

    void reinsert_vertices(Vertex_handle v)
    {
      Cell_handle hc = v->cell();
      for(typename std::vector<Vertex_handle>::iterator
           vi = vertices.local().begin(); vi != vertices.local().end(); ++vi)
      {
        if((*vi)->cell() != Cell_handle())
          continue;

        hc = t->locate((*vi)->point(), hc);
        hide_point(hc, (*vi)->point());
        t->tds().delete_vertex(*vi);
      }

      vertices.local().clear();
      for(typename std::vector<Weighted_point>::iterator
           hp = hidden_points.local().begin(); hp != hidden_points.local().end(); ++hp)
      {
        hc = t->locate(*hp, hc);
        hide_point (hc, *hp);
      }
      hidden_points.local().clear();
    }

    Vertex_handle replace_vertex(Cell_handle c, int index, const Weighted_point& p)
    {
      Vertex_handle v = c->vertex(index);
      hide_point(c, v->point());
      v->set_point(p);
      return v;
    }

    void hide_point(Cell_handle c, const Weighted_point& p)
    {
      c->hide_point(p);
    }
  };

  // Functor for parallel insert(begin, end) function
  template <typename RT>
  class Insert_point
  {
    typedef typename RT::Weighted_point                 Weighted_point;
    typedef typename RT::Vertex_handle                  Vertex_handle;

    RT& m_rt;
    const std::vector<Weighted_point>& m_points;
    tbb::enumerable_thread_specific<Vertex_handle>& m_tls_hint;

  public:
    // Constructor
    Insert_point(RT& rt,
                 const std::vector<Weighted_point>& points,
                 tbb::enumerable_thread_specific<Vertex_handle>& tls_hint)
      : m_rt(rt), m_points(points), m_tls_hint(tls_hint)
    {}

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
      static Profile_branch_counter_3 bcounter(
            "early withdrawals / late withdrawals / successes [Delaunay_tri_3::insert]");
#endif

      Vertex_handle& hint = m_tls_hint.local();
      Vertex_validity_checker<typename RT::Triangulation_data_structure> vertex_validity_check;

      for(size_t i_point = r.begin() ; i_point != r.end() ; ++i_point)
      {
        bool success = false;
        const Weighted_point& p = m_points[i_point];
        while(!success)
        {
          // The 'hint' is unsafe to use immediately because we are in a regular triangulation,
          // and the insertion of a (weighted) point in another thread might have hidden (deleted)
          // the hint.
          if(!vertex_validity_check(hint, m_rt.tds()))
          {
            hint = m_rt.finite_vertices_begin();
            continue;
          }

          // We need to make sure that while are locking the position P1 := hint->point(), 'hint'
          // does not get its position changed to P2 != P1.
          const Weighted_point hint_point_mem = hint->point();

          if(m_rt.try_lock_point(hint_point_mem) && m_rt.try_lock_point(p))
          {
            // Make sure that the hint is still valid (so that we can safely take hint->cell()) and
            // that its position hasn't changed to ensure that we will start the locate from where
            // we have locked.
            if(!vertex_validity_check(hint, m_rt.tds()) ||
               hint->point() != hint_point_mem)
            {
              hint = m_rt.finite_vertices_begin();
              m_rt.unlock_all_elements();
              continue;
            }

            bool could_lock_zone;
            Locate_type lt;
            int li, lj;

            Cell_handle c = m_rt.locate(p, lt, li, lj, hint->cell(), &could_lock_zone);
            Vertex_handle v;
            if(could_lock_zone)
              v = m_rt.insert(p, lt, c, li, lj, &could_lock_zone);

            if(could_lock_zone)
            {
              hint = (v == Vertex_handle() ? c->vertex(0) : v);
              m_rt.unlock_all_elements();
              success = true;
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
              ++bcounter;
#endif
            }
            else
            {
              m_rt.unlock_all_elements();
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
              bcounter.increment_branch_1(); // THIS is a late withdrawal!
#endif
            }
          }
          else
          {
            m_rt.unlock_all_elements();
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
            bcounter.increment_branch_2(); // THIS is an early withdrawal!
#endif
          }
        }
      }
    }
  };

  // Functor for parallel insert_with_info(begin, end) function
  template <typename RT>
  class Insert_point_with_info
  {
    typedef typename RT::Weighted_point                         Weighted_point;
    typedef typename RT::Vertex_handle                          Vertex_handle;
    typedef typename RT::Triangulation_data_structure::Vertex::Info Info;

    RT& m_rt;
    const std::vector<Weighted_point>& m_points;
    const std::vector<Info>& m_infos;
    const std::vector<std::size_t>& m_indices;
    tbb::enumerable_thread_specific<Vertex_handle>& m_tls_hint;

  public:
    // Constructor
    Insert_point_with_info(RT& rt,
                           const std::vector<Weighted_point>& points,
                           const std::vector<Info>& infos,
                           const std::vector<std::size_t>& indices,
                           tbb::enumerable_thread_specific<Vertex_handle>& tls_hint)
      : m_rt(rt), m_points(points), m_infos(infos), m_indices(indices),
        m_tls_hint(tls_hint)
    {}

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
      static Profile_branch_counter_3 bcounter(
            "early withdrawals / late withdrawals / successes [Delaunay_tri_3::insert]");
#endif

      Vertex_handle& hint = m_tls_hint.local();
      Vertex_validity_checker<typename RT::Triangulation_data_structure> vertex_validity_check;

      for(size_t i_idx = r.begin() ; i_idx != r.end() ; ++i_idx)
      {
        bool success = false;
        std::ptrdiff_t i_point = m_indices[i_idx];
        const Weighted_point& p = m_points[i_point];
        while(!success)
        {
          // The 'hint' is unsafe to use immediately because we are in a regular triangulation,
          // and the insertion of a (weighted) point in another thread might have hidden (deleted)
          // the hint.
          if(!vertex_validity_check(hint, m_rt.tds()))
          {
            hint = m_rt.finite_vertices_begin();
            continue;
          }

          // We need to make sure that while are locking the position P1 := hint->point(), 'hint'
          // does not get its position changed to P2 != P1.
          const Weighted_point hint_point_mem = hint->point();

          if(m_rt.try_lock_point(hint_point_mem) && m_rt.try_lock_point(p))
          {
            // Make sure that the hint is still valid (so that we can safely take hint->cell()) and
            // that its position hasn't changed to ensure that we will start the locate from where
            // we have locked.
            if(!vertex_validity_check(hint, m_rt.tds()) ||
               hint->point() != hint_point_mem)
            {
              hint = m_rt.finite_vertices_begin();
              m_rt.unlock_all_elements();
              continue;
            }

            bool could_lock_zone;
            Locate_type lt;
            int li, lj;

            Cell_handle c = m_rt.locate(p, lt, li, lj, hint->cell(), &could_lock_zone);
            Vertex_handle v;
            if(could_lock_zone)
              v = m_rt.insert(p, lt, c, li, lj, &could_lock_zone);

            if(could_lock_zone)
            {
              if(v == Vertex_handle())
              {
                hint = c->vertex(0);
              }
              else
              {
                v->info() = m_infos[i_point];
                hint = v;
              }

              m_rt.unlock_all_elements();
              success = true;
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
              ++bcounter;
#endif
            }
            else
            {
              m_rt.unlock_all_elements();
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
              bcounter.increment_branch_1(); // THIS is a late withdrawal!
#endif
            }
          }
          else
          {
            m_rt.unlock_all_elements();
#ifdef CGAL_CONCURRENT_TRIANGULATION_3_PROFILING
            bcounter.increment_branch_2(); // THIS is an early withdrawal!
#endif
          }
        }
      }
    }
  };

  // Functor for parallel remove(begin, end) function
  template <typename RT>
  class Remove_point
  {
    typedef typename RT::Weighted_point                 Weighted_point;
    typedef typename RT::Vertex_handle                  Vertex_handle;

    RT& m_rt;
    const std::vector<Vertex_handle>& m_vertices;
    tbb::concurrent_vector<Vertex_handle>& m_vertices_to_remove_sequentially;

  public:
    // Constructor
    Remove_point(RT& rt,
                 const std::vector<Vertex_handle>& vertices,
                 tbb::concurrent_vector<Vertex_handle>& vertices_to_remove_sequentially)
      : m_rt(rt), m_vertices(vertices),
        m_vertices_to_remove_sequentially(vertices_to_remove_sequentially)
    {}

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
      for(size_t i_vertex = r.begin() ; i_vertex != r.end() ; ++i_vertex)
      {
        Vertex_handle v = m_vertices[i_vertex];
        bool could_lock_zone, needs_to_be_done_sequentially;
        do
        {
          needs_to_be_done_sequentially = !m_rt.remove(v, &could_lock_zone);
          m_rt.unlock_all_elements();
        }
        while(!could_lock_zone);

        if(needs_to_be_done_sequentially)
          m_vertices_to_remove_sequentially.push_back(v);
      }
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  Hidden_point_visitor<Concurrency_tag>& get_hidden_point_visitor()
  {
    return hidden_point_visitor;
  }

  template < class RegularTriangulation_3 >
  class Vertex_remover;

  template < class RegularTriangulation_3 >
  class Vertex_inserter;

  Hidden_point_visitor<Concurrency_tag> hidden_point_visitor;
};


template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
nearest_power_vertex_in_cell(const Bare_point& p, Cell_handle c) const
// Returns the finite vertex of the cell c with smaller
// power distance  to p.
{
  CGAL_triangulation_precondition(dimension() >= 1);
  Vertex_handle nearest = nearest_power_vertex(p, c->vertex(0), c->vertex(1));
  if(dimension() >= 2)
  {
    nearest = nearest_power_vertex(p, nearest, c->vertex(2));
    if(dimension() == 3)
      nearest = nearest_power_vertex(p, nearest, c->vertex(3));
  }
  return nearest;
}


template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
nearest_power_vertex(const Bare_point& p, Cell_handle start) const
{
  if(number_of_vertices() == 0)
    return Vertex_handle();

  // Use a brute-force algorithm if dimension < 3.
  if(dimension() < 3)
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    Vertex_handle res = vit;
    ++vit;
    for(Finite_vertices_iterator end = finite_vertices_end(); vit != end; ++vit)
      res = nearest_power_vertex(p, res, vit);

    return res;
  }

  Locate_type lt;
  int li, lj;

  typename Gt::Construct_weighted_point_3 p2wp = geom_traits().construct_weighted_point_3_object();
  Cell_handle c = locate(p2wp(p), lt, li, lj, start);

  // - start with the closest vertex from the located cell.
  // - repeatedly take the nearest of its incident vertices if any
  // - if not, we're done.
  Vertex_handle nearest = nearest_power_vertex_in_cell(p, c);
  std::vector<Vertex_handle> vs;
  vs.reserve(32);
  while(true)
  {
    Vertex_handle tmp = nearest;
    adjacent_vertices_threadsafe(nearest, std::back_inserter(vs));
    for(typename std::vector<Vertex_handle>::const_iterator
         vsit = vs.begin(); vsit != vs.end(); ++vsit)
      tmp = nearest_power_vertex(p, tmp, *vsit);

    if(tmp == nearest)
      break;

    vs.clear();
    nearest = tmp;
  }
  return nearest;
}

template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Bare_point
Regular_triangulation_3<Gt,Tds,Lds>::
dual(Cell_handle c) const
{
  CGAL_triangulation_precondition(dimension()==3);
  CGAL_triangulation_precondition(! is_infinite(c));

  return c->weighted_circumcenter(geom_traits());
}

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_segment(Cell_handle c, int i, Bare_point& p, Bare_point&q) const
{
  Cell_handle n = c->neighbor(i);
  CGAL_assertion(! is_infinite(c) && ! is_infinite(n));

  p = dual(c);
  q = dual(n);
}

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_segment(const Facet& facet, Bare_point& p, Bare_point&q) const
{
  return dual_segment(facet.first, facet.second, p, q);
}

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_ray(Cell_handle c, int i, Ray& ray) const
{
  Cell_handle n = c->neighbor(i);
  CGAL_triangulation_precondition((!is_infinite(c) != !is_infinite(n))); // xor
  // either n or c is infinite
  int in;
  if(is_infinite(c))
  {
    in = n->index(c);
  }
  else
  {
    n = c;
    in = i;
  }

  // n now denotes a finite cell, either c or c->neighbor(i)
  int ind[3] = {(in+1)&3,(in+2)&3,(in+3)&3};
  if((in&1) == 1)
    std::swap(ind[0], ind[1]);

  const Weighted_point& p = n->vertex(ind[0])->point();
  const Weighted_point& q = n->vertex(ind[1])->point();
  const Weighted_point& r = n->vertex(ind[2])->point();
  const Bare_point& bp = construct_point(p);
  const Bare_point& bq = construct_point(q);
  const Bare_point& br = construct_point(r);

  Line l = construct_perpendicular_line(construct_plane(bp, bq, br),
                                        construct_weighted_circumcenter(p,q,r));

  ray = construct_ray(dual(n), l);
}

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_ray(const Facet& facet, Ray& ray) const
{
  return dual_ray(facet.first, facet.second, ray);
}

// Exact versions of dual_segment() and dual_ray() for Mesh_3.
// These functions are really dirty: they assume that the point type is nice enough
// such that EPECK can manipulate it (e.g. convert it to EPECK::Point_3) AND
// that the result of these manipulations will make sense.
template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_segment_exact(const Facet& facet, Bare_point& p, Bare_point&q) const
{
  typedef typename Kernel_traits<Bare_point>::Kernel           K;
  typedef Exact_predicates_exact_constructions_kernel          EK;
  typedef Cartesian_converter<K, EK>                           To_exact;
  typedef Cartesian_converter<EK,K>                            Back_from_exact;

  typedef EK                                                   Exact_Rt;

  To_exact to_exact;
  Back_from_exact back_from_exact;
  Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();

  Cell_handle c = facet.first;
  int i = facet.second;
  Cell_handle n = c->neighbor(i);

  const typename Exact_Rt::Weighted_point_3& cp = to_exact(c->vertex(0)->point());
  const typename Exact_Rt::Weighted_point_3& cq = to_exact(c->vertex(1)->point());
  const typename Exact_Rt::Weighted_point_3& cr = to_exact(c->vertex(2)->point());
  const typename Exact_Rt::Weighted_point_3& cs = to_exact(c->vertex(3)->point());

  const typename Exact_Rt::Weighted_point_3& np = to_exact(n->vertex(0)->point());
  const typename Exact_Rt::Weighted_point_3& nq = to_exact(n->vertex(1)->point());
  const typename Exact_Rt::Weighted_point_3& nr = to_exact(n->vertex(2)->point());
  const typename Exact_Rt::Weighted_point_3& ns = to_exact(n->vertex(3)->point());

  p = back_from_exact(exact_weighted_circumcenter(cp, cq, cr, cs));
  q = back_from_exact(exact_weighted_circumcenter(np, nq, nr, ns));
}

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
dual_ray_exact(const Facet& facet, Ray& ray) const
{
  Cell_handle c = facet.first;
  int i = facet.second;
  Cell_handle n = c->neighbor(i);
  CGAL_triangulation_precondition(!is_infinite(c) != !is_infinite(n)); // xor
  // either n or c is infinite
  int in;
  if(is_infinite(c))
  {
    in = n->index(c);
  }
  else
  {
    n = c;
    in = i;
  }

  // n now denotes a finite cell, either c or c->neighbor(i)
  int ind[3] = {(in+1)&3,(in+2)&3,(in+3)&3};
  if((in&1) == 1)
    std::swap(ind[0], ind[1]);

  // exact part
  typedef typename Kernel_traits<Bare_point>::Kernel           K;
  typedef Exact_predicates_exact_constructions_kernel          EK;
  typedef Cartesian_converter<K, EK>                           To_exact;
  typedef Cartesian_converter<EK,K>                            Back_from_exact;

  typedef EK                                                   Exact_Rt;

  To_exact to_exact;
  Back_from_exact back_from_exact;

  Exact_Rt::Construct_weighted_circumcenter_3 exact_weighted_circumcenter =
      Exact_Rt().construct_weighted_circumcenter_3_object();
  Exact_Rt::Construct_perpendicular_line_3 exact_perpendicular_line =
      Exact_Rt().construct_perpendicular_line_3_object();
  Exact_Rt::Construct_plane_3 exact_plane_3 = Exact_Rt().construct_plane_3_object();
  Exact_Rt::Construct_ray_3 exact_ray_3 = Exact_Rt().construct_ray_3_object();
  Exact_Rt::Construct_point_3 exact_point_3 = Exact_Rt().construct_point_3_object();

  const typename Exact_Rt::Weighted_point_3& p = to_exact(n->vertex(ind[0])->point());
  const typename Exact_Rt::Weighted_point_3& q = to_exact(n->vertex(ind[1])->point());
  const typename Exact_Rt::Weighted_point_3& r = to_exact(n->vertex(ind[2])->point());
  const typename Exact_Rt::Weighted_point_3& s = to_exact(n->vertex(in)->point());

  const typename Exact_Rt::Point_3& bp = exact_point_3(p);
  const typename Exact_Rt::Point_3& bq = exact_point_3(q);
  const typename Exact_Rt::Point_3& br = exact_point_3(r);

  typename Exact_Rt::Line_3 l = exact_perpendicular_line(
                                  exact_plane_3(bp, bq, br),
                                  exact_weighted_circumcenter(p, q, r));

  ray = back_from_exact(exact_ray_3(exact_weighted_circumcenter(p, q, r, s), l));
}

template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Object
Regular_triangulation_3<Gt,Tds,Lds>::
dual(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension()>=2);
  CGAL_triangulation_precondition(! is_infinite(c,i));

  if(dimension() == 2)
  {
    CGAL_triangulation_precondition(i == 3);
    return construct_object(construct_weighted_circumcenter(c->vertex(0)->point(),
                                                            c->vertex(1)->point(),
                                                            c->vertex(2)->point()));
  }

  // dimension() == 3
  Cell_handle n = c->neighbor(i);
  if(! is_infinite(c) && ! is_infinite(n))
  {
    // dual is a segment
    Bare_point bp = dual(c);
    Bare_point np = dual(n);
    return construct_object(construct_segment(bp, np));
  }

  // either n or c is infinite, dual is a ray
  Ray r;
  dual_ray(c, i, r);
  return construct_object(r);
}

template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Object
Regular_triangulation_3<Gt,Tds,Lds>::
dual(const Facet& facet) const
{
  return dual(facet.first, facet.second);
}

template < class Gt, class Tds, class Lds >
template < class Stream>
Stream&
Regular_triangulation_3<Gt,Tds,Lds>::
draw_dual(Stream& os) const
{
  for(Finite_facets_iterator fit = finite_facets_begin(), end = finite_facets_end();
       fit != end; ++fit)
  {
    Object o = dual(*fit);
    if(const Segment* s = object_cast<Segment>(&o))
      os << *s;
    else if(const Ray* r = object_cast<Ray>(&o))
      os << *r;
    else if(const Bare_point* p = object_cast<Bare_point>(&o))
      os << *p;
  }
  return os;
}

template < class Gt, class Tds, class Lds >
Oriented_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_oriented_power_sphere(const Weighted_point& p0,
                              const Weighted_point& p1,
                              const Weighted_point& p2,
                              const Weighted_point& p3,
                              const Weighted_point& p, bool perturb) const
{
  CGAL_triangulation_precondition(orientation(p0, p1, p2, p3) == POSITIVE);

  using namespace boost;

  Oriented_side os = power_test(p0, p1, p2, p3, p);

  if(os != ON_ORIENTED_BOUNDARY || !perturb)
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Weighted_point * points[5] = {&p0, &p1, &p2, &p3, &p};

  std::sort(points, points + 5, typename Tr_Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  for(int i=4; i>1; --i)
  {
    if(points[i] == &p)
      return ON_NEGATIVE_SIDE; // since p0 p1 p2 p3 are non coplanar and positively oriented
    Orientation o;
    if(points[i] == &p3 && (o = orientation(p0,p1,p2,p)) != COPLANAR)
      return o;
    if(points[i] == &p2 && (o = orientation(p0,p1,p,p3)) != COPLANAR)
      return o;
    if(points[i] == &p1 && (o = orientation(p0,p,p2,p3)) != COPLANAR)
      return o;
    if(points[i] == &p0 && (o = orientation(p,p1,p2,p3)) != COPLANAR)
      return o;
  }

  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template < class Gt, class Tds, class Lds >
Bounded_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_power_sphere(Cell_handle c, const Weighted_point& p, bool perturb) const
{
  CGAL_triangulation_precondition(dimension() == 3);
  int i3;
  if(! c->has_vertex(infinite_vertex(), i3))
  {
    return Bounded_side(side_of_oriented_power_sphere(c->vertex(0)->point(),
                                                       c->vertex(1)->point(),
                                                       c->vertex(2)->point(),
                                                       c->vertex(3)->point(),
                                                       p, perturb));
  }

  // else infinite cell :
  int i0,i1,i2;
  if((i3%2) == 1)
  {
    i0 = (i3+1)&3;
    i1 = (i3+2)&3;
    i2 = (i3+3)&3;
  }
  else
  {
    i0 = (i3+2)&3;
    i1 = (i3+1)&3;
    i2 = (i3+3)&3;
  }

  // general case
  Orientation o = orientation(c->vertex(i0)->point(),
                              c->vertex(i1)->point(),
                              c->vertex(i2)->point(), p);
  if(o != ZERO)
    return Bounded_side(o);

  // else p coplanar with i0,i1,i2
  return side_of_bounded_power_circle(c->vertex(i0)->point(),
                                      c->vertex(i1)->point(),
                                      c->vertex(i2)->point(),
                                      p, perturb);
}

template < class Gt, class Tds, class Lds >
Bounded_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_bounded_power_circle(const Weighted_point& p0,
                             const Weighted_point& p1,
                             const Weighted_point& p2,
                             const Weighted_point& p, bool perturb) const
{
  CGAL_triangulation_precondition(coplanar_orientation(p0, p1, p2) != 0);
  if(coplanar_orientation(p0, p1, p2) == POSITIVE)
    return Bounded_side (side_of_oriented_power_circle(p0, p1, p2, p, perturb));

  // Wrong because the low level power test already does a coplanar orientation test.
  // return Bounded_side (- side_of_oriented_power_circle (p0, p2, p1, p, perturb));

  return Bounded_side (side_of_oriented_power_circle(p0, p2, p1, p, perturb));
}


template < class Gt, class Tds, class Lds >
Oriented_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_oriented_power_circle(const Weighted_point& p0,
                              const Weighted_point& p1,
                              const Weighted_point& p2,
                              const Weighted_point& p, bool perturb) const
{
  CGAL_triangulation_precondition(coplanar_orientation(p0, p1, p2) == POSITIVE);

  using namespace boost;

  Oriented_side os = power_test(p0, p1, p2, p);

  if(os != ON_ORIENTED_BOUNDARY || !perturb)
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Weighted_point * points[4] = {&p0, &p1, &p2, &p};

  std::sort(points, points + 4, typename Tr_Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for(int i=3; i>1; --i)
  {
    if(points[i] == &p)
      return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear and positively oriented
    Orientation o;
    if(points[i] == &p2 && (o = coplanar_orientation(p0,p1,p)) != COPLANAR)
      return o;
    if(points[i] == &p1 && (o = coplanar_orientation(p0,p,p2)) != COPLANAR)
      return o;
    if(points[i] == &p0 && (o = coplanar_orientation(p,p1,p2)) != COPLANAR)
      return o;
  }

  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template < class Gt, class Tds, class Lds >
Bounded_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_power_circle(Cell_handle c, int i, const Weighted_point& p,
                     bool perturb) const
{
  CGAL_triangulation_precondition(dimension() >= 2);
  int i3 = 5;
  if(dimension() == 2)
  {
    CGAL_triangulation_precondition(i == 3);
    // the triangulation is supposed to be valid, ie the facet
    // with vertices 0 1 2 in this order is positively oriented
    if(! c->has_vertex(infinite_vertex(), i3))
      return Bounded_side(side_of_oriented_power_circle(c->vertex(0)->point(),
                                                         c->vertex(1)->point(),
                                                         c->vertex(2)->point(),
                                                         p, perturb));
    // else infinite facet
    // v1, v2 finite vertices of the facet such that v1,v2,infinite
    // is positively oriented
    Vertex_handle v1 = c->vertex(ccw(i3)),
                  v2 = c->vertex(cw(i3));
    CGAL_triangulation_assertion(
      coplanar_orientation(v1->point(), v2->point(), mirror_vertex(c, i3)->point()) == NEGATIVE);

    Orientation o = coplanar_orientation(v1->point(), v2->point(), p);
    if(o != ZERO)
      return Bounded_side(o);

    // case when p collinear with v1v2
    return side_of_bounded_power_segment(v1->point(),
                                         v2->point(),
                                         p, perturb);
  } // dim 2

  // else dimension == 3
  CGAL_triangulation_precondition((i >= 0) && (i < 4));
  if((! c->has_vertex(infinite_vertex(),i3)) || (i3 != i))
  {
    // finite facet
    // initialization of i0 i1 i2, vertices of the facet positively
    // oriented (if the triangulation is valid)
    int i0 = (i>0) ? 0 : 1;
    int i1 = (i>1) ? 1 : 2;
    int i2 = (i>2) ? 2 : 3;
    CGAL_triangulation_precondition(this->coplanar(c->vertex(i0)->point(),
                                                   c->vertex(i1)->point(),
                                                   c->vertex(i2)->point(), p));
    return side_of_bounded_power_circle(c->vertex(i0)->point(),
                                        c->vertex(i1)->point(),
                                        c->vertex(i2)->point(),
                                        p, perturb);
  }
  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle v1 = c->vertex(next_around_edge(i3,i)),
                v2 = c->vertex(next_around_edge(i,i3));
  Orientation o = (Orientation)
                    (coplanar_orientation(v1->point(), v2->point(),
                                          c->vertex(i)->point()) *
                     coplanar_orientation(v1->point(), v2->point(), p));
  // then the code is duplicated from 2d case
  if(o != ZERO)
    return Bounded_side(-o);
  // because p is in f iff
  // it is not on the same side of v1v2 as c->vertex(i)
  // case when p collinear with v1v2 :
  return side_of_bounded_power_segment(v1->point(),
                                       v2->point(),
                                       p, perturb);
}

template < class Gt, class Tds, class Lds >
Bounded_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_bounded_power_segment(const Weighted_point& p0,
                              const Weighted_point& p1,
                              const Weighted_point& p, bool perturb) const
{
  Oriented_side os = power_test(p0, p1, p);

  if(os != ON_ORIENTED_BOUNDARY || !perturb)
    return Bounded_side(os);

  // We are now in a degenerate case => we do a symbolic perturbation.
  switch (this->collinear_position(p0, p, p1))
  {
    case Tr_Base::BEFORE:
    case Tr_Base::AFTER:
      return ON_UNBOUNDED_SIDE;
    case Tr_Base::MIDDLE:
      return ON_BOUNDED_SIDE;
    default:
      ;
  }

  CGAL_triangulation_assertion(false);
  return ON_UNBOUNDED_SIDE;
}

template < class Gt, class Tds, class Lds >
Bounded_side
Regular_triangulation_3<Gt,Tds,Lds>::
side_of_power_segment(Cell_handle c, const Weighted_point& p, bool perturb) const
{
  CGAL_triangulation_precondition(dimension() == 1);
  if(! is_infinite(c,0,1))
    return side_of_bounded_power_segment(c->vertex(0)->point(),
                                         c->vertex(1)->point(),
                                         p, perturb);

  Locate_type lt; int i;
  Bounded_side soe = side_of_edge(p, c, lt, i);
  if(soe != ON_BOUNDARY)
    return soe;

  // Either we compare weights, or we use the finite neighboring edge
  Cell_handle finite_neighbor = c->neighbor(c->index(infinite_vertex()));
  CGAL_triangulation_assertion(!is_infinite(finite_neighbor,0,1));
  return side_of_bounded_power_segment(finite_neighbor->vertex(0)->point(),
                                       finite_neighbor->vertex(1)->point(),
                                       p, perturb);
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
greater_or_equal_power_distance(const Bare_point& p,
                                const Weighted_point& q,
                                const Weighted_point& r) const
{
  return ! less_power_distance(p, q, r);
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_Gabriel(const Facet& f) const
{
  return is_Gabriel(f.first, f.second);
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_Gabriel(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i));
  typename Geom_traits::Power_side_of_bounded_power_sphere_3
      side_of_bounded_orthogonal_sphere =
      geom_traits().power_side_of_bounded_power_sphere_3_object();

  if((!is_infinite(c->vertex(i))) &&
     side_of_bounded_orthogonal_sphere(c->vertex(vertex_triple_index(i,0))->point(),
                                       c->vertex(vertex_triple_index(i,1))->point(),
                                       c->vertex(vertex_triple_index(i,2))->point(),
                                       c->vertex(i)->point()) == ON_BOUNDED_SIDE)
    return false;

  Cell_handle neighbor = c->neighbor(i);
  int in = neighbor->index(c);

  if((!is_infinite(neighbor->vertex(in))) &&
     side_of_bounded_orthogonal_sphere(c->vertex(vertex_triple_index(i,0))->point(),
                                       c->vertex(vertex_triple_index(i,1))->point(),
                                       c->vertex(vertex_triple_index(i,2))->point(),
                                       neighbor->vertex(in)->point()) == ON_BOUNDED_SIDE)
    return false;

  return true;
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_Gabriel(const Edge& e) const
{
  return is_Gabriel(e.first, e.second, e.third);
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_Gabriel(Cell_handle c, int i, int j) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i,j));
  typename Geom_traits::Power_side_of_bounded_power_sphere_3
      side_of_bounded_orthogonal_sphere =
      geom_traits().power_side_of_bounded_power_sphere_3_object();

  Facet_circulator fcirc = incident_facets(c,i,j), fdone(fcirc);
  Vertex_handle v1 = c->vertex(i);
  Vertex_handle v2 = c->vertex(j);
  do
  {
    // test whether the vertex of cc opposite to *fcirc
    // is inside the sphere defined by the edge e = (s, i,j)
    Cell_handle cc = (*fcirc).first;
    int ii = (*fcirc).second;
    if(!is_infinite(cc->vertex(ii)) &&
        side_of_bounded_orthogonal_sphere(v1->point(),
                                           v2->point(),
                                           cc->vertex(ii)->point()) == ON_BOUNDED_SIDE)
      return false;
  }
  while(++fcirc != fdone);

  return true;
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_Gabriel(Vertex_handle v) const
{
  typename Geom_traits::Power_side_of_bounded_power_sphere_3
    side_of_bounded_orthogonal_sphere = geom_traits().power_side_of_bounded_power_sphere_3_object();

  Vertex_handle nearest_v = nearest_power_vertex(geom_traits().construct_point_3_object()(v->point()),
                                                                                          v->cell());

  return (side_of_bounded_orthogonal_sphere(v->point(), nearest_v->point()) != CGAL::ON_BOUNDED_SIDE);
}

// Returns
template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
insert(const Weighted_point& p, Cell_handle start, bool *could_lock_zone)
{
  Locate_type lt;
  int li, lj;

  // Parallel
  if(could_lock_zone)
  {
    Cell_handle c = locate(p, lt, li, lj, start, could_lock_zone);
    if(*could_lock_zone)
      return insert(p, lt, c, li, lj, could_lock_zone);
    else
      return Vertex_handle();
  }
  // Sequential
  else
  {
    Cell_handle c = locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
  }
}

template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
insert(const Weighted_point& p, Locate_type lt, Cell_handle c,
       int li, int lj, bool *could_lock_zone)
{
  switch(dimension())
  {
    case 3:
    {
      Conflict_tester_3 tester(p, this);
      return insert_in_conflict(p, lt,c,li,lj, tester,
                                get_hidden_point_visitor(),
                                could_lock_zone);
    }
    case 2:
    {
      Conflict_tester_2 tester(p, this);
      return insert_in_conflict(p, lt,c,li,lj, tester,
                                get_hidden_point_visitor(),
                                could_lock_zone);
    }
    case 1:
    {
      Conflict_tester_1 tester(p, this);
      return insert_in_conflict(p, lt,c,li,lj, tester,
                                get_hidden_point_visitor(),
                                could_lock_zone);
    }
  }

  Conflict_tester_0 tester(p, this);
  return insert_in_conflict(p, lt,c,li,lj, tester,
                            get_hidden_point_visitor(),
                            could_lock_zone);
}

template < class Gt, class Tds, class Lds >
template <class CellIt>
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
insert_in_hole(const Weighted_point& p, CellIt cell_begin, CellIt cell_end,
               Cell_handle begin, int i)
{
  CGAL_triangulation_precondition(cell_begin != cell_end);

  get_hidden_point_visitor().process_cells_in_conflict(cell_begin,cell_end);

  Vertex_handle v = Tr_Base::insert_in_hole(p, cell_begin, cell_end, begin, i);

  // Store the hidden points in their new cells and hide vertices that
  // have to be hidden
  get_hidden_point_visitor().reinsert_vertices(v);
  return v;
}

template < class Gt, class Tds, class Lds >
template <class CellIt>
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
insert_in_hole(const Weighted_point& p, CellIt cell_begin, CellIt cell_end,
               Cell_handle begin, int i, Vertex_handle newv)
{
  CGAL_triangulation_precondition(cell_begin != cell_end);

  get_hidden_point_visitor().process_cells_in_conflict(cell_begin,cell_end);

  Vertex_handle v = Tr_Base::insert_in_hole(p, cell_begin, cell_end, begin, i, newv);

  // Store the hidden points in their new cells and hide vertices that
  // have to be hidden
  get_hidden_point_visitor().reinsert_vertices(v);
  return v;
}

template <class Gt, class Tds, class Lds >
template <class RegularTriangulation_3>
class Regular_triangulation_3<Gt, Tds, Lds>::Vertex_remover
{
  typedef RegularTriangulation_3                          Regular;
  typedef typename Regular::Weighted_point                Weighted_point;

public:
  typedef typename std::vector<Weighted_point>::iterator  Hidden_points_iterator;

  Vertex_remover(Regular &tmp_) : tmp(tmp_) {}

  Regular& tmp;

  void add_hidden_points(Cell_handle ch)
  {
    std::copy(ch->hidden_points_begin(), ch->hidden_points_end(),
              std::back_inserter(hidden));
  }

  Hidden_points_iterator hidden_points_begin() { return hidden.begin(); }
  Hidden_points_iterator hidden_points_end() { return hidden.end(); }

  Bounded_side side_of_bounded_circle(const Weighted_point& p, const Weighted_point& q,
                                      const Weighted_point& r, const Weighted_point& s, bool perturb = false) const
  {
    return tmp.side_of_bounded_power_circle(p,q,r,s,perturb);
  }

private:
  // The removal of v may un-hide some points,
  // Space functions output them.
  std::vector<Weighted_point> hidden;
};

// The displacement method works only
// on regular triangulation without hidden points at any time
// the vertex inserter is used only
// for the purpose of displacements
template <class Gt, class Tds, class Lds >
template <class RegularTriangulation_3>
class Regular_triangulation_3<Gt, Tds, Lds>::Vertex_inserter
{
  typedef RegularTriangulation_3              Regular;

public:
  typedef std::nullptr_t                           Hidden_points_iterator;

  Vertex_inserter(Regular &tmp_) : tmp(tmp_) {}

  Regular& tmp;

  void add_hidden_points(Cell_handle) {}
  Hidden_points_iterator hidden_points_begin() { return nullptr; }
  Hidden_points_iterator hidden_points_end() { return nullptr; }

  Vertex_handle insert(const Weighted_point& p,
                       Locate_type lt, Cell_handle c, int li, int lj)
  {
    return tmp.insert(p, lt, c, li, lj);
  }

  Vertex_handle insert(const Weighted_point& p, Cell_handle c)
  {
    return tmp.insert(p, c);
  }

  Vertex_handle insert(const Weighted_point& p)
  {
    return tmp.insert(p);
  }
};

template < class Gt, class Tds, class Lds >
void
Regular_triangulation_3<Gt,Tds,Lds>::
remove(Vertex_handle v)
{
  Cell_handle c;
  if(dimension() > 0)
    c = v->cell()->neighbor(v->cell()->index(v));

  Self tmp;
  Vertex_remover<Self> remover(tmp);
  Tr_Base::remove(v,remover);

  // Re-insert the points that v was hiding.
  for(typename Vertex_remover<Self>::Hidden_points_iterator
       hi = remover.hidden_points_begin();
       hi != remover.hidden_points_end(); ++hi)
  {
    Vertex_handle hv = insert(*hi, c);
    if(hv != Vertex_handle())
      c = hv->cell();
  }
  CGAL_triangulation_expensive_postcondition(is_valid());
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
remove(Vertex_handle v, bool *could_lock_zone)
{
  bool removed = true;

  // Locking vertex v...
  if(!this->try_lock_vertex(v))
  {
    *could_lock_zone = false;
  }
  else
  {
    Vertex_validity_checker<Tds> vertex_validity_check;

    // Check that the vertex hasn't be deleted from the TDS while we were locking it
    if(!vertex_validity_check(v, tds()))
      return true; // vertex is already gone from the TDS, nothing to do

    Vertex_handle hint = v->cell()->vertex(0) == v ? v->cell()->vertex(1) : v->cell()->vertex(0);

    Self tmp;
    Vertex_remover<Self> remover(tmp);
    removed = Tr_Base::remove(v, remover, could_lock_zone);

    if(*could_lock_zone && removed)
    {
      // The vertex has been removed, re-insert the points that 'v' was hiding

      // Start by unlocking the area of the removed vertex to avoid deadlocks
      this->unlock_all_elements();

      for(typename Vertex_remover<Self>::Hidden_points_iterator
            hi = remover.hidden_points_begin();
            hi != remover.hidden_points_end(); ++hi)
      {
        const Weighted_point& wp = *hi;

        // try to lock the positions of the hint and the hidden point
        bool success = false;
        while(!success)
        {
          // The 'hint' is unsafe to use immediately because we are in a regular triangulation,
          // and the insertion of a (weighted) point in another thread might have hidden (deleted)
          // the hint.
          if(!vertex_validity_check(hint, tds()))
          {
            hint = finite_vertices_begin();
            continue;
          }

          // We need to make sure that while are locking the position P1 := hint->point(), 'hint'
          // does not get its position changed to P2 != P1.
          const Weighted_point hint_point_mem = hint->point();

          if(this->try_lock_point(hint_point_mem) && this->try_lock_point(wp))
          {
            // Make sure that the hint is still valid (so that we can safely take hint->cell()) and
            // that its position hasn't changed to ensure that we will start the locate from where
            // we have locked.
            if(!vertex_validity_check(hint, tds()) ||
               hint->point() != hint_point_mem)
            {
              hint = finite_vertices_begin();
              this->unlock_all_elements();
              continue;
            }

            Vertex_handle hv = insert(wp, hint, could_lock_zone);

            if(*could_lock_zone)
            {
              success = true;
              if(hv != Vertex_handle())
                hint = hv;
            }
          }

          // This unlocks everything in all cases: partial lock failure, failed insertion, successful insertion
          this->unlock_all_elements();
        }
      }

      CGAL_triangulation_expensive_postcondition(is_valid());
    }
  }

  return removed;
}


// Displacement works only for regular triangulation
// without hidden points at any time
template < class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
move_if_no_collision(Vertex_handle v, const Weighted_point& p)
{
  Self tmp;
  Vertex_remover<Self> remover(tmp);
  Vertex_inserter<Self> inserter(*this);
  Vertex_handle res = Tr_Base::move_if_no_collision(v,p,remover,inserter);

  CGAL_triangulation_expensive_postcondition(is_valid());
  return res;
}

template <class Gt, class Tds, class Lds >
typename Regular_triangulation_3<Gt,Tds,Lds>::Vertex_handle
Regular_triangulation_3<Gt,Tds,Lds>::
move(Vertex_handle v, const Weighted_point& p)
{
  CGAL_triangulation_precondition(!is_infinite(v));
  if(v->point() == p)
    return v;

  Self tmp;
  Vertex_remover<Self> remover(tmp);
  Vertex_inserter<Self> inserter(*this);
  return Tr_Base::move(v,p,remover,inserter);
}

template < class Gt, class Tds, class Lds >
bool
Regular_triangulation_3<Gt,Tds,Lds>::
is_valid(bool verbose, int level) const
{
  if(! Tr_Base::is_valid(verbose,level))
  {
    if(verbose)
      std::cerr << "invalid base triangulation" << std::endl;

    CGAL_triangulation_assertion(false);
    return false;
  }

  switch(dimension())
  {
    case 3:
    {
      for(Finite_cells_iterator it = finite_cells_begin(),
                                end = finite_cells_end(); it != end; ++it)
      {
        is_valid_finite(it, verbose, level);
        for(int i=0; i<4; i++)
        {
          if(!is_infinite(it->neighbor(i)->vertex(it->neighbor(i)->index(it))))
          {
            if(side_of_power_sphere(it,
                                    it->neighbor(i)->vertex(
                                      it->neighbor(i)->index(it))->point()) == ON_BOUNDED_SIDE)
            {
              if(verbose)
                std::cerr << "non-empty sphere " << std::endl;

              CGAL_triangulation_assertion(false);
              return false;
            }
          }
        }
      }
      break;
    }
    case 2:
    {
      for(Finite_facets_iterator it = finite_facets_begin(), end = finite_facets_end(); it!= end; ++it)
      {
        is_valid_finite((*it).first, verbose, level);
        for(int i=0; i<3; i++)
        {
          if(!is_infinite((*it).first->neighbor(i)->vertex(
                            (((*it).first)->neighbor(i))->index((*it).first))))
          {
            if(side_of_power_circle((*it).first, 3,
                                    (*it).first->neighbor(i)->
                                    vertex((((*it).first)->neighbor(i))
                                           ->index((*it).first))->point()) == ON_BOUNDED_SIDE)
            {
              if(verbose)
                std::cerr << "non-empty circle " << std::endl;

              CGAL_triangulation_assertion(false);
              return false;
            }
          }
        }
      }
      break;
    }
    case 1:
    {
      for(Finite_edges_iterator it = finite_edges_begin(),
                                end = finite_edges_end(); it != end; ++it)
      {
        is_valid_finite((*it).first, verbose, level);
        for(int i=0; i<2; i++)
        {
          if(!is_infinite
              ((*it).first->neighbor(i)->vertex(
                 (((*it).first)->neighbor(i))->index((*it).first))))
          {
            if(side_of_power_segment((*it).first,
                                     (*it).first->neighbor(i)->vertex(
                                       (((*it).first)->neighbor(i))->index(
                                         (*it).first))->point()) == ON_BOUNDED_SIDE)
            {
              if(verbose)
                std::cerr << "non-empty edge " << std::endl;

              CGAL_triangulation_assertion(false);
              return false;
            }
          }
        }
      }
      break;
    }
  }

  if(verbose)
    std::cerr << "valid regular triangulation" << std::endl;

  return true;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_REGULAR_TRIANGULATION_3_H
