// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Léo Valque, Sylvain Lazard

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SNAP_ROUNDING_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SNAP_ROUNDING_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Fraction_traits.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/mutex.h>

namespace CGAL
{

namespace Polygon_mesh_processing
{

namespace internal
{

template <class NT> double double_ceil(Lazy_exact_nt< NT > x);
template <class NT> double double_ceil(NT x);

template <class NT>
double double_ceil(Lazy_exact_nt< NT > x){
  // If both sides are in the same ceil, return this ceil
  double ceil_left=std::ceil(to_interval(x).first);
  if(ceil_left==std::ceil(to_interval(x).second))
    return ceil_left;
  // If not refine interval by contracting the DAG and try again
  x.exact();
  ceil_left=std::ceil(to_interval(x).first);
  if(ceil_left==std::ceil(to_interval(x).second))
    return ceil_left;
  // If not return the ceil of the exact value
  return double_ceil( x.exact());
};

template <class NT>
double double_ceil(NT x){
  using FT = Fraction_traits<NT>;
  if constexpr(std::is_same<typename FT::Is_fraction, Tag_true>::value){
    // If NT is a fraction, the ceil value is the result of the euclidian division of the numerator and the denominator.
    typename FT::Numerator_type num, r, e;
    typename FT::Denominator_type denom;
    typename FT::Decompose()(x,num,denom);
    div_mod(num, denom, r, e);
    if((r>=0) && e!=0) //If the result is positive, the ceil value is one above
      return to_double(r+1);
    return to_double(r);
  } else {
    // Return the ceil of the approximation
    return std::ceil(to_double(x));
  }
};

template <typename Range>
class Indexes_range : public Range{
  typedef std::remove_cv_t<typename std::iterator_traits<typename Range::iterator>::value_type> Value_type;
public:
  typedef typename Range::const_iterator const_iterator;
  typedef typename Range::iterator iterator;

  Indexes_range(){}
  Indexes_range(std::initializer_list<size_t>  l): Range(l), m_id(0), modified(true){}
  Indexes_range(const Indexes_range<Range> &ir): Range(ir), m_id(ir.id()), modified(ir.was_modified()){}
  Indexes_range(Range &p): Range(p), modified(true){}
  Indexes_range(Range &p, size_t id): Range(p), m_id(id),modified(false){}

  void operator=(const Indexes_range<Range> &ir)
  {
    Range::operator=(ir);
    modified=ir.was_modified();
    m_id=ir.id();
  }

  inline size_t id() const { return m_id; }
  inline void set_id(size_t id){ m_id=id; }
  inline bool was_modified() const { return modified; }

private:
  size_t m_id;
  bool modified;
};

//repair_polygon_soup for triangles
template <typename PointRange, typename PolygonRange,
          typename Polygon = typename Polygon_types<PointRange, PolygonRange>::Polygon_3,
          typename NamedParameters>
void repair_triangle_soup(PointRange& points,
                          PolygonRange& polygons,
                          const NamedParameters& np)
{
  merge_duplicate_points_in_polygon_soup(points, polygons, np);
  remove_invalid_polygons_in_polygon_soup(points, polygons);
  merge_duplicate_polygons_in_polygon_soup(points, polygons, np);
  remove_isolated_points_in_polygon_soup(points, polygons);
}

struct Wrapp_id_visitor : public Autorefinement::Default_visitor
{
  template< typename Triangle>
  inline void internal_new_subtriangle(Triangle& new_t, const Triangle& old_t) {
    new_t.set_id(old_t.id());
  }
};

/**
*
* Rounds the coordinates of the points so that they fit in doubles while making and keeping the model intersection free by potentially subdividing the triangles.
* The input can be any triangle soup and the output is an intersection-free triangle soup with Hausdorff distance
* between the input and the output bounded by M*2^-gs*k where M is the maximum absolute coordinate in the model, gs the snap_grid_size (see below) and k the number of iteration
* performed by the algorithm.
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam TriangleRange a model of the concepts `RandomAccessContainer`, `BackInsertionSequence` and `Swappable`, whose
* value type is a model of the concept `RandomAccessContainer` whose value type is convertible to `std::size_t` and that
* is constructible from an `std::initializer_list<std::size_t>` of size 3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param soup_triangles each element in the range describes a triangle using the indexed position of the points in `soup_points`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
*     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
*     \cgalParamDefault{`CGAL::Sequential_tag`}
*   \cgalParamNEnd
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `soup_points`}
*     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
*   \cgalParamNEnd
*   \cgalParamNBegin{snap_grid_size}
*     \cgalParamDescription{Scale the points to [-2^gs, 2^gs] where gs is the snap_grid_size before to round them on integer.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{23}
*     \cgalParamExtra{Must be lower than 52.}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{Maximum number of iteration performed by the algorithm.}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{5}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \return `true` if the modified triangle soup is free from self-intersection, and `false` if the algorithm was not
* able to provide such a triangle soup within the number of iterations.
*/
template <typename PointRange, typename PolygonRange, class NamedParameters>
bool polygon_soup_snap_rounding_impl(PointRange &points,
                                     PolygonRange &triangles,
                                     const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // typedef typename GetPolygonSoupGeomTraits<PointRange, NamedParameters>::type GT;
  typedef typename GetPointMap<PointRange, NamedParameters>::const_type    Point_map;
  Point_map pm = choose_parameter<Point_map>(get_parameter(np, internal_np::point_map));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  // visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Autorefinement::Default_visitor//default
  >::type Visitor;
  Visitor visitor(choose_parameter<Visitor>(get_parameter(np, internal_np::visitor)));

  constexpr bool has_visitor = !std::is_same_v<Autorefinement::Default_visitor, Visitor>;
  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;
  const size_t number_of_input_triangles = triangles.size();

#ifndef CGAL_LINKED_WITH_TBB
  static_assert (!parallel_execution,
                 "Parallel_tag is enabled but TBB is unavailable.");
#endif

  using Point_3 = std::remove_cv_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Kernel = typename Kernel_traits<Point_3>::Kernel;

  // Get the grid size from the named parameter, the grid size could not be greater than 52
  const unsigned int grid_size = (std::min)(52,choose_parameter(get_parameter(np, internal_np::snap_grid_size), 23));
  const unsigned int max_nb_of_iteration = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 5);

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
  std::cout << "Compute the scaling of the coordinates" << std::endl;
  std::cout << grid_size << std::endl;
#endif

  auto exp = [](const double v)
  {
    int n;
    frexp(v, &n);
    return n;
  };
  auto pow_2 = [](const int k)
  {
      return (k>=0)?std::pow(2,k):1./std::pow(2,-k);
  };

  // Get max absolute value of each dimension
  Bbox_3 bb = bbox_3(points.begin(), points.end());
  std::array<double, 3> max_abs{(std::max)(-bb.xmin(), bb.xmax()),
                                (std::max)(-bb.ymin(), bb.ymax()),
                                (std::max)(-bb.zmin(), bb.zmax())};
  // Compute scale so that the exponent of max absolute value are 52-1.
  std::array<double, 3> scale{pow_2(grid_size - exp(max_abs[0]) - 1),
                              pow_2(grid_size - exp(max_abs[1]) - 1),
                              pow_2(grid_size - exp(max_abs[2]) - 1)};

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
  std::cout << "Scaling: " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
#endif

  auto snap = [](typename Kernel::FT x, double scale)
  {
    // Scale the coordinate, round to nearest integer and scale back
    // TODO replace this ceil by the one of Algebraic_fondation when it will be add to master
    return internal::double_ceil((x * scale) - 0.5) / scale;
  };
  auto snap_p = [scale, snap](const Point_3 &p)
  {
    return Point_3( snap(p.x(),scale[0]),
                    snap(p.y(),scale[1]),
                    snap(p.z(),scale[2]) );
  };

  for (std::size_t i = 0; i <= max_nb_of_iteration; ++i)
  {
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Start Iteration " << i << std::endl;
    std::cout << "Model size: " << points.size() << " " << triangles.size() << std::endl;
#endif

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Round all coordinates on doubles" << std::endl;
#endif

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr(parallel_execution)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()),
                        [&](const tbb::blocked_range<size_t>& r){
                          for(size_t pi = r.begin(); pi != r.end(); ++pi)
                            points[pi] = Point_3(to_double(points[pi].x()), to_double(points[pi].y()), to_double(points[pi].z()));
                        }
                        );
    } else
#endif
    for (Point_3 &p : points)
      p = Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));

    repair_triangle_soup(points, triangles, np);

    // Get all intersecting triangles
    std::vector<std::pair<std::size_t, std::size_t>> pairs_of_intersecting_triangles;
    triangle_soup_self_intersections<Concurrency_tag>(points, triangles, std::back_inserter(pairs_of_intersecting_triangles), np);

    if (pairs_of_intersecting_triangles.empty())
    {
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "End of the snapping" << std::endl;
#endif
      CGAL_assertion(!does_triangle_soup_self_intersect<Concurrency_tag>(points, triangles));
      if constexpr(has_visitor)
      {
        std::vector<std::vector<size_t> > map_io(number_of_input_triangles);
        size_t id=0;
        for(auto &t: triangles)
          map_io[t.id()].push_back(id++);

        visitor.number_of_output_triangles(triangles.size());
        for(size_t src_id=0; src_id!=map_io.size(); ++src_id){
          if(map_io[src_id].size()==0)
            visitor.delete_triangle(src_id);
          else if(map_io[src_id].size()==1 && !triangles[map_io[src_id][0]].was_modified())
            visitor.verbatim_triangle_copy(map_io[src_id][0],src_id);
          else
            for(size_t new_id: map_io[src_id])
              visitor.new_subtriangle(new_id,src_id);
        }
      }
      return true;
    }

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Number of pairs of intersecting triangles: " << pairs_of_intersecting_triangles.size() << std::endl;
#endif

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Snap the coordinates of the vertices of the intersecting triangles" << std::endl;
#endif

#if defined(PMP_ROUNDING_VERTICES_NAIVE_SNAP)
    // Nothing

#elif defined(PMP_ROUNDING_VERTICES_CLOSEST_POINT_SNAP)
    // Version where points in a voxel are rounded to the closest point

    // Group the points of the vertices of the intersecting triangles by their voxel
    std::map<Point_3, size_t> snap_points;
    std::size_t index=0;
    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
      for (size_t pi : triangles[pair.second])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
    }

    std::vector<std::vector<size_t>> identical_points(index);
    for (size_t i=0; i!=points.size(); ++i)
    {
      Point_3 p_snap = snap_p(points[i]);
      auto it=snap_points.find(p_snap);
      if (it!=snap_points.end()){
        identical_points[it->second].push_back(i);
      }
    }

    // Replace all points in a voxel by the closest point
    for(const auto &v: identical_points){
      if(v.size()>1){
        Point_3 center = snap_p(points[v[0]]);
        size_t argmin(0);
        for(size_t i=1; i!=v.size(); ++i){
          if(Kernel().compare_squared_distance_3_object()(center, points[v[i]], center, points[v[argmin]])==SMALLER)
          {
            argmin=i;
          }
        }
        for(auto i: v){
          points[i]=points[v[argmin]];
        }
      }
    }

#elif defined(PMP_ROUNDING_VERTICES_BARYCENTER_SNAP)
    // Version where points in a voxel are rounded to their barycenter.

    // Group the points of the vertices of the intersecting triangles by their voxel
    std::map<Point_3, size_t> snap_points;
    std::size_t index=0;
    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
      for (size_t pi : triangles[pair.second])
      {
        auto res=snap_points.emplace(snap_p(points[pi]), index);
        if(res.second)
          ++index;
      }
    }

    std::vector<std::vector<size_t>> identical_points(index);
    for (size_t i=0; i!=points.size(); ++i)
    {
      Point_3 p_snap = snap_p(points[i]);
      auto it=snap_points.find(p_snap);
      if (it!=snap_points.end()){
        identical_points[it->second].push_back(i);
      }
    }

    // Replace all points in a voxel by their barycenter
    for(const auto &v: identical_points){
      if(v.size()>1){
        std::array<double, 3> a{0,0,0};
        for(auto i: v){
          a[0]+=to_double(points[i].x());
          a[1]+=to_double(points[i].y());
          a[2]+=to_double(points[i].z());
        }
        a[0]/=v.size();
        a[1]/=v.size();
        a[2]/=v.size();
        for(auto i: v){
          points[i]=Point_3(a[0],a[1],a[2]);
        }
      }
    }
#elif defined(PMP_ROUNDING_VERTICES_FLOAT_SNAP)
    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
        points[pi]=Point_3( (float) to_double(points[pi].x()), (float) to_double(points[pi].y()), (float) to_double(points[pi].z()));
      for (size_t pi : triangles[pair.second])
        points[pi]=Point_3( (float) to_double(points[pi].x()), (float) to_double(points[pi].y()), (float) to_double(points[pi].z()));
    }
#else // Default Version
    // Version where points are rounded to the center of their voxel.

    // Get all the snap version of the points of the vertices of the intersecting triangles
    // Note: points will not be modified here, they will be modified in the next for loop

    //TODO: TBB version of this for loop

    std::vector<Point_3> snap_points;
    snap_points.reserve(pairs_of_intersecting_triangles.size() * 3);

    for (auto &pair : pairs_of_intersecting_triangles)
    {
      for (size_t pi : triangles[pair.first])
        snap_points.emplace_back(snap_p(points[pi]));
      for (size_t pi : triangles[pair.second])
        snap_points.emplace_back(snap_p(points[pi]));
    }

#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Snap the coordinates of the vertices close-by the previous ones" << std::endl;
#endif

    std::sort(snap_points.begin(), snap_points.end());
    snap_points.erase(std::unique(snap_points.begin(), snap_points.end()), snap_points.end());

    // If the snapped version of a point correspond to one of the previous point, we snap it
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr(parallel_execution)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()),
                        [&](const tbb::blocked_range<size_t>& r){
                          for(size_t pi = r.begin(); pi != r.end(); ++pi){
                            Point_3 p_snap=snap_p(points[pi]);
                            if (std::binary_search(snap_points.begin(), snap_points.end(), p_snap))
                              points[pi] = p_snap;
                          }
                        }
                        );
    } else
#endif
    for (Point_3 &p : points)
    {
      Point_3 p_snap = snap_p(p);
      if (std::binary_search(snap_points.begin(), snap_points.end(), p_snap))
        p = p_snap;
    }
#endif

    repair_triangle_soup(points, triangles, np);
#ifdef PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
    std::cout << "Model size: " << points.size() << " " << triangles.size() << std::endl;
    std::cout << "Autorefine the soup" << std::endl;
#endif
    if constexpr(has_visitor)
    {
      Wrapp_id_visitor visitor;
      autorefine_triangle_soup(points, triangles, parameters::point_map(pm).concurrency_tag(Concurrency_tag()).visitor(visitor));
    }
    else
    {
      autorefine_triangle_soup(points, triangles, parameters::point_map(pm).concurrency_tag(Concurrency_tag()));
    }
  }
  return false;
}

template <typename PointRange, typename PolygonRange, class NamedParameters>
bool polygon_soup_snap_rounding(PointRange &soup_points,
                                PolygonRange &soup_triangles,
                                const NamedParameters& np)
{
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Autorefinement::Default_visitor//default
  > ::type Visitor;

  constexpr bool has_visitor = !std::is_same_v<Autorefinement::Default_visitor, Visitor>;
  if constexpr(has_visitor)
  {
    //If a visitor is provided, we color the triangles with the index of their source to correctly track the modification
    using Triangle = std::remove_cv_t<typename std::iterator_traits<typename PolygonRange::iterator>::value_type>;
    std::vector<Indexes_range<Triangle> > indexes_soup_triangles;
    size_t id=0;
    for(typename PolygonRange::iterator it=soup_triangles.begin(); it!=soup_triangles.end(); ++it)
      indexes_soup_triangles.emplace_back((*it), id++);

    bool res=polygon_soup_snap_rounding_impl(soup_points, indexes_soup_triangles, np);

    soup_triangles.clear();
    for(const Indexes_range<Triangle> &t: indexes_soup_triangles)
      soup_triangles.push_back({t[0],t[1],t[2]});
    return res;
  }
  else
  {
    return polygon_soup_snap_rounding_impl(soup_points, soup_triangles, np);
  }
}

} } } //end of CGAL::Polygon_mesh_processing::internal namespace

#endif //CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SNAP_ROUNDING_H
