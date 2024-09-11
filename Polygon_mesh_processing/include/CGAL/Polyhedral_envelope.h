// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: ( GPL-3.0-or-later OR LicenseRef-Commercial ) AND MIT
//
// Author(s)     : Andreas Fabri
//
// This file incorporates work covered by the following copyright and permission notice:
//
//     MIT License
//
//     Copyright (c) 2019 Bolun Wang, Teseo Schneider, Yixin Hu, Marco Attene, and Daniele Panozzo
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
//
// @article{Wang:2020:FE,
//    title={Exact and Efficient Polyhedral Envelope Containment Check},
//    author={Bolun Wang and Teseo Schneider and Yixin Hu and Marco Attene and Daniele Panozzo},
//    journal = {ACM Trans. Graph.},
//     volume = {39},
//     number = {4},
//     month = jul,
//     year = {2020},
//     publisher = {ACM}
// }
//
// The code below uses the version of
// https://github.com/wangbolun300/fast-envelope available on 7th of October 2020.
//
// The code below only uses the high-level algorithms checking that a query
// is covered by a set of prisms, where each prism is the offset of an input triangle.
// That is, we do not use indirect predicates.

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYHEDRAL_ENVELOPE_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYHEDRAL_ENVELOPE_H

#include <CGAL/license/Polygon_mesh_processing/Polyhedral_envelope.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_primitive.h>

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/assertions.h>

#ifdef CGAL_ENVELOPE_DEBUG
// This is for computing the surface mesh of a prism
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#endif

#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/has_range_iterator.hpp>

#include <string>
#include <fstream>
#include <type_traits>
#include <unordered_set>
#include <unordered_map>

namespace CGAL {

/**
 * \ingroup PkgPolygonMeshProcessingRef
 * This class can be used to check if a query point, segment, or triangle
 * is inside or outside a polyhedral envelope of a set of triangles, constructed for a given  \f$ \epsilon \f$ distance tolerance.
 * The polyhedral envelope is the union of <em>prisms</em> obtained. See Section \ref PMPEnvelope for more details.
 *
 * @tparam GeomTraits a geometric traits class, model of `Kernel`
 */


template <typename GeomTraits>
struct Polyhedral_envelope {

public:
  typedef typename GeomTraits::Point_3 Point_3;
#ifndef DOXYGEN_RUNNING
  typedef std::array<int, 3> Vector3i;
#endif

private:
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::Segment_3 Segment_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename GeomTraits::Plane_3 Plane_3;
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

  typedef Bbox_3 Bbox;

  typedef Exact_predicates_exact_constructions_kernel EK;
  typedef typename EK::Point_3 ePoint_3;
  typedef typename EK::Line_3 eLine_3;
  typedef typename EK::Plane_3 ePlane_3;
  typedef typename EK::Intersect_3 eIntersect_3;
  typedef typename EK::Oriented_side_3 eOriented_side_3;
  typedef typename EK::Are_parallel_3 eAre_parallel_3;


  template <typename K2>
  static int obtuse_angle(const CGAL::Point_3<K2>& p, const CGAL::Point_3<K2>& q, const CGAL::Point_3<K2>& r)
  {
    if(angle(r,p,q) == OBTUSE){
      return 0;
    }
    if(angle(p,q,r) == OBTUSE){
      return 1;
    }
    if(angle(q,r,p) == OBTUSE){
      return 2;
    }
    return -1;
  }


  template <typename K2>
  static CGAL::Vector_3<K2> normalize(const CGAL::Vector_3<K2>& v)
  {
    return v / approximate_sqrt(v*v);
  }


  // The class Triangle represents a query triangle
  struct Triangle {
    Triangle(const Triangle&) = delete;  // no need for making copies

    Triangle(const Point_3& t0, const Point_3& t1, const Point_3& t2)
    {
      triangle = { t0, t1, t2 };
      etriangle = { ePoint_3(t0.x(), t0.y(), t0.z()),
                    ePoint_3(t1.x(), t1.y(), t1.z()),
                    ePoint_3(t2.x(), t2.y(), t2.z()) };
      n = etriangle[0] + cross_product((etriangle[0] - etriangle[1]), (etriangle[0] - etriangle[2]));
    }

    void init_elines()
    {
      elines = { eLine_3(etriangle[1], etriangle[2]),
                 eLine_3(etriangle[0], etriangle[2]),
                 eLine_3(etriangle[0], etriangle[1]) };
      eplane = ePlane_3(etriangle[0], etriangle[1], etriangle[2]);
    }


    std::array<Point_3,3> triangle;
    std::array<ePoint_3,3> etriangle;
    ePlane_3 eplane;
    std::array<eLine_3,3> elines;  // the i'th line is opposite to vertex i
    ePoint_3 n; // triangle[0] offsetted by the triangle normal
  };


  // The class  `Plane` is used for the 7-8 walls of a prism.
  // We store at the same time three points and a plane.
  // That is easier than retrieving the 3 points of a lazy plane.
  struct Plane {
    Plane()
    {}

    Plane(const ePoint_3& ep, const ePoint_3& eq, const ePoint_3& er)
      : ep(ep), eq(eq), er(er), eplane(ep,eq,er)
    {}

    template <class Point>
    Plane(const Point& p, const Point& q, const Point& r,
          std::enable_if_t<!std::is_same<Point,ePoint_3>::value>* = 0)
      : ep(p.x(),p.y(),p.z()), eq(q.x(),q.y(),q.z()), er(r.x(),r.y(),r.z()), eplane(ep,eq,er)
    {}
    ePoint_3 ep, eq, er;
    ePlane_3 eplane;
  };

  struct Prism {

    std::size_t size() const
    {
      return planes.size();
    }

    void reserve(std::size_t n)
    {
      planes.reserve(n);
    }

    void emplace_back(const Plane& p)
    {
      planes.emplace_back(p);
    }


    Plane& operator[](std::size_t i)
    {
      return planes[i];
    }

    const Plane& operator[](std::size_t i) const
    {
      return planes[i];
    }

    std::vector<Plane> planes;
    int obtuse;
  };

  static const bool OUT_PRISM = 1;
  static const bool IN_PRISM = 0;
  static const int CUT_COPLANAR = 4;
  static const int CUT_EMPTY = -1;
  static const int CUT_FACE = 3;

  // For a query object the envelope test uses an AABB tree to find the relevant prisms

  // property maps for the primitive
  template <class Kernel>
  struct Datum_map
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t key_type;
    typedef typename Kernel::Iso_cuboid_3 value_type;
    typedef value_type reference;

    const std::vector<Iso_cuboid_3>* boxes_ptr;

    Datum_map() : boxes_ptr(nullptr) {}
    Datum_map(const std::vector<Iso_cuboid_3>& boxes) : boxes_ptr(&boxes) {}

    friend value_type get(const Datum_map& m, key_type k)
    {
      CGAL_assertion( m.boxes_ptr->size()>k );
      return (*m.boxes_ptr)[k];
    }
  };

  template <class Kernel>
  struct Point_map
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t key_type;
    typedef typename Kernel::Point_3 value_type;
    typedef value_type reference;

    const std::vector<Iso_cuboid_3>* boxes_ptr;

    Point_map() : boxes_ptr(nullptr) {}
    Point_map(const std::vector<Iso_cuboid_3>& boxes) : boxes_ptr(&boxes) {}

    friend value_type get(const Point_map& m, key_type k)
    {
      CGAL_assertion( m.boxes_ptr->size()>k );
      return ((*m.boxes_ptr)[k].min)();
    }
  };

  typedef AABB_primitive<unsigned int, Datum_map<GeomTraits>, Point_map<GeomTraits>, Tag_true /*UseSharedData*/, Tag_false /*CacheDatum*/> Primitive;
  typedef AABB_traits_3<GeomTraits, Primitive> Tree_traits;
  typedef AABB_tree<Tree_traits> Tree;

// Data members
  std::vector<Prism> halfspace; // should be renamed to "prisms"
  std::vector<Iso_cuboid_3> bounding_boxes;
  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;

  Tree tree;

  eOriented_side_3 oriented_side;

public:

  /// \name Initialization
  /// @{

  /**
   * Default constructor, envelope is empty
   */
  Polyhedral_envelope()
  {}

  // Disabled copy constructor & assignment operator
  Polyhedral_envelope(const Polyhedral_envelope<GeomTraits>&) = delete;
  Polyhedral_envelope<GeomTraits>& operator=(const Polyhedral_envelope<GeomTraits>&) = delete;

  Polyhedral_envelope<GeomTraits>& operator=(Polyhedral_envelope<GeomTraits>&& other) noexcept
  {
    halfspace = std::move(other.halfspace);
    bounding_boxes = std::move(other.bounding_boxes);
    env_vertices = std::move(other.env_vertices);
    env_faces = std::move(other.env_faces);
    tree = std::move(other.tree);
    oriented_side = std::move(other.oriented_side);

    const_cast<Tree_traits&>(tree.traits())
      .set_shared_data(Datum_map<GeomTraits>(bounding_boxes),
                       Point_map<GeomTraits>(bounding_boxes));

    return *this;
  }

  Polyhedral_envelope(Polyhedral_envelope<GeomTraits>&& other)
  {
    *this = std::move(other);
  }

  /**
   * Constructor with a triangulated surface mesh.
   * @tparam TriangleMesh a model of `FaceListGraph`
   * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * @param tmesh a triangle mesh
   * @param epsilon the distance of the Minkowski sum hull
   * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{vertex_point_map}
   *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
   *                    as key type and `%Point_3` as value type}
   *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
   *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
   *                     must be available in `TriangleMesh`.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{face_epsilon_map}
   *     \cgalParamDescription{a property map associating to each face of `tm` an epsilon value}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
   *                    as key type and `double` as value type}
   *     \cgalParamDefault{Use `epsilon` for all faces}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \note The triangle mesh gets copied internally, that is it can be modified after having passed as argument,
   *       while the queries are performed
   */
  template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
  Polyhedral_envelope(const TriangleMesh& tmesh,
                      double epsilon,
                      const NamedParameters& np = parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typedef boost::graph_traits<TriangleMesh> Graph_traits;
    typedef typename Graph_traits::face_descriptor face_descriptor;

    typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

    env_vertices.reserve(vertices(tmesh).size());
    env_faces.reserve(faces(tmesh).size());

    typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<int> >::const_type VIM;
    VIM vim  = get(CGAL::dynamic_vertex_property_t<int>(), tmesh);

    int id=0;
    for(typename Graph_traits::vertex_descriptor v : vertices(tmesh)){
      put(vim, v, id++);
      env_vertices.emplace_back(get(vpm, v));
    }

    GeomTraits gt;
    std::unordered_set<face_descriptor> deg_faces;
    for(face_descriptor f : faces(tmesh)){
      if(! Polygon_mesh_processing::is_degenerate_triangle_face(f, tmesh, parameters::geom_traits(gt).vertex_point_map(vpm))){
        typename Graph_traits::halfedge_descriptor h = halfedge(f, tmesh);
        int i = get(vim, source(h, tmesh));
        int j = get(vim, target(h, tmesh));
        int k = get(vim, target(next(h, tmesh), tmesh));

        Vector3i face = { i, j, k };
        env_faces.push_back(face);
      }
      else
        deg_faces.insert(f);
    }
    if (is_default_parameter<NamedParameters, internal_np::face_epsilon_map_t>::value)
      init(epsilon);
    else
    {
      std::vector<double> epsilon_values;
      epsilon_values.reserve(env_faces.size());

      typedef typename internal_np::Lookup_named_param_def<
        internal_np::face_epsilon_map_t,
        NamedParameters,
        Constant_property_map<face_descriptor, double>
      > ::type  Epsilon_map;

      Epsilon_map epsilon_map = choose_parameter(get_parameter(np, internal_np::face_epsilon_map),
                                                 Constant_property_map<face_descriptor, double>(epsilon));

      for(face_descriptor f : faces(tmesh))
        if(deg_faces.count(f)==0)
          epsilon_values.push_back( get(epsilon_map, f) );
      init(epsilon_values);
    }
  }

  /**
   * Constructor using a subset of faces of a triangulated surface mesh
   *
   * @tparam FaceRange a model of `ConstRange` with `ConstRange::const_iterator`
   *                   being a model of `InputIterator` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
   *                   as value type
   * @tparam TriangleMesh a model of `FaceListGraph`
   * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * @param tmesh a triangle mesh
   * @param face_range the subset of faces to be considered when computing the polyhedron envelope
   * @param epsilon the distance of the Minkowski sum hull
   * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{vertex_point_map}
   *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
   *                    as key type and `%Point_3` as value type}
   *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
   *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
   *                     must be available in `TriangleMesh`.}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{face_epsilon_map}
   *     \cgalParamDescription{a property map associating to each face of `tm` an epsilon value}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
   *                    as key type and `double` as value type}
   *     \cgalParamDefault{Use `epsilon` for all faces}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \note The triangle mesh gets copied internally, that is it can be modified after having passed as argument,
   *       while the queries are performed
   */
  template <typename FaceRange, typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
  Polyhedral_envelope(const FaceRange& face_range,
                      const TriangleMesh& tmesh,
                      double epsilon,
                      const NamedParameters& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                      , std::enable_if_t<!boost::has_range_const_iterator<TriangleMesh>::value>* = 0
#endif
  )
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    std::unordered_map<vertex_descriptor, int> local_vid;
    env_faces.reserve(face_range.size());
    env_vertices.reserve(3*face_range.size()/2); // only a guess

    GeomTraits gt;
    int curr_id=0;
    auto get_vid = [&local_vid,&curr_id, &vpm, this](vertex_descriptor v)
    {
      auto insert_res = local_vid.insert( std::make_pair(v, curr_id) );
      if (insert_res.second){
        ++curr_id;
        env_vertices.emplace_back(get(vpm, v));
      }
      return insert_res.first->second;
    };

    std::unordered_set<face_descriptor> deg_faces;
    for(face_descriptor f : face_range){
      if(! Polygon_mesh_processing::is_degenerate_triangle_face(f, tmesh, parameters::geom_traits(gt).vertex_point_map(vpm))){
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h = halfedge(f, tmesh);
        int i = get_vid(source(h, tmesh));
        int j = get_vid(target(h, tmesh));
        int k = get_vid(target(next(h, tmesh), tmesh));

        Vector3i face = { i, j, k };
        env_faces.push_back(face);
      }
      else
        deg_faces.insert(f);
    }

    if (is_default_parameter<NamedParameters, internal_np::face_epsilon_map_t>::value)
      init(epsilon);
    else
    {
      std::vector<double> epsilon_values;
      epsilon_values.reserve(env_faces.size());

      typedef typename internal_np::Lookup_named_param_def<
        internal_np::face_epsilon_map_t,
        NamedParameters,
        Constant_property_map<face_descriptor, double>
      > ::type  Epsilon_map;

      Epsilon_map epsilon_map = choose_parameter(get_parameter(np, internal_np::face_epsilon_map),
                                                 Constant_property_map<face_descriptor, double>(epsilon));

      for(face_descriptor f : face_range)
        if(deg_faces.count(f)==0)
          epsilon_values.push_back( get(epsilon_map, f) );
      init(epsilon_values);
    }
  }

  /**
    * Constructor with a triangle soup.
    *
    * @tparam PointRange a model of the concept `ConstRange` with `PointRange::const_iterator`
    *                    being a model of `InputIterator` with a point as value type
    * @tparam TriangleRange a model of the concept `ConstRange` with `TriangleRange::const_iterator`
    *                       being a model of `InputIterator` whose value type is model of
    *                       the concept `RandomAccessContainer` whose value_type is convertible to `std::size_t`.
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param points points of the soup of triangles
    * @param triangles each element in the range describes a triangle as a triple of indices of the points in `points`
    * @param epsilon the distance of the Minkowski sum hull
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{point_map}
    *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
    *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is `Point_3` and whose key
    *                    is the value type of `PointRange::const_iterator`}
    *     \cgalParamDefault{`CGAL::Identity_property_map`}
    *   \cgalParamNEnd
    *   \cgalParamNBegin{face_epsilon_map}
    *     \cgalParamDescription{a property map associating to each triangle an epsilon value}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with `std::size_t` as key type and `double` as value type}
    *     \cgalParamDefault{Use `epsilon` for all triangles}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    */
  template <typename PointRange, typename TriangleRange, typename NamedParameters = parameters::Default_named_parameters>
  Polyhedral_envelope(const PointRange& points,
                      const TriangleRange& triangles,
                      double epsilon,
                      const NamedParameters& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                      , std::enable_if_t<boost::has_range_const_iterator<TriangleRange>::value>* = 0
#endif
                      )
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    using parameters::is_default_parameter;

    typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type Point_map;
    Point_map pm = choose_parameter<Point_map>(get_parameter(np, internal_np::point_map));

    env_vertices.reserve(points.size());
    env_faces.reserve(triangles.size());

    typedef typename boost::range_value<PointRange>::type  Point_key;
    for (const Point_key& p : points)
      env_vertices.emplace_back(get(pm,p));

    typedef typename boost::range_value<TriangleRange>::type  Triangle;
    for (const Triangle& t : triangles)
    {
      Vector3i face = { int(t[0]), int(t[1]), int(t[2]) };
      env_faces.emplace_back(face);
    }

    if (is_default_parameter<NamedParameters, internal_np::face_epsilon_map_t>::value)
    {
      init(epsilon);
    }
    else
    {
      std::vector<double> epsilon_values;
      epsilon_values.reserve(env_faces.size());

      typedef typename internal_np::Lookup_named_param_def<
        internal_np::face_epsilon_map_t,
        NamedParameters,
        Constant_property_map<std::size_t, double>
      > ::type  Epsilon_map;

      Epsilon_map epsilon_map = choose_parameter(get_parameter(np, internal_np::face_epsilon_map),
                                                 Constant_property_map<std::size_t, double>(epsilon));

      for(std::size_t i=0; i<triangles.size(); ++i)
        epsilon_values.push_back( get(epsilon_map, i) );
      init(epsilon_values);
    }
  }

  /// @}

private:

  template <class Epsilons>
  void init(const Epsilons& epsilon_values)
  {
    halfspace_generation(env_vertices, env_faces, halfspace, bounding_boxes, epsilon_values);

    Datum_map<GeomTraits> datum_map(bounding_boxes);
    Point_map<GeomTraits> point_map(bounding_boxes);

    // constructs AABB tree
    tree.insert(boost::counting_iterator<unsigned int>(0),
                boost::counting_iterator<unsigned int>((unsigned int)bounding_boxes.size()),
                datum_map,
                point_map);
    tree.build();

#ifdef CGAL_ENVELOPE_DEBUG
    Surface_mesh<typename Exact_predicates_inexact_constructions_kernel::Point_3> sm;
    for(unsigned int i = 0; i < halfspace.size(); ++i){
      prism_to_mesh(i, sm);
    }
    std::ofstream("all_prisms.off") << std::setprecision(17) << sm;
#endif
  }


  struct INDEX
  {
    int Pi;
    std::vector<int> FACES;
  };


  // was marked todo???
  bool box_box_intersection(const Iso_cuboid_3 ic1, const Iso_cuboid_3 ic2) const
  {
    const Point_3& min1 = min_vertex(ic1);
    const Point_3& min2 = min_vertex(ic2);
    const Point_3& max1 = max_vertex(ic1);
    const Point_3& max2 = max_vertex(ic2);

    if ((compare_x(max1, min2) == SMALLER) ||(compare_y(max1, min2) == SMALLER) ||(compare_z(max1, min2) == SMALLER)){
      return 0;
    }
    if ((compare_x(max2, min1) == SMALLER) ||(compare_y(max2, min1) == SMALLER) ||(compare_z(max2, min1) == SMALLER)){
      return 0;
    }
    return 1;
  }


  bool
  point_out_prism(const ePoint_3 &point, const std::vector<unsigned int> &prismindex, unsigned int jump) const
  {
    Orientation ori;

    for (unsigned int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }

      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        const Plane& plane = halfspace[prismindex[i]][j];
        // As all points have coordinates with intervals with inf==sup the orientation test is faster
        // as it can exploit the static filters
        ori = orientation(plane.ep, plane.eq, plane.er, point);
        // ori = oriented_side(plane.eplane, point);
        if (ori != ON_NEGATIVE_SIDE){
          // if for a prism we are on the wrong side of one halfspace we are outside this prism
          // so no need to look at the other halfspaces
          break;
        }
        if (j == halfspace[prismindex[i]].size() - 1){
          // As we are in all halfspaces of one prism we are in the union of the prisms
          return false;
        }
      }
    }

    return true;
  }


  // \param jump is the index of the prism that shall be ignored
  // \param id is a return parameter for the prism with `point` inside
  bool
  point_out_prism_return_local_id(const Point_3 &point, const ePoint_3 &epoint, const std::vector<unsigned int> &prismindex, const unsigned int jump, int &id) const
  {
    Vector_3 bmin, bmax;

    Orientation ori;

    for (unsigned int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }
      if(bounding_boxes[prismindex[i]].has_on_unbounded_side(point)){
        continue;
      }
      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++){
        const Plane& plane = halfspace[prismindex[i]][j];
        ori = orientation(plane.ep, plane.eq, plane.er, epoint);
        if (ori != ON_NEGATIVE_SIDE){
          break;
        }

      }
      if (ori == ON_NEGATIVE_SIDE){
        id = i;
        return false;
      }

    }

    return true;
  }


  // \param cindex  the index of a prism
  // \param cid     a return parameter where the indices of the faces that intersect the segment `(source,target)`get inserted
  bool
  is_seg_cut_polyhedra(const int cindex,
                       const ePoint_3& source,
                       const ePoint_3& target,
                       const eLine_3& line,
                       std::vector<int>& cid) const
  {
    const Prism& prism = halfspace[cindex];
    cid.clear();
    std::array<bool,8> cut = { false, false,  false, false,  false, false,  false, false };

    std::array<Orientation,8> o1, o2;

    Oriented_side ori = ON_ORIENTED_BOUNDARY;
    int ct1 = 0, ct2 = 0;  //ori=0 to avoid the case that there is only one cut plane
    std::vector<int> cutp;
    cutp.reserve(8);

    for (unsigned int i = 0; i < prism.size(); i++){
      const Plane& plane = prism[i];
      // POSITIVE is outside the prism
      o1[i] = orientation(plane.ep, plane.eq, plane.er, source); // orientation exploits static filter as inf()==sup()
      o2[i] = orientation(plane.ep, plane.eq, plane.er, target);

      if (int(o1[i]) + int(o2[i]) >= 1)
        {
          return false;
        }

      if (o1[i] == ON_ORIENTED_BOUNDARY && o2[i] == ON_ORIENTED_BOUNDARY)
        {
          return false;
        }

      if (int(o1[i]) * int(o2[i]) == -1){
        cutp.emplace_back(i);
      }
      if (o1[i] == ON_POSITIVE_SIDE) ct1++;
      if (o2[i] == ON_POSITIVE_SIDE) ct2++; // if ct1 or ct2 >0, then NOT totally inside
    }

    if (cutp.size() == 0 && ct1 == 0 && ct2 == 0){
      // no intersected planes, and each point is either inside of poly,
      //or on one facet, since vertices are checked, then totally inside
      return true;
    }
    if (cutp.size() == 0) {
      return false;
    }

    /* The segment can have an intersection with several planes,
       but they may be outside the prism.
       So we have to test for the intersection points i if they are outside the prism.



                                    |
                                    |
                                    |      t
                                    |     /
       -------------*****************----i-------------
                    *               *   /
                    *               *  /
                    *    prism      * /
                    *               ./
                    *               i
                    *              /.
                    *             / *
                    *            /  *
                    *           s   *
                    *               *
       -------------*****************-----------------
                    |               |
                    |               |
                    |               |
    */


    for (unsigned int i = 0; i < cutp.size(); i++){
      const Plane& plane_i = prism[cutp[i]];

      std::optional<ePoint_3> op = intersection_point_for_polyhedral_envelope(line, plane_i.eplane);
      if(! op){
        std::cout <<  "there must be an intersection 2" << std::endl;
      }

      const ePoint_3& ip = *op;

      for(unsigned int j = 0; j < cutp.size(); j++) {
        if (i == j){
          continue;
        }
        const Plane& plane_j = prism[cutp[j]];
        ori = oriented_side(plane_j.eplane, ip);

        if(ori == ON_POSITIVE_SIDE){
          break;
        }
      }
      if (ori != ON_POSITIVE_SIDE) {
        cut[cutp[i]] = true;
      }
    }

    for (unsigned int i = 0; i < prism.size(); i++) {
      if (cut[i] == true){
        cid.emplace_back(i);
      }
    }
    if (cid.size() == 0){
      return false;// if no intersection points, and segment not totally inside, then not intersected
    }
    return true;
  }


  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(const ePoint_3 &ip,
                                                          const std::vector<unsigned int> &prismindex, const unsigned int &jump, int &id) const
  {
    Oriented_side ori = ON_POSITIVE_SIDE; // The compiler sees the
                                          // possibility that the
                                          // nested for loop body is
                                          // not executed and warns that
                                          // ori may not be initialized

    for (unsigned int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }

      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++){
        const Plane& plane = halfspace[prismindex[i]][j];
        ori = oriented_side(plane.eplane, ip);

        if (ori != ON_NEGATIVE_SIDE){
          break;
        }
      }
      if (ori == ON_NEGATIVE_SIDE){
        id = i;
        return IN_PRISM;
      }
    }

    return OUT_PRISM;
  }


  bool
  segment_out_of_envelope(const Point_3& source, const Point_3& target,
                          const std::vector<unsigned int>& prismindex) const
  {
    if (prismindex.size() == 0) {
      return true;
    }

    ePoint_3 esource(source.x(),source.y(),source.z());
    ePoint_3 etarget(target.x(),target.y(),target.z());

    unsigned int jump1 = -1;
    bool out, cut;
    int inter;

    int check_id, check_id1;

    // First check if the endpoints are outside the envelope
    out = point_out_prism_return_local_id(source, esource, prismindex, jump1, check_id);

    if (out) {
      return true;
    }
    out = point_out_prism_return_local_id(target, etarget, prismindex, jump1, check_id1);

    if (out) {
      return true;
    }

    // If both endpoints are in the same prism it is in the envelope
    if (check_id == check_id1){
      return false;
    }
    if (prismindex.size() == 1){
      return false;
    }
    eLine_3 line(esource,etarget);
    std::vector<unsigned int> queue, idlist;
    queue.emplace_back(check_id); //queue contains the id in prismindex
    idlist.emplace_back(prismindex[check_id]);

    std::vector<int> cidl;
    cidl.reserve(8);
    for (unsigned int i = 0; i < queue.size(); i++) {

      jump1 = prismindex[queue[i]];

      cut = is_seg_cut_polyhedra(jump1, esource, etarget, line, cidl);
      // cidl now contains the faces of the prism jump1
      if (cut&&cidl.size() == 0){
        return false;
      }
      if (!cut){
        continue;
      }

      for (unsigned int j = 0; j < cidl.size(); j++) {
        std::optional<ePoint_3> op = intersection_point_for_polyhedral_envelope(line,
                                                                                  halfspace[prismindex[queue[i]]][cidl[j]].eplane);
        const ePoint_3& ip = *op;
        inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id
          (ip, idlist, jump1, check_id);

        if (inter == 1){
          inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id
            (ip, prismindex, jump1, check_id);

          if (inter == 1) {
            return true; // outside envelope
          }
          if (inter == 0) {
            queue.emplace_back(check_id);
            idlist.emplace_back(prismindex[check_id]);
          }
        }
      }
    }
    return false; // fully inside the envelope
  }


  // This predicate checks if the intersection point of t/facet1/facet2   lies in the triangle
  int
  is_3_triangle_cut_float_fast(const ePoint_3& tp,
                               const ePoint_3& tq,
                               const ePoint_3& tr,
                               const ePoint_3& n,
                               const ePoint_3& ip
                               ) const
  {
    int o1 = int(orientation(n,tp,tq, ip));
    int o2 = int(orientation(n,tq,tr, ip));

    if (o1 * o2 == -1){

      return 0;
    }
    int o3 = int(orientation(n, tr, tp, ip));

    if (o1 * o3 == -1 || o2 * o3 == -1){

      return 0;
    }
    if (o1 * o2 * o3 == 0){
      return 2; // means we dont know
    }
    return 1;
  }


  bool
  is_3_triangle_cut(const ePoint_3& tp,
                    const ePoint_3& tq,
                    const ePoint_3& tr,
                    const ePoint_3& n,
                    const ePoint_3& ip) const
  {
    int o1 = int(orientation(n,tp,tq, ip));
    if (o1 == 0){
      return false;
    }

    int o2 = int(orientation(n, tq, tr, ip));

    if (o2 == 0 || o1 + o2 == 0){
      return false;
    }

    int o3 = int(orientation(n, tr, tp, ip));

    if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0){
      return false;
    }

    return true;
  }


  bool
  is_two_facets_neighboring(const unsigned int & pid, const unsigned int &i, const unsigned int &j)const
  {
    std::size_t facesize = halfspace[pid].size();
    if (i == j) return false;
    if( ((i == 0) && ( j==1))  || ((i == 1) && ( j==0)) ) return false;
    if (i == 0 && j != 1) return true;
    if (i == 1 && j != 0) return true;
    if (j == 0 && i != 1) return true;
    if (j == 1 && i != 0) return true;
    if (i - j == 1 || j - i == 1) return true;
    if (i == 2 && j == facesize - 1) return true;
    if (j == 2 && i == facesize - 1) return true;
    return false;
  }
  CGAL_DEPRECATED bool
  is_two_facets_neighbouring(const unsigned int & pid, const unsigned int &i, const unsigned int &j)const
  { return is_two_facets_neighboring(pid, i, j); }


  int
  is_triangle_cut_envelope_polyhedra(const int &cindex, //the triangle is not degenerated
                                     const Triangle& query,
                                     std::vector<unsigned int> &cid) const
  {
    const ePoint_3 &etri0 = query.etriangle[0];
    const ePoint_3 &etri1 = query.etriangle[1];
    const ePoint_3 &etri2 = query.etriangle[2];
    const ePoint_3& n = query.n;

    const Point_3 &tri0 = query.triangle[0];
    const Point_3 &tri1 = query.triangle[1];
    const Point_3 &tri2 = query.triangle[2];

    const Prism& prism = halfspace[cindex];
    cid.clear();
    cid.reserve(3);

    std::array<bool,8> cut = { false, false,  false, false,  false, false,  false, false };
    std::array<int,8> o1, o2, o3;
    std::vector<int>  cutp;
    cutp.reserve(8);

    Oriented_side ori = ON_ORIENTED_BOUNDARY;
    int ct1 = 0, ct2 = 0, ct3 = 0;

    for (unsigned int i = 0; i < prism.size(); i++)
      {
        const Plane& plane = prism[i];

        // As the orientation test operates on point with inf==sup of the coordinate intervals
        // we can profit from the static filters.
        // The oriented side test being made with interval arithmetic is more expensive
        o1[i] =  orientation(plane.ep, plane.eq, plane.er, etri0);
        o2[i] =  orientation(plane.ep, plane.eq, plane.er, etri1);
        o3[i] =  orientation(plane.ep, plane.eq, plane.er, etri2);

        if (o1[i] + o2[i] + o3[i] >= 2){ //1,1,0 case
          return 0;
        }
        if (o1[i] == 1) ct1++;
        if (o2[i] == 1) ct2++;
        if (o3[i] == 1) ct3++; // if ct1 or ct2 or ct3 >0, then NOT totally inside, otherwise, totally inside
        if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1){
            return 0;
          }
        if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0){
            return 0;
          }
        if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0){
            return 0;
          }

        if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1){
          cutp.emplace_back(i);
        }else if (o1[i] + o2[i] + o3[i] == -1 && o1[i] * o2[i] == 0) { //0,0,-1 case, we also want this face,really rare to happen
          cutp.emplace_back(i);
        }
      }
    if (cutp.size() == 0) {
      if (ct1 == 0 && ct2 == 0 && ct3 == 0) {
        return 2; // totally inside, or not any edge is on the facet
      }
    }

    if (cutp.size() == 0){
        return 0;
    }

    std::array<eLine_3*,2> seg;

    for (unsigned int i = 0; i < cutp.size(); i++)
      {
        int tmp = 0;
        if (o1[cutp[i]] * o2[cutp[i]] == -1|| o1[cutp[i]] + o2[cutp[i]] == -1) {
          seg[tmp] = const_cast<eLine_3*>(&(query.elines[2]));
          tmp++;
        }
        if (o1[cutp[i]] * o3[cutp[i]] == -1|| o1[cutp[i]] + o3[cutp[i]] == -1) {
          seg[tmp] = const_cast<eLine_3 *>(&(query.elines[1]));
          tmp++;
        }
        if (o2[cutp[i]] * o3[cutp[i]] == -1|| o2[cutp[i]] + o3[cutp[i]] == -1) {
          seg[tmp] = const_cast<eLine_3 *>(&(query.elines[0]));
          tmp++;
        }

        for (int k = 0; k < 2; k++){
          const Plane& plane_i = prism[cutp[i]];
          const eLine_3& eline = *(seg[k]);

          std::optional<ePoint_3> op = intersection_point_for_polyhedral_envelope(eline, plane_i.eplane);
          if(! op){
#ifdef CGAL_ENVELOPE_DEBUG
            std::cout <<  "there must be an intersection 6" << std::endl;
#endif
          }

          const ePoint_3& ip = *op;

          for (unsigned int j = 0; j < cutp.size(); j++){
            if (i == j){
                  continue;
            }
            const Plane& plane_j = prism[cutp[j]];
            ori = oriented_side(plane_j.eplane, ip);
            if (ori == ON_POSITIVE_SIDE){
              break;
            }
          }
          if (ori != 1){
            cut[cutp[i]] = true;
            break;
          }
        }

        ori = ON_ORIENTED_BOUNDARY;// initialize the orientation to avoid the j loop doesn't happen because cutp.size()==1
      }

    if (cutp.size() <= 2){
        for (unsigned int i = 0; i < prism.size(); i++){
          if (cut[i] == true){
              cid.emplace_back(i);
          }
        }
        return 1;
    }

    // triangle-facet-facet intersection

    const ePlane_3& tri_eplane = query.eplane;
    for (unsigned int i = 0; i < cutp.size(); i++)
      {
        for (unsigned int j = i + 1; j < cutp.size(); j++)// two facets and the triangle generate a point
          {
            if (cut[cutp[i]] == true && cut[cutp[j]] == true)
              continue;

            if (true /* USE_ADJACENT_INFORMATION*/ ) {
              bool neib = is_two_facets_neighboring(cindex, cutp[i], cutp[j]);
              if (neib == false) continue;
            }

            int inter = 0;

            static const int  edges[3][3] = { {2, 4, 6 }, {2, 3, 5}, {2, 4, 5} };

            int ob = prism.obtuse;
            if (ob == -1) ob = 0;

            if(((cutp[i] == 0)||(cutp[i] == 1)) && ( (cutp[j] == edges[ob][0])||(cutp[j] == edges[ob][1]) ||(cutp[j] == edges[ob][2]) )){ // ATTENTION Only correct together with CGAL_INITIAL
              int j0 = (cutp[j] == edges[ob][0])? 0 : (cutp[j]==edges[ob][1])? 1 : 2;
              int j1 = (cutp[j] == edges[ob][0])? 1 : (cutp[j]==edges[ob][1])? 2 : 0;
              const Vector3i & v3i = env_faces[cindex];
              const Point_3& pj0 = env_vertices[v3i[j0]];
              const Point_3& pj1 = env_vertices[v3i[j1]];

              if(pj0 == tri0){
                if((pj1 == tri1)||(pj1 == tri2)){
                  continue;
                }
              }
              if(pj0 == tri1){
                if((pj1 == tri2)||(pj1 == tri0)){
                  continue;
                }
              }
              if(pj0 == tri2){
                if((pj1 == tri0)||(pj1 == tri1)){
                  continue;
                }
              }


            }

            std::optional<ePoint_3>  ipp = intersection_point_for_polyhedral_envelope(tri_eplane, prism[cutp[i]].eplane, prism[cutp[j]].eplane);
            if(ipp){
                inter = is_3_triangle_cut_float_fast(etri0, etri1, etri2,
                                                     n,
                                                     *ipp);
            }
            // this was for a fast float check
            if (inter == 2)
              { //we dont know if point exist or if inside of triangle
                cut[cutp[i]] = true;
                cut[cutp[j]] = true;
                continue;
              }

            if (inter == 0){
              continue; // sure not inside
            }

            for (unsigned int k = 0; k < cutp.size(); k++){

              if (k == i || k == j){
                  continue;
              }

              ori = oriented_side(prism[cutp[k]].eplane, *ipp);

              if (ori == ON_POSITIVE_SIDE){
                  break;
              }
            }

            if (ori != 1) {
              cut[cutp[i]] = true;
              cut[cutp[j]] = true;
            }
          }
      }

    for (unsigned int i = 0; i < prism.size(); i++){
      if (cut[i] == true){
        cid.emplace_back(i);
      }
    }
    if (cid.size() == 0){
      return 0;// not cut and facets, and not totally inside, then not intersected
    }
    return 1;
  }


  int
  seg_cut_plane(const ePoint_3 &seg0, const ePoint_3 &seg1,
                const ePoint_3 &t0, const ePoint_3 &t1, const ePoint_3 &t2) const
  {
    int o1, o2;
    o1 = int(orientation(seg0, t0, t1, t2));
    o2 = int(orientation(seg1, t0, t1, t2));
    int op = o1 * o2;
    if (op >= 0)
      {
        return CUT_COPLANAR; //in fact, coplanar and not cut this plane
      }
    return CUT_FACE;
  }


  bool
  is_tpp_on_polyhedra(const ePoint_3& ip,
                      const int &prismid, const unsigned int &faceid)const
  {
    for (unsigned int i = 0; i < halfspace[prismid].size(); i++) {
      /*bool neib = is_two_facets_neighboring(prismid, i, faceid);// this works only when the polyhedron is convex and no two neighbor facets are coplanar
        if (neib == false) continue;*/
      if (i == faceid) continue;
      if(oriented_side(halfspace[prismid][i].eplane, ip) == ON_POSITIVE_SIDE){
        return false;
      }
    }
    return true;
  }



  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
                const ePoint_3& ip,
                const std::vector<unsigned int> &prismindex,
                const std::vector<std::vector<unsigned int>>& intersect_face, const unsigned int &jump, int &id) const
  {
    Oriented_side ori;
    unsigned int tot, fid;
    for (unsigned int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
          continue;
        }
      tot = 0; fid = 0;
      ori = ON_NEGATIVE_SIDE;
      const Prism& prism = halfspace[prismindex[i]];
      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {

        if (intersect_face[i][fid] == j)
          {
            if (fid + 1 < intersect_face[i].size()) fid++;
          }
        else continue;
        ori = oriented_side(prism[j].eplane, ip);

        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY)
          {
            break;
          }

        if (ori == ON_NEGATIVE_SIDE)
          {
            tot++;
          }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      fid = 0;
      ori = ON_NEGATIVE_SIDE;
      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j)
          {
            if (fid + 1 < intersect_face[i].size()) fid++;
            continue;
          }

        ori = oriented_side(prism[j].eplane, ip);
        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY)
          {
            break;
          }

        if (ori == ON_NEGATIVE_SIDE)
          {
            tot++;
          }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      if (tot == halfspace[prismindex[i]].size())
        {
          id = i;
          return IN_PRISM;
        }
    }
    return OUT_PRISM;
  }


  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
    const ePoint_3& ip,
    const std::vector<unsigned int> &prismindex,
    const std::vector<std::vector<unsigned int>>& intersect_face,
    const std::vector<int>& coverlist,
    const unsigned int &jump,
    int &id) const
  {
    Oriented_side ori;
    unsigned int tot, fid;
    for (unsigned int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }
      if (coverlist[i] == 1){
        continue;
      }
      tot = 0; fid = 0;
      ori = ON_NEGATIVE_SIDE;
      const Prism& prism = halfspace[prismindex[i]];

      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j)   {
          if (fid + 1 < intersect_face[i].size()) fid++;
        }else{
          continue;
        }
        ori = oriented_side(prism[j].eplane,ip);
        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
        continue;
      }
      fid = 0;
      ori = ON_NEGATIVE_SIDE;

      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j){
          if (fid + 1 < intersect_face[i].size()) fid++;
          continue;
        }

        ori = oriented_side(prism[j].eplane,ip);
        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
        continue;
      }
      if (tot == halfspace[prismindex[i]].size()) {
        id = i;
        return IN_PRISM;
      }
    }

    return OUT_PRISM;
  }


  int
  Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
    const ePoint_3& ip,
    const std::vector<unsigned int> &prismindex,
    const std::vector<std::vector<unsigned int>>&intersect_face,
    const unsigned int &jump1,
    const unsigned int &jump2,
    int &id) const
  {
    Oriented_side ori;
    unsigned int tot, fid;
    for (unsigned int i = 0; i < prismindex.size(); i++){

      if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
      if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump1])) continue;
      if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump2])) continue;

      tot = 0;
      fid = 0;
      ori = ON_NEGATIVE_SIDE;
      const Prism& prism = halfspace[prismindex[i]];
      for (unsigned int j = 0; j < prism.size(); j++) {
        if (intersect_face[i][fid] == j){
          if (fid + 1 < intersect_face[i].size()) fid++;
        }
        else continue;

        ori = oriented_side(prism[j].eplane, ip);

        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      fid = 0;
      ori = ON_NEGATIVE_SIDE;
      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j){
          if (fid + 1 < intersect_face[i].size()){
            fid++;
          }
          continue;
        }

        ori = oriented_side(prism[j].eplane, ip);

        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      if (tot == prism.size()){

        id = i;
        return IN_PRISM;
      }

    }

    return OUT_PRISM;
  }



  int
  Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
    const ePoint_3& ip,
    const std::vector<unsigned int>& prismindex,
    const std::vector<std::vector<unsigned int>*>& intersect_face,
    const std::vector<int>& coverlist,
    const unsigned int &jump1,
    const unsigned int &jump2,
    int &id) const
  {
    Oriented_side ori;
    unsigned int tot, fid;
    for (unsigned int i = 0; i < prismindex.size(); i++){


      if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
      if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump1])) continue;
      if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump2])) continue;
      if (coverlist[i] == 1) continue;
      tot = 0;
      fid = 0;
      ori = ON_NEGATIVE_SIDE;
      const Prism& prism = halfspace[prismindex[i]];
      for (unsigned int j = 0; j < prism.size(); j++) {
        if ((*intersect_face[i])[fid] == j){
          if (fid + 1 < intersect_face[i]->size()) fid++;
        }
        else continue;

        ori = oriented_side(prism[j].eplane, ip);

        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){
          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      fid = 0;
      ori = ON_NEGATIVE_SIDE;
      for (unsigned int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if ((*intersect_face[i])[fid] == j)
          {
            if (fid + 1 < intersect_face[i]->size()){
              fid++;
            }
            continue;
          }

        ori = oriented_side(prism[j].eplane, ip);

        if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY){

          break;
        }

        if (ori == ON_NEGATIVE_SIDE){
          tot++;
        }
      }
      if (ori == ON_POSITIVE_SIDE || ori == ON_ORIENTED_BOUNDARY) continue;
      if (tot == halfspace[prismindex[i]].size()){
        id = i;
        return IN_PRISM;
      }

    }

    return OUT_PRISM;
  }



  bool
  triangle_out_of_envelope(const Point_3 & t0,
                           const Point_3 & t1,
                           const Point_3 & t2,
                           const std::vector<unsigned int> &prismindex) const
  {

    if (prismindex.size() == 0) {
      return true;
    }

    Triangle query(t0, t1, t2);

    unsigned int jump1, jump2;
    static const std::array<std::array<int, 2>, 3> triseg = {
      {{{1, 2}}, {{2, 0}}, {{0, 1}}}
    };


    std::vector<unsigned int> filtered_intersection;
    filtered_intersection.reserve(prismindex.size() / 3);
    std::vector<std::vector<unsigned int>>intersect_face;
    intersect_face.reserve(prismindex.size() / 3);
    bool out, cut;

    int inter, tti; //triangle-triangle intersection

    jump1 = -1;

    int check_id;

    for (int i = 0; i < 3; i++) {
      out = point_out_prism_return_local_id(query.triangle[i], query.etriangle[i], prismindex, jump1, check_id);

      if (out) {
        return true;
      }
    }

    if (prismindex.size() == 1)
      return false;

    query.init_elines();

#ifdef DEGENERATION_FIX

    int degeneration = algorithms::is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

    if (degeneration == DEGENERATED_POINT){ //case 1 degenerate to a point
      return false;
    }

    if (degeneration == DEGENERATED_SEGMENT){
      std::vector<unsigned int > queue, idlist;
      queue.emplace_back(check_id);//queue contains the id in prismindex
      idlist.emplace_back(prismindex[check_id]);

      std::vector<int> cidl; cidl.reserve(8);
      for (unsigned int i = 0; i < queue.size(); i++) {

        jump1 = prismindex[queue[i]];
        int seg_inside = 0;
        for (int k = 0; k < 3; k++) {
          const eLine_3& eline = query.elines[k]; // (etriangle[triseg[k][0]], etriangle[triseg[k][1]]);
          cut = is_seg_cut_polyhedra(jump1, etriangle[triseg[k][0]], etriangle[triseg[k][1]], eline, cidl);
          if (cut&&cidl.size() == 0) {
            seg_inside++;
            if (seg_inside == 3) return false;// 3 segs are all totally inside of some polyhedrons
            continue;// means this seg is inside, check next seg
          }
          if (!cut) continue;

          for (unsigned int j = 0; j < cidl.size(); j++) {
            std::optional<ePoint_3> op = intersection_point_for_polyhedral_envelope(eline,
                                                                                      halfspace[prismindex[queue[i]]][cidl[j]].eplane);
            const ePoint_3& ip = *op;
            inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(ip, idlist, jump1, check_id);


            if (inter == 1){
              inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(ip, prismindex, jump1, check_id);

              if (inter == 1) {
                return true;
              }
              if (inter == 0) {
                queue.emplace_back(check_id);
                idlist.emplace_back(prismindex[check_id]);
              }
            }
          }
        }
      }

      return false;
    }
    //
#endif // DEGENERATION_FIX


    std::vector<unsigned int> cidl; cidl.reserve(8); // todo: std::array??
    for (unsigned int i = 0; i < prismindex.size(); i++) {
      tti = is_triangle_cut_envelope_polyhedra(prismindex[i],
                                               query,
                                               cidl);
      if (tti == 2) {
        return false; //totally inside of this polyhedron
      }
      else if (tti == 1 && cidl.size() > 0) {
        filtered_intersection.emplace_back(prismindex[i]);
        intersect_face.emplace_back(cidl);
      }
    }

    if (filtered_intersection.size() == 0) {
      return false;
    }

#ifdef CGAL_ENVELOPE_DEBUG
    for(std::size_t i = 0; i < filtered_intersection.size(); i++){
      prism_to_off(filtered_intersection[i], "filtered");
    }
#endif
    std::vector<unsigned int > queue, idlist;
    // coverlist shows if the element in filtered_intersection is one of the current covers
    std::vector<int> coverlist(filtered_intersection.size());

    queue.emplace_back(0);//queue contains the id in filtered_intersection
    idlist.emplace_back(filtered_intersection[queue[0]]);// idlist contains the id in prismid//it is fine maybe it is not really intersected
    coverlist[queue[0]] = 1 ;//when filtered_intersection[i] is already in the cover list, coverlist[i]=true

    std::vector<unsigned int> neighbors;//local id
    std::vector<unsigned int > list;
    std::vector<std::vector<unsigned int>*> neighbor_facets;
    std::vector<std::vector<unsigned int>>  idlistorder;
    std::vector<int> neighbor_cover;
    idlistorder.emplace_back(intersect_face[queue[0]]);

    for (unsigned int i = 0; i < queue.size(); i++) {

      jump1 = filtered_intersection[queue[i]];

      for (int k = 0; k < 3; k++) {
        const eLine_3& eline = query.elines[k];

        for (unsigned int j = 0; j < intersect_face[queue[i]].size(); j++) {
          tti = seg_cut_plane(query.etriangle[triseg[k][0]],
                              query.etriangle[triseg[k][1]],
                              halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].ep,
                              halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].eq,
                              halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].er);

          if (tti != CUT_FACE){
            continue;
          }
          // now we know that there exists an intersection point

          std::optional<ePoint_3> op = intersection_point_for_polyhedral_envelope(eline,
                                                                                    halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].eplane);
          const ePoint_3& ip = *op;

          inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(ip, idlist, idlistorder, jump1, check_id);



          if (inter == 1)
            {

              inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(ip, filtered_intersection, intersect_face, coverlist, jump1, check_id);


              CGAL_assertion(inter != 2);// the point must exist because it is a seg-halfplane intersection
              if (inter == 1) {

                return true;
              }
              if (inter == 0) {
                idlistorder.emplace_back(intersect_face[check_id]);
                queue.emplace_back(check_id);
                idlist.emplace_back(filtered_intersection[check_id]);
                coverlist[check_id] = 1;
              }
            }
        }
      }
    }

    //tpi part

    //tree

    Tree localtree;

    std::vector<Iso_cuboid_3> local_bounding_boxes;
    local_bounding_boxes.resize(filtered_intersection.size());

    for (unsigned int i = 0; i < filtered_intersection.size(); i++) {
      local_bounding_boxes[i] = bounding_boxes[filtered_intersection[i]];
    }

    Datum_map<GeomTraits> datum_map(local_bounding_boxes);
    Point_map<GeomTraits> point_map(local_bounding_boxes);

    // constructs AABB tree
    localtree.insert(boost::counting_iterator<unsigned int>(0),
                     boost::counting_iterator<unsigned int>((unsigned int)local_bounding_boxes.size()),
                     datum_map,
                     point_map);
    localtree.build();

    //tree end

    const ePlane_3& etriangle_eplane = query.eplane;

    for (unsigned int i = 1; i < queue.size(); i++){
      jump1 = filtered_intersection[queue[i]];

      localtree.all_intersected_primitives(bounding_boxes[jump1], std::back_inserter(list));

      neighbors.resize(list.size());
      neighbor_facets.resize(list.size());
      neighbor_cover.resize(list.size());
      for (unsigned int j = 0; j < list.size(); j++) {
        neighbors[j] = filtered_intersection[list[j]];
        neighbor_facets[j] = &(intersect_face[list[j]]);
        if (coverlist[list[j]] == 1) neighbor_cover[j] = 1;
        else neighbor_cover[j] = 0;
      }

      for (unsigned int j = 0; j < i; j++) {
        jump2 = filtered_intersection[queue[j]];

        if (! box_box_intersection(bounding_boxes[jump1], bounding_boxes[jump2])){
          continue;
        }
        for (unsigned int k = 0; k < intersect_face[queue[i]].size(); k++) {
          for (unsigned int h = 0; h < intersect_face[queue[j]].size(); h++) {

            // We moved the intersection here
            // In case there is no intersection point we continue
            std::optional<ePoint_3>
              op = intersection_point_for_polyhedral_envelope(etriangle_eplane,
                                                              halfspace[jump1][intersect_face[queue[i]][k]].eplane,
                                                              halfspace[jump2][intersect_face[queue[j]][h]].eplane);
            if(! op){
              continue;
            }

            const ePoint_3& ip = *op;

            cut = is_3_triangle_cut(query.etriangle[0], query.etriangle[1],
                                    query.etriangle[2], query.n, ip);

            if (!cut){
              continue;
            }

            cut = is_tpp_on_polyhedra(ip, jump1, intersect_face[queue[i]][k]);

            if (!cut){
              continue;
            }

            cut = is_tpp_on_polyhedra(ip, jump2, intersect_face[queue[j]][h]);

            if (!cut){
              continue;
            }

            inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(ip, idlist, idlistorder, jump1, jump2, check_id);
            if (inter == 1) {


              inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(ip, neighbors, neighbor_facets, neighbor_cover, jump1, jump2, check_id);


              if (inter == 1) {
                return true;
              }
              if (inter == 0) {
                idlistorder.emplace_back(intersect_face[list[check_id]]);
                queue.emplace_back(list[check_id]);
                idlist.emplace_back(filtered_intersection[list[check_id]]);
                coverlist[list[check_id]] = 1;
              }
            }
          }
        }
      }
    }

    return false;
  }


  Plane
  get_corner_plane(const Point_3& p0,
                   const Point_3& midp,
                   const Vector_3 &normal,
                   const double distance,
                   const bool /* robust */) const
  {
    Point_3 plane0, plane1, plane2;
    double distance_small = distance * 1;// to be conservative to reduce numerical error, can set the Scalar as 0.999
    Vector_3 direction = normalize(p0 - midp);
    plane0 = p0 + direction * distance_small;
    plane1 = plane0 + normal;
    Vector_3 axis =
      //(robust) ? robust_cross_product_direction(midp, p0, ORIGIN, normal) :
      normalize(cross_product(direction, normal));
    plane2 = plane0 + axis;

    return Plane(plane0, plane1,plane2);
  }

  double get_epsilon(double epsilon, std::size_t)
  {
    return epsilon;
  }

  double get_epsilon(const std::vector<double>& epsilon_values, std::size_t i)
  {
    return epsilon_values[i];
  }

  // build prisms for a list of triangles. each prism is represented by 7-8 planes, which are represented by 3 points
  template <class Epsilons>
  void
  halfspace_generation(const std::vector<Point_3> &ver, const std::vector<Vector3i> &faces,
                       std::vector<Prism>& halfspace,
                       std::vector<Iso_cuboid_3>& bounding_boxes, const Epsilons& epsilon_values)
  {
    Vector_3 AB, AC, BC, normal;
    Plane plane;
    std::array<Vector_3, 8> box;

#if 0
    static const std::array<Vector_3, 8> boxorder = {
      {
        {1, 1, 1},
        {-1, 1, 1},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, -1},
        {-1, 1, -1},
        {-1, -1, -1},
        {1, -1, -1},
      } };

  static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };
#endif
    bool use_accurate_cross = false;

    halfspace.resize(faces.size());
    bounding_boxes.resize(faces.size());
    for (unsigned int i = 0; i < faces.size(); ++i)
      {
        const double epsilon = get_epsilon(epsilon_values,i);
        double tolerance = epsilon / std::sqrt(3);// the envelope thickness, to be conservative
        double bbox_tolerance = epsilon *(1 + 1e-6);

        Bbox bb = ver[faces[i][0]].bbox () + ver[faces[i][1]].bbox() + ver[faces[i][2]].bbox();
        // todo: Add a grow() function to Bbox
        bounding_boxes[i] = Iso_cuboid_3(Point_3(bb.xmin()-bbox_tolerance, bb.ymin()-bbox_tolerance, bb.zmin()-bbox_tolerance),
                                         Point_3(bb.xmax()+bbox_tolerance, bb.ymax()+bbox_tolerance, bb.zmax()+bbox_tolerance));

        AB = ver[faces[i][1]] - ver[faces[i][0]];
        AC = ver[faces[i][2]] - ver[faces[i][0]];
        BC = ver[faces[i][2]] - ver[faces[i][1]];

#if 0
        int de = algorithms::is_triangle_degenerated(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);

        if (de == DEGENERATED_POINT)
          {
            for (int j = 0; j < 8; j++)
              {
                box[j] = ver[faces[i][0]] + boxorder[j] * tolerance;
              }
            halfspace[i].resize(6);
            for (int j = 0; j < 6; j++) {
              halfspace[i][j][0] = box[c_face[j][0]];
              halfspace[i][j][1] = box[c_face[j][1]];
              halfspace[i][j][2] = box[c_face[j][2]];
            }


            continue;
          }
        if (de == DEGENERATED_SEGMENT)
          {
            //logger().debug("Envelope Triangle Degeneration- Segment");
            Scalar length1 = AB.dot(AB), length2 = AC.dot(AC), length3 = BC.dot(BC);
            if (length1 >= length2 && length1 >= length3)
              {
                algorithms::seg_cube(ver[faces[i][0]], ver[faces[i][1]], tolerance, box);

              }
            if (length2 >= length1 && length2 >= length3)
              {
                algorithms::seg_cube(ver[faces[i][0]], ver[faces[i][2]], tolerance, box);

              }
            if (length3 >= length1 && length3 >= length2)
              {
                algorithms::seg_cube(ver[faces[i][1]], ver[faces[i][2]], tolerance, box);
              }
            halfspace[i].resize(6);
            for (int j = 0; j < 6; j++) {
              halfspace[i][j][0] = box[c_face[j][0]];
              halfspace[i][j][1] = box[c_face[j][1]];
              halfspace[i][j][2] = box[c_face[j][2]];
            }


            continue;
          }
        if (de == NERLY_DEGENERATED)
          {
            //logger().debug("Envelope Triangle Degeneration- Nearly");
            use_accurate_cross = true;

            normal = algorithms::accurate_normal_vector(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);

          }
        else
          {
            normal = normalize(cross_product(AB, AC));
          }
#endif
        normal = normalize(cross_product(AB, AC)); // remove as soon as #if 1 above

        halfspace[i].reserve(8);
        Vector_3 normaldist = normal * tolerance;
        Vector_3 edgedire, edgenormaldist;
        plane = Plane(ver[faces[i][0]] + normaldist,
                      ver[faces[i][1]] + normaldist,
                      ver[faces[i][2]] + normaldist);
        halfspace[i].emplace_back(plane);// number 0

        plane = Plane(ver[faces[i][0]] - normaldist,
                      ver[faces[i][2]] - normaldist,
                      ver[faces[i][1]] - normaldist);// order: 0, 2, 1
        halfspace[i].emplace_back(plane);// number 1

        int obtuse = obtuse_angle(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);
        halfspace[i].obtuse = obtuse;

        edgedire = normalize(AB);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(ORIGIN, edgedire, ORIGIN, normal)*tolerance;
        // else
        edgenormaldist = normalize(cross_product(edgedire,normal))*tolerance;
        plane = Plane(ver[faces[i][0]] + edgenormaldist,
                      ver[faces[i][1]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 2

        if (obtuse != 1) {
          plane = get_corner_plane(ver[faces[i][1]], midpoint(ver[faces[i][0]], ver[faces[i][2]]) , normal,
                                   tolerance, use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 3;
        }

        edgedire = normalize(BC);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(ORIGIN, edgedire, ORIGIN, normal)*tolerance;
        // else
        edgenormaldist = normalize(cross_product(edgedire, normal))*tolerance;

        plane = Plane(ver[faces[i][1]] + edgenormaldist,
                      ver[faces[i][2]] + edgenormaldist,
                      ver[faces[i][1]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 4

        if (obtuse != 2) {
          plane = get_corner_plane(ver[faces[i][2]], midpoint(ver[faces[i][0]], ver[faces[i][1]]), normal,
                                   tolerance,use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 5;
        }

        edgedire = -normalize(AC);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(ORIGIN, edgedire, ORIGIN , normal)*tolerance;
        // else
        edgenormaldist = normalize(cross_product(edgedire, normal))*tolerance;

        plane = Plane(ver[faces[i][2]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 6

        if (obtuse != 0) {
          plane = get_corner_plane(ver[faces[i][0]], midpoint(ver[faces[i][1]], ver[faces[i][2]]) , normal,
                                   tolerance,use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 7;
        }

#ifdef CGAL_ENVELOPE_DEBUG
        std::cout << "face "<< i << std::endl;
        for(unsigned int j = 0; j < halfspace[i].size(); j++){
          const Plane& p =  halfspace[i][j];
          std::cout << p.ep << " | "  << p.eq << " | "  << p.er << std::endl;
          ePoint_3 pv(ver[faces[i][0]].x(), ver[faces[i][0]].y(),ver[faces[i][0]].z());
          Orientation ori = orientation(p.ep, p.eq, p.er, pv);
          CGAL_assertion(ori == NEGATIVE);
        }
#endif

      }
  }

#ifdef CGAL_ENVELOPE_DEBUG
  template <class Mesh>
  void prism_to_mesh(unsigned int i, Mesh& sm) const
  {
    std::vector<ePlane_3> eplanes;
    for(unsigned int j = 0; j < halfspace[i].size(); j++){
      eplanes.push_back(halfspace[i][j].eplane);
    }
    ePoint_3 origin(env_vertices[env_faces[i][0]].x(), env_vertices[env_faces[i][0]].y(),env_vertices[env_faces[i][0]].z());
    Surface_mesh<ePoint_3> esm;
    halfspace_intersection_3(eplanes.begin(),eplanes.end(),esm , std::make_optional(origin));

    copy_face_graph(esm,sm);
  }

  void prism_to_off(unsigned int i, std::string fname) const
  {
    Surface_mesh<typename Exact_predicates_inexact_constructions_kernel::Point_3> sm;
    prism_to_mesh(i, sm);
    fname += "_";
    fname += std::to_string(i);
    fname += ".off";
    std::ofstream out(fname.c_str());
    out << sm << std::endl << std::endl;
  }
#endif

public:

  /// \name Query Operators
  /// @{

  /*!
   * returns `true`, iff the query point is inside the polyhedral envelope.
   */

  bool
  operator()(const Point_3& query) const
  {
    std::vector<unsigned int> prismindex;
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));
    if(prismindex.empty()){
      return false;
    }
    ePoint_3 equery(query.x(),query.y(),query.z());
    if(point_out_prism(equery, prismindex, -1)){
      return false;
    }
    return true;
  }


  /*!
   * returns `true`, iff the query segment defined by the points `source` and `target` is inside the polyhedral envelope.
   */
  bool
  operator()(const Point_3& source, const Point_3& target) const
  {
    if(source == target){
      return (*this)(source);
    }
    std::vector<unsigned int> prismindex;
    Segment_3 query(source,target);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));

    if(segment_out_of_envelope(source, target, prismindex)){
      return false;
    }
    return true;
  }


  /*!
   * returns `true`, iff the query triangle formed by the points `t0`, `t1`, and `t2` is inside the polyhedral envelope.
   */
  bool
  operator()(const Point_3& t0, const Point_3& t1, const Point_3& t2) const
  {
    if (collinear(t0,t1,t2)){
      Point_3 p0(t0), p1(t1), p2(t2);
      if(lexicographically_xyz_smaller(p1,p0)){
        std::swap(p1,p0);
      }
      if(lexicographically_xyz_smaller(p2,p1)){
        std::swap(p2,p1);
      }
      if(lexicographically_xyz_smaller(p1,p0)){
        std::swap(p1,p0);
      }
      return (*this)(p0,p2);
    }
    std::vector<unsigned int> prismindex;
    Triangle_3 query(t0, t1, t2);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));

    // only needed when we want to compare runs with FastEnvelope
    // std::sort(prismindex.begin(), prismindex.end());

    if(triangle_out_of_envelope(t0, t1, t2, prismindex)){
      return false;
    }
    return true;
  }

  /**
   * returns `true`, iff all the triangles of `tmesh` are inside the polyhedral envelope.
   * @tparam TriangleMesh a model of `FaceListGraph`
   * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * @param tmesh a triangle mesh
   * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{vertex_point_map}
   *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
   *                    as key type and `%Point_3` as value type}
   *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
   *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
   *                     must be available in `TriangleMesh`.}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   * \todo Add ConcurrencyTag as template parameter + use TBB parallel for
   * \todo Find a way to test the containment of the vertices first and then
   *       the triangles. It requires to have a map vertex->prism id so that
   *       we can test if the 3 vertices of a face are in the same face + have
   *       the initial list of prisms.
   * \todo apply that to the soup versions
   */
  template <typename TriangleMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
  bool
  operator()(const TriangleMesh& tmesh,
             const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
             , std::enable_if_t<!boost::has_range_const_iterator<TriangleMesh>::value>* = 0
#endif
    ) const
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

    for (typename boost::graph_traits<TriangleMesh>::face_descriptor f : faces(tmesh))
    {
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h =
        halfedge(f, tmesh);
      if (! this->operator()(get(vpm, source(h, tmesh)),
                             get(vpm, target(h, tmesh)),
                             get(vpm, target(next(h, tmesh), tmesh))) )
      {
        return false;
      }
    }

    return true;
  }

  /**
    * returns `true`, iff all the triangles in `triangles` are inside the polyhedral envelope.
    *
    * @tparam PointRange a model of the concept `ConstRange` with `PointRange::const_iterator`
    *                    being a model of `InputIterator` with a point as value type
    * @tparam TriangleRange a model of the concept `ConstRange` with `TriangleRange::const_iterator`
    *                       being a model of `InputIterator` whose value type is model of
    *                       the concept `RandomAccessContainer` whose value_type is convertible to `std::size_t`.
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param points points of the soup of triangles
    * @param triangles each element in the range describes a triangle as a triple of indices of the points in `points`
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{point_map}
    *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
    *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is `Point_3` and whose key
    *                    is the value type of `PointRange::const_iterator`}
    *     \cgalParamDefault{`CGAL::Identity_property_map`}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    */
  template <typename PointRange, typename TriangleRange, typename NamedParameters = parameters::Default_named_parameters>
  bool operator()(const PointRange& points,
                  const TriangleRange& triangles,
                  const NamedParameters& np = parameters::default_values()) const
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef typename std::iterator_traits<typename TriangleRange::const_iterator>::value_type Triangle;

    typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type Point_map;
    Point_map pm = choose_parameter<Point_map>(get_parameter(np, internal_np::point_map));

    std::array<Point_3, 3> pts;
    for (const Triangle& f : triangles)
    {
      typename Triangle::const_iterator t_it = f.begin();
      pts[0]=get(pm, points[*t_it]);
      pts[1]=get(pm, points[*(++t_it)]);
      pts[2]=get(pm, points[*(++t_it)]);
      if (! this->operator()(pts[0], pts[1], pts[2]) )
      {
        return false;
      }
    }

    return true;
  }

  /**
   * returns `true`, iff all the triangles in `triangle_range` are inside the polyhedral envelope.
   * @tparam TriangleRange a model of `ConstRange` with `ConstRange::const_iterator`
   *                       being a model of `InputIterator` with a value type being itself a model of
   *                       `ConstRange` with `ConstRange::const_iterator` being a model of `InputIterator`
   *                       with `Point_3` as value type.
   *
   * @param triangle_range a range of triangles
   */
  template <typename TriangleRange>
  bool
  operator()(const TriangleRange& triangle_range
#ifndef DOXYGEN_RUNNING
             , std::enable_if_t<boost::has_range_const_iterator<TriangleRange>::value>* = 0
#endif
    ) const
  {
    std::vector<Point_3> triangle;
    triangle.reserve(3);
    for (const auto& t : triangle_range)
    {
      triangle.clear();
      triangle.assign(t.begin(), t.end());
      CGAL_assertion(triangle.size()==3);
      if (! this->operator()(t[0], t[1], t[2]) )
      {
        return false;
      }
    }

    return true;
  }

  /// @}

  /// returns `true` if the polyhedral envelope is empty and `false` otherwise.
  bool is_empty() const
  {
    return env_faces.empty();
  }

}; // class Polyhedral_envelope

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_POLYHEDRAL_ENVELOPE_H
