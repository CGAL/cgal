// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

#ifndef P2T2_INTERFACE_TEST_H
#define P2T2_INTERFACE_TEST_H

#include <CGAL/point_generators_2.h>
#include <CGAL/use.h>

template <class T>
void test_constructor()
{
  typedef typename T::Iso_rectangle Iso_rectangle;
  typedef typename T::Geom_traits   Geom_traits;
  typedef typename T::Point         Point;

  T t;
  T t2(t);
  CGAL_assertion(t == t2);
  T t3 = t2;
  CGAL_assertion(t == t3);
  T t4(Iso_rectangle(0, 0, 2, 2));
  T t5(Iso_rectangle(0, 0, 2, 2), Geom_traits());
  t5.clear();

  t.insert(Point(0.5, 0.5));
  CGAL_assertion(t != t2);
  CGAL_assertion(t != t3);

  T t6(t);
  CGAL_assertion(t == t6);
  T t7 = t6;
  CGAL_assertion(t == t7);

  t.clear();
  CGAL_assertion(t != t6);
  CGAL_assertion(t != t7);

  t.swap(t7);
  CGAL_assertion(t7.empty());
  CGAL_assertion(!t.empty());
  CGAL_assertion(t != t7);
}

template <class T>
void test_global_access()
{
  T t;
  const T &t_const = t;

  typename T::Geom_traits gt;
  gt = t.geom_traits();

  t_const.tds();
  t.tds();
  t_const.domain();
  t_const.number_of_sheets();
  t_const.dimension();

  size_t number_of_vertices = t_const.number_of_vertices();
  CGAL_USE(number_of_vertices);
  CGAL_assertion(number_of_vertices == t.number_of_vertices());
  size_t number_of_faces = t_const.number_of_faces();
  CGAL_USE(number_of_faces);
  CGAL_assertion(number_of_faces == t.number_of_faces());
  size_t number_of_stored_vertices = t_const.number_of_stored_vertices();
  CGAL_USE(number_of_stored_vertices);
  CGAL_assertion(number_of_stored_vertices == t.number_of_stored_vertices());
  size_t number_of_stored_faces = t_const.number_of_stored_faces();
  CGAL_USE(number_of_stored_faces);
  CGAL_assertion(number_of_stored_faces == t.number_of_stored_faces());

  size_t number_of_edges = t_const.number_of_edges();
  CGAL_USE(number_of_edges);
  CGAL_assertion(number_of_edges == t.number_of_edges());
  size_t number_of_stored_edges = t_const.number_of_stored_edges();
  CGAL_USE(number_of_stored_edges);
  CGAL_assertion(number_of_stored_edges == t.number_of_stored_edges());

  bool ext1 = t_const.is_extensible_triangulation_in_1_sheet_h1();
  CGAL_USE(ext1);
  CGAL_assertion(ext1 == t.is_extensible_triangulation_in_1_sheet_h1());
  bool ext2 = t_const.is_extensible_triangulation_in_1_sheet_h2();
  CGAL_USE(ext2);
  CGAL_assertion(ext2 == t.is_extensible_triangulation_in_1_sheet_h2());
  bool is_triang1 = t_const.is_triangulation_in_1_sheet();
  CGAL_USE(is_triang1);
  CGAL_assertion(is_triang1 == t.is_triangulation_in_1_sheet());
  t.convert_to_1_sheeted_covering();
  t.convert_to_9_sheeted_covering();
}

template <class T>
void test_geometric_access()
{
  typedef typename T::Point             Point;
  typedef typename T::Segment           Segment;
  typedef typename T::Triangle          Triangle;
  typedef typename T::Periodic_point    Periodic_point;
  typedef typename T::Periodic_segment  Periodic_segment;
  typedef typename T::Periodic_triangle Periodic_triangle;
  typedef typename T::Vertex_iterator   Vertex_iterator;
  typedef typename T::Face_iterator     Face_iterator;
  typedef typename T::Face_handle       Face_handle;

  T t;
  const T &t_const = t;

  t.insert(Point(0.5, 0.5), Face_handle());
  t.insert(Point(0.7, 0.5));
  t.insert(Point(0.7, 0.7));

  Vertex_iterator vit = t.vertices_begin();
  Face_iterator fit = t.faces_begin();

  Periodic_point pp0 = t_const.periodic_point(vit);
  Periodic_point pp1 = t_const.periodic_point(fit, 0);
  CGAL_USE(pp1);
  Periodic_segment ps0 = t_const.periodic_segment(fit, 0);
  CGAL_USE(ps0);
  Periodic_segment ps1 = t_const.periodic_segment(*t.edges_begin());
  CGAL_USE(ps1);
  Periodic_triangle pt0 = t_const.periodic_triangle(fit);
  CGAL_USE(pt0);

  Point p0 = t_const.point(pp0);
  CGAL_USE(p0);
  Segment s0 = t_const.segment(ps0);
  CGAL_USE(s0);
  Triangle t0 = t_const.triangle(pt0);
  CGAL_USE(t0);
}

template <class T>
void test_predicates()
{
  typedef typename T::Vertex_handle   Vertex_handle;
  typedef typename T::Face_handle     Face_handle;
  typedef typename T::Point             Point;

  T t;

  Vertex_handle vh0 = t.insert(Point(0.5, 0.5));
  Vertex_handle vh1 = t.insert(Point(0.7, 0.5));
  Vertex_handle vh2 = t.insert(Point(0.7, 0.7));

  t.is_edge(vh0, vh1);
  Face_handle fh;
  int i;
  t.is_edge(vh0, vh1, fh, i);

  t.is_face(vh0, vh1, vh2);
  t.is_face(vh0, vh1, vh2, fh);
}

template <class T>
void test_queries()
{
  typedef typename T::Vertex_handle   Vertex_handle;
  typedef typename T::Face_handle     Face_handle;
  typedef typename T::Point             Point;

  T t;
  const T &t_const = t;
  CGAL_USE(t_const);

  Point p0(0.5, 0.5);
  Vertex_handle vh0 = t.insert(p0);
  CGAL_USE(vh0);
  Vertex_handle vh1 = t.insert(Point(0.7, 0.5));
  CGAL_USE(vh1);
  Vertex_handle vh2 = t.insert(Point(0.7, 0.7));
  CGAL_USE(vh2);

  Face_handle fh = t_const.locate(p0);
  fh = t_const.locate(Point(0.5, 0.5), fh);

  typename T::Locate_type lt;
  int li;
  fh = t_const.locate(p0, lt, li);
  fh = t_const.locate(p0, lt, li, fh);

  t_const.oriented_side(fh, p0);
}

template <class T>
void test_iterators()
{
  typedef typename T::Vertex_handle     Vertex_handle;
  typedef typename T::Point             Point;

  T t;
  const T &t_const = t;

  Point p0(0.5, 0.5);
  Vertex_handle vh0 = t.insert(p0);
  CGAL_USE(vh0);
  Vertex_handle vh1 = t.insert(Point(0.7, 0.5));
  CGAL_USE(vh1);
  Vertex_handle vh2 = t.insert(Point(0.7, 0.7));
  CGAL_USE(vh2);

  // vertices
  size_t size = 0;
  for (typename T::Vertex_iterator vit = t_const.vertices_begin();
       vit != t_const.vertices_end(); ++vit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_vertices());
  size = 0;
  for (typename T::Unique_vertex_iterator uvit = t_const.unique_vertices_begin();
       uvit != t_const.unique_vertices_end(); ++uvit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_vertices());
  size = 0;
  for (typename T::Vertex_iterator vit = t_const.all_vertices_begin();
       vit != t_const.all_vertices_end(); ++vit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_vertices());

  // edges
  size = 0;
  for (typename T::Edge_iterator eit = t_const.edges_begin();
       eit != t_const.edges_end(); ++eit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_edges());
  size = 0;
  for (typename T::Edge_iterator eit = t_const.all_edges_begin();
       eit != t_const.all_edges_end(); ++eit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_edges());

  // faces
  size = 0;
  for (typename T::Face_iterator fit = t_const.faces_begin();
       fit != t_const.faces_end(); ++fit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_faces());
  size = 0;
  for (typename T::All_faces_iterator fit = t_const.all_faces_begin();
       fit != t_const.all_faces_end(); ++fit)
    {
      ++size;
    }
  CGAL_assertion(size == t_const.number_of_stored_faces());

  /// Geometric iterators
  for (typename T::Periodic_point_iterator ppit = t_const.periodic_points_begin();
       ppit != t_const.periodic_points_end(); ++ppit)
    {
    }
  for (typename T::Periodic_point_iterator ppit =
         t_const.periodic_points_begin(T::STORED);
       ppit != t_const.periodic_points_end(T::STORED); ++ppit)
    {
    }
  for (typename T::Periodic_point_iterator ppit =
         t_const.periodic_points_begin(T::UNIQUE);
       ppit != t_const.periodic_points_end(T::UNIQUE); ++ppit)
    {
    }
  for (typename T::Periodic_point_iterator ppit =
         t_const.periodic_points_begin(T::STORED_COVER_DOMAIN);
       ppit != t_const.periodic_points_end(T::STORED_COVER_DOMAIN); ++ppit)
    {
    }
  for (typename T::Periodic_point_iterator ppit =
         t_const.periodic_points_begin(T::UNIQUE_COVER_DOMAIN);
       ppit != t_const.periodic_points_end(T::UNIQUE_COVER_DOMAIN); ++ppit)
    {
    }

  for (typename T::Periodic_segment_iterator psit = t_const.periodic_segments_begin();
       psit != t_const.periodic_segments_end(); ++psit)
    {
    }
  for (typename T::Periodic_segment_iterator psit =
         t_const.periodic_segments_begin(T::STORED);
       psit != t_const.periodic_segments_end(T::STORED); ++psit)
    {
    }
  for (typename T::Periodic_segment_iterator psit =
         t_const.periodic_segments_begin(T::UNIQUE);
       psit != t_const.periodic_segments_end(T::UNIQUE); ++psit)
    {
    }
  for (typename T::Periodic_segment_iterator psit =
         t_const.periodic_segments_begin(T::STORED_COVER_DOMAIN);
       psit != t_const.periodic_segments_end(T::STORED_COVER_DOMAIN); ++psit)
    {
    }
  for (typename T::Periodic_segment_iterator psit =
         t_const.periodic_segments_begin(T::UNIQUE_COVER_DOMAIN);
       psit != t_const.periodic_segments_end(T::UNIQUE_COVER_DOMAIN); ++psit)
    {
    }

  for (typename T::Periodic_triangle_iterator ptit = t_const.periodic_triangles_begin();
       ptit != t_const.periodic_triangles_end(); ++ptit)
    {
    }
  for (typename T::Periodic_triangle_iterator ptit =
         t_const.periodic_triangles_begin(T::STORED);
       ptit != t_const.periodic_triangles_end(T::STORED); ++ptit)
    {
    }
  for (typename T::Periodic_triangle_iterator ptit =
         t_const.periodic_triangles_begin(T::UNIQUE);
       ptit != t_const.periodic_triangles_end(T::UNIQUE); ++ptit)
    {
    }
  for (typename T::Periodic_triangle_iterator ptit =
         t_const.periodic_triangles_begin(T::STORED_COVER_DOMAIN);
       ptit != t_const.periodic_triangles_end(T::STORED_COVER_DOMAIN); ++ptit)
    {
    }
  for (typename T::Periodic_triangle_iterator ptit =
         t_const.periodic_triangles_begin(T::UNIQUE_COVER_DOMAIN);
       ptit != t_const.periodic_triangles_end(T::UNIQUE_COVER_DOMAIN); ++ptit)
    {
    }
}


template <class T>
void test_circulators()
{
  typedef typename T::Vertex_handle     Vertex_handle;
  typedef typename T::Point             Point;

  T t;
  const T &t_const = t;

  Point p0(0.5, 0.5);
  Vertex_handle vh0 = t.insert(p0);
  CGAL_USE(vh0);
  Vertex_handle vh1 = t.insert(Point(0.7, 0.5));
  CGAL_USE(vh1);
  Vertex_handle vh2 = t.insert(Point(0.7, 0.7));
  CGAL_USE(vh2);

  typename T::Face_circulator fcir = t_const.incident_faces(vh0);
  fcir = t_const.incident_faces(vh0, fcir);

  typename T::Edge_circulator ecir = t_const.incident_edges(vh0);
  ecir = t_const.incident_edges(vh0, fcir);

  typename T::Vertex_circulator vcir = t_const.adjacent_vertices(vh0);
  vcir = t_const.adjacent_vertices(vh0, fcir);

  Vertex_handle v_mirror = t_const.mirror_vertex(fcir, 0);
  CGAL_USE(v_mirror);
  t_const.mirror_index(fcir, 0);
}

template <class T>
void test_modifiers()
{
  typedef typename T::Vertex_handle     Vertex_handle;
  typedef typename T::Face_handle       Face_handle;
  typedef typename T::Point             Point;

  T t;
  const T &t_const = t;

  Point p0(0.5, 0.5);
  Point p1(0.8, 0.6);
  Point p2(0.7, 0.7);

  Vertex_handle vh0 = t.insert(p0);
  Face_handle fh = t.faces_begin();
  Vertex_handle vh1 = t.insert(Point(0.7, 0.5), fh);
  Vertex_handle vh2 = t.insert(Point(0.7, 0.7));


  t.flip(fh, 0);

  t.clear();
  vh0 = t.insert(p0);
  vh1 = t.insert(p1);
  vh2 = t.insert(p2);

  typename T::Locate_type lt;
  int li;
  fh = t_const.locate(p0, lt, li);
  t.insert(p0, lt, fh, li);
  t.push_back(p0);

  t.clear();
  vh0 = t.insert(p0);
  t.remove_first(vh0);

  t.clear();
  vh0 = t.insert_first(p0);
  vh1 = t.insert(p1);
  if (t.degree(vh1) == 3)
    {
      // The vertex has degree 6 in case we are testing the Delaunay triangulation
      t.remove_degree_3(vh1);
    }

  t.clear();
  vh0 = t.insert_first(p0);

  fh = t_const.locate(p1, lt, li);
  CGAL_assertion(lt == T::FACE);
  t.insert_in_face(p1, fh);
  fh = t_const.locate(p2, lt, li);
  CGAL_assertion(lt == T::EDGE);
  t.insert_in_edge(p2, fh, li);

  for (typename T::Vertex_iterator vit = t_const.vertices_begin();
       vit != t_const.vertices_end(); ++vit)
    {
      if (t_const.degree(vit) == 3)
        {
          t.remove_degree_3(vit);
          vit = t_const.vertices_begin();
        }
    }

  t.clear();
  vh0 = t.insert(p0);
  t.remove_first(vh0);

  // star_hole is not tested
}

template <class T>
void test_miscellaneous()
{
  typedef typename T::Vertex_handle     Vertex_handle;
  typedef typename T::Face_handle       Face_handle;
  typedef typename T::Point             Point;

  T t;

  Point p0(0.5, 0.5);
  Point p1(0.8, 0.6);
  Point p2(0.7, 0.7);
  Vertex_handle vh0, vh1, vh2;
  vh0 = t.insert(p0);
  vh1 = t.insert(p1);
  vh2 = t.insert(p2);

  t.set_domain(typename T::Iso_rectangle(0, 0, 2, 2));
  int i = t.ccw(0);
  int j = t.cw(0);
  CGAL_USE(i);
  CGAL_USE(j);
  CGAL_assertion(i + j == 3);

  t = T();
  vh0 = t.insert(p0);
  vh1 = t.insert(p1);
  vh2 = t.insert(p2);
  Face_handle fh = t.faces_begin();
  t.flippable(fh, 0);
  t.degree(vh0);

  t.is_valid();
  t.is_valid(true);
  t.is_valid(false);
  t.is_valid(true, 0);
  t.is_valid(false, 0);
}

template <class T>
void test_io(T &pt1, bool ex)
{
  std::cout << "I/O" << std::endl;
  std::cout << "  ascii" << std::endl;

  std::stringstream ss1;
  ss1 << pt1;

  T pt1r;
  ss1 >> pt1r;

  assert(CGAL::is_ascii(ss1));
  if (!ex)
  {
    assert(pt1 == pt1r);
  }

  std::cout << "  binary" << std::endl;
  pt1r.clear();
  // There are problems with the IO of exact number types in binary mode.
  if (!ex)
  {
    std::stringstream ss1b;
    CGAL::set_binary_mode(ss1b);
    ss1b << pt1;
    
    ss1b >> pt1r;
    assert(CGAL::is_binary(ss1b));
    
    assert(pt1 == pt1r);
  }

  std::cout << "  pretty" << std::endl;

  pt1r.clear();
  std::stringstream ss1p;
  CGAL::set_pretty_mode(ss1p);
  ss1p << pt1;

  assert(CGAL::is_pretty(ss1p));
}

template <class T>
void test_io(bool exact)
{
  typedef typename T::Point             Point;

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  T t;
  test_io(t, exact);

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  t.insert(Point(0.5, 0.5));
  test_io(t, exact);

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  t.insert(Point(0.6, 0.8));
  t.insert(Point(0.2, 0.6));
  test_io(t, exact);

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  // One cover for the Delaunay triangulation
  t.clear();
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y)
      t.insert(Point(x / 5.0, y / 5.0));
  CGAL_assertion(t.number_of_vertices() == 25);
  test_io(t, exact);

  std::cout << __FILE__ << ", " << __LINE__ << std::endl;
}

template <class T>
void test(bool exact)
{
  test_constructor<T>();
  test_global_access<T>();
  test_geometric_access<T>();
  test_predicates<T>();
  test_queries<T>();
  test_iterators<T>();
  test_circulators<T>();
  test_modifiers<T>();
  test_miscellaneous<T>();
  test_io<T>(exact);
}



template <class T>
void test_batch_insertion()
{
  typedef typename T::Point             Point;

  Point p0(0.5, 0.5);
  Point p1(0.8, 0.6);
  Point p2(0.7, 0.7);

  std::vector<Point> pts;
  pts.push_back(p0);
  pts.push_back(p1);
  pts.push_back(p0);
  pts.push_back(p2);

  T t(pts.begin(), pts.end());


  t.insert(pts.begin(), pts.end());

  t.clear();
  t.insert(pts.begin(), pts.end(), true);

  t.clear();
  t.insert(pts.begin(), pts.end(), false);
}

template <class T>
void test_nearest()
{
  typedef typename T::Vertex_handle     Vertex_handle;
  typedef typename T::Point             Point;

  Point p0(0.5, 0.5);
  Point p1(0.8, 0.6);
  Point p2(0.7, 0.7);
  Vertex_handle vh0, vh1, vh2;

  T t;
  CGAL_assertion(t.nearest_vertex(p0) == Vertex_handle());

  vh0 = t.insert(p0);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p0)) == vh0);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p1)) == vh0);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p2)) == vh0);

  vh1 = t.insert(p1);
  vh2 = t.insert(p2);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p0)) == vh0);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p1)) == vh1);
  CGAL_assertion(t.get_original_vertex(t.nearest_vertex(p2)) == vh2);
}

template <class T>
void test_locally_delaunay()
{
  typedef typename T::Geom_traits       Gt;
  typedef typename Gt::Vector_2         Vector;
  typedef typename T::Point             Point;
  typedef typename T::Face_iterator     Face_iterator;

  typedef CGAL::Creator_uniform_2<typename Gt::FT, Point>    Creator;
  typedef CGAL::Random_points_in_square_2<Point, Creator>   Random_points_in_square;

  T t;

  CGAL::Random random(1284141159);
  Random_points_in_square g(0.499, random);
  Vector midpoint(0.5, 0.5);

  for (int i = 0; i < 10; ++i)
    t.insert(*(++g) + midpoint);

  for (Face_iterator fit = t.faces_begin(); fit != t.faces_end(); ++fit)
    {
      for (int i = 0; i < 3; ++i)
        {
          CGAL_assertion(t.locally_Delaunay(fit, i, fit->neighbor(i)));
        }
    }

  while (!t.is_1_cover())
    t.insert(*(++g) + midpoint);

  for (Face_iterator fit = t.faces_begin(); fit != t.faces_end(); ++fit)
    {
      for (int i = 0; i < 3; ++i)
        {
          CGAL_assertion(t.locally_Delaunay(fit, i, fit->neighbor(i)));
        }
    }


}

template <class T>
void test_delaunay()
{
  test_batch_insertion<T>();
  test_nearest<T>();
  test_locally_delaunay<T>();
}

#endif  // P2T2_INTERFACE_TEST_H
