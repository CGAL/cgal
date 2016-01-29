// #define CGAL_IDENTIFICATION_XY  2
// #define CGAL_SL_VERBOSE 1

#include <string>
#include <cstring>
#include <vector>
#include <fstream>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/IO/Arr_text_formatter.h>

typedef CGAL::Exact_rational                              Number_type;
typedef CGAL::Cartesian<Number_type>                      Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel> Base_geom_traits;
typedef Base_geom_traits::Point_2                         Base_point_2;
typedef Base_geom_traits::Curve_2                         Base_curve_2;
typedef Base_geom_traits::X_monotone_curve_2
  Base_x_monotone_curve_2;
typedef CGAL::Arr_curve_data_traits_2<Base_geom_traits,
                                      unsigned int,
                                      std::plus<unsigned int> >
                                                          Geom_traits;
typedef Geom_traits::Point_2                              Point_2;
typedef Geom_traits::Curve_2                              Curve_2;
typedef Geom_traits::X_monotone_curve_2                   X_monotone_curve_2;

typedef CGAL::Arr_extended_dcel<Geom_traits,
                                unsigned int, unsigned int, unsigned int>
                                                          Dcel;

typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits, Dcel>
                                                          Topol_traits;
typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                          Arrangement;

// Iterators
typedef Arrangement::Vertex_iterator                    Vertex_iterator;
typedef Arrangement::Halfedge_iterator                  Halfedge_iterator;
typedef Arrangement::Face_iterator                      Face_iterator;
typedef Arrangement::Outer_ccb_iterator                 Outer_ccb_iterator;
typedef Arrangement::Inner_ccb_iterator                 Inner_ccb_iterator;

// Const iterators
typedef Arrangement::Vertex_const_iterator              Vertex_const_iterator;
typedef Arrangement::Halfedge_const_iterator            Halfedge_const_iterator;
typedef Arrangement::Face_const_iterator                Face_const_iterator;
typedef Arrangement::Outer_ccb_const_iterator           Outer_ccb_const_iterator;
typedef Arrangement::Inner_ccb_const_iterator           Inner_ccb_const_iterator;

// Handles
typedef Arrangement::Vertex_handle                      Vertex_handle;
typedef Arrangement::Halfedge_handle                    Halfedge_handle;
typedef Arrangement::Face_handle                        Face_handle;

// Const handles
typedef Arrangement::Vertex_const_handle                Vertex_const_handle;
typedef Arrangement::Halfedge_const_handle              Halfedge_const_handle;
typedef Arrangement::Face_const_handle                  Face_const_handle;

// Circulators
typedef Arrangement::Ccb_halfedge_circulator            Ccb_halfedge_circulator;
typedef Arrangement::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;

// Const Circulators
typedef Arrangement::Ccb_halfedge_const_circulator
  Ccb_halfedge_const_circulator;
typedef Arrangement::Halfedge_around_vertex_const_circulator
  Halfedge_around_vertex_const_circulator;

//! Overlay traits
class Overlay_traits {
public:
  /*! Destructor. */
  virtual ~Overlay_traits() {}

  /*! Create a vertex v that corresponds to the coinciding vertices v1 and v2. */
  virtual void create_vertex(Vertex_const_handle v1, Vertex_const_handle v2,
                             Vertex_handle v) const
  { v->set_data(v1->data() + v2->data()); }

  /*! Create a vertex v that mathces v1, which lies of the edge e2. */
  virtual void create_vertex(Vertex_const_handle  v1, Halfedge_const_handle e2,
                             Vertex_handle v) const
  { v->set_data(v1->data() + e2->data()); }

  /*! Create a vertex v that mathces v1, contained in the face f2. */
  virtual void create_vertex(Vertex_const_handle v1, Face_const_handle f2,
                             Vertex_handle v) const
  { v->set_data(v1->data() + f2->data());}

  /*! Create a vertex v that mathces v2, which lies of the edge e1. */
  virtual void create_vertex(Halfedge_const_handle e1, Vertex_const_handle v2,
                             Vertex_handle v) const
  { v->set_data(e1->data() + v2->data()); }

  /*! Create a vertex v that mathces v2, contained in the face f1. */
  virtual void create_vertex(Face_const_handle f1, Vertex_const_handle v2,
                             Vertex_handle v) const
  { v->set_data(f1->data() + v2->data()); }

  /*! Create a vertex v that mathces the intersection of the edges e1 and e2. */
  virtual void create_vertex(Halfedge_const_handle e1, Halfedge_const_handle e2,
                             Vertex_handle v) const
  { v->set_data(e1->data() + e2->data()); }

  /*! Create an edge e that matches the overlap between e1 and e2. */
  virtual void create_edge(Halfedge_const_handle e1, Halfedge_const_handle e2,
                           Halfedge_handle e) const
  {
    e->set_data(e1->data() + e2->data());
    e->twin()->set_data(e1->data() + e2->data());
  }

  /*! Create an edge e that matches the edge e1, contained in the face f2. */
  virtual void create_edge(Halfedge_const_handle e1, Face_const_handle f2,
                           Halfedge_handle e) const
  {
    e->set_data(e1->data() + f2->data());
    e->twin()->set_data(e1->data() + f2->data());
  }

  /*! Create an edge e that matches the edge e2, contained in the face f1. */
  virtual void create_edge(Face_const_handle f1, Halfedge_const_handle e2,
                           Halfedge_handle e) const
  {
    e->set_data(f1->data() + e2->data());
    e->twin()->set_data(f1->data() + e2->data());
  }

  /*! Create a face f that matches the overlapping region between f1 and f2. */
  virtual void create_face(Face_const_handle f1, Face_const_handle f2,
                           Face_handle f) const
  { f->set_data(f1->data() + f2->data()); }
};

//! Initialize an arrangement
void init_arr(Arrangement& arr, int verbose_level)
{
  // Initialize the data of the halfedges of the arrangement.
  Halfedge_iterator heit;
  for (heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit)
    heit->set_data(heit->curve().data() *
                   ((heit->direction() == CGAL::ARR_LEFT_TO_RIGHT) ? 1 : 2));

  // Initialize the data of the faces of the arrangement.
  Face_iterator fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    unsigned int count = 0;

    // Outer ccb
    Outer_ccb_iterator ocit;
    for (ocit = fit->outer_ccbs_begin(); ocit != fit->outer_ccbs_end(); ++ocit) {
      Ccb_halfedge_circulator curr = *ocit;
      do count += curr->data() * 2;
      while (++curr != *ocit);
    }

    // Inner ccbs
    Inner_ccb_iterator icit;
    for (icit = fit->inner_ccbs_begin(); icit != fit->inner_ccbs_end(); ++icit) {
      Ccb_halfedge_circulator curr = *icit;
      do count += curr->data();
      while (++curr != *icit);
    }

    fit->set_data(count);
  }

  // Initialize the data of the vertices of the arrangement.
  Vertex_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    unsigned int count = 0;
    if (vit->is_isolated()) count = vit->face()->data();
    else {
      Halfedge_around_vertex_const_circulator curr = vit->incident_halfedges();
      do count += curr->data();
      while (++curr != vit->incident_halfedges());
    }
    vit->set_data(count);
  }

  if (verbose_level > 0) std::cout << "Arrangement Input: " << std::endl;

  if (verbose_level > 1) {
    std::cout << "Halfedge Data: " << std::endl;
    Halfedge_iterator heit;
    for (heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit)
      std::cout << heit->source()->point() << " "
                << heit->target()->point() << " " << heit->data()
                << std::endl;
  }

  if (verbose_level > 0) {
    std::cout << "Face Data: " << std::endl;
    Face_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
      std::cout << fit->data() << std::endl;
  }
}

template <typename Xcurves_iterator, typename Point_iterator>
void read_arr(std::ifstream& in, Xcurves_iterator out_xcurves,
              Point_iterator out_isolated_points)
{
  unsigned int i;

  // Read the points:
  unsigned int num_of_points;
  in >> num_of_points;
  std::vector<Point_2> points(num_of_points);
  for (i = 0; i < num_of_points; ++i)
    in >> points[i];

  // Read the x-monotone curves:
  unsigned int num_of_curves;
  in >> num_of_curves;
  for (i = 0; i < num_of_curves; ++i) {
    unsigned int j, k;
    in >> j >> k;
    Base_x_monotone_curve_2 base_xcv(points[j], points[k]);
    X_monotone_curve_2 xcv(base_xcv, 1);
    *out_xcurves++ = xcv;
  }

  // Read the isolated points.
  unsigned int num_of_isolated_points;
  in >> num_of_isolated_points;
  for (i = 0; i < num_of_isolated_points; ++i) {
    unsigned int j;
    in >> j;
    *out_isolated_points++ = points[j];
  }

  points.clear();
}

template <typename Curve_iterator, typename Point_iterator>
void construct_arr(Arrangement& arr,
                   Curve_iterator xcurves_begin, Curve_iterator xcurves_end,
                   Point_iterator points_begin, Point_iterator points_end,
                   int verbose_level)
{
#if 1
  // Insert the curves incrementally.<
  Curve_iterator cit;
  for (cit = xcurves_begin; cit != xcurves_end; ++cit) {
    if (verbose_level > 2) std::cout << "inserting " << *cit << " ... ";
    std::cout.flush();
    CGAL::insert_non_intersecting_curve(arr, *cit);
    if (verbose_level > 2) std::cout << "inserted" << std::endl;
  }
#else
  // Insert the curves aggregately.
  if (verbose_level > 2) std::cout << "inserting x-monotone curves" << " ... ";
  std::cout.flush();
  CGAL::insert_non_intersecting_curves(arr, xcurves_begin, xcurves_end);
  if (verbose_level > 2) std::cout << "inserted" << std::endl;
#endif

  // Insert the isolated points.
  if (verbose_level > 2) std::cout << "inserting isolated vertices" << " ... ";
  Point_iterator pit;
  for (pit = points_begin; pit != points_end; ++pit) {
    Point_2 point(*pit);
    CGAL::insert_point(arr, point);
  }
  if (verbose_level > 2) std::cout << "inserted" << std::endl;
}

bool test_one_file(std::ifstream& in, int verbose_level)
{
  std::list<X_monotone_curve_2> xcurves;
  std::list<Point_2> isolated_points;
  read_arr(in, std::back_inserter(xcurves),
           std::back_inserter(isolated_points));
  Arrangement arr1;
  construct_arr(arr1, xcurves.begin(), xcurves.end(),
                isolated_points.begin(), isolated_points.end(), verbose_level);
  isolated_points.clear();
  xcurves.clear();
  init_arr(arr1, verbose_level);

  read_arr(in, std::back_inserter(xcurves),
           std::back_inserter(isolated_points));
  Arrangement arr2;
  construct_arr(arr2, xcurves.begin(), xcurves.end(),
                isolated_points.begin(), isolated_points.end(), verbose_level);
  isolated_points.clear();
  xcurves.clear();
  init_arr(arr2, verbose_level);

  // Read the number of cells left.
  Arrangement::Size num_vertices, num_edges, num_faces;
  in >> num_vertices >> num_edges >> num_faces;

  // Read the expected face data:
  std::vector<int> fdata(num_faces);
  std::copy(std::istream_iterator<int>(in), std::istream_iterator<int>(),
            fdata.begin());

  Arrangement arr;
  Overlay_traits overlay_traits;
  if (verbose_level > 2) std::cout << "overlaying" << " ... "; std::cout.flush();
  CGAL::overlay(arr1, arr2, arr, overlay_traits);
  if (verbose_level > 2) std::cout << "overlaid" << std::endl;

  // Verify the resulting arrangement.
  Arrangement::Size num_vertices_res = arr.number_of_vertices();
  Arrangement::Size num_edges_res = arr.number_of_edges();
  Arrangement::Size num_faces_res = arr.number_of_faces();

  if (verbose_level > 0) std::cout << "Arrangement Output: " << std::endl;

  if (verbose_level > 1) {
    std::cout << "Halfedge Data: " << std::endl;
    Halfedge_iterator heit;
    for (heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit)
      std::cout << heit->source()->point() << " "
                << heit->target()->point() << " " << heit->data()
                << std::endl;
  }

  if (verbose_level > 0) {
    std::cout << "Face Data: " << std::endl;
    Face_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
      std::cout << fit->data() << std::endl;
  }

  if ((num_vertices_res != num_vertices) ||
      (num_edges_res != num_edges) ||
      (num_faces_res != num_faces))
  {
    std::cerr << "ERROR: The number of arrangement cells is incorrect:"
              << std::endl
              << "   V = " << arr.number_of_vertices()
              << ", E = " << arr.number_of_edges()
              << ", F = " << arr.number_of_faces()
              << std::endl;
    arr.clear();
    return false;
  }

  std::vector<int> fdata_res(num_faces);
  std::vector<int>::iterator it = fdata_res.begin();
  Face_iterator fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
     *it++ = fit->data();
  std::sort(fdata_res.begin(), fdata_res.end());
  if (! std::equal(fdata_res.begin(), fdata_res.end(), fdata.begin())) {
    std::cerr << "ERROR: Incorrect face data:" << std::endl;
    std::copy(fdata_res.begin(), fdata_res.end(),
              std::ostream_iterator<int>(std::cout, "\n"));
    arr.clear();
    return false;
  }

  arr.clear();
  return true;
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    return -1;
  }

  // TBD: Replace with better parsing!
  int i = 1;
  int verbose_level = 0;
  if (argc > 2) {
    if ((argc > 3) && (std::strncmp(argv[1], "-v", 2) == 0)) {
      verbose_level = atoi(argv[2]);
      i += 2;
    }
  }

  int success = 0;
  for (; i < argc; ++i) {
    std::string str(argv[i]);
    if (str.empty()) continue;

    std::string::iterator itr = str.end();
    --itr;
    while (itr != str.begin()) {
      std::string::iterator tmp = itr;
      --tmp;
      if (*itr == 't')  break;

      str.erase(itr);
      itr = tmp;
    }
    if (str.size() <= 1) continue;
    std::ifstream inp(str.c_str());
    if (!inp.is_open()) {
      std::cerr << "Failed to open " << str << std::endl;
      return -1;
    }
    if (! test_one_file(inp, verbose_level)) {
      inp.close();
      std::cerr << str << ": ERROR" << std::endl;
      success = -1;
    }
    else std::cout <<str << ": succeeded" << std::endl;
    inp.close();
  }

  return success;
}
