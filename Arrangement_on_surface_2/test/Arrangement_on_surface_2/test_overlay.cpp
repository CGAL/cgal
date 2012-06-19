/* Test the overlay operation.
 * 1. We read a pair of sets---a set of curves followed by a set of
 *    (isolated) vertices.
 * 2. We read another pair of such sets.
 * 3. We read an extended arrangement and store the result in arr_result.
 *    Each type of cell, i.e., vertex, halfedge, and face, is extended with
 *    an unsigned int.
 * 4. We construct two extended arrangements, namely arr1 and arr2, induced
 *    by the elements in the two pairs above, respectively.
 * 5. We initialize the data field of each halfedge with the number of curves
 *    that induced that halfedge if the halfedge is directed left-to-right
 *    and twice the the number of curves that induced that halfedge if the
 *    halfedge is directed right-to-left. We initialize the data field of
 *    each face with the total sum of the data of the halfedges on the
 *    boundary of the face. We initialize the data field of each isolated
 *    vertex with the data of the containing face, and the data field of each
 *    non-isolated vertex with the total sum of the data of the incident
 *    halfedges (the vertex is their target.
 * 6. We compute the overlay of arr1 and arr2 and store the result in a third
 *    arrangement called arr. The data field of each cell of arr is assigned
 *    with the sum of the data of the 2 inducing cells from arr1 and arr2,
 *    respectively.
 * 7. Finally, we test whether arr and arr_result are equivalent (isomorphic).
 */

#include <iostream>
#include <vector>
#include <map>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/IO/Arr_text_formatter.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Base_traits_2;
typedef Base_traits_2::Point_2                          Base_point_2;
typedef Base_traits_2::Curve_2                          Base_curve_2;
typedef CGAL::Arr_curve_data_traits_2<Base_traits_2,
                                      unsigned int,
                                      std::plus<unsigned int> >  
                                                        Traits;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;

typedef CGAL::Arr_extended_dcel<Traits,unsigned int,unsigned int,unsigned int>
                                                        Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arrangement;

typedef CGAL::Arr_extended_dcel_text_formatter<Arrangement>  Formatter;

// Iterators
typedef Arrangement::Vertex_iterator                    Vertex_iterator;
typedef Arrangement::Halfedge_iterator                  Halfedge_iterator;
typedef Arrangement::Face_iterator                      Face_iterator;
typedef Arrangement::Hole_iterator                      Hole_iterator;

// Const iterators
typedef Arrangement::Vertex_const_iterator              Vertex_const_iterator;
typedef Arrangement::Halfedge_const_iterator            Halfedge_const_iterator;
typedef Arrangement::Face_const_iterator                Face_const_iterator;
typedef Arrangement::Hole_const_iterator                Hole_const_iterator;

// Handles
typedef Arrangement::Vertex_handle                      Vertex_handle;
typedef Arrangement::Halfedge_handle                    Halfedge_handle;
typedef Arrangement::Face_handle                        Face_handle;

// Const handles
typedef Arrangement::Vertex_const_handle                Vertex_const_handle;
typedef Arrangement::Halfedge_const_handle              Halfedge_const_handle;
typedef Arrangement::Face_const_handle                  Face_const_handle;

// Circulators
typedef Arrangement::Ccb_halfedge_const_circulator
  Ccb_halfedge_const_circulator;
typedef Arrangement::Halfedge_around_vertex_const_circulator
  Halfedge_around_vertex_const_circulator;

typedef std::vector<Curve_2>                            Curve_vector;
typedef std::vector<Point_2>                            Point_vector;

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

//! Construct and initialize an arrangement
template <typename Curve_iterator, typename Point_iterator>
void init_arr(Arrangement& arr,
              Curve_iterator curves_begin, Curve_iterator curves_end,
              Point_iterator points_begin, Point_iterator points_end)
{
  // Construct the arrangement.
  CGAL::insert(arr, curves_begin, curves_end);
  Point_iterator it;
  for (it = points_begin; it != points_end; ++it)
    CGAL::insert_point(arr, *it);

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
    if (!fit->is_unbounded()) {
      Ccb_halfedge_const_circulator curr = fit->outer_ccb();
      do count += curr->data();
      while (++curr != fit->outer_ccb());
    }
    
    // Inner ccbs
    Hole_iterator hit;
    for (hit = fit->holes_begin(); hit != fit->holes_end(); ++hit) {
      Ccb_halfedge_const_circulator curr = *hit;
      do count += curr->data();
      while (++curr != *hit);
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
}

// Maps
struct Less_than_handle {
  template <typename Type>
  bool operator()(Type s1, Type s2) const { return (&(*s1) < &(*s2)); }
};

typedef std::map<Vertex_const_handle, Vertex_const_handle,
                 Less_than_handle>                      Vertex_map;
typedef std::map<Halfedge_const_handle, Halfedge_const_handle,
                 Less_than_handle>                      Halfedge_map;
typedef std::map<Face_const_handle, Face_const_handle, Less_than_handle>
                                                        Face_map;

// Check whether two ccb's are equivalent
bool equivalent_ccb(Ccb_halfedge_const_circulator ccb1,
                    Ccb_halfedge_const_circulator ccb2,
                    const Halfedge_map& halfedge_map)
{
  
  // Find a matching starting point.
  Halfedge_map::const_iterator it = halfedge_map.find(ccb1);
  if (it == halfedge_map.end()) return false;
  Ccb_halfedge_const_circulator curr2 = ccb2;
  do if ((*it).second == curr2) break;
  while (++curr2 != ccb2);
  if ((*it).second != curr2) return false;
  ccb2 = curr2;
  
  // Match the rest
  Ccb_halfedge_const_circulator curr1 = ccb1;
  do {
    Halfedge_map::const_iterator it = halfedge_map.find(curr1);
    if ((it == halfedge_map.end()) || ((*it).second != curr2++)) return false;
  } while (++curr1 != ccb1);
  if (curr2 != ccb2) return false;
  
  return true;
}

// Check whether two faces are equivalent
bool equivalent_face(Face_const_handle fit1, Face_const_handle fit2,
                     const Halfedge_map& halfedge_map)
{
  if (fit1->is_unbounded() != fit2->is_unbounded()) return false;
  
  if (!fit1->is_unbounded()) {
    if (!equivalent_ccb(fit1->outer_ccb(), fit2->outer_ccb(), halfedge_map))
      return false;
  }
  
  Hole_const_iterator hit1;
  for (hit1 = fit1->holes_begin(); hit1 != fit1->holes_end(); ++hit1) {
    Hole_const_iterator hit2;
    bool found = false;
    for (hit2 = fit2->holes_begin(); hit2 != fit2->holes_end(); ++hit2) {
      if (equivalent_ccb(*hit1, *hit2, halfedge_map)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  
  return true;
}

// Check whether two arrangement are equivalent
bool equivalent_arr(const Arrangement& arr1, const Arrangement& arr2)
{
  if (arr1.number_of_vertices() != arr2.number_of_vertices()) return false;
  if (arr1.number_of_halfedges() != arr2.number_of_halfedges()) return false;
  if (arr1.number_of_faces() != arr2.number_of_faces()) return false;

  const Traits* traits = arr1.traits();
  Traits::Equal_2 equal = traits->equal_2_object();

  // Compare the vertices
  const Vertex_const_iterator invalid_vit;
  Vertex_map vertex_map;
  Vertex_const_iterator vit1;
  for (vit1 = arr1.vertices_begin(); vit1 != arr1.vertices_end(); ++vit1) {
    const Point_2& p1 = vit1->point();
    Vertex_const_iterator vit2;
    bool found = false;
    for (vit2 = arr2.vertices_begin(); vit2 != arr2.vertices_end(); ++vit2) {
      const Point_2& p2 = vit2->point();
      if (equal(p1, p2)) {
        if (vertex_map[vit1] != invalid_vit) {
          std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                    << ") has been mapped already!"
                    << std::endl;
          return false;
        }
        vertex_map[vit1] = vit2;
        found = true;

        if (vit1->data() != vit2->data()) {
          std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                    << ") data does not match (" << vit2->data()
                    << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                << ") was not found!" << std::endl;
       return false;
    }
  }

  // Compare the halfedges.
  const Halfedge_const_iterator invalid_heit;
  Halfedge_map halfedge_map;
  Halfedge_const_iterator heit1;
  for (heit1 = arr1.halfedges_begin(); heit1 != arr1.halfedges_end(); ++heit1) {
    const X_monotone_curve_2& xcv1 = heit1->curve();
    Halfedge_const_iterator heit2;
    bool found = false;
    for (heit2 = arr2.halfedges_begin(); heit2 != arr2.halfedges_end(); ++heit2)
    {
      const X_monotone_curve_2& xcv2 = heit2->curve();
      if ((vertex_map[heit1->source()] == heit2->source()) &&
          (vertex_map[heit1->target()] == heit2->target()))
      {
        if (halfedge_map[heit1] != invalid_heit) {
          std::cerr << "The halfedge ((" << heit1->source()->point()
                    << " => " << heit1->target()->point() << "), "
                    << heit1->data() << ") has been mapped already!"
                    << std::endl;
          return false;
        }
        halfedge_map[heit1] = heit2;
        found = true;

        if (!equal(xcv1, xcv2)) {
          std::cerr << "The halfedge ((" << heit1->source()->point()
                    << " => " << heit1->target()->point() << "), "
                    << heit1->data() << ") x-monotone curve does not match ("
                    << heit2->curve() << ")!" << std::endl;
          return false;
        }
        
        if (heit1->data() != heit2->data()) {
          std::cerr << "The halfedge ((" << heit1->source()->point()
                    << " => " << heit1->target()->point() << "), "
                    << heit1->data() << ") data does not match ("
                    << heit2->data() << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The halfedge ((" << heit1->source()->point()
                << " => " << heit1->target()->point() << "), "
                << heit1->data() << ") was not found!"
                << std::endl;
      return false;
    }
  }

  // Compare the faces.
  const Face_const_iterator invalid_fit;
  Face_map face_map;
  Face_const_iterator fit1;
  for (fit1 = arr1.faces_begin(); fit1 != arr1.faces_end(); ++fit1) {
    Face_const_iterator fit2;
    bool found = false;
    for (fit2 = arr2.faces_begin(); fit2 != arr2.faces_end(); ++fit2) { 
      if (equivalent_face(fit1, fit2, halfedge_map)) {
        if (face_map[fit1] != invalid_fit) {
          std::cerr << "The face (" << fit1->data()
                    << ") has been mapped already!" << std::endl;
          return false;
        }
        face_map[fit1] = fit2;
        found = true;

        if (fit1->data() != fit2->data()) {
          std::cerr << "The face (" << fit1->data()
                    << ") data does not match ("
                    << fit2->data() << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The face (" << fit1->data()
                << ") was not found!" << std::endl;
      return false;
    }
  }
  
  return true;
}

// Test one file
bool test_one_file(std::ifstream& in_file, bool verbose)
{
  Kernel kernel;
  
  unsigned int i;
  unsigned int num_of_elements;

  // Read the first set of segments.
  in_file >> num_of_elements;
  Curve_vector curves1(num_of_elements);
  for (i = 0; i < num_of_elements; ++i) { 
    Base_point_2 source, target;
    in_file >> source >> target;
    assert(kernel.compare_xy_2_object()(source, target) != CGAL::EQUAL);
    Base_curve_2 base_cv(source, target);
    curves1[i] = X_monotone_curve_2(base_cv, 1);
  }

  // Read the first set of points.
  in_file >> num_of_elements;
  Point_vector points1(num_of_elements);
  for (i = 0; i < num_of_elements; ++i) { 
    Point_2 point;
    in_file >> point;
    points1[i] = point;
  }

  // Read the second set of segments.
  in_file >> num_of_elements;
  Curve_vector curves2(num_of_elements);
  for (i = 0; i < num_of_elements; ++i) {
    Base_point_2 source, target;
    in_file >> source >> target;
    assert(kernel.compare_xy_2_object()(source, target) != CGAL::EQUAL);
    Base_curve_2 base_cv(source, target);
    curves2[i] = X_monotone_curve_2(base_cv, 1);
  }

  // Read the first second of points.
  in_file >> num_of_elements;
  Point_vector points2(num_of_elements);
  for (i = 0; i < num_of_elements; ++i) { 
    Point_2 point;
    in_file >> point;
    points2[i] = point;
  }

  // Consume eol
  int c;
  while ((c = in_file.get()) != '\n');
  
  // Read the resulting arrangement.
  Formatter formatter;
  Arrangement arr_result;
  CGAL::read(arr_result, in_file, formatter);
  
  // Initialize the arrangements.
  Arrangement arr1, arr2;
  init_arr(arr1, curves1.begin(), curves1.end(), points1.begin(), points1.end());
  init_arr(arr2, curves2.begin(), curves2.end(), points2.begin(), points2.end());

  // Overlay the arrangement.
  Arrangement arr;
  Overlay_traits overlay_traits;
  CGAL::overlay(arr1, arr2, arr, overlay_traits);

  // CGAL::write(arr, std::cout, formatter);
  // return true;
  return equivalent_arr(arr, arr_result);
}

//! main
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    std::exit(-1);
  }

  bool  verbose = false;
  
  if (argc > 2 && std::strncmp(argv[2], "-v", 2) == 0)
    verbose = true;

  int success = 0;
  for (int i = 1; i < argc; ++i) {
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
      return (-1);
    }
    if (! test_one_file(inp, verbose)) {
      inp.close();
      std::cerr << str << ": ERROR" << std::endl;
      success = -1;
    }
    else std::cout << str << ": succeeded" << std::endl;
    inp.close();
  }
  
  return success;
}
