//! \file examples/Arrangement_on_surface_2/aggregated_insertion.cpp
// Using the global aggregated insertion functions.

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_spherical_arc_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <list>

typedef CGAL::Gmpq                                            Number_type;
typedef CGAL::Cartesian<Number_type>                          Kernel;
typedef CGAL::Arr_spherical_arc_traits_2<Kernel>              Geom_traits_2;
typedef Geom_traits_2::Point_2                                Point_2;
typedef Geom_traits_2::X_monotone_curve_2                     Spherical_arc_3;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2>  Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                              Arrangement_2;
typedef Arrangement_2::Vertex_const_handle              Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle            Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle                Face_const_handle;

typedef std::pair<Vertex_const_handle, std::pair<CGAL::Object, CGAL::Object> >
                                                        Vert_decomp_entry;
typedef std::list<Vert_decomp_entry>                    Vert_decomp_list;

int main ()
{
  // Construct the arrangement of five intersecting segments.
  Arrangement_2 arr;
  std::list<Spherical_arc_3> arcs;

#if 1
  arcs.push_back(Spherical_arc_3(Point_2(1, 0, 0), Point_2(1, 1, 1)));
  arcs.push_back(Spherical_arc_3(Point_2(1, 1, 1), Point_2(0, 1, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 1, 0), Point_2(1, 0, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(1, 0, 0), Point_2(1, 1, -1)));
  arcs.push_back(Spherical_arc_3(Point_2(1, 1, -1), Point_2(0, 1, 0)));
#else
  arcs.push_back(Spherical_arc_3(Point_2(2, 1, 0), Point_2(0, 1, 2)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 1, 2), Point_2(0, 2, 0)));
  arcs.push_back(Spherical_arc_3(Point_2(0, 2, 0), Point_2(2, 1, 0)));
#endif
  
  insert_non_intersecting_curves (arr, arcs.begin(), arcs.end());

  Vert_decomp_list  vd_list;

  CGAL::decompose (arr, std::back_inserter(vd_list));

  // Print the results.
  Vert_decomp_list::const_iterator       vd_iter;
  std::pair<CGAL::Object, CGAL::Object>  curr;
  Vertex_const_handle                    vh;
  Halfedge_const_handle                  hh;
  Face_const_handle                      fh;

  for (vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
    curr = vd_iter->second;
    std::cout << "Vertex (" << vd_iter->first->point() << ") : ";

    std::cout << " feature below: ";
    if (CGAL::assign (hh, curr.first))
      std::cout << '[' << hh->curve() << ']';
    else if (CGAL::assign (vh, curr.first))
      std::cout << '(' << vh->point() << ')';
    else if (CGAL::assign (fh, curr.first))
      std::cout << "NONE";
    else
      std::cout << "EMPTY";

    std::cout << "   feature above: ";
    if (CGAL::assign (hh, curr.second))
      std::cout << '[' << hh->curve() << ']' << std::endl;
    else if (CGAL::assign (vh, curr.second))
      std::cout << '(' << vh->point() << ')' << std::endl;
    else if (CGAL::assign (fh, curr.second))
      std::cout << "NONE" << std::endl;
    else
      std::cout << "EMPTY" << std::endl;
  }

  return 0;
}
