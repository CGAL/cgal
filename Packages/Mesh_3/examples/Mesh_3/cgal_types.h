//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Constrained_regular_triangulation_3.h>
#include <CGAL/Constrained_triangulation_vertex_base_3.h>
#include <CGAL/Constrained_triangulation_cell_base_3.h>

#include <CGAL/Triangulation_mesher_3.h>

#include <CGAL/Triangulation_2_traits_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

typedef CGAL::Quotient<CGAL::MP_Float> FT;
//struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
struct K : public CGAL::Simple_cartesian<FT> {};

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Traits_3;
typedef CGAL::Constrained_triangulation_cell_base_3<Traits_3> Cb_3;
typedef CGAL::Constrained_triangulation_vertex_base_3<Traits_3> Vb_3;
typedef CGAL::Triangulation_data_structure_3<Vb_3, Cb_3> Tds_3;
typedef CGAL::Regular_triangulation_3<Traits_3, Tds_3>          Rt_3;
