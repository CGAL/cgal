// provide 3d kernel traits ...
#define CGAL_NO_DEPRECATED_CODE
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/basic.h>

#if (!defined(CGAL_USE_LEDA) || (__LEDA__ < 430)) 
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.3 or higher installed!\n";
 std::cout << "A LEDA version >= 4.3 is required !\n";
 return 0;
}
#else 

#include <CGAL/Kernel_special.h>
#include <CGAL/kernel_event_support.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/geowin_support.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

//typedef CGAL::leda_rat_kernel_traits                         K;

typedef CGAL::leda_rat_kernel_traits                         LEDA_KERNEL;
typedef CGAL::kernel_event<LEDA_KERNEL>                      KEV;
typedef CGAL::kernel_event<int>                              KRES;
typedef CGAL::Kernel_special<LEDA_KERNEL, KEV, KRES>         K;



typedef K::Point_3                                           Point_3;
typedef K::Segment_3                                         Segment_3;
typedef K::Triangle_3                                        Triangle_3;
typedef K::Tetrahedron_3                                     Tetrahedron_3;
typedef CGAL::Delaunay_triangulation_3<K>                    Delaunay_3;

// ------------------------------------------------------------------------
// event handling ...
// ------------------------------------------------------------------------

// triang ...
CGAL::event_item  ev_construct_segment_3;
CGAL::event_item  ev_construct_triangle_3;
CGAL::event_item  ev_construct_tetrahedron_3;
CGAL::event_item  ev_compare_xyz_3;
CGAL::event_item  ev_orientation_3;
CGAL::event_item  ev_coplanar_orientation_3;

// + Delaunay
CGAL::event_item  ev_side_of_oriented_sphere_3;
CGAL::event_item  ev_coplanar_side_of_bounded_circle_3;
CGAL::event_item  ev_compare_distance_3;

int ori_cnt = 0;
int sos_cnt = 0;

// functions ...
void  construct_segment(const LEDA_KERNEL::Construct_segment_3&, const Point_3& p1, const Point_3& p2, const Segment_3& s)
{
  std::cout << "construct_segment:" << p1 << " " << p2 << "\n";
}

void  construct_triangle(const LEDA_KERNEL::Construct_triangle_3&, const Point_3& p1, const Point_3& p2, const Point_3& p3, const Triangle_3&)
{
  std::cout << "construct_triangle:" << p1 << " " << p2 << " " << p3 << "\n";
}

void  construct_simplex(const LEDA_KERNEL::Construct_tetrahedron_3&, 
                        const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4,
			const Tetrahedron_3&)
{
  std::cout << "construct_simplex:" << p1 << " " << p2 << " " << p3 << " " << p4 << "\n";
}

void  cmp_xyz(const LEDA_KERNEL::Compare_xyz_3&, const Point_3& p1, const Point_3& p2, const CGAL::Comparison_result&)
{
  std::cout << "cmp_xy:" << p1 << " " << p2 << "\n";
}

void  orientation_3d(const LEDA_KERNEL::Orientation_3&, 
                     const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4,
		     const CGAL::Orientation&)
{
  std::cout << "orientation:" << p1 << " " << p2 << " " << p3 << " " << p4 << "\n";
  ori_cnt++;
}

void  coplanar_orientation(const LEDA_KERNEL::Coplanar_orientation_3&,
                           const Point_3& p1, const Point_3& p2, const Point_3& p3,
                           const CGAL::Orientation&)
{
  std::cout << "coplanar_orientation:" << p1 << " " << p2 << " " << p3 << "\n";
}

void  sos(const LEDA_KERNEL::Side_of_oriented_sphere_3&,
          const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4, const Point_3& p5,
	  const CGAL::Oriented_side&)
{
  std::cout << "side_of_sphere:" << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << "\n";
  sos_cnt++;
}

void  coplanar_soc(const LEDA_KERNEL::Coplanar_side_of_bounded_circle_3&,
                   const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4,
	           const CGAL::Bounded_side&)
{
  std::cout << "coplanar soc:" << p1 << " " << p2 << " " << p3 << " " << p4 << "\n";
}

void  cmp_dist(const LEDA_KERNEL::Compare_distance_3&,
               const Point_3& p1, const Point_3& p2, const Point_3& p3,
	       const CGAL::Comparison_result&)
{
  std::cout << "cmp_dist:" << p1 << " " << p2 << " " << p3 << "\n";
}

void attach_events()
{
  ev_construct_segment_3       = CGAL::attach(KRES::EVENT, construct_segment);
  ev_construct_triangle_3      = CGAL::attach(KRES::EVENT, construct_triangle);
  ev_construct_tetrahedron_3   = CGAL::attach(KRES::EVENT, construct_simplex);
  ev_compare_xyz_3             = CGAL::attach(KRES::EVENT, cmp_xyz);
  ev_orientation_3             = CGAL::attach(KRES::EVENT, orientation_3d);
  ev_coplanar_orientation_3    = CGAL::attach(KRES::EVENT, coplanar_orientation);
  ev_side_of_oriented_sphere_3 = CGAL::attach(KRES::EVENT, sos);
  ev_coplanar_side_of_bounded_circle_3 = CGAL::attach(KRES::EVENT, coplanar_soc);
  ev_compare_distance_3        = CGAL::attach(KRES::EVENT, cmp_dist);
}

void detach_events()
{
  CGAL::detach(ev_construct_segment_3);
  CGAL::detach(ev_construct_triangle_3);
  CGAL::detach(ev_construct_tetrahedron_3);
  CGAL::detach(ev_compare_xyz_3);
  CGAL::detach(ev_orientation_3);
  CGAL::detach(ev_coplanar_orientation_3);
  CGAL::detach(ev_side_of_oriented_sphere_3);
  CGAL::detach(ev_coplanar_side_of_bounded_circle_3);
  CGAL::detach(ev_compare_distance_3);
  
  std::cout << "orientation counter:   " << ori_cnt  << "\n";
  std::cout << "side_of_sphere counter:" << sos_cnt  << "\n";    
  
  ori_cnt = 0;
  sos_cnt = 0;  
}

void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene<std::list<Point_3> >* esc = (GeoEditScene<std::list<Point_3> > *) sc;
 std::list<Point_3> & L= esc->get_objref();

 Point_3 p;
 std::list<Point_3> ::const_iterator it = L.begin();

 // input points ...
 for(; it != L.end();++it) { p=*it; H.new_node(p.to_float()); }

 // compute triangulation...
 Delaunay_3  G_delaunay;

 attach_events();
 G_delaunay.insert(L.begin(),L.end()); 
 detach_events();

 // construct graph for visualization ...
 typedef Delaunay_3::Vertex_handle Vertex_handle;
 typedef Delaunay_3::Cell_handle Cell_handle;
 typedef Delaunay_3::Finite_edges_iterator Finite_edges_iterator;

 Vertex_handle v1, v2;
 Cell_handle f;
 int n1, n2;

 Finite_edges_iterator eit = G_delaunay.finite_edges_begin();
 Finite_edges_iterator beyond = G_delaunay.finite_edges_end();

 // edges ...
 for ( ;eit != beyond; ++eit) {
	f  = (*eit).first;
	n1 = (*eit).second;
       	n2 = (*eit).third;
	v1 = f->vertex(n1);
	v2 = f->vertex(n2);
        leda_node n1 = H.new_node((v1->point()).to_float() );
	leda_node n2 = H.new_node((v2->point()).to_float() );
	leda_edge e1 = H.new_edge(n1,n2);
	leda_edge e2 = H.new_edge(n2,n1);
	H.set_reversal(e1,e2);
 } 
 
 leda_node_array<leda_vector> pos(H);
 leda_node v;
 forall_nodes(v,H) pos[v] = H[v].to_vector();
 W.init(pos); 
}

int main()
{
  geowin_init_default_type((std::list<Point_3> *)0, leda_string("LEDA-d3_rat_point"));
 
  std::list<Point_3>  L;
  GeoWin GW("CGAL 3d Delaunay triangulation");
  GW.message("To show the 3d Delaunay triangulation use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
