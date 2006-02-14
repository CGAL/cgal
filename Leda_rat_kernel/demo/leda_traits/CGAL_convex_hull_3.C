// provide 3d kernel traits ...
#define CGAL_NO_POSTCONDITIONS  
#define CGAL_NO_PRECONDITIONS  
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#include <CGAL/basic.h>

#if (!defined(CGAL_USE_LEDA) || (__LEDA__ < 420)) 
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.2 or higher installed!\n";
 std::cout << "A LEDA version >= 4.2 is required !\n";
 return 0;
}
#else 

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/geowin_support.h>

// include spezializations for LEDA d3_rat_point's
// of the projective traits ...
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_xy_traits_leda_rat_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_xz_traits_leda_rat_2.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/Convex_hull_projective_yz_traits_leda_rat_2.h>

#include <CGAL/Convex_hull_traits_3.h>

namespace CGAL {
template <>
class Max_coordinate_3<leda_rat_vector> 
{
public:

    int operator()(const leda_rat_vector& v)
    {
      if (CGAL_NTS abs(v.xcoord()) >= CGAL_NTS abs(v.ycoord()))
      {
         if (CGAL_NTS abs(v.xcoord()) >= CGAL_NTS abs(v.zcoord())) return 0;
         return 2;
      }
      else
      {
         if (CGAL_NTS abs(v.ycoord()) >= CGAL_NTS abs(v.zcoord())) return 1;
         return 2;
      }
    }
};
}

#include <CGAL/convex_hull_3.h>
#include <CGAL/predicates_on_points_3.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits              K;
typedef K::Point_3                                Point_3;
typedef K::Segment_3                              Segment_3;
typedef K::Triangle_3                             Triangle_3;
typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef Polyhedron_3::Vertex_handle               Vertex_handle;
typedef Polyhedron_3::Halfedge_iterator           Halfedge_iterator;


void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene<std::list<Point_3> >* esc = (GeoEditScene<std::list<Point_3> > *) sc;
 std::list<Point_3>& L= esc->get_objref();

 Point_3 p;
 std::list<Point_3>::const_iterator it = L.begin();

 for(; it != L.end();++it) { p=*it; H.new_node(p.to_float()); }

 Traits       tr; 
 CGAL::Object ch_object;

 // compute convex hull 
 CGAL::convex_hull_3(L.begin(), L.end(), ch_object, tr);

 // visualize polyhedron...
 Polyhedron_3           Pol;
 Segment_3              segment;
 Triangle_3             triangle;
 Point_3                lpoint;
 
 if (CGAL::assign(Pol, ch_object)) {
  Halfedge_iterator hit = Pol.halfedges_begin();
  for (; hit != Pol.halfedges_end(); hit++) {
    Vertex_handle v1 = hit->vertex();
    Vertex_handle v2 = hit->opposite()->vertex();
    leda_d3_point ps = (v1->point()).to_float();
    leda_d3_point pt = (v2->point()).to_float();

    leda_node n1= H.new_node(ps);
    leda_node n2= H.new_node(pt);
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);
  }
 }
 
 else if (CGAL::assign(segment, ch_object) ){
    leda_d3_point ps = (segment.source()).to_float();
    leda_d3_point pt = (segment.target()).to_float(); 
 
    leda_node n1= H.new_node(ps);
    leda_node n2= H.new_node(pt);
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);    
 }
 else if (CGAL::assign(triangle, ch_object) ){
    leda_d3_point ps = triangle.point1().to_float();
    leda_d3_point pt = triangle.point2().to_float();
    leda_d3_point pu = triangle.point3().to_float(); 
    
    leda_node n1= H.new_node(ps);
    leda_node n2= H.new_node(pt);
    leda_node n3= H.new_node(pu);
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
    H.set_reversal(e1,e2);  
    leda_edge e3= H.new_edge(n2,n3), e4= H.new_edge(n3,n2);
    H.set_reversal(e3,e4);  
    leda_edge e5= H.new_edge(n3,n1), e6= H.new_edge(n1,n3);
    H.set_reversal(e5,e6);              
 }
 else if (CGAL::assign(lpoint, ch_object) ){
    leda_d3_point p = lpoint.to_float(); 
    H.new_node(p);
 } 
  
 leda_node_array<leda_vector> pos(H);
 leda_node v;
 forall_nodes(v,H) pos[v] = H[v].to_vector();
 W.init(pos); 
}

int main()
{
  geowin_init_default_type((std::list<Point_3>*)0, leda_string("LEDA-d3_rat_point"));
 
  std::list<Point_3> L;
  GeoWin GW("3d Convex hull");
  GW.message("To show the convex hull use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
