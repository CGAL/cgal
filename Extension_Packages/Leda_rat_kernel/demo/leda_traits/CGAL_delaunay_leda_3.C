// provide 3d kernel traits ...
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

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/geowin_support.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                         K;
typedef K::Point_3                                           Point_3;
typedef CGAL::Delaunay_triangulation_3<K>                    Delaunay_3;

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
 G_delaunay.clear();
 G_delaunay.insert(L.begin(),L.end());

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
