#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/geowin_support.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef double NT;
typedef CGAL::Cartesian<NT> rep3_t;
typedef CGAL::Point_3<rep3_t> point_3;
typedef CGAL::Triangulation_geom_traits_3<rep3_t>  traits_3;
typedef CGAL::Triangulation_vertex_base_3<traits_3>     Vb ;
typedef CGAL::Triangulation_cell_base_3<traits_3>       Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> TDS3 ;
typedef CGAL::Triangulation_3< traits_3 , TDS3> Triangulation_3;
typedef CGAL::Delaunay_triangulation_3<traits_3,TDS3> Delaunay_3;

static void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene< CGALPoint_3_list >* esc = (GeoEditScene< CGALPoint_3_list > *) sc;
 CGALPoint_3_list& L= esc->get_objref();

 CGALPoint_3 p;
 CGALPoint_3_list::const_iterator it = L.begin();

 // input points ...
 for(; it != L.end();++it) { p=*it; H.new_node(convert_to_leda(p)); }

 // compute triangulation...
 Delaunay_3  G_delaunay;
 G_delaunay.clear();
 G_delaunay.insert(L.begin(),L.end());

 // construct graph for visualization ...
 typedef Delaunay_3::Vertex_handle Vertex_handle;
 typedef Delaunay_3::Cell_handle Cell_handle;
 typedef Delaunay_3::Vertex Vertex;
 typedef Delaunay_3::Cell Cell;
 typedef Delaunay_3::Edge_iterator Edge_iterator;
 typedef Delaunay_3::Cell_iterator Cell_iterator;
 typedef Delaunay_3::Triangle  Triangle; 

 Vertex_handle v1, v2;
 Cell_handle f;
 int n1, n2;

 Edge_iterator eit = G_delaunay.finite_edges_begin();
 Edge_iterator beyond = G_delaunay.edges_end();

 // edges ...
 for ( ;eit != beyond; ++eit) {
	f = (*eit).first;
	n1 = (*eit).second;
       	n2 = (*eit).third;
	v1 = f->vertex(n1);
	v2 = f->vertex(n2);
        leda_node n1 = H.new_node(convert_to_leda(v1->point()));
	leda_node n2 = H.new_node(convert_to_leda(v2->point()));
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
  geowin_init_default_type((CGALPoint_3_list*)0, leda_string("CGALPoint_3_List"));
 
  CGALPoint_3_list L;
  GeoWin GW("CGAL 3d Delaunay triangulation");
  GW.message("To show the 3d Delaunay triangulation use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
