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
#include <CGAL/leda_rational.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Cartesian<leda_rational>                         R;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R>  HDS;
typedef CGAL::Polyhedron_default_traits_3<R>                   Traits;
typedef CGAL::Polyhedron_3<Traits,HDS>                         Polyhedron;
typedef CGAL::Polyhedron_3<Traits,HDS>::Vertex_handle          Vertex_handle;
typedef Polyhedron::Halfedge_iterator                          Halfedge_iterator;


CGAL::Point_3<CGAL::Cartesian<leda_rational> > conversion_to_rat(CGAL::Point_3<CGAL::Cartesian<double> > pin)
{
 double x = pin.x();
 double y = pin.y();
 double z = pin.z();
 leda_d3_point pf(x,y,z);
 leda_d3_rat_point prat(leda_integer(x*100000),leda_integer(y*100000),leda_integer(z*100000),leda_integer(100000));
 return CGAL::Point_3<CGAL::Cartesian<leda_rational> >(prat.xcoord(),prat.ycoord(),prat.zcoord());
}

CGAL::Point_3<CGAL::Cartesian<double> > conversion_to_float(CGAL::Point_3<CGAL::Cartesian<leda_rational> > pin)
{
 leda_rational x = pin.x();
 leda_rational y = pin.y();
 leda_rational z = pin.z();
 leda_d3_rat_point prat(x,y,z);
 leda_d3_point pf = prat.to_d3_point();
 return CGAL::Point_3<CGAL::Cartesian<double> >(pf.xcoord(),pf.ycoord(),pf.zcoord());
}

static void show_d3_points(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoEditScene< CGALPoint_3_list >* esc = (GeoEditScene< CGALPoint_3_list > *) sc;
 CGALPoint_3_list& L= esc->get_objref();

 CGALPoint_3 p;
 CGALPoint_3_list::const_iterator it = L.begin();

 for(; it != L.end();++it) { p=*it; H.new_node(convert_to_leda(p)); }

 // compute hull...
 Polyhedron Pol;

 std::list<CGAL::Point_3<CGAL::Cartesian<leda_rational> > > L2;
 // convert list with FP values to rational values ...
 it = L.begin();
 for(; it != L.end();++it) {
    CGAL::Point_3<CGAL::Cartesian<leda_rational> > rpt = conversion_to_rat(*it);
    L2.push_back(rpt);
 }

 CGAL::convex_hull_3( L2.begin(), L2.end(), Pol);

 // visualize polyhedron...
 Halfedge_iterator hit = Pol.halfedges_begin();
 for (; hit != Pol.halfedges_end(); hit++) {
    Vertex_handle v1 = hit->vertex();
    Vertex_handle v2 = hit->opposite()->vertex();
    CGALPoint_3 ps = conversion_to_float(v1->point());
    CGALPoint_3 pt = conversion_to_float(v2->point());

    leda_node n1= H.new_node(convert_to_leda(ps));
    leda_node n2= H.new_node(convert_to_leda(pt));
    leda_edge e1= H.new_edge(n1,n2), e2= H.new_edge(n2,n1);
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
  GeoWin GW("Convex hull in 3d");
  GW.message("To show the convex hull use 'Show d3 output' in Window menu.");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_d3_fcn(my_scene, show_d3_points);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
