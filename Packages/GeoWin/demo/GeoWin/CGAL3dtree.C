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
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>
#include  <CGAL/kdtree_d.h>

typedef CGAL::Kdtree_interface_3d<CGALPoint_3>  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>            kd_tree;
typedef kd_tree::Box                            box;

geo_scene rec_scene;

static void show_d3(geo_scene sc, leda_d3_window& W, GRAPH<leda_d3_point,int>& H)
{
 GeoWin* gw = get_geowin(sc);
 CGALPoint_3_list L;
 gw->get_objects(sc,L);
 CGALRectanglelist rects;
 gw->get_objects(rec_scene,rects);

 CGALPoint_3 p;
 CGALPoint_3_list::const_iterator itp = L.begin();

 CGAL::Kdtree_d<kd_interface>  tree(3);
 tree.build(L);    

 std::list<CGALPoint_3> Out;
 std::list<CGALRectangle>::const_iterator it = rects.begin();
 CGALPoint_3 left,right;
 CGALPoint lhelp,rhelp;

 for (; it != rects.end(); it++) {
      lhelp= it->min(); rhelp= it->max();
      left = CGALPoint_3(lhelp.x(),lhelp.y(),-1000.0);
      right= CGALPoint_3(rhelp.x(),rhelp.y(), 1000.0);
      box B(left,right, 3);
      std::list<CGALPoint_3> res;
      tree.search( std::back_inserter( res ), B );
      std::copy(res.begin(), res.end(), std::back_inserter(Out));    
 }
 
 for (itp= Out.begin(); itp != Out.end();++itp) { p=*itp; H.new_node(convert_to_leda(p)); }
 
 leda_node_array<leda_vector> pos(H);
 leda_node v;
 forall_nodes(v,H) pos[v] = H[v].to_vector();
 W.init(pos);
}

int main()
{
  geowin_init_default_type((CGALPoint_3_list*)0, leda_string("CGALPoint_3_List"));
  geowin_init_default_type((CGALRectanglelist*)0, leda_string("CGALRectangleList"));

  CGALPoint_3_list L;
  CGALRectanglelist RL;

  GeoWin GW("CGALTEST - 3d tree");
 
  geo_scene my_scene= GW.new_scene(L);  
  GW.set_color(my_scene,leda_black);
  GW.set_d3_fcn(my_scene, show_d3);

  rec_scene= GW.new_scene(RL);  
  GW.set_color(rec_scene,leda_red);  
  GW.set_fill_color(rec_scene,leda_invisible);

  GW.set_all_visible(true);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
