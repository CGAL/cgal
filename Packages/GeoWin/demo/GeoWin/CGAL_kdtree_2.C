// ======================================================================
//
// Copyright (c) 1999 The GALIA Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the GALIA Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the GALIA Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The GALIA Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Trier University
// (Germany), Max-Planck-Institute Saarbrucken (Germany),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : demo/GeoWin/CGAL_kdtree_2.C
//
// ======================================================================

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
#include <CGAL/kdtree_d.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>               K;
typedef K::Point_2                            Point;
typedef K::Circle_2                           Circle;
typedef K::Iso_rectangle_2                    Rectangle;

typedef CGAL::Kdtree_interface_2d<Point>      kd_interface;
typedef CGAL::Kdtree_d<kd_interface>          kd_tree;
typedef kd_tree::Box                          box;

geo_scene rec_scene;

class geo_rs : public geowin_update<std::list<Point>,std::list<Circle> >
{
 void update(const std::list<Point>& L, std::list<Circle>& Outcl)
 {
    GeoWin* gw = get_geowin(rec_scene);
    std::list<Rectangle> rects;
    gw->get_objects(rec_scene,rects);
     
    Outcl.clear();
    CGAL::Kdtree_d<kd_interface>  tree(2);
    tree.build((std::list<Point>&)L); 
    std::list<Point> Out;
    std::list<Rectangle>::const_iterator it = rects.begin();
    std::list<Point>::const_iterator pit;

    for (; it != rects.end(); it++) {
      box B(it->min(), it->max(), 2);
      std::list<Point> res;
      tree.search( std::back_inserter( res ), B );
      std::copy(res.begin(), res.end(), std::back_inserter(Out));    
    }
    pit=Out.begin();
    for (; pit != Out.end(); pit++) Outcl.push_back(Circle(*pit,3.0));
 }
};

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
  geowin_init_default_type((std::list<Rectangle>*)0, leda_string("CGALRectangleList"));
 
  std::list<Point>     L;
  std::list<Rectangle> RL;

  GeoWin GW("CGALTEST - 2d tree");
  GW.message("Activate rectangle scene to perform rectangular range searches on the input point set");
 
  geo_scene my_scene= GW.new_scene(L);  
  GW.set_color(my_scene,leda_black);

  rec_scene= GW.new_scene(RL);  
  GW.set_color(rec_scene,leda_red);  
  GW.set_fill_color(rec_scene,leda_invisible);
  
  geo_rs RSR;
  geo_scene result  = GW.new_scene(RSR, my_scene, leda_string("Range Search"));
  GW.set_color(result,leda_blue2);
  GW.set_fill_color(result, leda_blue);
  GW.set_line_width(result, 3);
  GW.set_point_style(result, leda_disc_point);

  GW.add_dependence(rec_scene,result);
  GW.set_all_visible(true);

  GW.edit(my_scene);
  
  return 0;  
}

#endif
