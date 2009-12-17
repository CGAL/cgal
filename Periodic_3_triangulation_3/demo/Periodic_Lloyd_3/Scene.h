#ifndef SCENE_H
#define SCENE_H

#include "typedefs.h"

struct Scene {

  std::list<Point_3> points;
  P3DT3 periodic_triangulation;

  bool eight_copies;
  bool two_dimensional;

  void lloyd_step() {
    std::vector<Vertex_handle> vts;
    for (Periodic_point_iterator ppit
	   = periodic_triangulation.periodic_points_begin(P3DT3::UNIQUE) ;
	 ppit != periodic_triangulation.periodic_points_end(P3DT3::UNIQUE) ;
	 ++ppit)
      vts.push_back(ppit.get_vertex());
    
    points.clear();
    std::vector<Point_3> dual_vertices;
    Point_3 new_point;
    for (std::vector<Vertex_handle>::iterator vit = vts.begin();
	 vit != vts.end(); ++vit) {
      if (two_dimensional) {
	dual_vertices.clear();
	periodic_triangulation.dual(*vit,std::back_inserter(dual_vertices));
	new_point = compute_barycenter_2D(dual_vertices);
      } else {
	new_point = periodic_triangulation.dual_centroid(*vit);
      }
      points.push_back(new_point);
    }
    periodic_triangulation.clear();
    periodic_triangulation.insert(points.begin(),points.end());
  }

  Point_3 compute_barycenter_2D(std::vector<Point_3> dual_pts) const {
    T2 t;
    for (unsigned int i=0 ; i<dual_pts.size() ; i++)
      t.insert(Point_2(dual_pts[i].x(),dual_pts[i].y()));
    FT area(0);
    FT x(0),y(0);
    for (T2::Finite_faces_iterator ffit = t.finite_faces_begin() ;
	 ffit != t.finite_faces_end() ; ++ffit) {
      Triangle_2 tri = t.triangle(ffit);
      FT triarea = tri.area();
      Point_2 tricentr = centroid(tri);
      area += triarea;
      x += triarea * tricentr.x();
      y += triarea * tricentr.y();
    }

    x /= area;
    y /= area;

    Iso_cuboid_3 d = periodic_triangulation.domain();
    x = (x < d.xmin() ? x+d.xmax()-d.xmin() 
	: (x >= d.xmax() ? x-d.xmax()+d.xmin() : x));
    y = (y < d.ymin() ? y+d.ymax()-d.ymin() 
	: (y >= d.ymax() ? y-d.ymax()+d.ymin() : y));

    CGAL_triangulation_postcondition((d.xmin()<=x)&&(x<d.xmax()));
    CGAL_triangulation_postcondition((d.ymin()<=y)&&(y<d.ymax()));
  
    return Point_3(x,y,0);
  }

};

#endif //SCENE_H
