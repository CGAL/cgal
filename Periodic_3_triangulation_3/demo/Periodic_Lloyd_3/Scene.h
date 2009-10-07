#ifndef SCENE_H
#define SCENE_H

#include "typedefs.h"

struct Scene {

  std::list<Point_3> points;
  P3DT3 periodic_triangulation;
  Vertex_handle const_vertex;

  bool eight_copies;
  bool two_dimensional;

  void lloyd_step() {
    Timer timer;
    timer.reset();

    std::vector<Vertex_handle> vts;
    for (Periodic_point_iterator ppit
	   = periodic_triangulation.periodic_points_begin(P3DT3::UNIQUE) ;
	 ppit != periodic_triangulation.periodic_points_end(P3DT3::UNIQUE) ;
	 ++ppit)
      vts.push_back(ppit.get_vertex());
    
    points.clear();
    for (std::vector<Vertex_handle>::iterator vit = vts.begin();
	 vit != vts.end(); ++vit) {
      if (*vit == const_vertex) continue;
      std::vector<Point_3> dual_vertices;
      periodic_triangulation.dual(*vit,std::back_inserter(dual_vertices));
      Point_3 new_point = (two_dimensional ?
	  compute_barycenter_2D(dual_vertices) :
	  compute_barycenter(dual_vertices) );
      dual_vertices.clear();
      points.push_back(new_point);
    }
    periodic_triangulation.clear();
    const_vertex = periodic_triangulation.insert(const_vertex->point());
    periodic_triangulation.insert(points.begin(),points.end());
  }

  Point_3 compute_barycenter(std::vector<Point_3> dual_pts) const {
    FT x(0), y(0), z(0);
    int i;
    for ( i=0 ; i<dual_pts.size() ; i++) {
      x += dual_pts[i].x();
      y += dual_pts[i].y();
      z += dual_pts[i].z();
    }
    x /= i;
    y /= i;
    z /= i;

    x = (x < -1 ? x+2 : (x >= 1 ? x-2 : x));
    y = (y < -1 ? y+2 : (y >= 1 ? y-2 : y));
    z = (z < -1 ? z+2 : (z >= 1 ? z-2 : z));
  
    return Point_3(x,y,z);
  }

  Point_3 compute_barycenter_2D(std::vector<Point_3> dual_pts) const {
    FT x(0), y(0);
    int i;
    for ( i=0 ; i<dual_pts.size() ; i++) {
      x += dual_pts[i].x();
      y += dual_pts[i].y();
    }
    x /= i;
    y /= i;

    x = (x < -1 ? x+2 : (x >= 1 ? x-2 : x));
    y = (y < -1 ? y+2 : (y >= 1 ? y-2 : y));
  
    return Point_3(x,y,0);
  }

};

#endif //SCENE_H
