#include "Scene.h"
#include "Polyhedron_type.h"
#include "Nef_type.h"

Scene::Bbox Scene::bbox()
{
  if(polyhedra.empty()) {
    Bbox bbox; // default constructor defined
    //bbox.xmin = bbox.ymin = bbox.zmin = 0.0;
    //bbox.xmax = bbox.ymax = bbox.zmax = 1.0;
    return bbox;
  }
  else
  {
    CGAL::Bbox_3 bbox;
    if(polyhedra.begin()->polyhedron_ptr.which() == POLYHEDRON_ENTRY) {
      Polyhedron* poly = boost::get<Polyhedron*>(polyhedra.begin()->polyhedron_ptr);
      const Point& p = *(poly->points_begin());
      bbox = CGAL::Bbox_3(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    }
    else {
      Nef_polyhedron* nef_p = boost::get<Nef_polyhedron*>(polyhedra.begin()->polyhedron_ptr);
      bbox = nef_p->vertices_begin()->point().bbox();
    }
    for(Polyhedra::iterator poly_it = polyhedra.begin();
        poly_it != polyhedra.end();
        ++poly_it)
    {
      if(poly_it->polyhedron_ptr.which() == POLYHEDRON_ENTRY)
      {
	Polyhedron* poly = boost::get<Polyhedron*>(poly_it->polyhedron_ptr);
	for(Polyhedron::Point_iterator it = poly->points_begin();
	    it != poly->points_end();
	    ++it)
	  bbox = bbox + it->bbox();
      }
      else {
	Nef_polyhedron* nef_p = boost::get<Nef_polyhedron*>(poly_it->polyhedron_ptr);
	for(Nef_polyhedron::Vertex_const_iterator it = nef_p->vertices_begin();
	    it != nef_p->vertices_end();
	    ++it)
	  bbox = bbox + it->point().bbox();
      }
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
  }
}
