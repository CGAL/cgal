#include "Scene.h"
#include "Polyhedron_type.h"
#include "Nef_type.h"

Scene::Bbox Scene::bbox()
{
  if(polyhedra.empty())
  {
    Bbox bbox; // default constructor defined
    return bbox;
  }
  else
  {
    CGAL::Bbox_3 bbox;
    switch(polyhedra.begin()->polyhedron_ptr.which())
    {
      case POLYHEDRON_ENTRY:
      {
	Polyhedron* poly = boost::get<Polyhedron*>(polyhedra.begin()->polyhedron_ptr);
	const Point& p = *(poly->points_begin());
	bbox = CGAL::Bbox_3(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
        break;
      }
      case NEF_ENTRY:
      {
	Nef_polyhedron* nef_p = boost::get<Nef_polyhedron*>(polyhedra.begin()->polyhedron_ptr);
	bbox = nef_p->vertices_begin()->point().bbox();
        break;
      }
      case TEX_POLYHEDRON_ENTRY:
      {
	Textured_polyhedron* poly = boost::get<Textured_polyhedron*>(polyhedra.begin()->polyhedron_ptr);
	const Point& p = *(poly->points_begin());
	bbox = CGAL::Bbox_3(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
      }
    }

    for(Polyhedra::iterator poly_it = polyhedra.begin();
      poly_it != polyhedra.end();
      ++poly_it)
    {
      switch(poly_it->polyhedron_ptr.which())
      {
      case POLYHEDRON_ENTRY:
	{
	  Polyhedron* poly = boost::get<Polyhedron*>(poly_it->polyhedron_ptr);
	  for(Polyhedron::Point_iterator it = poly->points_begin();
	    it != poly->points_end();
	    ++it)
	    bbox = bbox + it->bbox();
	  break;
	}
      case NEF_ENTRY:
	{
	  Nef_polyhedron* nef_p = boost::get<Nef_polyhedron*>(poly_it->polyhedron_ptr);
	  for(Nef_polyhedron::Vertex_const_iterator it = nef_p->vertices_begin();
	    it != nef_p->vertices_end();
	    ++it)
	    bbox = bbox + it->point().bbox();
	  break;
	}
      case TEX_POLYHEDRON_ENTRY:
	{
	  Textured_polyhedron* poly = boost::get<Textured_polyhedron*>(poly_it->polyhedron_ptr);
	  for(Textured_polyhedron::Point_iterator it = poly->points_begin();
	    it != poly->points_end();
	    ++it)
	    bbox = bbox + it->bbox();
	}
      }
    }

    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
      bbox.xmax(),bbox.ymax(),bbox.zmax());
  }
}
