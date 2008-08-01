#include "Scene.h"
#include "Polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

Polyhedron* Scene::new_polyhedron() 
{
  return new Polyhedron;
}

Polyhedron* Scene::copy_polyhedron(Polyhedron* poly)
{
  return new Polyhedron(*poly);
}

void Scene::destroy_polyhedron(Polyhedron* poly)
{
  delete poly;
}

bool Scene::load_polyhedron(Polyhedron* poly, std::istream& in)
{
  in >> *poly;
  return in;
}

bool Scene::save_polyhedron(Polyhedron* poly, std::ostream& out)
{
  out << *poly;
  return out;
}

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

QString Scene::polyhedronToolTip(int index) const
{
  Polyhedron* poly = polyhedron(index);
  if(!poly)
    return QString();

  return tr("<p><b>%1</b> (mode: %5, color: %6)</p>"
            "<p>Number of vertices: %2<br />"
            "Number of edges: %3<br />"
            "Number of facets: %4</p>")
    .arg(polyhedronName(index))
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(polyhedra[index].rendering_mode == Wireframe ? "wireframe" : "fill")
    .arg(polyhedra[index].color.name());
}
