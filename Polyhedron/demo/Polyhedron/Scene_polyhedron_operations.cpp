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

void Scene::destroy(Polyhedron* poly)
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

CGAL::Bbox_3 
Scene::bbox()
{
  if(polyhedra.empty()) {
    return CGAL::Bbox_3(0.0, 0.0, 0.0, 
			                  1.0, 1.0, 1.0);
  }
  else
  {
    Point p = polyhedra.begin()->polyhedron_ptr->vertices_begin()->point();
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedra::iterator 
          poly_it = polyhedra.begin(),
          poly_end = polyhedra.end();
        poly_it != poly_end; ++poly_it) {
      for(Polyhedron::Vertex_iterator
            v = poly_it->polyhedron_ptr->vertices_begin(),
            v_end = poly_it->polyhedron_ptr->vertices_end();
          v != v_end; ++v)
      {
        bbox = bbox + v->point().bbox();
      }
    }
    return bbox;
  }
}
