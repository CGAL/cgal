#ifndef CGAL_GOCAD_IO_H
#define CGAL_GOCAD_IO_H

#include <iostream>
#include <string>

template <typename Polyhedron>
bool
write_gocad(Polyhedron& polyhedron, std::ostream& os, const std::string& name)
{
  os << "GOCAD TSurf 1\n"
    "HEADER {\n"
    "name:";
  os << name << std::endl;
  os << "*border:on\n"
    "*border*bstone:on\n"
    "}\n"
    "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
    "NAME Default\n"
    "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
    "AXIS_UNIT \"m\" \"m\" \"m\"\n"
    "ZPOSITIVE Elevation\n"
    "END_ORIGINAL_COORDINATE_SYSTEM\n"
    "TFACE\n";

  os.precision(16);
  {
    typename Polyhedron::Vertex_iterator it, end;
    it = polyhedron.vertices_begin();
    end = polyhedron.vertices_end();
    int i = 0;
    for(; it != end; ++it){
      it->id() = i;
      os << "VRTX " << i << " " << it->point() << "\n";
      ++i;
    }
  }

  {
    typename Polyhedron::Facet_iterator it, end;
    it = polyhedron.facets_begin();
    end = polyhedron.facets_end();
    for(; it != end; ++it){
      os << "TRGL " << it->halfedge()->prev()->vertex()->id() << " " 
         << it->halfedge()->vertex()->id() << " "
         << it->halfedge()->next()->vertex()->id()<< "\n";
    }
  }

  os << "END" << std::endl;

  return true;
}

#endif // CGAL_GOCAD_IO_H
