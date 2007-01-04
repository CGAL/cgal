#ifndef CGAL_AFSR_VRML_H
#define CGAL_AFSR_VRML_H

namespace CGAL {

template < class Vb, class Fb>
void
afsr_vrml_output(const Triangulation_data_structure_2<Vb,Fb>& tds,
		 std::ostream& os, double r, double g, double b,
		 typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle v, bool skip_infinite)
{
  typedef typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure_2<Vb,Fb>::Vertex_iterator Vertex_iterator;
  typedef typename Triangulation_data_structure_2<Vb,Fb>::Face_iterator Face_iterator;

  // ouput to a vrml file style
  // Point are assumed to be 3d points with a stream operator <<
  // if non NULL, v is the vertex to be output first
  // if skip_inf is true, the point in the first vertex is not output
  // and the faces incident to v are not output
  // (it may be for instance the infinite vertex of the terrain)

  os << "#VRML V2.0 utf8" << std::endl;
  os << "Shape {\n"
     << "appearance Appearance {\n"
     << "material Material { diffuseColor " << r << " " << g << " " << b << "}}\n";
  os << "\tgeometry IndexedFaceSet {" << std::endl;
  os << "\t\tcoord Coordinate {" << std::endl;
  os << "\t\t\tpoint [" << std::endl;

  std::map<Vertex_handle,int> vmap;
  Vertex_iterator vit;
  Face_iterator fit;

  int inum = 0;
  for( vit= tds.vertices_begin(); vit != tds.vertices_end() ; ++vit) {
    if ( v != vit) {
      vmap[vit] = inum++;
      os << "\t\t\t\t" << *vit << ","<< std::endl;
    }
  }

   os << "\t\t\t]" << std::endl;
   os << "\t\t}" << std::endl;
   os << "\t\tsolid FALSE\n"
     "\t\tcoordIndex [" << std::endl;

   // faces
   for(fit= tds.faces_begin(); fit != tds.faces_end(); ++fit) {
     if (!skip_infinite || !fit->has_vertex(v)) {
   	os << "\t\t\t";
	os << vmap[(*fit).vertex(0)] << ", ";
	os << vmap[(*fit).vertex(1)] << ", ";
	os << vmap[(*fit).vertex(2)] << ", ";
	os << "-1, " << std::endl;  
     }
   }
   os << "\t\t]" << std::endl;
   os << "\t}" << std::endl;
   os << "}" << std::endl;
   return;
}


} // namespace CGAL

#endif

