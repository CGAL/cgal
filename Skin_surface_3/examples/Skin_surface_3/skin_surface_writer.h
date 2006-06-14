#ifndef CGAL_SKIN_SURFACE_WRITER_H
#define CGAL_SKIN_SURFACE_WRITER_H

#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

template <class Polyhedron, class SkinSurface>
/// Write polyhedron with normals:
void write_polyhedron_with_normals(Polyhedron &p, 
				   SkinSurface &skin, 
				   std::ostream &out)
{
  typedef typename Polyhedron::Vertex_iterator                  Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
  typedef typename Polyhedron::Vertex_handle                    Vertex_handle;
  typedef typename Polyhedron::Traits::Vector_3                 Vector_3;
  
  // Write header
  out << "NOFF " << p.size_of_vertices ()
      << " " << p.size_of_facets()
      << " " << p.size_of_halfedges()
      << std::endl;
  
  // Write vertices
  typedef CGAL::Skin_surface_subdivision_policy_base_3<Polyhedron, SkinSurface> 
    Subdivision_policy;
  Subdivision_policy *policy = get_subdivision_policy(p, skin);
  for (Vertex_iterator vit = p.vertices_begin();
       vit != p.vertices_end(); vit ++) {
    Vector_3 n = policy->normal(vit);
    n = n/sqrt(n*n);
    out << vit->point() << " "
 	<< n
	<< std::endl;
  }

  // Write faces
  CGAL::Inverse_index<Vertex_handle> index(p.vertices_begin(),
					   p.vertices_end());
  for(Facet_iterator fi = p.facets_begin();
      fi != p.facets_end(); ++fi) {
    HFC hc = fi->facet_begin();
    HFC hc_end = hc;
    std::size_t n = circulator_size( hc);
    out << n;
    do {
      Vertex_handle vh = (*hc).vertex();
      out << " " << index[vh];
    } while (++hc != hc_end);
    out << "\n";
  }
}  

template <class Polyhedron, class SkinSurface>
/// Write polyhedron with normals:
void write_union_of_balls_with_normals(Polyhedron &p, 
				       SkinSurface &skin, 
				       std::ostream &out_intersection_curves,
				       std::ostream &out_union_of_balls)
{
  
  write_skin_surface_with_normals(p, skin, out_union_of_balls);
//   typedef typename Polyhedron::Vertex_iterator                  Vertex_iterator;
//   typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
//   typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
//   typedef typename Polyhedron::Vertex_handle                    Vertex_handle;
//   typedef typename Polyhedron::Traits::Vector_3                 Vector_3;

//   // Write header
//   out << "NOFF " << p.size_of_vertices ()
//       << " " << p.size_of_facets()
//       << " " << p.size_of_halfedges()
//       << std::endl;

//   // Write vertices
//   typedef CGAL::Skin_surface_subdivision_policy_base_3<Polyhedron, SkinSurface> 
//     Subdivision_policy;
//   Subdivision_policy *policy = get_subdivision_policy(p, skin);
//   for (Vertex_iterator vit = p.vertices_begin();
//        vit != p.vertices_end(); vit ++) {
//     Vector_3 n = policy->normal(vit);
//     n = n/sqrt(n*n);
//     out << vit->point() << " "
//  	<< n
// 	<< std::endl;
//   }

//   // Write faces
//   CGAL::Inverse_index<Vertex_handle> index(p.vertices_begin(),
// 					   p.vertices_end());
//   for(Facet_iterator fi = p.facets_begin();
//       fi != p.facets_end(); ++fi) {
//     HFC hc = fi->facet_begin();
//     HFC hc_end = hc;
//     std::size_t n = circulator_size( hc);
//     out << n;
//     do {
//       Vertex_handle vh = (*hc).vertex();
//       out << " " << index[vh];
//     } while (++hc != hc_end);
//     out << "\n";
//   }
}  

#endif // CGAL_SKIN_SURFACE_WRITER_H
