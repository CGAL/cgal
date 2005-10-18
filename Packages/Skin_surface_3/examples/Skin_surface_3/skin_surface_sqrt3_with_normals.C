// examples/Skin_surface_3/skin_surface_sqrt3.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Skin_surface_sqrt3_3.h>
#include <CGAL/Skin_surface_refinement_traits_3.h>

#include <list>
#include <fstream>
#include <CGAL/Inverse_index.h>

typedef CGAL::Skin_surface_traits_3<>                    Skin_surface_traits;
typedef Skin_surface_traits::Regular_traits              Regular_traits;
typedef CGAL::Regular_triangulation_3<Regular_traits>    Regular;
typedef Regular_traits::Bare_point                       Reg_point;
typedef Regular_traits::Weighted_point                   Reg_weighted_point;
typedef CGAL::Triangulated_mixed_complex_3<Skin_surface_traits>
                                                     Triangulated_mixed_complex;
typedef Skin_surface_traits::Polyhedron_kernel       Polyhedron_kernel;
typedef CGAL::Polyhedron_3<Polyhedron_kernel>        Polyhedron;

typedef CGAL::Marching_tetrahedra_traits_skin_surface_3<
  Triangulated_mixed_complex,
  Polyhedron,
  Skin_surface_traits::T2P_converter>  Marching_tetrahedra_traits;
typedef CGAL::Marching_tetrahedra_observer_skin_surface_3<
  Triangulated_mixed_complex, Polyhedron>     Marching_tetrahedra_observer;

typedef CGAL::Skin_surface_refinement_traits_3<
          Triangulated_mixed_complex, 
          Polyhedron_kernel, 
          Skin_surface_traits::T2P_converter,
          Skin_surface_traits::P2T_converter>    Skin_surface_refinement_traits;

int main(int argc, char *argv[]) {
  std::list<Reg_weighted_point> l;

  Skin_surface_traits skin_surface_traits(.85);
  
  std::ifstream in("data/caffeine.cin");
  Reg_weighted_point wp;
  while (in >> wp)
    l.push_back(Reg_weighted_point(
		  wp.point(),
		  wp.weight()/skin_surface_traits.shrink_factor()));
    
//   l.push_front(Reg_weighted_point(Reg_point(0,0,0), 1.1));
//   l.push_front(Reg_weighted_point(Reg_point(0,1,0), 2));
//   l.push_front(Reg_weighted_point(Reg_point(0,0,2), 1));

  // Code
  Regular regular;
  Triangulated_mixed_complex triangulated_mixed_complex;
  Polyhedron polyhedron;
  
  // Construct regular triangulation ...
  CGAL::Bbox_3 bbox = (*l.begin()).bbox();
  double max_weight=1;
  for (std::list<Reg_weighted_point>::iterator it= l.begin();
       it != l.end(); it++) {
    max_weight = std::max(max_weight, (*it).weight());
    bbox = bbox + (*it).bbox();
    regular.insert((*it));
  }

  // add a bounding octahedron:
  Reg_point mid((bbox.xmin() + bbox.xmax())/2,
                (bbox.ymin() + bbox.ymax())/2,
                (bbox.zmin() + bbox.zmax())/2);
  double size = 1.5*((bbox.xmax() - bbox.xmin() +
                bbox.ymax() - bbox.ymin() +
	        bbox.zmax() - bbox.zmin())/2 + max_weight);
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x()+size,mid.y(),mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x()-size,mid.y(),mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y()+size,mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y()-size,mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()+size),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()-size),-1));

  // Construct the triangulated mixed complex:
  CGAL::triangulate_mixed_complex_3(
    regular, triangulated_mixed_complex, skin_surface_traits);

  // Extract the coarse mesh using marching_tetrahedra
  Marching_tetrahedra_traits marching_traits;
  Marching_tetrahedra_observer marching_observer;
  
  CGAL::marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits, marching_observer);

  { 
    std::ofstream out("coarse.off");
    out << polyhedron;
  }

  Skin_surface_refinement_traits refinement_traits(triangulated_mixed_complex);
  CGAL::skin_surface_sqrt3(polyhedron, refinement_traits);
  
  { 
    typedef Polyhedron::Vertex_const_iterator                  VCI;
    typedef Polyhedron::Facet_const_iterator                   FCI;
    typedef Polyhedron::Halfedge_around_facet_const_circulator HFCC;

    std::ofstream out("mesh.off");
    out << "NOFF " << polyhedron.size_of_vertices ()
	<< " " << polyhedron.size_of_facets()
	<< " " << polyhedron.size_of_halfedges()
	<< std::endl;

    for (VCI vit = polyhedron.vertices_begin();
	 vit != polyhedron.vertices_end(); vit ++) {
      out << vit->point() << " "
	  << refinement_traits.normal(vit->point())
	  << std::endl;
    }
    CGAL::Inverse_index<VCI> index(
      polyhedron.vertices_begin(),
      polyhedron.vertices_end());
    for(FCI fi = polyhedron.facets_begin();
	fi != polyhedron.facets_end(); ++fi) {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      std::size_t n = circulator_size( hc);
      CGAL_assertion( n >= 3);
      out << n;
      double x,y,z;
      x = y = z = 0;
      do {
	x += hc->vertex()->point().x();
	y += hc->vertex()->point().y();
	z += hc->vertex()->point().z();
	out << " " << index[ VCI(hc->vertex())];
	++hc;
      } while( hc != hc_end);
      int dim = refinement_traits.dimension(
	Polyhedron_kernel::Point_3(x/n,y/n,z/n));
      out << " " << 250-75*dim << " " << 250-75*dim << " " << 250-75*dim;
      
      out << std::endl;
    }
    
  }
  
  return 0;
}
