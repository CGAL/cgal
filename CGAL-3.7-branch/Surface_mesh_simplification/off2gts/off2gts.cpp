#include <iostream>
#include <fstream>
#include <map>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> Kernel ;

typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> OffSurface; 

typedef OffSurface::Vertex_handle   Vertex_handle ;
typedef OffSurface::Halfedge_handle Halfedge_handle ;
typedef OffSurface::Facet_handle    Facet_handle ;
typedef OffSurface::Vertex_iterator Vertex_iterator ;
typedef OffSurface::Edge_iterator   Edge_iterator ;
typedef OffSurface::Facet_iterator  Facet_iterator ;

void Convert ( OffSurface& off, char const* gts_name )
{
  std::ofstream gts(gts_name);  
  if ( gts )
  {
    std::cout << "Writting " << gts_name << std::endl ;
    
    gts << off.size_of_vertices() << " " << (off.size_of_halfedges()/2) << " " << off.size_of_facets() << std::endl ;
    
    int vid = 1 ;
    for ( Vertex_iterator vit = off.vertices_begin() ; vit != off.vertices_end() ; ++ vit )
    {
      Vertex_handle v = vit ;
      gts << v->point().x() << " " << v->point().y() << " " << v->point().z() << std::endl ;
      v->id() = vid ++ ;
    }
    
    int eid = 1 ;
    for ( Edge_iterator eit = off.edges_begin(); eit != off.edges_end() ; ++ eit )
    {
      Halfedge_handle e = eit ;
      Vertex_handle s = e->opposite()->vertex();
      Vertex_handle t = e->vertex();
      gts << s->id() << " " << t->id() << std::endl ;
      e            ->id() = eid ;
      e->opposite()->id() = eid ;
      ++ eid ;
    }
    
    for ( Facet_iterator fit = off.facets_begin(); fit != off.facets_end() ; ++ fit )
    {
      Facet_handle f = fit ;
      Halfedge_handle e0 = f->halfedge();
      Halfedge_handle e1 = e0->next();
      Halfedge_handle e2 = e1->next();
      gts << e0->id() << " " << e1->id() << " " << e2->id() << std::endl ;
    }
    
    
  }
  else std::cerr << "Unable to open output file: " << gts_name << std::endl ;
}

int main( int argc, char* argv[] )
{
  if ( argc > 2 )
  {
    std::ifstream in(argv[1]);
    if ( in )
    {
      OffSurface off ;
      in >> off ;
      std::cout << "Converting " << argv[1] << " with " << off.size_of_vertices() << " vertices, " << (off.size_of_halfedges()/2) << " edges and " << off.size_of_facets() << " faces" << std::endl ;
      if ( off.is_valid() )
      {
        if ( off.is_pure_triangle() )
        {
          Convert(off,argv[2]);  
        }
        else std::cerr << "Polyhedron is not a purely triangulated surface" << std::endl ;
      }
      else std::cerr << "Invalid polyhedron" << std::endl ;
    }
    else
    {
      std::cerr << "Unable to open " << argv[1] << std::endl ;
    }
  }
  else
  {
    std::cerr << "USAGE: off2gts input.off output.gts" << std::endl ;
  }
  return 0 ;
}
