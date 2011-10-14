#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

template <class LCC>
class Build_triangle123
{
public:
    Build_triangle123(LCC& lcc)
    {
        CGAL::Linear_cell_complex_incremental_builder_3<LCC> B( lcc );
        B.begin_surface( 5, 3, 9);
        typedef typename LCC::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.add_vertex( Point(-1, 0, 0));
        B.add_vertex( Point( 0,-1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 2);
        B.add_vertex_to_facet( 3);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 4);
        B.add_vertex_to_facet( 1);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 3);
        B.add_vertex_to_facet( 4);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 4);
        B.add_vertex_to_facet( 3);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};

template <class HDS>
class Build_triangle123_Polyhedron : public CGAL::Modifier_base<HDS>
{
public:
  Build_triangle123_Polyhedron()
  {};
  void operator()( HDS& p)
    {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( p );
        B.begin_surface( 5, 3, 9);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.add_vertex( Point(-1, 0, 0));
        B.add_vertex( Point( 0,-1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 2);
        B.add_vertex_to_facet( 3);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 4);
        B.add_vertex_to_facet( 1);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 3);
        B.add_vertex_to_facet( 4);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 4);
        B.add_vertex_to_facet( 3);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};

typedef CGAL::Combinatorial_map_with_points<2,3> LCC_3;
typedef LCC_3::Dart_handle                     Dart_handle;
typedef LCC_3::Point                           Point;

typedef CGAL::Polyhedron_3<LCC_3::Traits>  Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;


int main() {
    LCC_3 lcc;
    Build_triangle123<LCC_3> triangle_order123(lcc);

    // Version Linear_cell_complex
    std::cout<<"LCC characteristics: ";
    lcc.display_characteristics(std::cout) << ", valid="
					   << lcc.is_valid()<< std::endl;

    int amark=lcc.get_new_mark();
    for ( LCC_3::Dart_range::iterator it=lcc.darts().begin(),
	    itend=lcc.darts().end(); it!=itend; ++it )
      {
	if ( !lcc.is_marked(it, amark) )
	  {
	    for ( LCC_3::Dart_of_orbit_range<1>::iterator it2=lcc.darts_of_orbit<1>(it).begin(),
		    itend2=lcc.darts_of_orbit<1>(it).end(); it2!=itend2; ++it2 )
	      {
		lcc.mark(it2, amark);
		std::cout<<"("<<LCC_3::point(it2)<<")  ";
	      }
	    std::cout<<std::endl;
	  }
      }
    
		
    // Version Polyhedron_3
    std::cout<<"Polyhedron_3"<<std::endl;
    Polyhedron P;
    Build_triangle123_Polyhedron<HalfedgeDS> triangle;
    P.delegate( triangle);

    for (Polyhedron::Facet_iterator it = P.facets_begin(); it != P.facets_end(); ++it)
   {
      Polyhedron::Halfedge_around_facet_const_circulator it2 = it->facet_begin();
      do
      {
	std::cout<<"("<<it2->vertex()->point()<<")  ";
      }
      while (++it2 != it->facet_begin());
      std::cout<<std::endl;   
   }
    
}
