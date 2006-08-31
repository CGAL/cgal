#include <iostream>
#include <iomanip>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

// Target surface type. (this include Polyhedron_3.h itself)
#include <CGAL/Surface_mesh_simplification/Polyhedron.h>

// Policies
#include <CGAL/Surface_mesh_simplification/Vertex_is_fixed_map_stored.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

// Simplification method
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

#include <CGAL/IO/Polyhedron_iostream.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel;
typedef Kernel::Vector_3         Vector;
typedef Kernel::Point_3          Point;

//
// Setup an enriched polyhedron type which stores in a halfedge each extra pointer needed by the algorithm
// and in each vertex whether it is fixed or not.
//

template <class Refs, class Traits>
struct My_vertex : public HalfedgeDS_vertex_base<Refs,Tag_true,Point>
{ 
  typedef HalfedgeDS_vertex_base<Refs,Tag_true,Point> Base ;
  
  My_vertex() : is_fixed_(false) {}
  
  My_vertex( Point const& p ) : Base(p), is_fixed_(false) {}
 
  bool is_fixed() const { return is_fixed_ ; }
    
  bool is_fixed_ ;
};

template <class Refs, class Traits>
struct My_halfedge : public HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() : extra_ptr_(0) {}
 
  void*& extra_pointer() { return extra_ptr_ ; }
    
  void* extra_ptr_ ;
};

struct My_items : public Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef My_vertex<Refs,Traits>  Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper { 
        typedef My_halfedge<Refs,Traits>  Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true> Face;
    };
};

typedef Polyhedron_3<Kernel,My_items> Surface; 

typedef Surface::Halfedge_handle Halfedge_handle ;
typedef Surface::Vertex_handle   Vertex_handle ;

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

struct Visitor
{
  Visitor() : collected(0), processed(0), collapsed(0), non_collapsable(0), cost_uncomputable(0) {} 
  
  void OnStarted( Surface& aSurface ) {} 
  
  void OnFinished ( Surface& aSurface ) { cerr << "\n" << flush ; } 
  
  void OnStopConditionReached( Surface& aSurface ) {} 
  
  void OnCollected( Halfedge_handle const& aEdge
                  , bool                   aIsFixed
                  , Surface&            aSurface
                  )
  {
    ++ collected ;
    cerr << "\rEdges collected: " << collected << flush ;
 }                
  
  void OnProcessed(Halfedge_handle const& aEdge
                  ,Surface&            aSurface
                  ,optional<double>       aCost
                  ,Vertex_handle const&   aVertex
                  )
  {
    ++ processed ;
    if ( aVertex == Vertex_handle() )
    {
      if ( !aCost )
           ++ cost_uncomputable ;
      else ++ non_collapsable ;
    }
    else 
    {
      ++ collapsed;
    }
  }                
  
  void OnStep(Halfedge_handle const& aEdge, Surface& aSurface, size_t aInitial, size_t aCurrent)
  {
    if ( aCurrent == aInitial )
      cerr << "\n" << flush ;
    cerr << "\r" << aCurrent << flush ;
  }                
  
  size_t collected, processed, collapsed, non_collapsable, cost_uncomputable ; 
} ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  //
  // Simplify surface
  //  Down to a 10% of the number of edges.
  //  Per-edge extra pointer embeeded in each halfedge.
  //  "is-fixed" flag embeeded in each vertex.
  //  Caching cost with not placement (partial collapse data)
  //  Using LindstromTurk cost strategy from the cached collapse data.
  //  Using LindstromTurk placement strategy computed on demand.
  //  Using default LindstromTurk params
  //  Using a visitor to track progress and collect stats.
  //

  Visitor visitor ;
  
  LindstromTurk_params params ;
    
  int r = edge_collapse(surface
                       ,Count_ratio_stop_condition<Surface>(0.10)
                       ,Edge_extra_pointer_map_stored<Surface>()
                       ,Vertex_is_fixed_map_stored<Surface>()
                       ,Set_partial_collapse_data_LindstromTurk<Surface>()
                       ,Cached_cost<Surface>()
                       ,LindstromTurk_placement<Surface>()
                       ,&params
                       ,&params
                       ,&visitor
                       );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges." ;
  
  cout << "\nEdges collected: " << visitor.collected
       << "\nEdges proccessed: " << visitor.processed
       << "\nEdges collapsed: " << visitor.collapsed
       << endl
       << "\nEdges not collapsed due to topological constrians: " << visitor.non_collapsable
       << "\nEdge not collapsed due to computational constrians: " << visitor.cost_uncomputable 
       << endl ; 
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
