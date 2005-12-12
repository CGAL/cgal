#include<vector>
#include<iterator>
#include<iostream>
#include<iomanip>
#include<string>

#include<boost/shared_ptr.hpp>

#include<CGAL/basic.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Straight_skeleton_builder_2.h>
#include<CGAL/Polygon_offset_builder_2.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL_HALFEDGEDS_DEFAULT<Kernel,Straight_skeleton_items_2>  Ssds;

typedef Kernel::Point_2 Point_2;

typedef Polygon_2<Kernel>   Contour;
typedef shared_ptr<Contour> ContourPtr;
typedef vector<ContourPtr>  ContourSequence ;

typedef Ssds::Halfedge_iterator Halfedge_iterator;
typedef Ssds::Vertex_handle     Vertex_handle;

typedef Straight_skeleton_builder_traits_2<Kernel>      SsBuilderTraits;
typedef Straight_skeleton_builder_2<BuilderTraits,Ssds> SsBuilder;

typedef Polygon_offset_builder_traits_2<Kernel>                    OffsetBuilderTraits;
typedef Polygon_offset_builder_2<Ssds,OffsetBuilderTraits,Contour> OffsetBuilder;


int main()
{
  Point_2 outer[] = { Point_2(0,0)
                    , Point_2(20,0)
                    , Point_2(20,10)
                    , Point_2(0,10)
                    } ;

  Point_2 hole[] = { Point_2(20,70)
                   , Point_2(50,90)
                   , Point_2(70,50)
                   } ;

  // Instantiate the skeleton builder
  SsBuilder ssb ;

  // Enter the outer contour
  ssb.enter_countour(outer,outer+4);

  // Enter the hole
  ssb.enter_countour(hole,hole+3);

  // Construct the skeleton
  Ssds ss = ssb.construct_skeleton();

  // Dump the edges of the skeleton
  for ( Halfedge_iterator i = ss.halfedges_begin()
      ; i != ss.halfedges_end()
      ; ++i
      )
  {
    string edge_type = (i->is_bisector())? " bisector" : " contour";
    Vertex_handle s = i->source();
    Vertex_handle t = i->target();
    cout << s->point() << "->" << t->point() << edge_type << endl;
  }

  // Instantiate the container of offset polygons
  ContourSequence offset_polygons ;

  // Instantiate the offset polygon builder with the skeleton
  OffsetBuilder ob(ss);

  // Obtain offset polygons at distance 2
  ob.create_offset_polygons(2, std::back_inserter(offset_polygons));

  // Obtain offset polygons at distance 5
  ob.create_offset_polygons(5, std::back_inserter(offset_polygons));

  // Dump the generated offset polygons

  cout << offset_polygons.size() << " offset polygons obtained\n" ;

  for ( typedef ContourSequence::const_iterator i = offset_polygons.begin()
      ; i != offset_polygons.end()
      ; ++ i
      )
  {
    // Each element in the offset_polygons sequence is a shared pointer to a Polygon_2 instance.

    cout << i->size() << " vertices in offset polygon\n" ;

    for ( typedef Contour::const_vertex_iterator j = i->vertices_begin()
        ; j != i->vertices_end()
        ; ++ j
        )
       cout << j->x() << ',' << j->y() << endl ;
  }

  return 0;
}
