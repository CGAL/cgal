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


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Straight_skeleton_2<Kernel> Ss;
typedef Kernel::Point_2 Point_2;

typedef CGAL::Polygon_2<Kernel>    Contour;
typedef boost::shared_ptr<Contour> ContourPtr;
typedef std::vector<ContourPtr>    ContourSequence ;

typedef Ss::Halfedge_iterator Halfedge_iterator;
typedef Ss::Halfedge_handle   Halfedge_handle;
typedef Ss::Vertex_handle     Vertex_handle;

typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>      SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<Kernel>                  OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<Ss,OffsetBuilderTraits,Contour> OffsetBuilder;


int main()
{
  Point_2 outer[] = { Point_2(0,0)
                    , Point_2(20,0)
                    , Point_2(20,10)
                    , Point_2(0,10)
                    } ;

  Point_2 hole[] = { Point_2(5,5)
                   , Point_2(10,7)
                   , Point_2(15,5)
                   } ;

  // Instantiate the skeleton builder
  SsBuilder ssb ;

  // Enter the outer contour
  ssb.enter_contour(outer,outer+4);

  // Enter the hole
  ssb.enter_contour(hole,hole+3);

  // Construct the skeleton
  Ss ss = ssb.construct_skeleton();

  Halfedge_handle null_halfedge ;
  Vertex_handle   null_vertex ;

  // Dump the edges of the skeleton
  for ( Halfedge_iterator i = ss.halfedges_begin()
      ; i != ss.halfedges_end()
      ; ++i
      )
  {
    std::string edge_type = (i->is_bisector())? " bisector" : " contour";
    Vertex_handle s = i->opposite()->vertex(); //i->source();
    Vertex_handle t = i->vertex();
    std::cout << s->point() << "->" << t->point() << edge_type << std::endl;
  }

  // Instantiate the container of offset polygons
  ContourSequence offset_contours ;

  // Instantiate the offset polygon builder with the skeleton
  OffsetBuilder ob(ss);

  // Obtain offset contours at distance 2
  ob.construct_offset_contours(2, std::back_inserter(offset_contours));

  // Obtain offset contours at distance 5
  ob.construct_offset_contours(5, std::back_inserter(offset_contours));

  // Dump the generated offset polygons

  std::cout << offset_contours.size() << " offset contours obtained\n" ;

  for (ContourSequence::const_iterator i = offset_contours.begin()
      ; i != offset_contours.end()
      ; ++ i
      )
  {
    // Each element in the offset_contours sequence is a shared pointer to a Polygon_2 instance.

    std::cout << (*i)->size() << " vertices in offset contour\n" ;

    for (Contour::Vertex_const_iterator j = (*i)->vertices_begin()
        ; j != (*i)->vertices_end()
        ; ++ j
        )
       std::cout << j->x() << ',' << j->y() << std::endl ;
  }

  return 0;
}
