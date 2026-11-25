#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT                         FT ;
typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Ss ;

typedef std::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_2 outer ;
  outer.push_back( Point(-1,-1) ) ;
  outer.push_back( Point(0,-12) ) ;
  outer.push_back( Point(1,-1) ) ;
  outer.push_back( Point(12,0) ) ;
  outer.push_back( Point(1,1) ) ;
  outer.push_back( Point(0,12) ) ;
  outer.push_back( Point(-1,1) ) ;
  outer.push_back( Point(-12,0) ) ;

  Polygon_2 hole ;
  hole.push_back( Point(-1,0) ) ;
  hole.push_back( Point(0,1 ) ) ;
  hole.push_back( Point(1,0 ) ) ;
  hole.push_back( Point(0,-1) ) ;

  assert(outer.is_counterclockwise_oriented());
  assert(hole.is_clockwise_oriented());

  Polygon_with_holes poly( outer ) ;
  poly.add_hole( hole ) ;

  // One range of weights per contour: first the outer boundary, then the holes (in the same order
  // as in the polygon). Within a range, the i-th weight is associated to the contour edge between
  // the (i-1)-th and the i-th vertices (the first weight is thus for the edge from the last vertex
  // to the first vertex). The offset of a contour edge with weight w at time t is at distance t/w
  // from the edge: a larger weight gives a slower front.
  std::vector<std::vector<FT> > weights =
  {
    { 2, 1, 1, 2, 2, 1, 1, 2 }, // the horizontal spikes move slower than the vertical ones
    { 0.5, 0.5, 0.5, 0.5 }      // the hole grows faster
  };

  SsPtr iss = CGAL::create_interior_weighted_straight_skeleton_2(poly, weights);
  if(!iss)
  {
    std::cerr << "Error: failed to create the weighted straight skeleton" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);

  // The weight of a contour edge is stored in the corresponding halfedge of the straight skeleton
  for(Ss::Halfedge_const_iterator h = iss->halfedges_begin(); h != iss->halfedges_end(); ++h)
  {
    if(!h->is_bisector() && !h->is_border())
      std::cout << "Contour edge " << h->opposite()->vertex()->point() << " -> " << h->vertex()->point()
                << " has weight " << h->weight() << std::endl;
  }

  return EXIT_SUCCESS;
}
