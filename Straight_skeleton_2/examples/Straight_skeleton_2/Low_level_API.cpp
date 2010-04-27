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
#include<CGAL/compute_outer_frame_margin.h>

#include "print.h"

//
// This example illustrates how to use the CGAL Straight Skeleton package
// to construct an offset contour on the outside of a polygon
//

// This is the recommended kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel>    Contour;
typedef boost::shared_ptr<Contour> ContourPtr;
typedef std::vector<ContourPtr>    ContourSequence ;

typedef CGAL::Straight_skeleton_2<Kernel> Ss;

typedef Ss::Halfedge_iterator Halfedge_iterator;
typedef Ss::Halfedge_handle   Halfedge_handle;
typedef Ss::Vertex_handle     Vertex_handle;

typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>      SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<Kernel>                  OffsetBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<Ss,OffsetBuilderTraits,Contour> OffsetBuilder;

int main()
{
  // A start-shaped polygon, oriented counter-clockwise as required for outer contours.
  Point_2 pts[] = { Point_2(-1,-1)
                  , Point_2(0,-12)
                  , Point_2(1,-1)
                  , Point_2(12,0)
                  , Point_2(1,1)
                  , Point_2(0,12)
                  , Point_2(-1,1)
                  , Point_2(-12,0)
                  } ;

  std::vector<Point_2> star   (pts,pts+8);

  // Instantiate the skeleton builder
  SsBuilder ssb ;

  // Pass -1.0 as uniform weight to obtain an exterior skeleton
  ssb.enter_contour(star.begin(), star.end(), -1.0 );

  // Construct the skeleton
  boost::shared_ptr<Ss> ss = ssb.construct_skeleton();

  // Proceed only if the skeleton was correctly constructed.
  if ( ss )
  {
    print_straight_skeleton(*ss);
    
    // Instantiate the container of offset contours
    ContourSequence offset_contours ;

    // Instantiate the offset builder with the skeleton
    OffsetBuilder ob(*ss);

    // Obtain the offset contours
    double offset = 3 ; // The offset distance
    ob.construct_offset_contours(offset, std::back_inserter(offset_contours));
    
    print_polygons(offset_contours);
  }

  return 0;
}
