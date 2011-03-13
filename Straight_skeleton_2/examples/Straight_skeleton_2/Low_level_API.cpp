#include<vector>
#include<iterator>
#include<iostream>
#include<iomanip>
#include<string>

#include<boost/shared_ptr.hpp>

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

  std::vector<Point_2> star(pts,pts+8);

  // We want an offset contour in the outside.
  // Since the package doesn't support that operation directly, we use the following trick:
  // (1) Place the polygon as a hole of a big outer frame.
  // (2) Construct the skeleton on the interior of that frame (with the polygon as a hole)
  // (3) Construc the offset contours
  // (4) Identify the offset contour that corresponds to the frame and remove it from the result


  double offset = 3 ; // The offset distance

  // First we need to determine the proper separation between the polygon and the frame.
  // We use this helper function provided in the package.
  boost::optional<double> margin = CGAL::compute_outer_frame_margin(star.begin(),star.end(),offset);

  // Proceed only if the margin was computed (an extremely sharp corner might cause overflow)
  if ( margin )
  {
    // Get the bbox of the polygon
    CGAL::Bbox_2 bbox = CGAL::bbox_2(star.begin(),star.end());

    // Compute the boundaries of the frame
    double fxmin = bbox.xmin() - *margin ;
    double fxmax = bbox.xmax() + *margin ;
    double fymin = bbox.ymin() - *margin ;
    double fymax = bbox.ymax() + *margin ;

    // Create the rectangular frame
    Point_2 frame[4]= { Point_2(fxmin,fymin)
                      , Point_2(fxmax,fymin)
                      , Point_2(fxmax,fymax)
                      , Point_2(fxmin,fymax)
                      } ;

    // Instantiate the skeleton builder
    SsBuilder ssb ;

    // Enter the frame
    ssb.enter_contour(frame,frame+4);

    // Enter the polygon as a hole of the frame (NOTE: as it is a hole we insert it in the opposite orientation)
    ssb.enter_contour(star.rbegin(),star.rend());

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
      ob.construct_offset_contours(offset, std::back_inserter(offset_contours));

      // Locate the offset contour that corresponds to the frame
      // That must be the outmost offset contour, which in turn must be the one
      // with the largetst unsigned area.
      ContourSequence::iterator f = offset_contours.end();
      double lLargestArea = 0.0 ;
      for (ContourSequence::iterator i = offset_contours.begin(); i != offset_contours.end(); ++ i  )
      {
        double lArea = CGAL_NTS abs( (*i)->area() ) ; //Take abs() as  Polygon_2::area() is signed.
        if ( lArea > lLargestArea )
        {
          f = i ;
          lLargestArea = lArea ;
        }
      }

      // Remove the offset contour that corresponds to the frame.
      offset_contours.erase(f);
      
      print_polygons(offset_contours);
    }
  }

  return 0;
}
