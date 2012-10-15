#include<CGAL/Cartesian.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Straight_skeleton_builder_2.h>
#include<CGAL/Polygon_offset_builder_2.h>
#include<CGAL/compute_outer_frame_margin.h>

#include<iostream>
#include<boost/shared_ptr.hpp>


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
  Point_2 pts[] = { Point_2(0,0),
		    Point_2(700,0),
		    Point_2(700,600),
		    Point_2(0,600),
		    Point_2(0,300),
		    Point_2(300,300),
		    Point_2(300,400),
		    Point_2(600,400),
		    Point_2(600,100),
		    Point_2(300,100),
		    Point_2(300,200),
		    Point_2(0,200)
                  } ;

  std::vector<Point_2> star(pts,pts+12);

  double offset = 50; // The offset distance
  boost::optional<double> margin = CGAL::compute_outer_frame_margin(star.begin(),star.end(),50.);

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

    SsBuilder ssb ;

    ssb.enter_contour(frame,frame+4);
    ssb.enter_contour(star.rbegin(),star.rend());

    boost::shared_ptr<Ss> ss = ssb.construct_skeleton();

    if ( ss )
    {
      ContourSequence offset_contours ;
      OffsetBuilder ob(*ss);

      ob.construct_offset_contours(offset, std::back_inserter(offset_contours));
      
      if ( offset_contours.size()!=3 )
      {
        std::cerr << "ERROR: invalid number of componant!\n";
        return 1;
      }
    }
  }

  return 0;
}
