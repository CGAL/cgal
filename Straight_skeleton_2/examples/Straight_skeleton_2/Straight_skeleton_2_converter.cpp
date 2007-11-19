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
#include<CGAL/Exact_predicates_exact_constructions_kernel.h>
#include<CGAL/Straight_skeleton_builder_2.h>
#include<CGAL/Straight_skeleton_2_converter.h>
#include<CGAL/Polygon_offset_builder_2.h>
#include<CGAL/compute_outer_frame_margin.h>

//
// This example illustrates how to use the CGAL Straight Skeleton package
// to construct an offset contour on the outside of a polygon
//

// This is the recommended kernel for the skeleton builder
typedef CGAL::Exact_predicates_inexact_constructions_kernel SsKernel;

// This is a the recommended kernel for the offset builder if
// simple offset polygons are required
typedef CGAL::Exact_predicates_exact_constructions_kernel OfKernel;

typedef SsKernel::Point_2 SsPoint_2;
typedef OfKernel::Point_2 OfPoint_2;

typedef CGAL::Polygon_2<OfKernel>    OfContour;
typedef boost::shared_ptr<OfContour> OfContourPtr;
typedef std::vector<OfContourPtr>    OfContourSequence ;

typedef CGAL::Straight_skeleton_2<SsKernel> Ss;
typedef CGAL::Straight_skeleton_2<OfKernel> Ss2;

typedef Ss::Halfedge_iterator Halfedge_iterator;
typedef Ss::Halfedge_handle   Halfedge_handle;
typedef Ss::Vertex_handle     Vertex_handle;

typedef CGAL::Straight_skeleton_builder_traits_2<SsKernel>    SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

typedef CGAL::Polygon_offset_builder_traits_2<OfKernel>               OfBuilderTraits;
typedef CGAL::Polygon_offset_builder_2<Ss2,OfBuilderTraits,OfContour> OfBuilder;

typedef CGAL::Straight_skeleton_2_items_converter<Ss2,Ss> SsItemsConverter ;

typedef CGAL::Straight_skeleton_2_converter<Ss2,Ss,SsItemsConverter> SsConverter ;

int main()
{
  // A start-shaped polygon, oriented counter-clockwise as required for outer contours.
  SsPoint_2 pts[] = { SsPoint_2(-1,-1)
                    , SsPoint_2(0,-12)
                    , SsPoint_2(1,-1)
                    , SsPoint_2(12,0)
                    , SsPoint_2(1,1)
                    , SsPoint_2(0,12)
                    , SsPoint_2(-1,1)
                    , SsPoint_2(-12,0)
                    } ;

  std::vector<SsPoint_2> star(pts,pts+8);

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
    SsPoint_2 frame[4]= { SsPoint_2(fxmin,fymin)
                        , SsPoint_2(fxmax,fymin)
                        , SsPoint_2(fxmax,fymax)
                        , SsPoint_2(fxmin,fymax)
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
      // Instantiate the container of offset contours
      OfContourSequence offset_contours ;

      SsConverter ConvertSkeleton;

      boost::shared_ptr<Ss2> ss2 = ConvertSkeleton(*ss);

      // Instantiate the offset builder with the skeleton
      OfBuilder ob(*ss2);

      // Obtain the offset contours
      ob.construct_offset_contours(offset, std::back_inserter(offset_contours));

      // Locate the offset contour that corresponds to the frame
      // That must be the outmost offset contour, which in turn must be the one
      // with the largetst unsigned area.
      OfContourSequence::iterator f = offset_contours.end();
      double lLargestArea = 0.0 ;
      for ( OfContourSequence::iterator i = offset_contours.begin(); i != offset_contours.end(); ++ i  )
      {
        //Take abs() as  Polygon_2::area() is signed.
        double lArea = CGAL_NTS to_double( CGAL_NTS abs( (*i)->area() ) ) ; 
        if ( lArea > lLargestArea )
        {
          f = i ;
          lLargestArea = lArea ;
        }
      }

      // Remove the offset contour that corresponds to the frame.
      offset_contours.erase(f);


      // Print out the skeleton
      Halfedge_handle null_halfedge ;
      Vertex_handle   null_vertex ;

      // Dump the edges of the skeleton
      for ( Halfedge_iterator i = ss->halfedges_begin(); i != ss->halfedges_end(); ++i )
      {
        std::string edge_type = (i->is_bisector())? "bisector" : "contour";
        Vertex_handle s = i->opposite()->vertex();
        Vertex_handle t = i->vertex();
        std::cout << "(" << s->point() << ")->(" << t->point() << ") " << edge_type << std::endl;
      }

      // Dump the generated offset polygons

      std::cout << offset_contours.size() << " offset contours obtained\n" ;

      for (OfContourSequence::const_iterator i = offset_contours.begin(); i != offset_contours.end(); ++ i )
      {
        // Each element in the offset_contours sequence is a shared pointer to a Polygon_2 instance.

        std::cout << (*i)->size() << " vertices in offset contour\n" ;

        for (OfContour::Vertex_const_iterator j = (*i)->vertices_begin(); j != (*i)->vertices_end(); ++ j )
        {
          // The coordinates of offset vertices are not double since exact constructions are used
          std::cout << "(" 
                    << CGAL_NTS to_double(j->x()) 
                    << "," 
                    << CGAL_NTS to_double(j->y()) 
                    << ")" 
                    << std::endl ;
        }
      }
    }
  }

  return 0;
}
