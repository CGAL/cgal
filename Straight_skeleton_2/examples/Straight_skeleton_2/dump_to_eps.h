#include <fstream>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL
{

template<class Poly> Bbox_2 bbox_2 ( Poly const& aPoly ) { return bbox_2(aPoly.begin(), aPoly.end()); }

template<class K> Bbox_2 bbox_2 ( Polygon_2<K> const& aPoly ) { return bbox_2(aPoly.vertices_begin(), aPoly.vertices_end()); }

template<class K>
Bbox_2 bbox_2 ( Polygon_with_holes_2<K> const& aPolyWH )
{
  Bbox_2 rBbox = bbox_2(aPolyWH.outer_boundary());
  
  for ( typename Polygon_with_holes_2<K>::Hole_const_iterator hit = aPolyWH.holes_begin()
      ; hit != aPolyWH.holes_end() 
      ; ++ hit 
      )
    rBbox = rBbox + bbox_2(*hit);
    
  return rBbox ;
}

}


template<class InputIterator>
void dump_poly_to_eps( InputIterator aBegin, InputIterator aEnd, char const* aType, double aScale, std::ostream& rOut )
{
  InputIterator lLast = aEnd - 1 ;
  
  for( InputIterator curr = aBegin ; curr != aEnd ; ++ curr )
  {
    InputIterator next = curr == lLast ? aBegin : curr + 1 ;
    
    rOut << aType << std::endl
         << aScale * curr->x() 
         << " " 
         << aScale * curr->y()
         << " "
         << aScale * next->x()
         << " "
         << aScale * next->y() 
         << " E\n";
  }
}

template<class K>
void dump_poly_to_eps( CGAL::Polygon_2<K> const& aPoly, char const* aType, double aScale, std::ostream& rOut )
{
  dump_poly_to_eps(aPoly.vertices_begin(), aPoly.vertices_end(), aType, aScale, rOut);
}

template<class K>
void dump_poly_to_eps( CGAL::Polygon_with_holes_2<K> const& aPWH, char const* aType, double aScale, std::ostream& rOut )
{
  dump_poly_to_eps(aPWH.outer_boundary(), aType, aScale, rOut ) ;
      
  for ( typename CGAL::Polygon_with_holes_2<K>::Hole_const_iterator hit = aPWH.holes_begin()
      ; hit != aPWH.holes_end()
      ; ++ hit
      )
    dump_poly_to_eps(*hit, aType, aScale, rOut ) ;
}

template<class K>
void dump_ss_to_eps( CGAL::Straight_skeleton_2<K> const& aSkeleton, char const* aType, double aScale, std::ostream& rOut )
{
  typedef typename CGAL::Straight_skeleton_2<K>::Halfedge_const_iterator Halfedge_const_iterator ;
  typedef typename CGAL::Straight_skeleton_2<K>::Halfedge_const_handle   Halfedge_const_handle ;
  typedef typename CGAL::Straight_skeleton_2<K>::Traits                  Traits ;

  typedef typename Traits::Point_2   Point_2 ;
  typedef typename Traits::Vector_2  Vector_2 ;
  typedef typename Traits::Segment_2 Segment_2 ;
  
  for(Halfedge_const_iterator hit = aSkeleton.halfedges_begin(); hit != aSkeleton.halfedges_end(); ++hit)
  {
    Halfedge_const_handle h = hit ;

    if( h->is_bisector() && ((h->id()%2)==0) )
    { 
      Point_2 s,t;

      if ( !h->has_infinite_time() && !h->opposite()->has_infinite_time() )
      {
        s = h->vertex()->point() ;
        t = h->opposite()->vertex()->point();
      }
      else
      {
        Halfedge_const_handle outh = h->has_infinite_time() ? h : h->opposite();

        Halfedge_const_handle contour_edge_0 = outh            ->defining_contour_edge();
        Halfedge_const_handle contour_edge_1 = outh->opposite()->defining_contour_edge();

        Point_2 const& p0 = contour_edge_0->opposite()->vertex()->point();  
        Point_2 const& p1 = contour_edge_0            ->vertex()->point();  
        Point_2 const& p2 = contour_edge_1            ->vertex()->point();

        Vector_2 bisect = CGAL::ccw_angular_bisector_2(p0, p1, p2, contour_edge_0->weight(), contour_edge_1->weight() );

        s = outh->opposite()->vertex()->point() ;

        t = s - bisect ;
      }

      rOut << aType 
           << std::endl 
           << aScale * CGAL::to_double(s.x()) 
           << " " 
           << aScale * CGAL::to_double(s.y()) 
           << " " 
           << aScale * CGAL::to_double(t.x()) 
           << " " 
           << aScale * CGAL::to_double(t.y()) 
           << " E\n";
    }
  }
}

template<class Poly>
void dump_offset_to_eps( Poly const&                                   aInput
                       , std::vector< boost::shared_ptr<Poly> > const& aOutput
                       , std::ostream&                                 rOut
                       ) 
{
  typedef std::vector< boost::shared_ptr<Poly> > Poly_vector ;
    
  CGAL::Bbox_2 lBbox = CGAL::bbox_2(aInput);
  
  for( typename Poly_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
    lBbox = lBbox + CGAL::bbox_2(**it);
    
  double lWidth  = (lBbox.xmax() - lBbox.xmin()) ;
  double lHeight = (lBbox.ymax() - lBbox.ymin()) ;

  double lSize = (std::max)(lWidth, lHeight);

  double lScale = 1000 / lSize  ;

  double lMargin = lSize * 0.05 ;

  if ( lScale < 1 )
    lScale = 1 ;
    
  CGAL::Bbox_2 lViewport(std::floor(lScale * (lBbox.xmin() - lMargin) )
                        ,std::floor(lScale * (lBbox.ymin() - lMargin) )
                        ,std::ceil (lScale * (lBbox.xmax() + lMargin) )
                        ,std::ceil (lScale * (lBbox.ymax() + lMargin) )
                        );

  rOut << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" 
       << static_cast<int>(lViewport.xmin()) 
       << " " 
       << static_cast<int>(lViewport.ymin()) 
       << " "
       << static_cast<int>(lViewport.xmax()) 
       << " " 
       << static_cast<int>(lViewport.ymax())
       << std::endl;

  rOut << "%%EndComments\n"
          "gsave\n"
          "1.0 setlinewidth\n"
          "/input { 1 0 0 setrgbcolor } bind def\n"
          "/input_w { 0.1 setlinewidth } bind def\n"
          "/output { 0 1 0 setrgbcolor } bind def\n"
          "/output_w { 0.1 setlinewidth } bind def\n"
          "% stroke - x1 y1 x2 y2 E\n"
          "/E {newpath moveto lineto stroke} bind def\n\n"  ;
       
  dump_poly_to_eps(aInput,"input",lScale,rOut);
   
  for( typename Poly_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
    dump_poly_to_eps(**it,"output",lScale,rOut);
   
  rOut << "grestore\nshowpage" << std::endl;
  
}

template<class Input, class K>
void dump_ss_to_eps ( Input                                             const& aInput
                    , boost::shared_ptr< CGAL::Straight_skeleton_2<K> > const& aInSkeleton
                    , boost::shared_ptr< CGAL::Straight_skeleton_2<K> > const& aOutSkeleton
                    , std::ostream&                                             rOut
                    ) 
{
  CGAL::Bbox_2 lBbox = CGAL::bbox_2(aInput);
    
  double lWidth  = (lBbox.xmax() - lBbox.xmin()) ;
  double lHeight = (lBbox.ymax() - lBbox.ymin()) ;

  double lSize = (std::max)(lWidth, lHeight);

  double lScale = 1000 / lSize  ;

  double lMargin = lSize * 0.25 ;

  if ( lScale < 1 )
    lScale = 1 ;
    
  CGAL::Bbox_2 lViewport(std::floor(lScale * (lBbox.xmin() - lMargin) )
                        ,std::floor(lScale * (lBbox.ymin() - lMargin) )
                        ,std::ceil (lScale * (lBbox.xmax() + lMargin) )
                        ,std::ceil (lScale * (lBbox.ymax() + lMargin) )
                        );
    
  rOut << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" 
       << static_cast<int>(lViewport.xmin()) 
       << " " 
       << static_cast<int>(lViewport.ymin()) 
       << " "
       << static_cast<int>(lViewport.xmax()) 
       << " " 
       << static_cast<int>(lViewport.ymax())
       << std::endl;

  rOut << "%%EndComments\n"
          "gsave\n"
          "1.0 setlinewidth\n"
          "/input { 1 0 0 setrgbcolor } bind def\n"
          "/input_w { 0.1 setlinewidth } bind def\n"
          "/in_skeleton { 0 1 0 setrgbcolor } bind def\n"
          "/in_skeleton_w { 0.1 setlinewidth } bind def\n"
          "/out_skeleton { 0 1 0 setrgbcolor } bind def\n"
          "/out_skeleton_w { 0.1 setlinewidth } bind def\n"
          "% stroke - x1 y1 x2 y2 E\n"
          "/E {newpath moveto lineto stroke} bind def\n\n"  ;
       
  dump_poly_to_eps(aInput,"input",lScale,rOut);
   
  if ( aInSkeleton )
    dump_ss_to_eps(*aInSkeleton,"in_skeleton",lScale,rOut);

  if ( aOutSkeleton )
    dump_ss_to_eps(*aOutSkeleton,"out_skeleton",lScale,rOut);
   
  rOut << "grestore\nshowpage" << std::endl;
  
}

template<class Input, class K>
void dump_ss_to_eps ( Input                                              const& aInput
                     , boost::shared_ptr< CGAL::Straight_skeleton_2<K> > const& aInSkeleton
                     , boost::shared_ptr< CGAL::Straight_skeleton_2<K> > const& aOutSkeleton
                    , std::string                                               aFilename
                    ) 
{
  std::ofstream eps(aFilename.c_str()) ;
  if ( eps )  
  {
    std::cerr << "Result: " << aFilename << std::endl ;
    dump_ss_to_eps(aInput,aInSkeleton,aOutSkeleton,eps);
  }
  else
  {
    std::cerr << "Could not open result file: " << aFilename << std::endl ;
  }  
}

template<class Input, class K>
void dump_ss_to_eps ( Input                                             const& aInput
                    , boost::shared_ptr< CGAL::Straight_skeleton_2<K> > const& aSkeleton
                    , std::string                                              aFilename
                    ) 
{
  dump_ss_to_eps(aInput, aSkeleton, boost::shared_ptr< CGAL::Straight_skeleton_2<K> >() , aFilename);
}

template<class Input, class Output>
void dump_offset_to_eps ( Input  const& aInput
                        , Output const& aOutput
                        , std::string   aFilename
                        ) 
{
  std::ofstream eps(aFilename.c_str()) ;
  if ( eps )  
  {
    std::cerr << "Result: " << aFilename << std::endl ;
    dump_offset_to_eps(aInput,aOutput,eps);
  }
  else
  {
    std::cerr << "Could not open result file: " << aFilename << std::endl ;
  }  
}
