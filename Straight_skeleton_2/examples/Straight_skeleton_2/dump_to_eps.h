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
void dump_to_eps( InputIterator aBegin, InputIterator aEnd, char const* aType, double aScale, std::ostream& rOut )
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

template<class Poly>
void dump_to_eps( Poly const& aPoly, char const* aType, double aScale, std::ostream& rOut )
{
  dump_to_eps(aPoly.begin(), aPoly.end(), aType, aScale, rOut);
}

template<class K>
void dump_to_eps( CGAL::Polygon_2<K> const& aPoly, char const* aType, double aScale, std::ostream& rOut )
{
  dump_to_eps(aPoly.vertices_begin(), aPoly.vertices_end(), aType, aScale, rOut);
}

template<class K>
void dump_to_eps( CGAL::Polygon_with_holes_2<K> const& aPWH, char const* aType, double aScale, std::ostream& rOut )
{
  dump_to_eps(aPWH.outer_boundary(), aType, aScale, rOut ) ;
      
  for ( typename CGAL::Polygon_with_holes_2<K>::Hole_const_iterator hit = aPWH.holes_begin()
      ; hit != aPWH.holes_end()
      ; ++ hit
      )
    dump_to_eps(*hit, aType, aScale, rOut ) ;
}

template<class K>
void dump_to_eps( CGAL::Straight_skeleton_2<K> const& aSkeleton, char const* aType, double aScale, std::ostream& rOut )
{
  typedef typename CGAL::Straight_skeleton_2<K>::Halfedge_const_iterator Halfedge_const_iterator ;
  typedef typename CGAL::Straight_skeleton_2<K>::Halfedge_const_handle   Halfedge_const_handle ;
  
  for(Halfedge_const_iterator hit = aSkeleton.halfedges_begin(); hit != aSkeleton.halfedges_end(); ++hit)
  {
    Halfedge_const_handle h = hit ;

    if( h->is_bisector() && ((h->id()%2)==0) && !h->has_infinite_time() && !h->opposite()->has_infinite_time() )
    { 
      rOut << aType << std::endl 
           << aScale * h->vertex()->point().x() 
           << " " 
           << aScale * h->vertex()->point().y()
           << " "
           << aScale * h->opposite()->vertex()->point().x()
           << " "
           << aScale * h->opposite()->vertex()->point().y() 
           << " E\n";
    }
  }
}

template<class Poly>
void dump_to_eps ( Poly const&                                   aInput
                 , std::vector< boost::shared_ptr<Poly> > const& aOutput
                 , std::ostream&                                 rOut
                 ) 
{
  typedef std::vector< boost::shared_ptr<Poly> > Poly_vector ;
    
  CGAL::Bbox_2 lBbox = CGAL::bbox_2(aInput);
  
  for( typename Poly_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
    lBbox = lBbox + CGAL::bbox_2(**it);
    
  double lScale = 1000 / (lBbox.xmax() - lBbox.xmin()) ;

  if ( lScale < 1 )
    lScale = 1 ;
    
  rOut << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" 
       << static_cast<int>(std::floor(lScale* lBbox.xmin()-1)) 
       << " " 
       << static_cast<int>(std::floor(lScale* lBbox.ymin()-1)) 
       << " "
       << static_cast<int>(std::ceil(lScale*lBbox.xmax()+1)) 
       << " " 
       << static_cast<int>(std::ceil(lScale*lBbox.ymax()+1))
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
       
  dump_to_eps(aInput,"input",lScale,rOut);
   
  for( typename Poly_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
    dump_to_eps(**it,"output",lScale,rOut);
   
  rOut << "grestore\nshowpage" << std::endl;
  
}

template<class Input, class K>
void dump_to_eps ( Input                         const& aInput
                 , CGAL::Straight_skeleton_2<K>  const& aSkeleton
                 , std::ostream&                        rOut
                 ) 
{
  CGAL::Bbox_2 lBbox = CGAL::bbox_2(aInput);
    
  double lScale = 1000 / (lBbox.xmax() - lBbox.xmin()) ;

  if ( lScale < 1 )
    lScale = 1 ;
    
  rOut << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" 
       << static_cast<int>(std::floor(lScale* lBbox.xmin()-1)) 
       << " " 
       << static_cast<int>(std::floor(lScale* lBbox.ymin()-1)) 
       << " "
       << static_cast<int>(std::ceil(lScale*lBbox.xmax()+1)) 
       << " " 
       << static_cast<int>(std::ceil(lScale*lBbox.ymax()+1))
       << std::endl;

  rOut << "%%EndComments\n"
          "gsave\n"
          "1.0 setlinewidth\n"
          "/input { 1 0 0 setrgbcolor } bind def\n"
          "/input_w { 0.1 setlinewidth } bind def\n"
          "/skeleton { 0 1 0 setrgbcolor } bind def\n"
          "/skeleton_w { 0.1 setlinewidth } bind def\n"
          "% stroke - x1 y1 x2 y2 E\n"
          "/E {newpath moveto lineto stroke} bind def\n\n"  ;
       
  dump_to_eps(aInput,"input",lScale,rOut);
   
  dump_to_eps(aSkeleton,"skeleton",lScale,rOut);
   
  rOut << "grestore\nshowpage" << std::endl;
  
}

template<class Input, class K>
void dump_to_eps ( Input                         const& aInput
                 , CGAL::Straight_skeleton_2<K>  const& aSkeleton
                 , std::string                          aFilename
                 ) 
{
  std::ofstream eps(aFilename.c_str()) ;
  if ( eps )  
  {
    std::cerr << "Result: " << aFilename << std::endl ;
    dump_to_eps(aInput,aSkeleton,eps);
  }
  else
  {
    std::cerr << "Could not open result file: " << aFilename << std::endl ;
  }  
}

template<class Input, class Output>
void dump_to_eps ( Input  const& aInput
                 , Output const& aOutput
                 , std::string   aFilename
                 ) 
{
  std::ofstream eps(aFilename.c_str()) ;
  if ( eps )  
  {
    std::cerr << "Result: " << aFilename << std::endl ;
    dump_to_eps(aInput,aOutput,eps);
  }
  else
  {
    std::cerr << "Could not open result file: " << aFilename << std::endl ;
  }  
}