#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL
{

template<class K, class C>
Bbox_2 bbox_2 ( Polygon_with_holes_2<K,C> const& aPolyWH )
{
  Bbox_2 rBbox = bbox_2(aPolyWH.outer_boundary().vertices_begin(), aPolyWH.outer_boundary().vertices_end());
  
  for ( typename Polygon_with_holes_2<K,C>::Hole_const_iterator hit = aPolyWH.holes_begin()
      ; hit != aPolyWH.holes_end() 
      ; ++ hit 
      )
    rBbox = rBbox + bbox_2(hit->vertices_begin(), hit->vertices_end());
    
  return rBbox ;
}

}

template<class K, class C>
void dump_to_eps( CGAL::Polygon_2<K,C> const& aPoly, char const* aType, double aScale, std::ostream& rOut )
{
  typedef typename CGAL::Polygon_2<K,C>::const_iterator vertex_const_iterator ;
    
  vertex_const_iterator begin = aPoly.vertices_begin() ;
  vertex_const_iterator end   = aPoly.vertices_end  () ;
  vertex_const_iterator last  = end - 1 ;
  
  for( vertex_const_iterator curr = begin ; curr != end ; ++ curr )
  {
    vertex_const_iterator next = curr == last ? begin : curr + 1 ;
    
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

template<class K, class C>
void dump_to_eps( CGAL::Polygon_with_holes_2<K,C> const& aPWH, char const* aType, double aScale, std::ostream& rOut )
{
  dump_to_eps(aPWH.outer_boundary(), aType, aScale, rOut ) ;
      
  for ( typename CGAL::Polygon_with_holes_2<K,C>::Hole_const_iterator hit = aPWH.holes_begin()
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

template<class K, class C>
void dump_to_eps ( CGAL::Polygon_with_holes_2<K,C> const&                                     aInput
                 , std::vector< boost::shared_ptr< CGAL::Polygon_with_holes_2<K,C> > > const& aOutput
                 , std::ostream&                                                            rOut
                 ) 
{
  typedef std::vector< boost::shared_ptr< CGAL::Polygon_with_holes_2<K,C> > > PolyWH_vector ;
    
  CGAL::Bbox_2 lBbox = CGAL::bbox_2(aInput);
  
  for( typename PolyWH_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
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
   
  for( typename PolyWH_vector::const_iterator it = aOutput.begin() ; it != aOutput.end(); ++ it )
    dump_to_eps(**it,"output",lScale,rOut);
   
  rOut << "grestore\nshowpage" << std::endl;
  
}

template<class K, class C>
void dump_to_eps ( CGAL::Polygon_with_holes_2<K,C> const& aInput
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
