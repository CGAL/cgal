// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Polygon_offset_traits_2.c
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_TRAITS_2_C
#define CGAL_POLYGON_OFFSET_TRAITS_2_C 1

CGAL_BEGIN_NAMESPACE

template<class R>
typename Polygon_offset_traits_2<R>::Optional_point_2
Polygon_offset_traits_2<R>::Construct_offset_point_2( FT             aTime
                                                    , Point_2 const& aBorderS
                                                    , Point_2 const& aBorderT
                                                    , Point_2 const& aBisectorS
                                                    , Point_2 const& aBisectorT
                                                    ) const
{  
  Line_2 lBorderL   = construct_line(aBorderS,aBorderT);
  Line_2 lBisectorL = construct_line(aBisectorS,aBisectorT);

  Line_2 lOffsetBorderL = Construct_offset_line(lBorderL,aTime);
  
  return Intersect(lOffsetBorderL,lBisectorL);   
}
  
template<class R>
typename Polygon_offset_traits_2<R>::Line_2
Polygon_offset_traits_2<R>::Construct_offset_line( Line_2 const& aLine, FT aSquaredOffsetDist ) const
{
  typedef typename Rep::Vector_2 Vector_2 ;
  
  typename Rep::Compute_squared_length_2   squared_length ;
    
  Vector_2 lDir = aLine.to_vector();
  
  Vector_2 lN = lDir.perpendicular( LEFT_TURN ) ;
                                  
  Vector_2 lShift = lN * CGAL::sqrt ( aSquaredOffsetDist / squared_length(lN) )  ;
  
  Point_2 lQ = aLine.point(0);
  
  Point_2 lP = lQ + lShift ;
    
  return Line_2(lP,lDir);
}

template<class R>
typename Polygon_offset_traits_2<R>::Optional_point_2
Polygon_offset_traits_2<R>::Intersect( Line_2 const& aA, Line_2 const& aB ) const
{
  typename Rep::Compute_squared_distance_2 squared_distance ;
  
  Optional_point_2 rResult ; 
  Object lIsec = Rep().intersect_2_object()(aA,aB);
  if ( !lIsec.is_empty() )
  {
    Point_2 lP ; 
    if ( assign(lP,lIsec) )
    {
      //
      // Hack for unfiltered floating-point kernels!!!!!!!
      //
      FT const lDiskR = 1e-5 ;
      FT lDevL = squared_distance(aA,lP);
      FT lDevR = squared_distance(aB,lP);
      if ( lDevL <= lDiskR && lDevR <= lDiskR )
        rResult.reset(lP);
    }  
  }
  return rResult ;

}
CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_TRAITS_2_C //
// EOF //
