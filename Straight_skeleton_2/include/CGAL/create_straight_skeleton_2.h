// Copyright (c) 2006-2010 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//

// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>
//
#ifndef CGAL_CREATE_STRAIGHT_SKELETON_2_H
#define CGAL_CREATE_STRAIGHT_SKELETON_2_H

#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_2.h>

namespace CGAL {

namespace CGAL_SS_i
{

template<class Poly>
inline typename Poly::const_iterator vertices_begin ( Poly const& aPoly ) { return aPoly.begin() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_end ( Poly const& aPoly ) { return aPoly.end() ; }


template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_begin ( Polygon_2<K,C> const& aPoly ) 
{ return aPoly.vertices_begin() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_end( Polygon_2<K,C> const& aPoly ) 
{ return aPoly.vertices_end() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_begin ( boost::shared_ptr<Poly> const& aPoly ) { return aPoly->begin() ; }

template<class Poly>
inline typename Poly::const_iterator vertices_end ( boost::shared_ptr<Poly> const& aPoly ) { return aPoly->end() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_begin ( boost::shared_ptr< Polygon_2<K,C> > const& aPoly ) 
{ return aPoly->vertices_begin() ; }

template<class K, class C>
inline typename Polygon_2<K,C>::Vertex_const_iterator vertices_end( boost::shared_ptr< Polygon_2<K,C> > const& aPoly ) 
{ return aPoly->vertices_end() ; }

}

template<class PointIterator, class HoleIterator, class WeightIterator, class WeightSequenceIterator, class NT, class K >
boost::shared_ptr< Straight_skeleton_2<K> >
create_straight_skeleton_2 ( PointIterator              aOuterContour_VerticesBegin
                           , PointIterator              aOuterContour_VerticesEnd
                           , WeightIterator             aOuterContour_WeightsBegin
                           , WeightIterator             aOuterContour_WeightsEnd
                           , HoleIterator               aHolesBegin
                           , HoleIterator               aHolesEnd
                           , WeightSequenceIterator     aHolesWeightBegin
                           , WeightSequenceIterator     aHolesWeightEnd
                           , boost::optional<NT> const& aMaxTime  
                           , K const&                             
                          )
{
  typedef Straight_skeleton_2<K> Ss ;
  typedef boost::shared_ptr<Ss>  SsPtr ;

  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
  
  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;
  
  Cartesian_converter<InputKernel, K> Point_converter ;
  
  SsBuilder ssb(aMaxTime) ;
  
  ssb.enter_contour( aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aOuterContour_WeightsBegin, aOuterContour_WeightsEnd, true, Point_converter ) ;
  
  WeightSequenceIterator whi = aHolesWeightBegin   ;
  for ( HoleIterator hi = aHolesBegin ; hi != aHolesEnd ; ++ hi, ++ whi )
    ssb.enter_contour( CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi), (*whi)->begin(), (*whi)->end(), true, Point_converter ) ;
  
  return ssb.construct_skeleton();
}

template<class PointIterator, class HoleIterator, class WeightIterator, class WeightSequenceIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_straight_skeleton_2 ( PointIterator                  aOuterContour_VerticesBegin
                           , PointIterator                  aOuterContour_VerticesEnd
                           , WeightIterator                 aOuterContour_WeightsBegin
                           , WeightIterator                 aOuterContour_WeightsEnd
                           , HoleIterator                   aHolesBegin
                           , HoleIterator                   aHolesEnd
                           , WeightSequenceIterator         aHolesWeightBegin
                           , WeightSequenceIterator         aHolesWeightEnd
                           , boost::optional<double> const& aMaxTime = boost::optional<double>()
                           )
{
  return create_straight_skeleton_2(aOuterContour_VerticesBegin
                                   ,aOuterContour_VerticesEnd
                                   ,aOuterContour_WeightsBegin
                                   ,aOuterContour_WeightsEnd
                                   ,aHolesBegin
                                   ,aHolesEnd
                                   ,aHolesWeightBegin
                                   ,aHolesWeightEnd
                                   ,aMaxTime
                                   ,Exact_predicates_inexact_constructions_kernel()
                                   );
}


template<class PointIterator, class HoleIterator, class NT, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_straight_skeleton_2 ( PointIterator              aOuterContour_VerticesBegin
                           , PointIterator              aOuterContour_VerticesEnd
                           , HoleIterator               aHolesBegin
                           , HoleIterator               aHolesEnd
                           , NT                         aWeight   
                           , boost::optional<NT> const& aMaxTime  
                           , K const&                             
                           )
{
  typedef Straight_skeleton_2<K> Ss ;
  typedef boost::shared_ptr<Ss>  SsPtr ;

  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
  
  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;
  
  Cartesian_converter<InputKernel, K> Point_converter ;
  
  SsBuilder ssb(aMaxTime) ;
  
  ssb.enter_contour( aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aWeight, true, Point_converter ) ;
  
  for ( HoleIterator hi = aHolesBegin ; hi != aHolesEnd ; ++ hi )
    ssb.enter_contour( CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi), aWeight, true, Point_converter ) ;
  
  return ssb.construct_skeleton();
}

template<class PointIterator, class HoleIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_straight_skeleton_2 ( PointIterator                  aOuterContour_VerticesBegin
                           , PointIterator                  aOuterContour_VerticesEnd
                           , HoleIterator                   aHolesBegin
                           , HoleIterator                   aHolesEnd
                           , double                         aWeight  
                           , boost::optional<double> const& aMaxTime = boost::optional<double>()
                           )
{
  return create_straight_skeleton_2(aOuterContour_VerticesBegin
                                   ,aOuterContour_VerticesEnd
                                   ,aHolesBegin
                                   ,aHolesEnd
                                   ,aWeight
                                   ,aMaxTime
                                   ,Exact_predicates_inexact_constructions_kernel()
                                   );
}

template<class PointIterator, class WeightIterator, class NT, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_straight_skeleton_2 ( PointIterator              aOuterContour_VerticesBegin
                           , PointIterator              aOuterContour_VerticesEnd
                           , WeightIterator             aOuterContour_WeightsBegin
                           , WeightIterator             aOuterContour_WeightsEnd
                           , boost::optional<NT> const& aMaxTime  
                           , K const&                             
                           )
{
  typedef Straight_skeleton_2<K> Ss ;
  typedef boost::shared_ptr<Ss>  SsPtr ;

  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
  
  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;
  
  Cartesian_converter<InputKernel, K> Point_converter ;
  
  SsBuilder ssb(aMaxTime) ;
  
  ssb.enter_contour( aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aOuterContour_WeightsBegin, aOuterContour_WeightsEnd, true, Point_converter ) ;
  
  return ssb.construct_skeleton();
}

template<class PointIterator, class WeightIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_straight_skeleton_2 ( PointIterator                  aOuterContour_VerticesBegin
                           , PointIterator                  aOuterContour_VerticesEnd
                           , WeightIterator                 aOuterContour_WeightsBegin
                           , WeightIterator                 aOuterContour_WeightsEnd
                           , boost::optional<double> const& aMaxTime  = boost::optional<double>()
                           )
{
  return create_straight_skeleton_2(aOuterContour_VerticesBegin
                                   ,aOuterContour_VerticesEnd
                                   ,aOuterContour_WeightsBegin
                                   ,aOuterContour_WeightsEnd
                                   ,aMaxTime
                                   ,Exact_predicates_inexact_constructions_kernel()
                                   );
}

template<class PointIterator, class NT, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_straight_skeleton_2 ( PointIterator              aVerticesBegin
                           , PointIterator              aVerticesEnd
                           , NT                         aWeight   
                           , boost::optional<NT> const& aMaxTime  
                           , K const&                            
                           )
{
  typedef Straight_skeleton_2<K> Ss ;
  typedef boost::shared_ptr<Ss>  SsPtr ;

  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
  
  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;
  
  Cartesian_converter<InputKernel, K> Point_converter ;
  
  SsBuilder ssb(aMaxTime) ;
  
  ssb.enter_contour( aVerticesBegin, aVerticesEnd, aWeight, true, Point_converter ) ;
  
  return ssb.construct_skeleton();
}

template<class PointIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_straight_skeleton_2 ( PointIterator                  aVerticesBegin
                           , PointIterator                  aVerticesEnd
                           , double                         aWeight   = 1.0
                           , boost::optional<double> const& aMaxTime  = boost::optional<double>()
                           )
{
  return create_straight_skeleton_2( aVerticesBegin
                                   , aVerticesEnd
                                   , aWeight
                                   , aMaxTime
                                   , Exact_predicates_inexact_constructions_kernel()
                                   );
}

template<class Polygon, class NT, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_straight_skeleton_2 ( Polygon const&             aContour
                           , NT                         aWeight 
                           , boost::optional<NT> const& aMaxTime
                           , K const&                   aK 
                           )
{
  return create_straight_skeleton_2(CGAL_SS_i::vertices_begin(aContour)
                                   ,CGAL_SS_i::vertices_end  (aContour)
                                   ,aWeight
                                   ,aMaxTime
                                   ,aK
                                   );
}

template<class Polygon>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_straight_skeleton_2 ( Polygon const&                 aContour
                           , double                         aWeight   = 1.0
                           , boost::optional<double> const& aMaxTime  = boost::optional<double>()
                           )
{
  return create_straight_skeleton_2(aContour
                                   ,aWeight
                                   ,aMaxTime
                                   ,Exact_predicates_inexact_constructions_kernel()
                                   );
}


//
//
//


template<class PointIterator, class HoleIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , HoleIterator  aHolesBegin
                                    , HoleIterator  aHolesEnd
                                    , K const&      aK                       
                                    )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aHolesBegin, aHolesEnd, FT(1.0), boost::optional<FT>(), aK ) ;
}


template<class PointIterator, class HoleIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_interior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , HoleIterator  aHolesBegin
                                    , HoleIterator  aHolesEnd
                                    )
{
  return create_interior_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aHolesBegin, aHolesEnd, Exact_predicates_inexact_constructions_kernel() ) ;
}

template<class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , K const&      aK                       
                                    )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, FT(1.0), boost::optional<FT>(), aK ) ;
}

template<class PointIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_interior_straight_skeleton_2 ( PointIterator aVerticesBegin
                                    , PointIterator aVerticesEnd 
                                    )
{ 
  return create_interior_straight_skeleton_2(aVerticesBegin, aVerticesEnd, Exact_predicates_inexact_constructions_kernel()  ) ;
}

template<class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
inline
create_interior_straight_skeleton_2 ( Polygon const& aContour, K const& aK )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aContour, FT(1.0), boost::optional<FT>(), aK );
}

template<class Polygon>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_interior_straight_skeleton_2 ( Polygon const& aContour )
{
  return create_interior_straight_skeleton_2(aContour, Exact_predicates_inexact_constructions_kernel() );
}



//
//
//


template<class PointIterator, class HoleIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_exterior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , HoleIterator  aHolesBegin
                                    , HoleIterator  aHolesEnd
                                    , K const&      aK                       
                                    )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aHolesBegin, aHolesEnd, FT(-1.0), boost::optional<FT>(), aK ) ;
}


template<class PointIterator, class HoleIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
create_exterior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , HoleIterator  aHolesBegin
                                    , HoleIterator  aHolesEnd
                                    )
{
  return create_exterior_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, aHolesBegin, aHolesEnd, Exact_predicates_inexact_constructions_kernel() ) ;
}

template<class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_exterior_straight_skeleton_2 ( PointIterator aOuterContour_VerticesBegin
                                    , PointIterator aOuterContour_VerticesEnd
                                    , K const&      aK                       
                                    )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, FT(-1.0), boost::optional<FT>(), aK ) ;
}

template<class PointIterator>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_exterior_straight_skeleton_2 ( PointIterator aVerticesBegin
                                    , PointIterator aVerticesEnd 
                                    )
{ 
  return create_exterior_straight_skeleton_2(aVerticesBegin, aVerticesEnd, Exact_predicates_inexact_constructions_kernel()  ) ;
}

template<class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
inline
create_exterior_straight_skeleton_2 ( Polygon const& aContour, K const& aK )
{
  typedef typename K::FT FT ;
  return create_straight_skeleton_2(aContour, FT(-1.0), boost::optional<FT>(), aK );
}

template<class Polygon>
boost::shared_ptr< Straight_skeleton_2<Exact_predicates_inexact_constructions_kernel> >
inline
create_exterior_straight_skeleton_2 ( Polygon const& aContour )
{
  return create_exterior_straight_skeleton_2(aContour, Exact_predicates_inexact_constructions_kernel() );
}

} //namespace CGAL


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
