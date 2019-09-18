// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_CONVERTER_2_H
#define CGAL_STRAIGHT_SKELETON_CONVERTER_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/disable_warnings.h>

#include <boost/shared_ptr.hpp>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Straight_skeleton_2.h>

namespace CGAL {

template<class Source_skeleton_
        ,class Target_skeleton_
        ,class NT_converter = typename internal::Default_converter<typename Source_skeleton_::Traits
                                                               ,typename Target_skeleton_::Traits
                                                               >::Type 
        >
struct Straight_skeleton_items_converter_2: Cartesian_converter< typename Source_skeleton_::Traits
                                                               , typename Target_skeleton_::Traits
                                                               , NT_converter
                                                               >
{
  typedef Source_skeleton_ Source_skeleton ;
  typedef Target_skeleton_ Target_skeleton ;
  
  typedef typename Source_skeleton::Traits Source_traits ;
  typedef typename Target_skeleton::Traits Target_traits ;
  
  typedef Cartesian_converter<Source_traits,Target_traits> Base ;
  
  typedef typename Source_skeleton::Vertex_const_handle   Source_vertex_const_handle ;
  typedef typename Source_skeleton::Halfedge_const_handle Source_halfedge_const_handle ;
  typedef typename Source_skeleton::Face_const_handle     Source_face_const_handle  ;
  
  typedef typename Target_skeleton::Vertex    Target_vertex ;
  typedef typename Target_skeleton::Halfedge  Target_halfedge ;
  typedef typename Target_skeleton::Face      Target_face ;
  
  Target_vertex operator() ( Source_vertex_const_handle aV ) const 
  {
    CGAL_assertion( handle_assigned(aV) ) ;
    
    return Target_vertex( aV->id()
                        , this->Base::operator()(aV->point())
                        , this->Base::operator()(aV->time ())
                        , aV->is_split()
                        , aV->has_infinite_time()
                        ) ;
  }
  
  Target_halfedge operator() ( Source_halfedge_const_handle aH ) const 
  {
    CGAL_assertion( handle_assigned(aH) ) ;
    
    return Target_halfedge( aH->id(), aH->slope() ) ;
  }
  
  Target_face operator() ( Source_face_const_handle aF ) const 
  {
    CGAL_assertion( handle_assigned(aF) ) ;
    
    return Target_face( aF->id() );
  }
} ;

template<class Source_skeleton_, class Target_skeleton_, class Items_converter_>
struct Straight_skeleton_converter_2
{
  typedef Source_skeleton_ Source_skeleton ;
  typedef Target_skeleton_ Target_skeleton ;
  typedef Items_converter_ Items_converter ;
  
  typedef typename Source_skeleton::Traits Source_traits ;
  typedef typename Target_skeleton::Traits Target_traits ;
  
  typedef boost::shared_ptr<Target_skeleton> Target_skeleton_ptr ;
  
  typedef typename Source_skeleton::Vertex_const_iterator   Source_vertex_const_iterator ;
  typedef typename Source_skeleton::Halfedge_const_iterator Source_halfedge_const_iterator ;
  typedef typename Source_skeleton::Face_const_iterator     Source_face_const_iterator  ;
  
  typedef typename Source_skeleton::Halfedge_handle Source_halfedge_handle ;
  
  typedef typename Target_skeleton::Vertex   Target_vertex ;
  typedef typename Target_skeleton::Halfedge Target_halfedge ;
  typedef typename Target_skeleton::Face     Target_face ;
  
  typedef typename Target_skeleton::Vertex_handle   Target_vertex_handle ;
  typedef typename Target_skeleton::Halfedge_handle Target_halfedge_handle ;
  typedef typename Target_skeleton::Face_handle     Target_face_handle  ;
  
  typedef typename Target_skeleton::Vertex_iterator   Target_vertex_iterator ;
  typedef typename Target_skeleton::Halfedge_iterator Target_halfedge_iterator ;
  typedef typename Target_skeleton::Face_iterator     Target_face_iterator  ;

  typedef typename Target_skeleton::Base      SlsBase ;
  typedef typename Target_halfedge::Base_base HBase_base ;
  typedef typename Target_halfedge::Base      HBase ;
  typedef typename Target_vertex::Base        VBase ;
  typedef typename Target_face::Base          FBase ;
  
  typedef CGAL_SS_i::Triedge<Source_halfedge_handle> Source_triedge ;
  typedef CGAL_SS_i::Triedge<Target_halfedge_handle> Target_triedge ;
  
  Straight_skeleton_converter_2 ( Items_converter const& aCvt = Items_converter() )
    :
    cvt(aCvt)
  {}
  
  Target_skeleton_ptr operator() ( Source_skeleton const& aSkeleton ) const
  {
    CGAL_assertion(aSkeleton.is_valid());
    Target_skeleton_ptr rResult = create_unconnected_copy(aSkeleton);
    connect_items(aSkeleton,*rResult);
    CGAL_assertion(rResult->is_valid());
    return rResult ;
  }
  
private :
  
  Target_skeleton_ptr create_unconnected_copy ( Source_skeleton const& aSource ) const
  {
    Target_skeleton_ptr rCopy ( new Target_skeleton ) ;
    
    int lMaxVID =-1, lMaxHID = -1, lMaxFID = -1 ;
    
    for ( Source_vertex_const_iterator vit = aSource.vertices_begin(); vit != aSource.vertices_end(); ++ vit )
      if ( vit->id() > lMaxVID )
        lMaxVID = vit->id();
        
    for ( Source_halfedge_const_iterator hit = aSource.halfedges_begin(); hit != aSource.halfedges_end(); ++ hit )
      if ( hit->id() > lMaxHID )
        lMaxHID = hit->id();
        
    for ( Source_face_const_iterator fit = aSource.faces_begin(); fit != aSource.faces_end(); ++ fit )
      if ( fit->id() > lMaxFID )
        lMaxFID = fit->id();
    
    Target_vertices .clear();
    Target_halfedges.clear();
    Target_faces    .clear();
    Target_vertices .resize(lMaxVID+1);
    Target_halfedges.resize(lMaxHID+1);
    Target_faces    .resize(lMaxFID+1);
    
    for ( Source_vertex_const_iterator vit = aSource.vertices_begin(); vit != aSource.vertices_end(); ++ vit )
    {
      Target_vertex_handle tv = rCopy->SlsBase::vertices_push_back( cvt(vit) ) ; 
       
      Target_vertices.at(tv->id()) = tv ;
    }  
      
    for ( Source_halfedge_const_iterator hit = aSource.halfedges_begin(); hit != aSource.halfedges_end(); ++++ hit )
    {
      // In this loop, `hit` is incremented twice, to iterate on edges. We
      // could have used `edges_begin()` and `edges_end()` instead.

      Target_halfedge_handle    th = rCopy->SlsBase::edges_push_back( cvt(hit), cvt(hit->opposite()) ) ; 
      Target_halfedge_handle oppth = th->opposite();
      
      Target_halfedges.at(   th->id()) = th ;
      Target_halfedges.at(oppth->id()) = oppth ;
      
    } 
    
    for ( Source_face_const_iterator fit = aSource.faces_begin(); fit != aSource.faces_end(); ++ fit )
    {
      Target_face_handle tf = rCopy->SlsBase::faces_push_back( cvt(fit) ) ; 
      
      Target_faces.at(tf->id()) = tf ;
    }  
    
    return rCopy ;
  }
  
  void connect_items ( Source_skeleton const& aSource, Target_skeleton& aTarget ) const
  {
    Target_vertex_iterator tvit = aTarget.vertices_begin();
    for ( Source_vertex_const_iterator svit = aSource.vertices_begin(); svit != aSource.vertices_end(); ++ svit, ++ tvit )
    {
      CGAL_assertion( handle_assigned(svit) ) ;
      CGAL_assertion( handle_assigned(svit->halfedge()) ) ;
      
      Target_halfedge_handle tgt_halfedge = Target_halfedges.at(svit->halfedge()->id());
          
      CGAL_assertion( handle_assigned(tgt_halfedge) ) ;
      tvit->VBase::set_halfedge(tgt_halfedge);
      
      Target_halfedge_handle tgt_striedge_e0, tgt_striedge_e1, tgt_striedge_e2 ;
      
      Source_triedge const& stri = svit->event_triedge() ;
      
      if ( handle_assigned(stri.e0()) )
        tgt_striedge_e0 = Target_halfedges.at(stri.e0()->id());
        
      if ( handle_assigned(stri.e1()) )
        tgt_striedge_e1 = Target_halfedges.at(stri.e1()->id());
      
      if ( handle_assigned(stri.e2()) )
        tgt_striedge_e2 = Target_halfedges.at(stri.e2()->id());
        
      tvit->VBase::set_event_triedge( Target_triedge(tgt_striedge_e0, tgt_striedge_e1, tgt_striedge_e2) ) ;
    }
    
    Target_halfedge_iterator thit = aTarget.halfedges_begin();
    for ( Source_halfedge_const_iterator shit = aSource.halfedges_begin(); shit != aSource.halfedges_end(); ++ shit, ++ thit )
    {
      CGAL_assertion( handle_assigned(shit->opposite()) ) ;
      CGAL_assertion( handle_assigned(shit->next    ()) ) ;
      CGAL_assertion( handle_assigned(shit->prev    ()) ) ;
      CGAL_assertion( handle_assigned(shit->vertex  ()) ) ;
      
      Target_halfedge_handle tgt_opposite = Target_halfedges.at(shit->opposite()->id());
      Target_halfedge_handle tgt_next     = Target_halfedges.at(shit->next    ()->id());
      Target_halfedge_handle tgt_prev     = Target_halfedges.at(shit->prev    ()->id());
      Target_vertex_handle   tgt_vertex   = Target_vertices .at(shit->vertex  ()->id());
      
      CGAL_assertion( handle_assigned(tgt_opposite) ) ;
      CGAL_assertion( handle_assigned(tgt_next)     ) ;
      CGAL_assertion( handle_assigned(tgt_prev)     ) ;
      CGAL_assertion( handle_assigned(tgt_vertex)   ) ;
      
      thit->HBase_base::set_opposite (tgt_opposite);
      thit->HBase_base::set_next     (tgt_next);
      thit->HBase_base::set_prev     (tgt_prev);
      thit->HBase_base::set_vertex   (tgt_vertex);
      
      if ( handle_assigned(shit->face()) )
      {
        Target_face_handle tgt_face = Target_faces.at(shit->face()->id());
        CGAL_assertion( handle_assigned(tgt_face) ) ;
        thit->HBase_base::set_face(tgt_face);
      }
     }  
    
    Target_face_iterator tfit = aTarget.faces_begin();
    for ( Source_face_const_iterator sfit = aSource.faces_begin(); sfit != aSource.faces_end(); ++ sfit, ++ tfit )
    {
      CGAL_assertion( handle_assigned(sfit->halfedge()) ) ;
      
      Target_halfedge_handle tgt_halfedge = Target_halfedges.at(sfit->halfedge()->id());
          
      CGAL_assertion( handle_assigned(tgt_halfedge) ) ;
      
      tfit->FBase::set_halfedge(tgt_halfedge);
    }  
  }
  
  Items_converter const& cvt ;
  
  mutable std::vector<Target_vertex_handle>   Target_vertices ;
  mutable std::vector<Target_halfedge_handle> Target_halfedges;
  mutable std::vector<Target_face_handle>     Target_faces    ;
  
} ;

template<class Target_skeleton, class Source_skeleton, class Items_converter>
boost::shared_ptr<Target_skeleton> 
convert_straight_skeleton_2 ( Source_skeleton const& aSrc, Items_converter const& ic )
{
  typedef Straight_skeleton_converter_2<Source_skeleton,Target_skeleton,Items_converter> Skeleton_converter ;
  
  Skeleton_converter c(ic) ;
  
  return c(aSrc);
    
}

template<class Target_skeleton, class Source_skeleton>
boost::shared_ptr<Target_skeleton> 
convert_straight_skeleton_2 ( Source_skeleton const& aSrc )
{
  typedef Straight_skeleton_items_converter_2<Source_skeleton,Target_skeleton> Items_converter ;
  
  typedef Straight_skeleton_converter_2<Source_skeleton,Target_skeleton,Items_converter> Skeleton_converter ;
  
  Skeleton_converter c ;
  
  return c(aSrc);
    
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_STRAIGHT_SKELETON_2_CONVERTER_H //
// EOF //

