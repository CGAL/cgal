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
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>

#ifndef CGAL_STRAIGHT_SKELETON_2_H
#define CGAL_STRAIGHT_SKELETON_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_items_2.h>
#include <CGAL/HalfedgeDS_default.h>

namespace CGAL {

template<  class Traits_
         , class Items_ = Straight_skeleton_items_2
         , class Alloc_ = CGAL_ALLOCATOR(int)
        >
class Straight_skeleton_2 : public CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_>
{
public :

  typedef Traits_ Traits ;

  typedef Straight_skeleton_2<Traits_,Items_,Alloc_> Self ;
  
  typedef CGAL_HALFEDGEDS_DEFAULT <Traits_,Items_,Alloc_> Base ;
  
  typedef typename Base::Vertex_base     Vertex ;
  typedef typename Base::Halfedge_base   Halfedge ;
  typedef typename Base::Face_base       Face ;
  
  typedef typename Base::Vertex_handle   Vertex_handle ;
  typedef typename Base::Halfedge_handle Halfedge_handle  ;
  typedef typename Base::Face_handle     Face_handle  ;
  
  typedef typename Base::Vertex_const_handle   Vertex_const_handle ;
  typedef typename Base::Halfedge_const_handle Halfedge_const_handle  ;
  typedef typename Base::Face_const_handle     Face_const_handle  ;
  
  typedef typename Base::Vertex_iterator   Vertex_iterator ;
  typedef typename Base::Halfedge_iterator Halfedge_iterator  ;
  typedef typename Base::Face_iterator     Face_iterator  ;

  typedef typename Base::Vertex_const_iterator   Vertex_const_iterator ;
  typedef typename Base::Halfedge_const_iterator Halfedge_const_iterator  ;
  typedef typename Base::Face_const_iterator     Face_const_iterator  ;

  typedef typename Base::size_type size_type ;
    
  Straight_skeleton_2() {}
  
private :

    Vertex_handle   vertices_push_back( const Vertex& v)                   { return Base::vertices_push_back(v); }
    Halfedge_handle edges_push_back( const Halfedge& h, const Halfedge& g) { return Base::edges_push_back(h,g); }
    Halfedge_handle edges_push_back( const Halfedge& h)                    { return Base::edges_push_back(h); }
    Face_handle     faces_push_back( const Face& f)                        { return Base::faces_push_back(f); }
    
    void vertices_pop_front()                                          { Base::vertifces_pop_front(); }
    void vertices_pop_back()                                           { Base::vertifces_pop_back(); }
    void vertices_erase( Vertex_handle v)                              { Base::vertices_erase(v); }
    void vertices_erase( Vertex_iterator first, Vertex_iterator last)  { Base::vertices_erase(first,last); }
    void edges_erase( Halfedge_handle h)                               { Base::edges_erase(h) ; }
    void edges_pop_front()                                             { Base::edges_pop_front(); }
    void edges_pop_back()                                              { Base::edges_pop_back(); }
    void edges_erase( Halfedge_iterator first, Halfedge_iterator last) { Base::edges_erase(first,last); }
    void faces_pop_front()                                             { Base::faces_pop_front(); }
    void faces_pop_back()                                              { Base::faces_pop_back(); }
    void faces_erase( Face_handle f)                                   { Base::faces_erase(f); }
    void faces_erase( Face_iterator first, Face_iterator last)         { Base::faces_erase(first,last); }
    void vertices_clear()                                              { Base::vertices_clear(); }
    void edges_clear()                                                 { Base::edeges_clear(); }
    void faces_clear()                                                 { Base::faces_clear(); }
    void clear()                                                       { Base::clear();}
    
    void vertices_splice( Vertex_iterator target, Self &source, Vertex_iterator begin, Vertex_iterator end) 
      { Base::vertices_splice(target,source,begin,end); }

    void halfedges_splice( Halfedge_iterator target, Self &source, Halfedge_iterator begin, Halfedge_iterator end)
      { Base::halfedges_splice(target,source,begin,end); }

    void faces_splice( Face_iterator target, Self &source,  Face_iterator begin, Face_iterator end)
      { Base::faces_splice(target,source,begin,end); }
    
    void normalize_border() { Base::normalize_border(); }
    
public :

    static int id ( Vertex_const_handle h )
    {
      Vertex_const_handle null ;
      return h != null ? h->id() : -1 ; 
    }
    static int id ( Halfedge_const_handle h )
    {
      Halfedge_const_handle null ;
      return h != null ? h->id() : -1 ; 
    }
    static int id ( Face_const_handle h )
    {
      Face_const_handle null ;
      return h != null ? 0 : -1 ; 
    }

    bool is_valid() const
    {
      //
      // This is a copy of the validity code in Halfedge_const_decorator with a different reporting mechanism
      //
      CGAL_STSKEL_VALIDITY_TRACE("begin Straight_skeleton::is_valid()" );
  
      bool valid = ( 1 != (this->size_of_halfedges() & 1));
      
      CGAL_STSKEL_VALIDITY_TRACE_IF(!valid,"number of halfedges: " << this->size_of_halfedges() << " is odd." ) ;
      
      // All halfedges.
      Halfedge_const_iterator begin = this->halfedges_begin();
      Halfedge_const_iterator end   = this->halfedges_end();
      size_type  n = 0;
      size_type nb = 0;
      for( ; valid && (begin != end); begin++)
      {
          CGAL_STSKEL_VALIDITY_TRACE("he["<< id(begin) << "]" << ( begin->is_border() ?  " [border]" : "" ) );
             
          // Pointer integrity.
          valid = valid && ( begin->next() != Halfedge_const_handle());
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->next() == NULL!");
            break;
          }
          valid = valid && ( begin->opposite() != Halfedge_const_handle());
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->opposite() == NULL!");
            break;
          }
          // opposite integrity.
          valid = valid && ( begin->opposite() != begin);
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->opposite() == he!");
            break;
          }
          valid = valid && ( begin->opposite()->opposite() == begin);
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->opposite()["<< id(begin->opposite())
                                       <<"]->opposite()["<< id(begin->opposite()->opposite()) <<"] != he!"
                                      );
            break;
          }
          // previous integrity.
          valid = valid && begin->next()->prev() == begin;
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<< id(begin) <<"]->next()["<< id(begin->next())
                                      <<"]->prev()["<< id(begin->next()->prev()) <<"] != he."
                                      );
            break;
          }
          // vertex integrity.
          valid = valid && begin->vertex() != Vertex_const_handle();
          if ( ! valid) 
          {
              CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->vertex() == NULL!");
              break;
          }
          if ( ! begin->vertex()->has_infinite_time() )
          {
            valid = valid && ( begin->vertex() == begin->next()->opposite()->vertex());
            if ( ! valid) 
            {
                CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<< id(begin) <<"]->vertex()["<< id(begin->vertex())
                                          <<"] != he->next()["<< id(begin->next())
                                          <<"]->opposite()["<< id(begin->next()->opposite()) 
                                          <<"]->vertex()["<< id(begin->next()->opposite()->vertex())<<"]"
                                          );
                break;
            }
          }
          // face integrity.
          valid = valid && ( begin->is_border() || begin->face() != Face_const_handle() );
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<<id(begin)<<"]->face() == NULL.");
            break;
          }
          valid = valid && ( begin->face() == begin->next()->face());
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: he["<< id(begin) <<"]->face()["<< id(begin->face())
                                      <<"] != he->next()["<< id(begin->next()) <<"]->face()["<< id(begin->next()->face())<<"]."
                                      );
            break;
          }
          ++n;
          if ( begin->is_border())
              ++nb;
      }
      CGAL_STSKEL_VALIDITY_TRACE("summe border halfedges (2*nb) = " << 2 * nb );
      
      bool nvalid = ( n == this->size_of_halfedges());
      
      CGAL_STSKEL_VALIDITY_TRACE_IF(valid && !nvalid
                                   ,"ERROR: counted number of halfedges:" << n 
                                   << " mismatch with this->size_of_halfedges():" << this->size_of_halfedges() 
                                   );
          
      valid = valid && nvalid ;
      
      // All vertices.
      Vertex_const_iterator vbegin = this->vertices_begin();
      Vertex_const_iterator vend   = this->vertices_end();
      
      size_type v = 0;
      n = 0;
      bool is_partial_skeleton = false ;
      
      for( ; valid && (vbegin != vend); ++vbegin) 
      {
          // Pointer integrity.
          valid = valid && vbegin->halfedge() != Halfedge_const_handle()  ;
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: v["<< id(vbegin) <<"]->halfedge() == NULL.");
            break;
          }
          
          // cycle-around-vertex test.
          if ( !vbegin->has_infinite_time() )
          {
            valid = valid && vbegin->halfedge()->vertex() == vbegin;
            if ( ! valid) 
            {
              CGAL_STSKEL_VALIDITY_TRACE("ERROR: v["<< id(vbegin) <<"]->halfedge()["<< id(vbegin->halfedge()) 
                                        <<"]->vertex()["<< id(vbegin->halfedge()->vertex()) <<"] != v."
                                        );
              break;
            }
            
            CGAL_STSKEL_VALIDITY_TRACE("Circulating halfedges around v["<<id(vbegin)<<"]");
            
            Halfedge_const_handle h =  vbegin->halfedge();
            if ( h != Halfedge_const_handle()) 
            {
              Halfedge_const_handle g = h;
              do 
              {
                CGAL_STSKEL_VALIDITY_TRACE("  v->halfedge(): " << id(h) << ", ->next(): " << id(h->next()) 
                                          << ", ->next()->opposite(): " << id(h->next()->opposite())
                                          );
                ++n;
                h = h->next()->opposite();
                valid = valid && ( n <= this->size_of_halfedges() && n!=0);
                CGAL_STSKEL_VALIDITY_TRACE_IF(!valid,"ERROR: more than " << this->size_of_halfedges() 
                                             << " halfedges around v["<< id(vbegin)<<"]"
                                             );
              } while ( valid && (h != g));
            }
          }
          else is_partial_skeleton = true ;
          
          ++v;
      }
      
      if ( ! is_partial_skeleton )
      {
        bool vvalid = (v == this->size_of_vertices());
        
        CGAL_STSKEL_VALIDITY_TRACE_IF(valid && !vvalid
                                     ,"ERROR: counted number of vertices:" << v 
                                     << " mismatch with this->size_of_vertices():" << this->size_of_vertices()
                                     ); 
            
        bool vnvalid = n == this->size_of_halfedges() ;
        CGAL_STSKEL_VALIDITY_TRACE_IF(valid && !vnvalid
                                     ,"ERROR: counted number of halfedges via vertices:" << n 
                                     << " mismatch with this->size_of_halfedges():" << this->size_of_halfedges() 
                                     );
        
        valid = valid && vvalid && vnvalid ;
      }
      
      // All faces.
      Face_const_iterator fbegin = this->faces_begin();
      Face_const_iterator fend   = this->faces_end();
      size_type f = 0;
      n = 0;
      for( ; valid && (fbegin != fend); ++fbegin) 
      {
      
          valid = valid && ( begin->is_border() || fbegin->halfedge() != Halfedge_const_handle()  );
          if ( ! valid)
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: f["<<id(fbegin)<<"]->halfedge() == NULL." );
            break;
          }
          
          valid = valid && fbegin->halfedge()->face() == fbegin ;
          if ( ! valid) 
          {
            CGAL_STSKEL_VALIDITY_TRACE("ERROR: f["<<id(fbegin)<<"]->halfedge()["<< id(fbegin->halfedge()) 
                                       <<"]->face()["<< id(fbegin->halfedge()->face()) <<"] != f."
                                      );
            break;
          }
          // cycle-around-face test.
          CGAL_STSKEL_VALIDITY_TRACE("Circulating halfedges around f["<<id(fbegin)<<"]" );
          Halfedge_const_handle h = fbegin->halfedge();
          if ( h != Halfedge_const_handle()) 
          {
            Halfedge_const_handle g = h;
            do 
            {
              CGAL_STSKEL_VALIDITY_TRACE("  f->halfedge():" << id(h) << ", ->next(): " << id(h->next()));
              ++n;
              h = h->next();
              valid = valid && ( n <= this->size_of_halfedges() && n!=0);
              CGAL_STSKEL_VALIDITY_TRACE_IF(!valid,"ERROR: more than " << this->size_of_halfedges() 
                                           << " halfedges around f["<< id(fbegin)<<"]"
                                           );
            } while ( valid && (h != g));
          }
          ++f;
      }
      
      bool fvalid = ( f == this->size_of_faces());
      
      CGAL_STSKEL_VALIDITY_TRACE_IF(valid && !fvalid
                                   ,"ERROR: counted number of faces:" << f 
                                   << " mismatch with this->size_of_faces():" << this->size_of_faces() 
                                   );
          
      bool fnvalid = ( n + nb  == this->size_of_halfedges() );
                     
      CGAL_STSKEL_VALIDITY_TRACE_IF(valid && !fnvalid
                                   ,"ERROR: counted number of halfedges via faces:" << n
                                   << " plus counted number of border halfedges: " << nb  
                                   << " mismatch with this->size_of_halfedges():" << this->size_of_halfedges() 
                                   );
      
      valid = valid && fvalid && fnvalid ;
      
      CGAL_STSKEL_VALIDITY_TRACE ("end of Straight_skeleton_2>::is_valid(): " << ( valid ? "valid." : "NOT VALID.") );
      
      return valid;
    }    
};


} // end namespace CGAL


#endif // CGAL_STRAIGHT_SKELETON_2_H //
// EOF //

