// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#include <CGAL/IO/Qt_widget_layer.h>

template <class Ssds>
class Qt_layer_show_skeleton : public CGAL::Qt_widget_layer
{
public:

  Qt_layer_show_skeleton(Ssds const& aSsds) : mSsds(aSsds)
  {};
  
  void draw()
  {
    typename Ssds::Rep::Construct_segment_2 construct_segment ;
    
    widget->lock();
    
    for ( Face_const_iterator fit = mSsds.faces_begin(), efit = mSsds.faces_end()
          ; fit != efit
          ; ++ fit
        )
    {
      Halfedge_const_handle hstart = fit->halfedge();
      Halfedge_const_handle he     = hstart ;
      do
      {
        if ( he == null_halfedge )
          break ;
        if ( he->is_bisector() )
        {  
          bool lOK =    he->vertex() != null_vertex 
	             && he->opposite() != null_halfedge 
		     && he->opposite()->vertex() != null_vertex ;
		     
          *widget << ( lOK ? (he->is_contour_bisector()? CGAL::GREEN : CGAL::BLUE ) : CGAL::YELLOW ) ;
          *widget << construct_segment(he->opposite()->vertex()->point(),he->vertex()->point()) ;
        }
	he = he->next();
      }
      while ( he != hstart ) ;
    }
    
    widget->unlock();
  }
private:

  Ssds const& mSsds;
  const Halfedge_const_handle null_halfedge ;
  const Vertex_const_handle   null_vertex ;

}
;//end class


template <class PolygonalRegion>
class Qt_layer_show_polygon : public CGAL::Qt_widget_layer
{
public:

  Qt_layer_show_polygon(PolygonalRegion const& aRegion, CGAL::Color aColor ) : mRegion(aRegion),  mColor(aColor) {};
  
  void draw()
  {
    widget->lock();
    
    *widget << mColor;

    for ( typename PolygonalRegion::const_iterator bit = mRegion.begin(), ebit = mRegion.end()
        ; bit != ebit 
	; ++ bit 
	)
    {
      Polygon::const_iterator first = (*bit)->vertices_begin();
      Polygon::const_iterator end   = (*bit)->vertices_end  ();
      Polygon::const_iterator last  = end - 1 ;
      for ( Polygon::const_iterator it = first ; it != end ; ++ it )
      {
        Polygon::const_iterator nx = ( it != last ? it + 1 : first ) ;
        *widget << Segment(*it,*nx) ;  
      }  
    }
    
    widget->unlock();

  }

private:

  PolygonalRegion const& mRegion;
  CGAL::Color mColor ;
}
;//end class



