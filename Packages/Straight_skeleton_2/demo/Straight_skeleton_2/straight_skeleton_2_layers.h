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

template <class Sls>
class Qt_layer_show_skeleton : public CGAL::Qt_widget_layer
{
public:

  Qt_layer_show_skeleton(Sls const& aSls) : mSls(aSls)
  {};

  void draw()
  {
    typename Sls::Rep::Construct_segment_2 construct_segment ;

    widget->lock();

    for ( Face_const_iterator fit = mSls.faces_begin(), efit = mSls.faces_end()
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

  Sls const& mSls;
  const Halfedge_const_handle null_halfedge ;
  const Vertex_const_handle   null_vertex ;

}
;//end class


template <class PolygonalRegion>
class Qt_layer_show_polygon : public CGAL::Qt_widget_layer
{
  typedef typename PolygonalRegion::value_type BoundaryPtr ;

  typedef typename BoundaryPtr::element_type Boundary ;

  typedef typename PolygonalRegion::const_iterator const_boundary_iterator ;

  typedef typename Boundary::const_iterator const_vertex_iterator ;

  typename Boundary::Traits::Construct_segment_2 construct_segment ;

public:

  Qt_layer_show_polygon(PolygonalRegion const& aRegion, CGAL::Color aColor ) : mRegion(aRegion),  mColor(aColor) {};

  void draw()
  {
    widget->lock();

    *widget << mColor;

    for ( const_boundary_iterator bit = mRegion.begin(), ebit = mRegion.end(); bit != ebit; ++ bit )
    {
      const_vertex_iterator first = (*bit)->vertices_begin();
      const_vertex_iterator end   = (*bit)->vertices_end  ();
      const_vertex_iterator last  = end - 1 ;
      for ( const_vertex_iterator it = first ; it != end ; ++ it )
      {
        const_vertex_iterator nx = ( it != last ? it + 1 : first ) ;
        *widget << construct_segment(*it,*nx) ;
      }
    }

    widget->unlock();

  }

private:

  PolygonalRegion const& mRegion;
  CGAL::Color mColor ;
}
;//end class



