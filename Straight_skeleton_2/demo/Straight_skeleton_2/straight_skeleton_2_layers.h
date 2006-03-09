// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <CGAL/IO/Qt_widget_layer.h>

template <class Sls>
class Qt_layer_show_skeleton : public CGAL::Qt_widget_layer
{
  typedef typename Sls::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename Sls::Vertex_const_handle   Vertex_const_handle ;
  typedef typename Sls::Face_const_iterator   Face_const_iterator ;

public:

  Qt_layer_show_skeleton(Sls const& aSls) : mSls(aSls)
  {};

  void draw()
  {
    typename Sls::Traits::Construct_segment_2 construct_segment ;

    widget->lock();

    int watchdog_limit = mSls.size_of_halfedges();
    
    for ( Face_const_iterator fit = mSls.faces_begin(), efit = mSls.faces_end()
          ; fit != efit
          ; ++ fit
        )
    {
      Halfedge_const_handle hstart = fit->halfedge();
      Halfedge_const_handle he     = hstart ;
      
      int watchdog = watchdog_limit ;
      
      do
      {
        if ( he == null_halfedge )
          break ;

        if ( he->is_bisector() )
        {
          bool lVertexOK      = he->vertex() != null_vertex ;
          bool lOppositeOK    = he->opposite() != null_halfedge ;
          bool lOppVertexOK   = lOppositeOK && he->opposite()->vertex() != null_vertex ;
          bool lVertexHeOK    = lVertexOK && he->vertex()->halfedge() != null_halfedge ;
          bool lOppVertexHeOK = lOppVertexOK && he->opposite()->vertex()->halfedge() != null_halfedge ;

          if ( lVertexOK && lOppVertexOK && lVertexHeOK && lOppVertexHeOK )
          {
            *widget << ( he->is_inner_bisector()? CGAL::BLUE : CGAL::GREEN ) ;
            *widget << construct_segment(he->opposite()->vertex()->point(),he->vertex()->point()) ;
          }
          /*
          else
          {
            int id       = he->id();
            int oppid    = lOppositeOK ? he->opposite()->id() : -1 ;
            int vid      = lVertexOK ? he->vertex()->id() : -1 ;
            int oppvid   = lOppVertexOK ? he->opposite()->vertex()->id() : -1 ;
            int vheid    = lVertexHeOK ? he->vertex()->halfedge()->id() : -1 ;
            int oppvheid = lOppVertexHeOK ? he->opposite()->vertex()->halfedge()->id() : -1 ;

            std::printf("B%d ill-connected:\n"
                        "  Opposite: %d\n"
                        "  Vertex: %d\n"
                        "  Opposite Vertex: %d\n"
                        "  Vertex Halfedge: %d\n"
                        "  Opposite Vertex Halfedge: %d\n"
                        , id, oppid, vid, oppvid, vheid, oppvheid
                        ) ;
          }
          */
        }
        he = he->next();
      }
      while ( -- watchdog > 0 && he != hstart ) ;
    }

    widget->unlock();
  }
private:

  Sls const& mSls;
  const Halfedge_const_handle null_halfedge ;
  const Vertex_const_handle   null_vertex ;

}
;//end class


template <class RegionList>
class Qt_layer_show_regions : public CGAL::Qt_widget_layer
{
  typedef typename RegionList::value_type RegionPtr ;

  typedef typename RegionPtr::element_type Region ;

  typedef typename Region::value_type BoundaryPtr ;

  typedef typename BoundaryPtr::element_type Boundary ;

  typedef typename RegionList::const_iterator const_region_iterator ;

  typedef typename Region::const_iterator const_boundary_iterator ;

  typedef typename Boundary::const_iterator const_vertex_iterator ;

  typename Boundary::Traits::Construct_segment_2 construct_segment ;

public:

  Qt_layer_show_regions(RegionList const& aRegions, CGAL::Color aColor ) : mRegions(aRegions),  mColor(aColor) {};

  void draw()
  {
    widget->lock();

    *widget << mColor;

    for ( const_region_iterator rit = mRegions.begin(), erit = mRegions.end(); rit != erit; ++ rit )
    {
      for ( const_boundary_iterator bit = (*rit)->begin(), ebit = (*rit)->end(); bit != ebit; ++ bit )
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
    }

    widget->unlock();

  }

private:

  RegionList const& mRegions;
  CGAL::Color mColor ;
}
;//end class



