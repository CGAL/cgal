// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 21
//
// file          : 
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : chull_traits.lw
// revision      : 2.3  
// revision_date : 01 Feb 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CHULL_SUPPORT_3_H
#define CGAL_CHULL_SUPPORT_3_H

#include <CGAL/dd_geo/chull.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

CGAL_BEGIN_NAMESPACE
template <class _HDS>
class Build_polyhedron_from_GRAPH : public Modifier_base< _HDS>
{
  public:
    typedef _HDS                           HDS;
    typedef typename HDS::Vertex           HdsVertex;
    typedef typename HdsVertex::Point      Point;

    Build_polyhedron_from_GRAPH( GRAPH<Point, int>& G ) 
      : gr(G) 
      { CGAL_assertion( G.is_map() ); }

    void operator()( HDS& hds )
         {
           Polyhedron_incremental_builder_3<HDS> B( hds, true);
           leda_node_array<int> index(gr);
           leda_node v;
           CGAL_assertion( gr.is_map() );
           if (gr.number_of_faces() == 0) gr.compute_faces();
           CGAL_assertion( gr.number_of_faces() > 0 );
           B.begin_surface( gr.number_of_nodes(),
                            gr.number_of_faces(),
                            gr.number_of_edges() );
           int i = 0;
#if ( __LEDA__ > 361 )
           forall_nodes(v, gr) 
#else
           for ( v = gr.first_node(); v != nil; v = gr.succ_node( v) )
#endif // __LEDA__
           {
               index[v] = i++;
               B.add_vertex( gr[v] );
           }
           leda_list<leda_face> Llf = gr.all_faces();
           leda_face  lf;
           leda_list<leda_node> Lln;
#if ( __LEDA__ > 361 )
           forall( lf, Llf )
           {
#else
           for ( list_item lli = Llf.first(); lli != nil; lli = Llf.succ( lli) )
           {
               lf = Llf.contents(lli);
#endif // __LEDA__
               B.begin_facet();
               Lln = gr.adj_nodes( lf );
#if ( __LEDA__ > 361 )
               forall( v, Lln)
               {
#else
               for ( list_item lli2 = Lln.first(); lli2 != nil; lli2 = Lln.succ( lli2) )
               {
                   v = Lln.contents(lli2);
#endif // __LEDA__
                   B.add_vertex_to_facet( index[v] );
               }
               B.end_facet();
           } 
           B.end_surface();
         }
           
  private:
    GRAPH<Point, int>& gr;
};

template <class _HDS, class _ChullType>
class Build_polyhedron_from_chull : public Modifier_base< _HDS>
{
  public:
    typedef _HDS                           HDS;
    typedef _ChullType                     ChullType;
    typedef typename HDS::Vertex           HdsVertex;
    typedef typename HdsVertex::Point      Point;
    typedef typename ChullType::ch_facet   Facet;
    typedef typename ChullType::ch_facet_iterator   
                                           FacetIterator;
    typedef typename ChullType::ch_vertex  chVertex;

    Build_polyhedron_from_chull( ChullType& CH) 
      : ch(CH) {}

    void operator()( HDS& hds )
         {
           CGAL_assertion( ch.dcurrent() == 3); 
           Polyhedron_incremental_builder_3<HDS> B( hds, true);
           B.begin_surface( 100, 300);  // would be nice to have statistical data on
                                        // Chull available other than print_statistics()
           leda_map< chVertex, int>  index( -1);
           FacetIterator fit;
           Facet f;
           chVertex v;
           int i = 0;
           for ( fit = ch.facets_begin(); fit != ch.facets_end(); ++fit)
           {
               f = *fit;
               for (int k=0; k < 3; ++k) 
               {
                   v = ch.vertex_of_facet(f,k);
                   if ( index[v] == -1 )
                   {
                       B.add_vertex( ch.associated_point( v));
                       index[v] = i++;
                   }
                }
           }
           Point center = ch.center();
           for ( fit = ch.facets_begin(); fit != ch.facets_end(); ++fit)
           {
               B.begin_facet();
               f = *fit;
               chVertex v0 = ch.vertex_of_facet(f,0);
               chVertex v1 = ch.vertex_of_facet(f,1);
               chVertex v2 = ch.vertex_of_facet(f,2);
               if ( orientation( 
                                            ch.associated_point(v0),
                                            ch.associated_point(v1),
                                            ch.associated_point(v2),
                                            center
                                          ) == POSITIVE )
               {
                   B.add_vertex_to_facet( index[v0] );
                   B.add_vertex_to_facet( index[v1] );
                   B.add_vertex_to_facet( index[v2] );
               }
               else
               {
                   B.add_vertex_to_facet( index[v0] );
                   B.add_vertex_to_facet( index[v2] );
                   B.add_vertex_to_facet( index[v1] );
               }
               B.end_facet();

           }
           B.end_surface();
         }
           
  private:
    ChullType& ch;
};

CGAL_END_NAMESPACE

#endif // CGAL_CHULL_SUPPORT_3_H
