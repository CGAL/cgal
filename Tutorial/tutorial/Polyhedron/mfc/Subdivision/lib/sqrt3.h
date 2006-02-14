/***************************************************************************
sqrt3.cpp  -  Sqrt3 subdivision (Leif Kobbelt)
c++ wrapper for Lutz Kettner CGAL code
----------------------------------------------------------------------------
begin                : march 2003
copyright            : (C) 2003 by GEOMETRICA / INRIA
email                : pierre.alliez@sophia.inria.fr
/**************************************************************************/

#ifndef SQRT3_H
#define SQRT3_H

#include "Enriched_polyhedron.h"

template <class Polyhedron,class kernel>
class CSubdivider_sqrt3
{
  typedef typename kernel::Point_3 Point;
  typedef typename kernel::Vector_3 Vector;
  typedef typename Polyhedron::Vertex                                   Vertex;
  typedef typename Polyhedron::Vertex_iterator                          Vertex_iterator;
  typedef typename Polyhedron::Halfedge_handle                          Halfedge_handle;
  typedef typename Polyhedron::Edge_iterator                            Edge_iterator;
  typedef typename Polyhedron::Facet_iterator                           Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator         HF_circulator;
  #define PI 3.1415926535897932

  public:
    CSubdivider_sqrt3() {}
    ~CSubdivider_sqrt3() {}   

  public:

    // Subdivision
    int subdivide(Polyhedron &P,
                  int iter)
    {
      // check for valid polygon mesh
      if(P.size_of_facets() == 0)
        return false;

      // normalize border
      P.normalize_border();

      for(int i=0;i<iter;i++)
      {
        // subdivision
        subdivide(P);

        // boundary subdivision   
        if(i & 1)
          subdivide_border(P);
      }
      return true;
    }

    // Flip edge
    void flip_edge(Polyhedron& P,
                   Halfedge_handle e)
    {
      if(e->is_border_edge())
        return;
      Halfedge_handle h = e->next();
      P.join_facet( e);
      P.split_facet( h, h->next()->next());
    }

    // Subdivide
    void subdivide(Polyhedron& P)
    {
      // We use that new vertices/halfedges/facets are appended at the end.
      std::size_t nv = P.size_of_vertices();
      Vertex_iterator last_v = P.vertices_end();
      -- last_v;  // the last of the old vertices
      Edge_iterator last_e = P.edges_end();
      -- last_e;  // the last of the old edges
      Facet_iterator last_f = P.facets_end();
      -- last_f;  // the last of the old facets

      Facet_iterator f = P.facets_begin();    // create new center vertices
      do {
          create_center_vertex( P, f);
      } while ( f++ != last_f);

      std::vector<Point> pts;                    // smooth the old vertices
      pts.reserve( nv);  // get intermediate space for the new points
      ++ last_v; // make it the past-the-end position again
      std::transform( P.vertices_begin(), last_v, std::back_inserter( pts),
                      Smooth_old_vertex());
      std::copy( pts.begin(), pts.end(), P.points_begin());

      Edge_iterator e = P.edges_begin();              // flip the old edges
      ++ last_e; // make it the past-the-end position again
      while ( e != last_e) {
          Halfedge_handle h = e;
          ++e; // careful, incr. before flip since flip destroys current edge
          flip_edge( P, h);
      };
      CGAL_postcondition( P.is_valid());
    }

    // Trisect border halfedge
    void trisect_border_halfedge(Polyhedron& P,
                                 Halfedge_handle e)
    {
        CGAL_precondition( e->is_border());
        // Create two new vertices on e.
        e = e->prev();
        P.split_vertex( e, e->next()->opposite());
        P.split_vertex( e, e->next()->opposite());
        e = e->next();
        // We use later for the smoothing step that e->next()->next()
        // is our original halfedge we started with, i.e., its vertex is
        // from the unrefined mesh.  Split the face twice.
        Halfedge_handle h = e->opposite()->next();
        P.split_facet( e->next()->next()->opposite(), h);
        P.split_facet( e->next()->opposite(), h);
    }

    // Subdivide border
    void subdivide_border(Polyhedron& P)
    {
      if ( P.size_of_facets() == 0)
          return;
      // We use that new halfedges are appended at the end.
      Edge_iterator last_e = P.edges_end();
      -- last_e;  // the last of the old edges
      Edge_iterator e = P.edges_begin(); // create trisected border edges
      do {
          if ( e->opposite()->is_border())
              trisect_border_halfedge( P, e->opposite());
          else if ( e->is_border())
              trisect_border_halfedge( P, e);
      } while ( e++ != last_e);
      e = P.edges_begin();               // smooth points on border edges
      std::vector<Point> pts;  // store new smoothed points temporarily
      do {
          if ( e->opposite()->is_border())
              smooth_border_vertices( e->opposite(), std::back_inserter(pts));
          else if ( e->is_border())
              smooth_border_vertices( e, std::back_inserter(pts));
      } while ( e++ != last_e);
      e = P.edges_begin(); // copy smoothed points back
      std::vector<Point>::iterator i = pts.begin();
      do {
          if ( e->opposite()->is_border()) {
              e->vertex()->point() = *i++;
              e->opposite()->vertex()->point() = *i++;
              e->opposite()->next()->vertex()->point() = *i++;
          } else if ( e->is_border()) {
              e->opposite()->vertex()->point() = *i++;
              e->vertex()->point() = *i++;
              e->next()->vertex()->point() = *i++;
          }
      } while ( e++ != last_e);
      CGAL_assertion( i == pts.end());
      CGAL_postcondition( P.is_valid());
    }

    // Create center vertex
    void create_center_vertex(Polyhedron& P,
                              Facet_iterator f)
    {
      Vector vec( 0.0, 0.0, 0.0);
      std::size_t order = 0;
      HF_circulator h = f->facet_begin();
      do {
          vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
          ++ order;
      } while ( ++h != f->facet_begin());
      CGAL_assertion( order >= 3); // guaranteed by definition of Polyhedron
      Point center =  CGAL::ORIGIN + (vec / (kernel::FT)order);
      Halfedge_handle new_center = P.create_center_vertex( f->halfedge());
      new_center->vertex()->point() = center;
    }
    
    struct Smooth_old_vertex
    {
        Point operator()( const Vertex& v) const
        {
            std::size_t degree = CGAL::circulator_size( v.vertex_begin());
            if ( degree & 1) // odd degree only at border vertices
                return v.point();
            degree = degree / 2;
            kernel::FT alpha = (4.0f - 2.0f * (kernel::FT)cos( 2.0f * PI / (kernel::FT)degree)) / 9.0f;
            Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0f - alpha);
            HV_circulator h = v.vertex_begin();
            do {
                // Even degree and border edges indicated non-manifold
                // vertex with two holes.
                if( h->is_border())
                {
                    std::cerr << "Error: non-manifold vertex. Erroneous smoothing."
                              << std::endl;
                    return v.point();
                }
                vec = vec + ( h->opposite()->vertex()->point() - CGAL::ORIGIN)
                  * alpha / (kernel::FT)degree;
                ++ h;
                CGAL_assertion( h != v.vertex_begin()); // even degree guaranteed
                ++ h;
            } while ( h != v.vertex_begin());
            return (CGAL::ORIGIN + vec);
        }
    };    

    template <class OutputIterator>
    void smooth_border_vertices(typename Halfedge_handle e,
                                OutputIterator out)
    {
        CGAL_precondition( e->is_border());
        // We know that the vertex at this edge is from the unrefined mesh.
        // Get the locus vectors of the unrefined vertices in the neighborhood.
        Vector v0 = e->prev()->prev()->opposite()->vertex()->point() -CGAL::ORIGIN;
        Vector v1 = e->vertex()->point() - CGAL::ORIGIN;
        Vector v2 = e->next()->next()->next()->vertex()->point() - CGAL::ORIGIN;
        *out++ = CGAL::ORIGIN + (10.0 * v0 + 16.0 * v1 +        v2) / 27.0;
        *out++ = CGAL::ORIGIN + ( 4.0 * v0 + 19.0 * v1 +  4.0 * v2) / 27.0;
        *out++ = CGAL::ORIGIN + (       v0 + 16.0 * v1 + 10.0 * v2) / 27.0;
    }

};


#endif // SQRT3_H