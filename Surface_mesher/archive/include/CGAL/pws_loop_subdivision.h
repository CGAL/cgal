// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez

#ifndef PWS_LOOP_SUBDIVISION
#define PWS_LOOP_SUBDIVISION

#include "enriched_polyhedron.h"
#include "builder.h"

template <class HDS,class Polyhedron,class Kernel>
class CModifierLoop : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename Kernel::FT           FT;
  typedef typename Kernel::Point_3      Point;
  typedef typename Kernel::Vector_3     Vector;
  typedef typename HDS::Face_handle     Face_handle;
  typedef typename HDS::Vertex_handle   Vertex_handle;
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> Builder;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator HV_circulator;
  typedef typename Polyhedron::Vertex_type Vertex_type;
  Polyhedron *m_pMesh;
  bool m_piecewise_smooth;

public:
  // life cycle
  CModifierLoop(Polyhedron *pMesh,
                const bool piecewise_smooth = true)
  {
    CGAL_assertion(pMesh != NULL);
    m_pMesh = pMesh;
    m_piecewise_smooth = piecewise_smooth;
  }
  ~CModifierLoop() {}

  // subdivision
  void operator()( HDS& hds)
  {
    Builder builder(hds,true);
    builder.begin_surface(3,1,6);
      add_vertices(builder);
      add_facets(builder);
    builder.end_surface();
  }

  // add vertices
  void add_vertices(Builder &builder)
  {
    // put original vertices
    int index = 0;
    typename Polyhedron::Vertex_iterator v;
    for(v = m_pMesh->vertices_begin();
        v != m_pMesh->vertices_end();
        v++)
    {
      v->tag(index);

      // border vertices
      if(Polyhedron::is_border(v))
      {
        const Point& curr = v->point();
        const typename Polyhedron::Halfedge_handle& he = m_pMesh->get_border_halfedge(v);
        CGAL_assertion(he != NULL);
        const Point& next = he->next()->vertex()->point();
        const Point& prev = he->prev()->vertex()->point();
        Point p = CGAL::ORIGIN + 0.25 *(prev - CGAL::ORIGIN) +
                                 0.50 *(curr - CGAL::ORIGIN) +
                                 0.25 *(next - CGAL::ORIGIN);
        builder.add_vertex(p);
      } // end is border
      else
      {
        Vertex_type vertex_type = m_pMesh->type(v);

        switch(vertex_type)
        {
          // corner - do not move
          case Polyhedron::CORNER:
            builder.add_vertex(v->point());
            break;

          // crease type (regular or irregular)
          case Polyhedron::CREASE_REGULAR:
          case Polyhedron::CREASE_IRREGULAR:
          {
            const Point& curr = v->point();
            // get two other points along crease
            // prev ==== curr ==== next
            Point prev,next;
            m_pMesh->incident_points_on_crease(v,prev,next);
            Point p = CGAL::ORIGIN + 0.125 *(prev - CGAL::ORIGIN) +
                                     0.750 *(curr - CGAL::ORIGIN) +
                                     0.125 *(next - CGAL::ORIGIN);
            builder.add_vertex(p);
            break;
          }

          // smooth or dart
          case Polyhedron::SMOOTH:
          case Polyhedron::DART:
          {
            unsigned int n = Polyhedron::valence(v);
            FT w = alpha((FT)n);
            FT sum_weights = w + (FT)n;
            // new position
            Vector pos = w/sum_weights * (v->point() - CGAL::ORIGIN);
            HV_circulator he = v->vertex_begin();
            HV_circulator end = he;
            CGAL_For_all(he,end)
            {
              // add edge-vertex contribution
              const Point& point = he->opposite()->vertex()->point();
              pos = pos + (point - CGAL::ORIGIN) / sum_weights;
            }
            builder.add_vertex(CGAL::ORIGIN + pos);
          } // end smooth or dart
        } // end switch
      } // end !is border
      index++;
    }

    // as many as #edges
    m_pMesh->tag_halfedges(-1);
    typename Polyhedron::Halfedge_iterator he;
    for(he = m_pMesh->halfedges_begin();
        he != m_pMesh->halfedges_end();
        he++)
    {
      if(he->tag() != -1)
        continue;

      if(he->is_border() || he->opposite()->is_border())
      {
        const Point& p1 = he->vertex()->point();
        const Point& p2 = he->opposite()->vertex()->point();
        Point point = CGAL::midpoint(p1, p2);
        builder.add_vertex(point);
      }
      else
      {
        int mask_type = 1; //  smooth by default
        Vertex_type type_v1;
        Vertex_type type_v2;

        // if sharp, get type of end vertices to deduce mask type
        if(he->sharp())
        {
          Vertex_handle v1 = he->vertex();
          Vertex_handle v2 = he->opposite()->vertex();
          type_v1 = m_pMesh->type(v1);
          type_v2 = m_pMesh->type(v2);

          if(type_v1 == Polyhedron::DART || type_v2 == Polyhedron::DART)
            mask_type = 1; // smooth
          else
          {
            if(type_v1 == type_v2) // all diagonals of minor in table 1
              mask_type = 2; // regular crease case
            else
              if(type_v1 == Polyhedron::CREASE_REGULAR ||
                 type_v2 == Polyhedron::CREASE_REGULAR)
                mask_type = 3; // irregular crease case
              else
              {
                mask_type = 2; // regular crease case
                // this is all must remain
                CGAL_assertion( (type_v1 == Polyhedron::CREASE_IRREGULAR && type_v2 == Polyhedron::CORNER) ||
                        (type_v2 == Polyhedron::CREASE_IRREGULAR && type_v1 == Polyhedron::CORNER));
              }
          }
        }

        const Point& p1 = he->vertex()->point();
        const Point& p2 = he->opposite()->vertex()->point();

        switch(mask_type)
        {
          case 1:
          {
            const Point& p3 = he->next()->vertex()->point();
            const Point& p4 = he->opposite()->next()->vertex()->point();
            Point p = CGAL::ORIGIN + 0.375 *(p1 - CGAL::ORIGIN) +
                                     0.375 *(p2 - CGAL::ORIGIN) +
                                     0.125 *(p3 - CGAL::ORIGIN) +
                                     0.125 *(p4 - CGAL::ORIGIN);
            builder.add_vertex(p);
            break;
          }
          case 2: // midpoint
          {
            builder.add_vertex(CGAL::midpoint(p1,p2));
            break;
          }
          case 3: // 5/3 (regular vertex gets weight 5)
          {
            // 0.625
            Point regular,irregular;
            if(type_v1 == Polyhedron::CREASE_REGULAR)
            {
              regular = p1;
              irregular = p2;
              CGAL_assertion(type_v2 == Polyhedron::CREASE_IRREGULAR);
            }
            else
            {
              regular = p2;
              irregular = p1;
              CGAL_assertion(type_v1 == Polyhedron::CREASE_IRREGULAR);
              CGAL_assertion(type_v2 == Polyhedron::CREASE_REGULAR);
            }
            Point p = CGAL::ORIGIN + 0.625 *(regular   - CGAL::ORIGIN) +
                                     0.375 *(irregular - CGAL::ORIGIN);
            builder.add_vertex(p);
          }
        }
      }

      // put vertex indices on both halfedges
      he->tag(index);
      he->opposite()->tag(index);
      index++;
    }
  }

  // add facets
  void add_facets(Builder &builder)
  {
    typename Polyhedron::Facet_iterator f;
    for(f = m_pMesh->facets_begin();
        f != m_pMesh->facets_end();
        f++)
    {
      unsigned int degree = Polyhedron::degree(f);
      CGAL_assertion(degree == 3);

      typename Polyhedron::Halfedge_handle he = f->halfedge();
      int i0 = he->tag();
      int i1 = he->vertex()->tag();
      int i2 = he->next()->tag();
      int i3 = he->next()->vertex()->tag();
      int i4 = he->next()->next()->tag();
      int i5 = he->next()->next()->vertex()->tag();

      // create 4 triangles
      builder.begin_facet();
        builder.add_vertex_to_facet(i1);
        builder.add_vertex_to_facet(i2);
        builder.add_vertex_to_facet(i0);
      const Halfedge_handle& h1 = builder.end_facet();
      builder.begin_facet();
        builder.add_vertex_to_facet(i3);
        builder.add_vertex_to_facet(i4);
        builder.add_vertex_to_facet(i2);
      const Halfedge_handle& h2 = builder.end_facet();
      builder.begin_facet();
        builder.add_vertex_to_facet(i5);
        builder.add_vertex_to_facet(i0);
        builder.add_vertex_to_facet(i4);
      const Halfedge_handle& h3 = builder.end_facet();

      // center face
      builder.begin_facet();
        builder.add_vertex_to_facet(i0);
        builder.add_vertex_to_facet(i2);
        builder.add_vertex_to_facet(i4);
      Halfedge_handle h4 = builder.end_facet();

      h1->next()->next()->control_edge() = false;
      h2->next()->next()->control_edge() = false;
      h3->next()->next()->control_edge() = false;
      h4->control_edge() = false;
      h4->next()->control_edge() = false;
      h4->next()->next()->control_edge() = false;

      h1->next()->next()->sharp() = false;
      h2->next()->next()->sharp() = false;
      h3->next()->next()->sharp() = false;
      h4->sharp() = false;
      h4->next()->sharp() = false;
      h4->next()->next()->sharp() = false;

      // propagate control edge tags
      h1->control_edge() = he->control_edge();
      h3->next()->control_edge() = he->control_edge();

      h1->next()->control_edge() = he->next()->control_edge();
      h2->control_edge() = he->next()->control_edge();

      h2->next()->control_edge() = he->prev()->control_edge();
      h3->control_edge() = he->prev()->control_edge();

      // propagate sharp tags
      h1->sharp() = he->sharp();
      h3->next()->sharp() = he->sharp();

      h1->next()->sharp() = he->next()->sharp();
      h2->sharp() = he->next()->sharp();

      h2->next()->sharp() = he->prev()->sharp();
      h3->sharp() = he->prev()->sharp();
    }
  }

  static FT alpha(const FT n)
  {
    FT t = (3.0 + 2.0 * cos(2 * 3.1415926535 / n));
    FT a = 5.0 / 8.0 - t*t / 64.0;
    return n * (1-a) / a;
  }

  // smooth vertex positions
  static void smooth(Polyhedron *pMesh)
  {
    CGAL_assertion(pMesh != NULL);

    // alloc position vectors
    std::vector<Vector> pos(pMesh->size_of_vertices());

    // compute new positions
    unsigned int index = 0;
    typename Polyhedron::Vertex_iterator v;
    for(v = pMesh->vertices_begin();
        v != pMesh->vertices_end();
        v++)
    {
      // border vertices
      if(Polyhedron::is_border(v))
      {
        const Point& curr = v->point();
        const typename Polyhedron::Halfedge_handle& he = pMesh->get_border_halfedge(v);
        CGAL_assertion(he != NULL);
        const Point& next = he->next()->vertex()->point();
        const Point& prev = he->prev()->vertex()->point();
        pos[index] = 0.25 *(prev - CGAL::ORIGIN) +
                     0.50 *(curr - CGAL::ORIGIN) +
                     0.25 *(next - CGAL::ORIGIN);
      } // end is border
      else
      {
        unsigned int n = Polyhedron::valence(v);
        FT w = alpha((FT)n);
        FT sum_weights = w + (FT)n;

        // new position
        pos[index] = w/sum_weights * (v->point() - CGAL::ORIGIN);

        // rotate around vertex to compute new position
        HV_circulator he = v->vertex_begin();
        HV_circulator end = he;
        CGAL_For_all(he,end)
        {
          // add edge-vertex contribution
          const Point& point = he->opposite()->vertex()->point();
          pos[index] = pos[index] + (point - CGAL::ORIGIN) / sum_weights;
        }
      } // end !is border
      index++;
    }

    // set new positions
    for(v  = pMesh->vertices_begin(), index = 0;
        v != pMesh->vertices_end();
        v++,index++)
      v->point() = CGAL::ORIGIN + pos[index];
  }
};

template <class Polyhedron,class Kernel>
class CSubdivider_loop
{
  public:
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
    CSubdivider_loop() {}
    ~CSubdivider_loop() {}

  public:
    void subdivide(Polyhedron &input,
                   Polyhedron &output)
    {
      // subdivide, then smooth
      CModifierLoop<HalfedgeDS,Polyhedron,Kernel> builder(&input);
      output.delegate(builder);
    }
};


#endif // PWS_LOOP_SUBDIVISION
