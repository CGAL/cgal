/***************************************************************************
quad-triangle.h  -  quad/triangle subdivision
----------------------------------------------------------------------------
begin                : march 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/

#ifndef QUAD_TRIANGLE
#define QUAD_TRIANGLE

#include "enriched_polyhedron.h"
#include "builder.h"

// "Quad/Triangle Subdivision"
// Levin-Levin 02 
// Stam-Loop 03

template <class HDS,class Polyhedron,class kernel>
class CModifierQuadTriangle : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename HDS::Vertex          Vertex;
  typedef typename Vertex::Point        Point;
  typedef typename HDS::Face_handle     Face_handle;
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename CGAL::Enriched_polyhedron_incremental_builder_3<HDS> builder;

  Polyhedron *m_pMesh;

public:

  // life cycle
  CModifierQuadTriangle(Polyhedron *pMesh)
  {
    CGAL_assertion(pMesh != NULL);
    m_pMesh = pMesh;
  }
	~CModifierQuadTriangle() {}

	// subdivision
  void operator()( HDS& hds)
  {
    builder B(hds,true);
    B.begin_surface(3,1,6);
      add_vertices(B);
      add_facets(B);
    B.end_surface();
  }

  //***************************************************
  // add vertices
  //***************************************************
  void add_vertices(builder &B)
  {
    // put original vertices
    int index = 0;
    Polyhedron::Vertex_iterator pVertex;
    for(pVertex = m_pMesh->vertices_begin();
        pVertex != m_pMesh->vertices_end();
        pVertex++)
    {
      pVertex->tag(index);
      B.add_vertex(pVertex->point());
      index++;
    }

    // as many as #edges
    m_pMesh->tag_halfedges(-1);
    Polyhedron::Halfedge_iterator pHalfedge;
    for(pHalfedge = m_pMesh->halfedges_begin();
        pHalfedge != m_pMesh->halfedges_end();
        pHalfedge++)
    {
      if(pHalfedge->tag() != -1)
        continue;

      // simple edge bissection
      const Point& p1 = pHalfedge->vertex()->point();
      const Point& p2 = pHalfedge->opposite()->vertex()->point();
      Point point = Point(0.5f*(p1.x()+p2.x()),
                          0.5f*(p1.y()+p2.y()),
                          0.5f*(p1.z()+p2.z()));
      B.add_vertex(point);

      // put vertex indices on both halfedges
      pHalfedge->tag(index);
      pHalfedge->opposite()->tag(index);
      index++;
    }

    // and as many as #facets with degree > 3
    m_pMesh->tag_facets(-1);
    Polyhedron::Facet_iterator pFacet;
    for(pFacet = m_pMesh->facets_begin();
        pFacet != m_pMesh->facets_end();
        pFacet++)
    {
      unsigned int degree = Polyhedron::degree(pFacet);
      CGAL_assertion(degree >= 3);
      // forget about triangles, they are
      // simply 1-to-4 subdivided by edge bisection
      if(degree == 3)
        continue; 

      // barycentric subdivision
      Point barycenter;
	    m_pMesh->compute_facet_center(pFacet,barycenter);
      B.add_vertex(barycenter);
      pFacet->tag(index);
      index++;
    }
  }

  //***************************************************
  // add facets
  //***************************************************
  void add_facets(builder &B)
  {
    Polyhedron::Facet_iterator pFacet;
    for(pFacet = m_pMesh->facets_begin();
        pFacet != m_pMesh->facets_end();
        pFacet++)
    {
      unsigned int degree = Polyhedron::degree(pFacet);
      CGAL_assertion(degree >= 3);

      if(degree == 3)
      {
        Polyhedron::Halfedge_handle pHalfedge = pFacet->halfedge();
        int i0 = pHalfedge->tag();
        int i1 = pHalfedge->vertex()->tag();
        int i2 = pHalfedge->next()->tag();
        int i3 = pHalfedge->next()->vertex()->tag();
        int i4 = pHalfedge->next()->next()->tag();
        int i5 = pHalfedge->next()->next()->vertex()->tag();
        CGAL_assertion(i0 >= 0);
        CGAL_assertion(i1 >= 0);
        CGAL_assertion(i2 >= 0);
        CGAL_assertion(i3 >= 0);
        CGAL_assertion(i4 >= 0);
        CGAL_assertion(i5 >= 0);

        // create 4 triangles
        B.begin_facet();
          B.add_vertex_to_facet(i1);
          B.add_vertex_to_facet(i2);
          B.add_vertex_to_facet(i0);
        const Halfedge_handle& h1 = B.end_facet();
        B.begin_facet();
          B.add_vertex_to_facet(i3);
          B.add_vertex_to_facet(i4);
          B.add_vertex_to_facet(i2);
        const Halfedge_handle& h2 = B.end_facet();
        B.begin_facet();
          B.add_vertex_to_facet(i5);
          B.add_vertex_to_facet(i0);
          B.add_vertex_to_facet(i4);
        const Halfedge_handle& h3 = B.end_facet();

        // center face
        B.begin_facet();
          B.add_vertex_to_facet(i0);
          B.add_vertex_to_facet(i2);
          B.add_vertex_to_facet(i4);
        Halfedge_handle h4 = B.end_facet();

        h1->control_edge(false);
        h1->next()->control_edge(false);
        h1->next()->next()->control_edge(false);

        h2->control_edge(false);
        h2->next()->control_edge(false);
        h2->next()->next()->control_edge(false);

        h3->control_edge(false);
        h3->next()->control_edge(false);
        h3->next()->next()->control_edge(false);

        h4->control_edge(false);
        h4->next()->control_edge(false);
        h4->next()->next()->control_edge(false);
        
        if(pHalfedge->control_edge())
        {
          h1->control_edge(true);
          h3->next()->control_edge(true);
        }
        if(pHalfedge->next()->control_edge())
        {
          h1->next()->control_edge(true);
          h2->control_edge(true);
        }
        if(pHalfedge->next()->next()->control_edge())
        {
          h2->next()->control_edge(true);
          h3->control_edge(true);
        }
      }
      else
      {
        // i1: index of barycenter vertex
        int i1 = pFacet->tag();
        CGAL_assertion(i1 >= 0);

        // for each halfedge
        Polyhedron::Halfedge_around_facet_circulator h;
        h = pFacet->facet_begin();
        do
        {
          // i2,3,4: indices of three consecutive
          // vertices on halfedges
          int i2 = h->tag();
          int i3 = h->vertex()->tag();
          int i4 = h->next()->tag();
          CGAL_assertion(i2 >= 0);
          CGAL_assertion(i3 >= 0);
          CGAL_assertion(i4 >= 0);

          // create a quad
          B.begin_facet();
            B.add_vertex_to_facet(i3);
            B.add_vertex_to_facet(i4);
            B.add_vertex_to_facet(i1);
            B.add_vertex_to_facet(i2);
          const Halfedge_handle& pNewHalfedge = B.end_facet();

          pNewHalfedge->control_edge(false);
          pNewHalfedge->next()->control_edge(false);
          pNewHalfedge->next()->next()->control_edge(false);
          pNewHalfedge->next()->next()->next()->control_edge(false);

          if(h->control_edge())
            pNewHalfedge->control_edge(true);
          if(h->next()->control_edge())
            pNewHalfedge->next()->control_edge(true);
        }
        while(++h != pFacet->facet_begin());
      }
    }
  }

  //***************************************************
  // correcting factor
  //***************************************************
	static float correcting_factor(unsigned int ne,
		                             unsigned int nq)
  {
    if(ne == 2 && nq == 1)
      return -0.20505f;
    if(ne == 3 && nq == 1)
      return 0.80597f;
    if(ne == 3 && nq == 2)
      return 0.61539f;
    if(ne == 4 && nq == 1)
      return 0.34792f;
    if(ne == 4 && nq == 2)
      return 0.21380f;
    if(ne == 4 && nq == 3)
      return 0.10550f;
    return 0.0f;
  }
  
  //***************************************************
  // smooth vertex positions
  //***************************************************
  static void smooth(Polyhedron *pMesh,
                     bool smooth_boundary = true)
  {
    CGAL_assertion(pMesh != NULL);

    // alloc position vectors
    unsigned int nb_vertices = pMesh->size_of_vertices();
    kernel::FT *pPos = new kernel::FT[3*nb_vertices];
    CGAL_assertion(pPos != NULL);

    // compute new positions
    unsigned int index = 0;
    Polyhedron::Vertex_iterator pVertex;
    for(pVertex = pMesh->vertices_begin();
        pVertex != pMesh->vertices_end();
        pVertex++)
    {
      // border vertices will not move
      if(Polyhedron::is_border(pVertex))
      {
        // do not smooth it
        const Point& curr = pVertex->point();
        if(!smooth_boundary)
        {
          pPos[3*index]   = curr.x();
          pPos[3*index+1] = curr.y();
          pPos[3*index+2] = curr.z();
        }
        // smooth using [1/4 1/2 1/4] cubic B-spline averaging mask
        else 
        {
          const Polyhedron::Halfedge_handle& pHalfedge =
            pMesh->get_border_halfedge(pVertex);
          CGAL_assertion(pHalfedge != NULL);
          const Point& next = pHalfedge->next()->vertex()->point();
          const Point& prev = pHalfedge->prev()->vertex()->point();
          pPos[3*index]   = 0.25f*prev.x() + 0.5f*curr.x() + 0.25f*next.x();
          pPos[3*index+1] = 0.25f*prev.y() + 0.5f*curr.y() + 0.25f*next.y();
          pPos[3*index+2] = 0.25f*prev.z() + 0.5f*curr.z() + 0.25f*next.z();
        }
      } // end is border
      else
      {

        unsigned int nb_quads = 0;
        unsigned int nb_edges = 0;

        // rotate around vertex to count #edges and #quads
        Polyhedron::Halfedge_around_vertex_circulator
          pHalfEdge = pVertex->vertex_begin();
        Polyhedron::Halfedge_around_vertex_circulator end = pHalfEdge;
        CGAL_For_all(pHalfEdge,end)
        {
          const Polyhedron::Facet_handle& pFacet = pHalfEdge->facet();
          CGAL_assertion(pFacet != NULL);
          unsigned int degree = Polyhedron::degree(pFacet);
          CGAL_assertion(degree == 4 || degree == 3);
          if(degree == 4)
            nb_quads++;
          nb_edges++;
        }

        // compute coefficients
        kernel::FT ne = (kernel::FT)nb_edges;
        kernel::FT nq = (kernel::FT)nb_quads;
        kernel::FT alpha = 1.0f / (1.0f + ne/2.0f + nq/4.0f);
        kernel::FT beta = alpha / 2.0f;  // edges
        kernel::FT gamma = alpha / 4.0f; // corners of incident quads
        kernel::FT eta = correcting_factor(nb_edges,nb_quads);

        // new position
        pPos[3*index]   = alpha * pVertex->point().x();
        pPos[3*index+1] = alpha * pVertex->point().y();
        pPos[3*index+2] = alpha * pVertex->point().z();

        // rotate around vertex to compute new position
        pHalfEdge = pVertex->vertex_begin();
        end = pHalfEdge;
        CGAL_For_all(pHalfEdge,end)
        {
          const Polyhedron::Facet_handle& pFacet = pHalfEdge->facet();
          CGAL_assertion(pFacet != NULL);
          unsigned int degree = Polyhedron::degree(pFacet);
          CGAL_assertion(degree == 4 || degree == 3);

          // add edge-vertex contribution
          const Point& point = pHalfEdge->prev()->vertex()->point();
          pPos[3*index]   += beta * point.x();
          pPos[3*index+1] += beta * point.y();
          pPos[3*index+2] += beta * point.z();

          // add corner vertex contribution
          if(degree == 4)
          {
            const Point& corner = pHalfEdge->next()->next()->vertex()->point();
            pPos[3*index]   += gamma * corner.x();
            pPos[3*index+1] += gamma * corner.y();
            pPos[3*index+2] += gamma * corner.z();
          }
        }

        // apply correction
        pPos[3*index] = pPos[3*index] +
          eta*(pPos[3*index]-pVertex->point().x());
        pPos[3*index+1] = pPos[3*index+1] +
          eta*(pPos[3*index+1]-pVertex->point().y());
        pPos[3*index+2] = pPos[3*index+2] +
          eta*(pPos[3*index+2]-pVertex->point().z());

      } // end !is border
      index++;
    }

    // set new positions
    index = 0;
    for(pVertex = pMesh->vertices_begin();
        pVertex != pMesh->vertices_end();
        pVertex++)
    {
      Point& point = pVertex->point();
      point = Point(pPos[3*index],
                    pPos[3*index+1],
                    pPos[3*index+2]);
      index++;
    }

    // cleanup
    delete [] pPos;
  }
};

template <class Polyhedron,class kernel>
class CSubdivider_quad_triangle
{
	public:
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
		CSubdivider_quad_triangle() {}
		~CSubdivider_quad_triangle() {}

	public:
		void subdivide(Polyhedron &OriginalMesh,
                   Polyhedron &NewMesh,
                   bool smooth_boundary = true)
    {
      // subdivide, then smooth
      CModifierQuadTriangle<HalfedgeDS,Polyhedron,kernel> builder(&OriginalMesh);
      NewMesh.delegate(builder);
      builder.smooth(&NewMesh,smooth_boundary);
    }
};


#endif // QUAD_TRIANGLE
