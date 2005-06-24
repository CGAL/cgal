/***************************************************************************
Polyhedron_ex.h  -  description
                             -------------------
begin                : jan 2002
copyright            : (C) 2002 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef POLYHEDRON_EX_H_INCLUDED
#define POLYHEDRON_EX_H_INCLUDED

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Cartesian.h>

#include <algorithm>
#include <vector>
#include <list>
#include <stdio.h>


// CGAL kernel
typedef CGAL::Cartesian<double> My_kernel;


// compute facet center
struct Facet_center
{
    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;

    template<class Facet>
    void operator()(Facet& f)
    {
        Vector_3 vec(0.0,0.0,0.0);
        int degree = 0;
        typedef typename Facet::Halfedge_around_facet_const_circulator circ;
        circ h = f.facet_begin();
        do
        {
            vec = vec + (h->vertex()->point()-CGAL::ORIGIN);
            degree++;
        }
        while (++h != f.facet_begin());
        f.center() = CGAL::ORIGIN + (vec/degree);
    }
};

template<class Refs, class T>
class My_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
public:

    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;

    // life cycle
    // no constructors to repeat, since only
    // default constructor mandatory
    My_facet()
    {
        m_tag = -1;             // uninitialized
    }

    // center
    Point_3& center() { return m_center; }
    const Point_3& center() const { return m_center; }

    // tag
    int tag() const { return m_tag; }
    void tag(int tag) { m_tag = tag; }

    // distance
    double distance(Point_3 *pPoint) const
    {
        Vector_3 vec = (*pPoint-m_center);
        return std::sqrt(vec*vec);
    }

// Fields
private:

    // facet data
    int m_tag;
    Point_3 m_center;
};

template<class Refs, class Tprev, class Tvertex, class Tface>
class My_halfedge
    : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
    int m_tag;

    // parameterization
    bool m_is_parameterized;
    int m_seaming;              // seaming status
    double m_u;                 // texture coordinates
    double m_v;
    int m_index;                // for parameterization

    // surface cutting
    float m_distance;

public:
    // life cycle
    // no constructors to repeat, since only
    // default constructor mandatory
    My_halfedge()
    {
        m_tag = -1;             // uninitialized
        m_u = 0.0;
        m_v = 0.0;
        m_index = -1;           // uninitialized
        m_seaming = -1;         // uninitialized
        m_is_parameterized = false;
    }

    // tag
    int tag() const { return m_tag; }
    void tag(int tag) { m_tag = tag; }

    // seaming status
    int seaming() const { return m_seaming; }
    void seaming(int seaming) { m_seaming = seaming; }

    // precomputed distance
    float distance() const { return m_distance; }
    void distance(float distance) { m_distance = distance; }

    // texture coordinates
    double u() const { return m_u; }
    double v() const { return m_v; }
    void uv(double u, double v) 
    { 
#ifdef DEBUG_TRACE
        std::cerr << "      H" << index() << "(" << opposite()->vertex()->index() << "->" << vertex()->index() << ") <- (u=" << u << ",v=" << v << ")\n";
#endif
        m_u = u; 
        m_v = v; 
    }

    // param.
    bool is_parameterized() const { return m_is_parameterized; }
    void is_parameterized(bool is)  { m_is_parameterized = is; }

    // index
    int index() const { return m_index; }
    void index(int i) { m_index = i; }
};


// A redefined vertex class for the Polyhedron_3
template<class Refs, class T, class P>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
    // index
    int m_index;

    // misc
    int m_tag;

    // seaming status
    int m_seaming;              

public:
    // life cycle
    My_vertex()  { init(); }
    // repeat mandatory constructors
    My_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
    {
        init();
    }

    void init()
    {
        m_index = -1;           // uninitialized
        m_tag = -1;             // uninitialized
        m_seaming = -1;         // uninitialized
    }

    // index
    int index() const { return m_index; }
    void index(int i) { m_index = i; }

    // tag
    int tag() const { return m_tag; }
    void tag(int tag) { m_tag = tag; }

    // seaming status
    int seaming() const { return m_seaming; }
    void seaming(int seaming) { m_seaming = seaming; }
};


// A redefined items class for the Polyhedron_3 with a refined vertex, facet and halfedge classes
struct My_items : public CGAL::Polyhedron_items_3
{
    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;

    // wrap vertex
    template<class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point_3;
        typedef My_vertex<Refs,
                          CGAL::Tag_true,
                          Point_3> Vertex;
    };

    // wrap facet
    template<class Refs, class Traits>
    struct Face_wrapper
    {
        typedef My_facet<Refs,
                         CGAL::Tag_true> Face;
    };

    // wrap halfedge
    template<class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef My_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true> Halfedge;
    };
};


class Polyhedron_ex : public CGAL::Polyhedron_3<My_kernel,My_items>
{
public:

    typedef My_kernel::Vector_3 Vector_3;
    typedef My_kernel::Point_3  Point_3;

  public:

    // life cycle
    Polyhedron_ex() {}
    virtual ~Polyhedron_ex() {}

    // facet centers
    void compute_facet_centers()
    {
        fprintf(stderr,"  compute facet centers...");
        std::for_each(facets_begin(),facets_end(),Facet_center());
        fprintf(stderr,"ok\n");
    }

    // tag all facets
    void tag_facets(const int tag)
    {
        Facet_iterator pFace;
        for(pFace = facets_begin();
            pFace != facets_end();
            pFace++)
        pFace->tag(tag);
    }

    // get closest inner facet
    Facet_handle get_closest_inner_facet(Point_3 *pPoint)
    {
        Facet_iterator pFace = facets_begin();
        Facet_handle pClosest = pFace;
        double min = pFace->distance(pPoint);
        for(;pFace != facets_end();
            pFace++)
        {
            if(is_inner(pFace))
            {
                double distance = pFace->distance(pPoint);
                if(distance < min)
                {
                    pClosest = pFace;
                    min = distance;
                }
            }
        }
        return pClosest;
    }

    bool is_inner(Facet_handle pFace)
    {
        typedef Halfedge_around_facet_const_circulator circ;
        circ h = pFace->facet_begin();
        do
        {
            if(h->opposite()->is_border())
                return false;
        }
        while(++h != pFace->facet_begin());
        return true;
    }

    // tag all vertices
    void tag_vertices(const int tag)
    {
        Vertex_iterator iter;
        for(iter = vertices_begin(); iter != vertices_end(); iter++)
            iter->tag(tag);
    }

    // tag all halfedges
    void tag_halfedges(const int tag)
    {
        Halfedge_iterator iter;
        for(iter = halfedges_begin(); iter != halfedges_end(); iter++)
            iter->tag(tag);
    }

    // compute bounding interval
    double min(int coord)
    {
        CGAL_assertion(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        double min = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
            min = std::min(min,pVertex->point()[coord]);
        return min;
    }
    double max(int coord)
    {
        CGAL_assertion(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        double max = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
            max = std::max(max,pVertex->point()[coord]);
        return max;
    }
    Vertex_handle vertex_min(int coord,
                             double &min)
    {
        CGAL_assertion(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        Vertex_handle pBest = pVertex;
        min = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
        {
            double value = pVertex->point()[coord];
            if(value < min)
            {
                min = std::min(min,value);
                pBest = pVertex;
            }
        }
        return pBest;
    }
    Vertex_handle vertex_max(int coord,
                             double &max)
    {
        CGAL_assertion(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        Vertex_handle pBest = pVertex;
        max = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
        {
            double value = pVertex->point()[coord];
            if(value > max)
            {
                max = std::max(max,value);
                pBest = pVertex;
            }
        }
        return pBest;
    }

    // Index all mesh vertices following the order of the vertices_begin() iterator
    void precompute_vertex_indices()
    {
        Vertex_iterator pVertex;
        unsigned int i = 0;
        for(pVertex = vertices_begin();
            pVertex != vertices_end();
            pVertex++)
            pVertex->index(i++);
    }

    // Index all mesh half edges following the order of the halfedges_begin() iterator
    void precompute_halfedge_indices()
    {
        Halfedge_iterator pHalfedge;
        unsigned int i = 0;
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        pHalfedge->index(i++);
    }


    // Dump parameterized mesh to an eps file
    //********************************
    bool dump_eps(const char *pFilename,
                  double scale = 500.0)
    {
        std::cerr << "  dump mesh to " << pFilename << "..." << std::endl;
        FILE *pFile = fopen(pFilename,"wt");
        if(pFile == NULL)
        {
            std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
            return false;
        }

        // compute bounding box
        double xmin,xmax,ymin,ymax;
        xmin = ymin = xmax = ymax = 0;
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        {
            double x1 = scale * pHalfedge->prev()->u();
            double y1 = scale * pHalfedge->prev()->v();
            double x2 = scale * pHalfedge->u();
            double y2 = scale * pHalfedge->v();
            xmin = std::min(xmin,x1);
            xmin = std::min(xmin,x2);
            xmax = std::max(xmax,x1);
            xmax = std::max(xmax,x2);
            ymax = std::max(ymax,y1);
            ymax = std::max(ymax,y2);
            ymin = std::min(ymin,y1);
            ymin = std::min(ymin,y2);
        }

        fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
        fprintf(pFile,"%%%%BoundingBox: %d %d %d %d\n", int(xmin+0.5), int(ymin+0.5), int(xmax+0.5), int(ymax+0.5));
        fprintf(pFile,"%%%%HiResBoundingBox: %g %g %g %g\n",xmin,ymin,xmax,ymax);
        fprintf(pFile,"%%%%EndComments\n");
        fprintf(pFile,"gsave\n");
        fprintf(pFile,"0.1 setlinewidth\n");

        // color macros
        fprintf(pFile,"\n%% RGB color command - r g b C\n");
        fprintf(pFile,"/C { setrgbcolor } bind def\n");
        fprintf(pFile,"/white { 1 1 1 C } bind def\n");
        fprintf(pFile,"/black { 0 0 0 C } bind def\n");

        // edge macro -> E
        fprintf(pFile,"\n%% Black stroke - x1 y1 x2 y2 E\n");
        fprintf(pFile,"/E {moveto lineto stroke} bind def\n");
        fprintf(pFile,"black\n\n");

        // for each halfedge
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        {
            double x1 = scale * pHalfedge->prev()->u();
            double y1 = scale * pHalfedge->prev()->v();
            double x2 = scale * pHalfedge->u();
            double y2 = scale * pHalfedge->v();
            fprintf(pFile,"%g %g %g %g E\n",x1,y1,x2,y2);
        }

        /* Emit EPS trailer. */
        fputs("grestore\n\n",pFile);
        fputs("showpage\n",pFile);

        fclose(pFile);
        return true;
    }

	// Dump parameterized mesh to a Wavefront OBJ file
    // v x y z
    // f 1 2 3 4 (1-based)
    //********************************
    bool write_file_obj(const char *pFilename)
    {
        std::cerr << "  dump mesh to " << pFilename << "..." << std::endl;
        FILE *pFile = fopen(pFilename,"wt");
        if(pFile == NULL)
        {
            std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
            return false;
        }

        // Index all mesh vertices following the order of the vertices_begin() iterator
        precompute_vertex_indices();
        // Index all mesh half edges following the order of the halfedges_begin() iterator
        precompute_halfedge_indices();

        // write the name of material file
        fprintf(pFile, "mtllib parameterization.mtl\n") ;

        // output coordinates
        fprintf(pFile, "# vertices\n") ;
        Vertex_iterator pVertex;
        for(pVertex = vertices_begin(); pVertex != vertices_end(); pVertex++)
            fprintf(pFile,"v %g %g %g\n", (double)pVertex->point().x(), (double)pVertex->point().y(), (double)pVertex->point().z());

        // Write UVs (1 UV / halfedge)
        // Implementation note: the UV is meaningless for NON parameterized halfedges
        fprintf(pFile, "# uv coordinates\n") ;
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin(); pHalfedge != halfedges_end(); pHalfedge++)
        {
            if (!pHalfedge->is_border())
                fprintf(pFile, "vt %f %f\n", pHalfedge->u(), pHalfedge->v());
            else
                fprintf(pFile, "vt %f %f\n", 0.0, 0.0);
        }

        // Write facets using the unique material # 1
        fprintf(pFile, "# facets\nusemtl Mat_1\n");
        Facet_iterator pFacet;
        for(pFacet = facets_begin(); pFacet != facets_end(); pFacet++)
        {
            Halfedge_around_facet_circulator h = pFacet->facet_begin();
            fprintf(pFile,"f");
            do {
                fprintf(pFile, " %d", (int)h->vertex()->index()+1);
                if (h->is_parameterized())
                    fprintf(pFile, "/%d", (int)h->index()+1);
            }
            while(++h != pFacet->facet_begin());
            fprintf(pFile,"\n");
        }

        fclose(pFile);
        return true;
    }

    // is vertex on border ?
    static bool is_border(Vertex_const_handle pVertex)
    {
        Halfedge_around_vertex_const_circulator pHalfedge = pVertex->vertex_begin();
        Halfedge_around_vertex_const_circulator end = pHalfedge;
        if(pHalfedge == NULL) // isolated vertex
            return true;
        CGAL_For_all(pHalfedge,end)
            if(pHalfedge->is_border())
            return true;
        return false;
    }

    // compute distance from facet center to halfedge center
    double distance(Facet_handle pFacet,
                    Halfedge_handle pHalfedge)
    {
        // we assume
        Point_3 center_facet = pFacet->center();

        Vector_3 v = (pHalfedge->opposite()->vertex()->point()
                 - pHalfedge->vertex()->point());
        Point_3 center_halfedge = pHalfedge->vertex()->point() + (v/2);
        Vector_3 d = center_facet-center_halfedge;
        return std::sqrt(d*d);
    }

    void farthest_point_aligned(Vertex_handle &pVertexMin,
                                Vertex_handle &pVertexMax)
    {
        double xmin,xmax,ymin,ymax,zmin,zmax;
        Vertex_handle pVertex_xMin = vertex_min(0,xmin);
        Vertex_handle pVertex_xMax = vertex_max(0,xmax);
        Vertex_handle pVertex_yMin = vertex_min(1,ymin);
        Vertex_handle pVertex_yMax = vertex_max(1,ymax);
        Vertex_handle pVertex_zMin = vertex_min(2,zmin);
        Vertex_handle pVertex_zMax = vertex_max(2,zmax);
        double xdiff = xmax-xmin;
        double ydiff = ymax-ymin;
        double zdiff = zmax-zmin;
        double max = std::max(std::max(xdiff,ydiff),zdiff);
        if(max == xdiff)
        {
            pVertexMin = pVertex_xMin;
            pVertexMax = pVertex_xMax;
        }
        else if(max == ydiff)
        {
            pVertexMin = pVertex_yMin;
            pVertexMax = pVertex_yMax;
        }
        else
        {
            pVertexMin = pVertex_zMin;
            pVertexMax = pVertex_zMax;
        }
    }

}; // end class PolyhedronEx


#endif // POLYHEDRON_EX_H_INCLUDED

