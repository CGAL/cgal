/***************************************************************************
cgal_types.h  -  description
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

#ifndef CGAL_TYPES_H_INCLUDED
#define CGAL_TYPES_H_INCLUDED


#include <CGAL/basic.h>

// CGAL stuff

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <algorithm>
#include <vector>
#include <list>
#include <stdio.h>

#include "my_kernel.h"
#include "Feature_skeleton.h"

typedef My_kernel::Vector_3 Vector;
typedef My_kernel::Vector_2 Vector_2;
typedef My_kernel::Point_3  Point;
typedef My_kernel::Point_2  Point_2;
typedef My_kernel::Segment_2 Segment_2;

// compute facet center
struct Facet_center
{
    template<class Facet>
    void operator()(Facet& f)
    {
        Vector vec(0.0,0.0,0.0);
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
    // facet data
    int m_tag;
    Point m_center;

public:

    // life cycle
    // no constructors to repeat, since only
    // default constructor mandatory
    My_facet()
    {
        m_tag = 0;
    }

    // center
    Point& center() { return m_center; }
    const Point& center() const { return m_center; }

    // tag
    int tag() const { return m_tag; }
    void tag(int tag) { m_tag = tag; }

    // distance
    double distance(Point *pPoint) const
    {
        Vector vec = (*pPoint-m_center);
        return My_kernel::len(vec);
    }
};

template<class Refs, class Tprev, class Tvertex, class Tface>
class My_halfedge
    : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
    int m_tag;

    // parameterization
    bool m_is_parameterized;
    int m_seaming;  // seaming status
    double m_u;             // texture coordinates
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
        m_tag = 0;
        m_u = 0.0;
        m_v = 0.0;
        m_index = 0;
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
    void u(double u) { m_u = u; }
    void v(double v) { m_v = v; }
    void uv(double u,double v) { m_u = u; m_v = v; }

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
    }

    // index
    int index() const { return m_index; }
    void index(int i) { m_index = i; }

    // tag
    int tag() const { return m_tag; }
    void tag(int tag) { m_tag = tag; }
};

// A redefined items class for the Polyhedron_3 with a refined vertex, facet and halfedge classes
struct My_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template<class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef My_vertex<Refs,
                          CGAL::Tag_true,
                          Point> Vertex;
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


typedef CGAL::Polyhedron_3<My_kernel,My_items> Polyhedron;

class Polyhedron_ex : public Polyhedron
{
public:

    // feature/boundary skeletons
    typedef Feature_backbone<Vertex_handle,Halfedge_handle> backbone;
    typedef Feature_skeleton<Vertex_handle,Halfedge_handle> skeleton;

 private:

    // feature/boundary skeletons
    skeleton m_skeleton;
    backbone m_seaming_backbone;

  public:

    // data access
    skeleton* get_skeleton() { return &m_skeleton; }
    const skeleton* get_skeleton() const { return &m_skeleton; }
    backbone* get_seaming_backbone() { return &m_seaming_backbone; }
    const backbone* get_seaming_backbone() const { return &m_seaming_backbone; }

    // life cycle
    //********************************
    Polyhedron_ex()
    {
    }
    virtual ~Polyhedron_ex()
    {
        free_skeleton();
    }

    void free_skeleton()
    {
        m_skeleton.free();
    }


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

    // get any facet with tag
    Facet_handle get_any_facet_tag(int tag)
    {
        Facet_iterator pFace;
        for(pFace = facets_begin();
            pFace != facets_end();
            pFace++)
        if(pFace->tag() == tag)
            return pFace;
        return NULL;
    }

    // get closest inner facet
    Facet_handle get_closest_inner_facet(Point *pPoint)
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
        for(iter = vertices_begin();
            iter != vertices_end();
            iter++)
        {
            Vertex_handle hVertex = iter;
            hVertex->tag(tag);
        }
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

    // tag all halfedges
    void tag_halfedges(const int tag)
    {
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        pHalfedge->tag(tag);
    }

    // Set seaming status of all halfedges
    void flag_halfedges_seaming(int flag)
    {
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        pHalfedge->seaming(flag);
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


    // dump the param to an eps file
    //********************************
    bool dump_param(const char *pFilename)
    {
        std::cerr << "  dump parameterization to " << pFilename << "..." << std::endl;
        FILE *pFile = fopen(pFilename,"wt");
        if(pFile == NULL)
        {
            std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
            return false;
        }
        dump_param(pFile);
        fclose(pFile);
        return true;
    }

    // dump the param to the stdout
    //********************************
    void dump_param()  { dump_param(stdout); }

    // dump the param to an eps file
    //********************************
    void dump_param(FILE *pFile,
                    double scale = 500.0)
    {
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
        fprintf(pFile,"%%%%BoundingBox: %g %g %g %g\n",xmin,ymin,xmax,ymax);
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
    }

    // output to a Wavefront OBJ file
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

        bool ok = write_file_obj(pFile);

        fclose(pFile);

        return ok;
    }

    float round(double x)
    {
        return x;                                       // Default implementation
        //return float(int(x*100.0 + 0.5))/100.0;       // Round number to ease files comparison
    }

    // output to a Wavefront OBJ file
    // v x y z
    // f 1 2 3 4 (1-based)
    //********************************
    bool write_file_obj(FILE *pFile)
    {
        fprintf(stderr,"  write_file_obj()...");

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
            fprintf(pFile, "vt %f %f\n", round(pHalfedge->u()), round(pHalfedge->v()));

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

        fprintf(stderr,"ok\n");

        return true;
    }

    // dump Wavefront OBJ file to stdout
    //********************************
    bool write_file_obj()  { return write_file_obj(stdout); }

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

    // halfedge len
    double len(Polyhedron::Halfedge_handle halfedge)
    {
        Vector v = (halfedge->vertex()->point()-
                    halfedge->prev()->vertex()->point());
        return std::sqrt(v*v);
    }

    // get any border halfedge with tag
    Halfedge_handle get_border_halfedge_tag(int tag)
    {
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        {
            if(pHalfedge->is_border() &&
                pHalfedge->tag() == tag)
            return pHalfedge;
        }
        return NULL;
    }

    // get index of the longest backbone
    int get_index_longest_backbone()
    {
        int index = 0;
        double max = 0.0;
        // #backbones
        int nb = (*m_skeleton.backbones()).size();
        for(int i=0;i<nb;i++)
        {
            backbone *pBackbone = (*m_skeleton.backbones())[i];
            double length = len(pBackbone);
            if(length>max)
            {
            index = i;
            max = length;
            }
        }
        return index;
    }

    // count #boundaries
    // return the number of boundary backbones
    int nb_boundaries()
    {
        int nb = 0;
        tag_halfedges(0);
        Halfedge_handle seed_halfedge = NULL;
        while((seed_halfedge = get_border_halfedge_tag(0)) != NULL)
        {
            nb++;
            seed_halfedge->tag(1);
            Vertex_handle seed_vertex = seed_halfedge->prev()->vertex();
            Halfedge_handle current_halfedge = seed_halfedge;
            Halfedge_handle next_halfedge;
            do
            {
                next_halfedge = current_halfedge->next();
                next_halfedge->tag(1);
                current_halfedge = next_halfedge;
            }
            while(next_halfedge->prev()->vertex() != seed_vertex);
        }
        return nb;
    }

    // tag component
    void tag_component(Facet_handle pSeedFacet,
                       const int tag_free,
                       const int tag_done)
    {
        pSeedFacet->tag(tag_done);
        std::list<Facet_handle> facets;
        facets.push_front(pSeedFacet);
        while(!facets.empty())
        {
            Facet_handle pFacet = facets.front();
            facets.pop_front();
            pFacet->tag(tag_done);
            Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
            Halfedge_around_facet_circulator end = pHalfedge;
            CGAL_For_all(pHalfedge,end)
            {
                Facet_handle pNFacet = pHalfedge->opposite()->facet();
                if(pNFacet != NULL && pNFacet->tag() == tag_free)
                {
                    facets.push_front(pNFacet);
                    pNFacet->tag(tag_done);
                }
            }
        }
    }

    // count # of connected components
    unsigned int nb_components()
    {
        unsigned int nb = 0;
        tag_facets(0);
        Facet_handle seed_facet = NULL;
        while((seed_facet = get_any_facet_tag(0)) != NULL)
        {
            nb++;
            tag_component(seed_facet,0,1);
        }
        return nb;
    }

    // compute the genus
    // G = (2*C + E - B - F - V)/2 with
    // G : genus
    // C : # of connected components
    // E : # of edges
    // B : # of boundaries
    // F : # of facets
    // V : # of vertices
    int genus()
    {
        int c = nb_components();
        int b = nb_boundaries();
        int v = size_of_vertices();
        int e = size_of_halfedges()/2;
        int f = size_of_facets();
        int genus = (2*c+e-b-f-v)/2;
        std::cerr << "  " << v << " vertices, " << f << " facets, ";
        std::cerr << e << " edges, " << b << " boundary(ies), genus " << genus << std::endl;
        return genus;
    }

    // compute  total len of a backbone
    double len(backbone *pBackbone)
    {
        std::list<Polyhedron_ex::Halfedge_handle> *pHalfedges = pBackbone->halfedges();
        std::list<Polyhedron_ex::Halfedge_handle>::iterator pHalfedge;
        double len = 0.0;
        for(pHalfedge = pHalfedges->begin();
            pHalfedge != pHalfedges->end();
            pHalfedge++)
        {
            Polyhedron_ex::Halfedge_handle he = (*pHalfedge);
            Vector v = (he->vertex()->point()-he->prev()->vertex()->point());
            len += My_kernel::len(v);
        }
        return len;
    }

    // compute distance from facet center to halfedge center
    double distance(Facet_handle pFacet,
                    Halfedge_handle pHalfedge)
    {
        // we assume
        Point center_facet = pFacet->center();

        Vector v = (pHalfedge->opposite()->vertex()->point()
                 - pHalfedge->vertex()->point());
        Point center_halfedge = pHalfedge->vertex()->point() + (v/2);
        Vector d = center_facet-center_halfedge;
        return My_kernel::len(d);
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


#endif // CGAL_TYPES_H_INCLUDED

