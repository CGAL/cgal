#ifndef POLYHEDRON_EX_H_INCLUDED
#define POLYHEDRON_EX_H_INCLUDED

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <algorithm>
#include <vector>
#include <list>
#include <fstream>
#include <cassert>


// CGAL kernel
typedef CGAL::Simple_cartesian<double> My_kernel;


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
    double distance(Point_3& point) const
    {
        Vector_3 vec = (point-m_center);
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
    double m_distance;

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
    double distance() const { return m_distance; }
    void distance(double distance) { m_distance = distance; }

    // texture coordinates
    double u() const { return m_u; }
    double v() const { return m_v; }
    void uv(double u, double v) { m_u = u; m_v = v; }

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
        std::for_each(facets_begin(),facets_end(),Facet_center());
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
    Facet_handle get_closest_inner_facet(Point_3& point)
    {
        Facet_iterator pFace = facets_begin();
        Facet_handle pClosest = pFace;
        double minimum = pFace->distance(point);
        for(;pFace != facets_end();
            pFace++)
        {
            if(is_inner(pFace))
            {
                double distance = pFace->distance(point);
                if(distance < minimum)
                {
                    pClosest = pFace;
                    minimum = distance;
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
    double minimum (int coord)
    {
        assert(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        double minimum = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
            minimum = (std::min)(minimum,pVertex->point()[coord]);
        return minimum;
    }
    double maximum (int coord)
    {
        assert(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        double maximum = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
            maximum = (std::max)(maximum,pVertex->point()[coord]);
        return maximum;
    }
    Vertex_handle vertex_min(int coord,
                             double &minimum)
    {
        assert(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        Vertex_handle pBest = pVertex;
        minimum = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
        {
            double value = pVertex->point()[coord];
            if(value < minimum)
            {
                minimum = (std::min)(minimum,value);
                pBest = pVertex;
            }
        }
        return pBest;
    }
    Vertex_handle vertex_max(int coord,
                             double &maximum)
    {
        assert(size_of_vertices() > 0);
        Vertex_iterator pVertex = vertices_begin();
        Vertex_handle pBest = pVertex;
        maximum = pVertex->point()[coord];
        for(;pVertex != vertices_end();pVertex++)
        {
            double value = pVertex->point()[coord];
            if(value > maximum)
            {
                maximum = (std::max)(maximum,value);
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


#ifdef DEBUG_TRUNCATE_OUTPUT
    // Debug: write coordinates with 2 digits precision
    #define FORMAT_EPS_COORD(x) (int(x/10.0+0.5)*10)
#else
    #define FORMAT_EPS_COORD(x) (x)
#endif

    // Dump parameterized mesh to an eps file
    bool write_file_eps(const char *pFilename,
                        double scale = 500.0)
    {
        assert(pFilename != NULL);

        std::ofstream out(pFilename);
        if(!out)
            return false;
        CGAL::set_ascii_mode(out);

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
            xmin = (std::min)(xmin,x1);
            xmin = (std::min)(xmin,x2);
            xmax = (std::max)(xmax,x1);
            xmax = (std::max)(xmax,x2);
            ymax = (std::max)(ymax,y1);
            ymax = (std::max)(ymax,y2);
            ymin = (std::min)(ymin,y1);
            ymin = (std::min)(ymin,y2);
        }

        out << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
        out << "%%BoundingBox: " << int(xmin+0.5) << " "
                                 << int(ymin+0.5) << " "
                                 << int(xmax+0.5) << " "
                                 << int(ymax+0.5) << std::endl;
        out << "%%HiResBoundingBox: " << xmin << " "
                                      << ymin << " "
                                      << xmax << " "
                                      << ymax << std::endl;
        out << "%%EndComments" << std::endl;
        out << "gsave" << std::endl;
        out << "0.1 setlinewidth" << std::endl;

        // color macros
        out << std::endl;
        out << "% RGB color command - r g b C" << std::endl;
        out << "/C { setrgbcolor } bind def" << std::endl;
        out << "/white { 1 1 1 C } bind def" << std::endl;
        out << "/black { 0 0 0 C } bind def" << std::endl;

        // edge macro -> E
        out << std::endl;
        out << "% Black stroke - x1 y1 x2 y2 E" << std::endl;
        out << "/E {moveto lineto stroke} bind def" << std::endl;
        out << "black" << std::endl << std::endl;

        // output edge coordinates
        for(pHalfedge = halfedges_begin();
            pHalfedge != halfedges_end();
            pHalfedge++)
        {
            double x1 = scale * pHalfedge->prev()->u();
            double y1 = scale * pHalfedge->prev()->v();
            double x2 = scale * pHalfedge->u();
            double y2 = scale * pHalfedge->v();
            out << FORMAT_EPS_COORD(x1) << " "
                << FORMAT_EPS_COORD(y1) << " "
                << FORMAT_EPS_COORD(x2) << " "
                << FORMAT_EPS_COORD(y2) << " E" << std::endl;
        }

        /* Emit EPS trailer. */
        out << "grestore" << std::endl;
        out << std::endl;
        out << "showpage" << std::endl;

        return true;
    }

#ifdef DEBUG_TRUNCATE_OUTPUT
    // Debug: write coordinates with 2 digits precision
    #define FORMAT_UV(x) (float(int(x*100.0+0.5))/100.0)
#else
    #define FORMAT_UV(x) (x)
#endif

    // Dump parameterized mesh to a Wavefront OBJ file
    // v x y z
    // f 1 2 3 4 (1-based)
    //
    // Implementation note: the UV is meaningless for a NON parameterized halfedge
    bool write_file_obj(const char *pFilename)
    {
        assert(pFilename != NULL);

        std::ofstream out(pFilename);
        if(!out)
            return false;
        CGAL::set_ascii_mode(out);

        // Index all mesh vertices following the order of vertices_begin() iterator
        precompute_vertex_indices();
        // Index all mesh half edges following the order of halfedges_begin() iterator
        precompute_halfedge_indices();

        // write the name of material file
        out <<  "mtllib parameterization.mtl" << std::endl ;

        // output coordinates
        out <<  "# vertices" << std::endl ;
        Vertex_iterator pVertex;
        for(pVertex = vertices_begin(); pVertex != vertices_end(); pVertex++)
            out << "v " << pVertex->point().x() << " "
                        << pVertex->point().y() << " "
                        << pVertex->point().z() << std::endl;

        // Write UVs (1 UV / halfedge)
        out <<  "# uv coordinates" << std::endl ;
        Halfedge_iterator pHalfedge;
        for(pHalfedge = halfedges_begin(); pHalfedge != halfedges_end(); pHalfedge++)
        {
            if (pHalfedge->is_parameterized())
                out << "vt " << FORMAT_UV(pHalfedge->u()) << " " << FORMAT_UV(pHalfedge->v()) << std::endl;
            else
                out << "vt " << 0.0 << " " << 0.0 << std::endl;
        }

        // Write facets using the unique material # 1
        out << "# facets" << std::endl;
        out << "usemtl Mat_1" << std::endl;
        Facet_const_iterator pFacet;
        for(pFacet = facets_begin(); pFacet != facets_end(); pFacet++)
        {
            Halfedge_around_facet_const_circulator h = pFacet->facet_begin();
            out << "f";
            do {
                out << " " << h->vertex()->index()+1;
                if (h->is_parameterized())
                    out <<  "/" << h->index()+1;
            }
            while(++h != pFacet->facet_begin());
            out << std::endl;
        }

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
        if (xdiff >= (std::max)(ydiff,zdiff))
        {
            pVertexMin = pVertex_xMin;
            pVertexMax = pVertex_xMax;
        }
        else if (ydiff >= (std::max)(xdiff,zdiff))
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

