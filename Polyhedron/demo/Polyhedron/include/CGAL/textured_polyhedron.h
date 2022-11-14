#ifndef _TEXTURED_MESH_
#define _TEXTURED_MESH_

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_default.h>


namespace CGAL

{

template <class Refs, class T, class P, class Norm>
class Textured_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
    // normal
    Norm m_normal;
    //id
    std::size_t m_id;

public:

    // life cycle
    // no constructors to repeat, since only
    // default constructor mandatory

    Textured_facet()
    {
    }

    // normal
    typedef Norm Normal_3;
    Normal_3& normal() { return m_normal; }
    const Normal_3& normal() const { return m_normal; }
    // id
    std::size_t& id() {        return m_id; }
    const std::size_t& id() const { return m_id; }
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Textured_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
    double m_u;
    double m_v;
    std::size_t m_id;
public:
    // life cycle
    Textured_halfedge()
    {
    }
    // u,v coordinates
    double& u() {        return m_u; }
    const double& u() const { return m_u;        }
    double& v() {        return m_v; }
    const double& v() const { return m_v;        }
    // id
    std::size_t& id() {        return m_id; }
    const std::size_t& id() const { return m_id;        }


};

template <class Refs, class T, class P, class Norm>
class Textured_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
    // normal
    Norm m_normal;
    double m_u;
    double m_v;
    std::size_t m_id;

public:
    // life cycle
    Textured_vertex()  {}

    // repeat mandatory constructors
    Textured_vertex(const P& pt)
        : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
    {
    }

    // normal
    typedef Norm Normal_3;
    Normal_3& normal() { return m_normal; }
    const Normal_3& normal() const { return m_normal; }

    // u,v coordinates
    double& u() {        return m_u; }
    const double& u() const { return m_u;        }
    double& v() {        return m_v; }
    const double& v() const { return m_v;        }
    // id
    std::size_t& id() {        return m_id; }
    const std::size_t& id() const { return m_id;        }
};

struct Textured_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template<class Refs, class Traits> struct Vertex_wrapper
    {
        typedef typename Traits::Point_3 Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Textured_vertex<Refs,
        CGAL::Tag_true,
        Point,
        Normal> Vertex;
    };

    // wrap face
    template<class Refs, class Traits> struct Face_wrapper
    {
        typedef typename Traits::Point_3 Point;
        typedef typename Traits::Vector_3 Normal;
        typedef Textured_facet<Refs,
        CGAL::Tag_true,
        Point,
        Normal> Face;
    };

    // wrap halfedge
    template<class Refs, class Traits> struct Halfedge_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef Textured_halfedge<Refs,
        CGAL::Tag_true,
        CGAL::Tag_true,
        CGAL::Tag_true,
        Normal> Halfedge;
    };
};

// compute facet normal
struct Facet_normal // (functor)
{
    template<class Facet> void operator()(Facet& f)
    {
        typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
        typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
        do
        {
            typename Facet::Normal_3 normal = CGAL::cross_product(h->next()->vertex()->point() - h->vertex()->point(), h->next()->next()->vertex()->point() - h->next()->vertex()->point());
            double sqnorm = normal * normal;
            if (sqnorm != 0)
                normal = normal / (float)std::sqrt(sqnorm);
            sum = sum + normal;
        } while (++h != f.facet_begin());
        float sqnorm = sum * sum;
        if (sqnorm != 0.0)
            f.normal() = sum / std::sqrt(sqnorm);
        else
            f.normal() = CGAL::NULL_VECTOR;
    }
};

// compute vertex normal
struct Vertex_normal // (functor)
{
    template<class Vertex> void operator()(Vertex& v)
    {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge =
                v.vertex_begin();
        typename Vertex::Halfedge_around_vertex_const_circulator begin =
                pHalfedge;
        CGAL_For_all(pHalfedge,begin)
                if(!pHalfedge->is_border())
                normal = normal + pHalfedge->facet()->normal();
        float sqnorm = normal * normal;
        if (sqnorm != 0.0f)
            v.normal() = normal / (float)std::sqrt(sqnorm);
        else
            v.normal() = CGAL::NULL_VECTOR;
    }
};

//*********************************************************
template <class Kernel, class Items>
class Textured_polyhedron : public CGAL::Polyhedron_3<Kernel,Items>
{
public :
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename CGAL::Polyhedron_3<Kernel,Items> Base;
    typedef typename CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_3> Basic_polyhedron;

    typedef typename Base::Vertex_handle Vertex_handle;
    typedef typename Base::Vertex_iterator Vertex_iterator;
    typedef typename Base::Halfedge_handle Halfedge_handle;
    typedef typename Base::Halfedge_iterator Halfedge_iterator;
    typedef typename Base::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    typedef typename Base::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
    typedef typename Base::Edge_iterator Edge_iterator;
    typedef typename Base::Facet Facet;
    typedef typename Base::Facet_iterator Facet_iterator;
    typedef typename Base::Facet_handle Facet_handle;

public :

    // life cycle
    Textured_polyhedron()
    {
    }

    virtual ~Textured_polyhedron()
    {
    }

    // normals (per facet, then per vertex)
    void compute_normals_per_facet()
    {
        std::for_each(this->facets_begin(),this->facets_end(),Facet_normal());
    }
    void compute_normals_per_vertex()
    {
        std::for_each(this->vertices_begin(),this->vertices_end(),Vertex_normal());
    }
    void compute_normals()
    {
        compute_normals_per_facet();
        compute_normals_per_vertex();
    }

}; // end class Textured_polyhedron

} // end namespace CGAL

#endif // _TEXTURED_MESH_
