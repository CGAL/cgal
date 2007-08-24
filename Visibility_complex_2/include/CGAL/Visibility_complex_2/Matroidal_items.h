#ifndef VISIBILITY_COMPLEX_2_MATROIDAL_ITEMS_H
#define VISIBILITY_COMPLEX_2_MATROIDAL_ITEMS_H

#ifndef VISIBILITY_COMPLEX_2_ITEMS_H
#include <CGAL/Visibility_complex_2/Items.h>
#endif

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template < class Vc_ >
class Matroidal_vertex
    : public Vertex_base<Vc_>
{
public:
    typedef Matroidal_vertex<Vc_>  Self;
    typedef Vertex_base<Vc_>     Vertex_base;
    typedef typename Vertex_base::Vertex_handle     Vertex_handle;
    typedef typename Vertex_base::Edge_handle       Edge_handle;
    typedef typename Vertex_base::Disk_handle    Disk_handle;
    typedef typename Vertex_base::Bitangent_2       Bitangent_2;
    typedef typename Vertex_base::Type              Type;

    Matroidal_vertex() {}
    Matroidal_vertex(Type t , Disk_handle start , 
				        Disk_handle finish) 
	: Vertex_base(t,start,finish) , phiR_(0) , phiL_(0) {}
    Matroidal_vertex(Edge_handle start , Edge_handle finish)
	: Vertex_base(start,finish) , phiR_(0) , phiL_(0)   {}
    Matroidal_vertex(const Bitangent_2& b) 
	: Vertex_base(b) , phiR_(0) , phiL_(0) {}
    Matroidal_vertex( const Vertex_base& v)   // down cast
        : Vertex_base(v) , phiR_(0) , phiL_(0) {}
/*     Self& operator=( const Self& v) { */
/*         this->Vertex_base::operator=(v); */
/*         phi_R_=v.phi_R_; */
/*         phi_L_=v.phi_L_; */
/*         return *this; */
/*     } */
    Edge_handle phiR() const { return phiR_; }
    void set_phiR(const Edge_handle& e) { phiR_ = e; }
    Edge_handle phiL() const { return phiL_; }
    void set_phiL(const Edge_handle& e) { phiL_ = e; }

    Vertex_handle phiR_vertex() const { return phiR_vertex_; }
    void set_phiR_vertex(const Vertex_handle& v) { phiR_vertex_ = v; }
    Vertex_handle phiL_vertex() const { return phiL_vertex_; }
    void set_phiL_vertex(const Vertex_handle& v) { phiL_vertex_ = v; }

private:
    Edge_handle phiR_,phiL_;
    Vertex_handle phiR_vertex_,phiL_vertex_;
};

// -----------------------------------------------------------------------------

class Matroidal_items 
: public Items {
public:
    template <class VCA>
    struct Vertex_wrapper {
	typedef Matroidal_vertex<VCA> Vertex;
    };
};

// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_MATROIDAL_ITEMS_H
