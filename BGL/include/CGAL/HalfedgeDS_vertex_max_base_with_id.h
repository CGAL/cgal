#ifndef CGAL_HALFEDGEDS_VERTEX_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_VERTEX_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_vertex_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class P, class ID>
class HalfedgeDS_vertex_max_base_with_id : public HalfedgeDS_vertex_base< Refs, Tag_true, P>
{
public:
    typedef HalfedgeDS_vertex_max_base_with_id< Refs, P, ID> Base;
    
    typedef HalfedgeDS_vertex_base< Refs, Tag_true, P> Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type id ;
    
public:

    HalfedgeDS_vertex_max_base_with_id() : id ( size_type(-1) )  {}
    HalfedgeDS_vertex_max_base_with_id( Point const& p) : Base_base(p), id ( size_type(-1) ) {}
    HalfedgeDS_vertex_max_base_with_id( Point const& p, size_type i ) : Base_base(p), id(i) {}
    
    size_type&       id()       { return id; }
    size_type const& id() const { return id; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_VERTEX_MAX_BASE_WITH_ID_H //
// EOF //
