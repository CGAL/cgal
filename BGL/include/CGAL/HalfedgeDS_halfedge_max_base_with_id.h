#ifndef CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_halfedge_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class ID>
class HalfedgeDS_halfedge_max_base_with_id : public HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true >
{
public:
    typedef HalfedgeDS_halfedge_max_base_with_id< Refs, P, ID> Base;
    
    typedef HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true> Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type id ;
    
public:

    HalfedgeDS_halfedge_max_base_with_id( size_type i = size_type(-1) ) : id(i) {}
    
    size_type&       id()       { return id; }
    size_type const& id() const { return id; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H //
// EOF //
