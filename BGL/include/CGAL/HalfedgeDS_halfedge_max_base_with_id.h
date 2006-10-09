#ifndef CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_halfedge_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class ID>
class HalfedgeDS_halfedge_max_base_with_id : public HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true >
{
public:
    typedef HalfedgeDS_halfedge_base< Refs, Tag_true, Tag_true, Tag_true> Base ;
    
    typedef typename Base::Base_base Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type mID ;
    
public:

    HalfedgeDS_halfedge_max_base_with_id( size_type i = size_type(-1) ) : mID(i) {}
    
    size_type&       id()       { return mID; }
    size_type const& id() const { return mID; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_HALFEDGE_MAX_BASE_WITH_ID_H //
// EOF //
