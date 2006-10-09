#ifndef CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_vertex_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class Pln, class ID>
class HalfedgeDS_face_max_base_with_id : public HalfedgeDS_face_base< Refs, Tag_true, Pln>
{
public:
    
    typedef HalfedgeDS_face_base< Refs, Tag_true, Pln> Base ;
    
    typedef ID size_type ;
    
private:

    size_type mID ;
    
public:

    HalfedgeDS_face_max_base_with_id() : mID ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p) : Base(p), mID ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p, size_type i ) : Base(p), mID (i) {}
    
    size_type&       id()       { return mID; }
    size_type const& id() const { return mID; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H //
// EOF //
