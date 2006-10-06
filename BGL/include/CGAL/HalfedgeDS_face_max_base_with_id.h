#ifndef CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H
#define CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H 1

#include <CGAL/HalfedgeDS_vertex_base.h>

CGAL_BEGIN_NAMESPACE

template < class Refs, class Pln, class ID>
class HalfedgeDS_face_max_base_with_id : public HalfedgeDS_face_base< Refs, Tag_true, Pln>
{
public:
    typedef HalfedgeDS_face_max_base_with_id< Refs, Pln, ID> Base;
    
    typedef HalfedgeDS_face_base< Refs, Tag_true, Pln> Base_base ;
    
    typedef ID size_type ;
    
private:

    size_type id ;
    
public:

    HalfedgeDS_face_max_base_with_id() : id ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p) : Base_base(p), id ( size_type(-1) ) {}
    HalfedgeDS_face_max_base_with_id( Pln const& p, size_type i ) : Base_base(p), id (i) {}
    
    size_type&       id()       { return id; }
    size_type const& id() const { return id; }
};

CGAL_END_NAMESPACE

#endif // CGAL_HALFEDGEDS_FACE_MAX_BASE_WITH_ID_H //
// EOF //
