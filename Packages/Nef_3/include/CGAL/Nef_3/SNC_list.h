#ifndef CGAL_SNC_LIST_H
#define CGAL_SNC_LIST_H

#include <CGAL/In_place_list.h>

CGAL_BEGIN_NAMESPACE

template < class Sphere_map>
class SNC_in_place_list_sm
    : public Sphere_map, 
      public In_place_list_base<SNC_in_place_list_sm<Sphere_map> > {
public:
    typedef SNC_in_place_list_sm<Sphere_map> Self;
    //    typedef typename Vertex::Vertex_handle       Vertex_handle;
    //    typedef typename Vertex::Vertex_const_handle Vertex_const_handle;
    SNC_in_place_list_sm() {}
    SNC_in_place_list_sm(const Sphere_map& sm)   // down cast
        : Sphere_map(sm) {}
    Self& operator=( const Self& sm) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Sphere_map*)this) = ((const Sphere_map&)sm);
        return *this;
    }
};

template < class Vertex>
class SNC_in_place_list_vertex
    : public Vertex, 
      public In_place_list_base<SNC_in_place_list_vertex<Vertex> > {
public:
    typedef SNC_in_place_list_vertex<Vertex> Self;
    //    typedef typename Vertex::Vertex_handle       Vertex_handle;
    //    typedef typename Vertex::Vertex_const_handle Vertex_const_handle;
    SNC_in_place_list_vertex() {}
    SNC_in_place_list_vertex(const Vertex& v)   // down cast
        : Vertex(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Vertex*)this) = ((const Vertex&)v);
        return *this;
    }
};

template < class Halffacet>
class SNC_in_place_list_halffacet
    : public Halffacet, 
      public In_place_list_base<SNC_in_place_list_halffacet<Halffacet> > {
public:
    typedef SNC_in_place_list_halffacet<Halffacet> Self;
    //    typedef typename Halffacet::Halffacet_handle       Halffacet_handle;
    //    typedef typename Halffacet::Halffacet_const_handle Halffacet_const_handle;
    SNC_in_place_list_halffacet() {}
    SNC_in_place_list_halffacet(const Halffacet& v)   // down cast
        : Halffacet(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Halffacet*)this) = ((const Halffacet&)v);
        return *this;
    }
};

template < class Volume>
class SNC_in_place_list_volume
    : public Volume, 
      public In_place_list_base<SNC_in_place_list_volume<Volume> > {
public:
    typedef SNC_in_place_list_volume<Volume> Self;
    //    typedef typename Volume::Volume_handle       Volume_handle;
    //    typedef typename Volume::Volume_const_handle Volume_const_handle;
    SNC_in_place_list_volume() {}
    SNC_in_place_list_volume(const Volume& v)   // down cast
        : Volume(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Volume*)this) = ((const Volume&)v);
        return *this;
    }
};

template < class SVertex>
class SNC_in_place_list_svertex
    : public SVertex, 
      public In_place_list_base<SNC_in_place_list_svertex<SVertex> > {
public:
    typedef SNC_in_place_list_svertex<SVertex> Self;
    //    typedef typename SVertex::SVertex_handle       SVertex_handle;
    //    typedef typename SVertex::SVertex_const_handle SVertex_const_handle;
    SNC_in_place_list_svertex() {}
    SNC_in_place_list_svertex(const SVertex& v)   // down cast
        : SVertex(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SVertex*)this) = ((const SVertex&)v);
        return *this;
    }
};

template < class SHalfedge>
class SNC_in_place_list_shalfedge
    : public SHalfedge, 
      public In_place_list_base<SNC_in_place_list_shalfedge<SHalfedge> > {
public:
    typedef SNC_in_place_list_shalfedge<SHalfedge> Self;
    //    typedef typename SHalfedge::SHalfedge_handle       SHalfedge_handle;
    //    typedef typename SHalfedge::SHalfedge_const_handle SHalfedge_const_handle;
    SNC_in_place_list_shalfedge() {}
    SNC_in_place_list_shalfedge(const SHalfedge& v)   // down cast
        : SHalfedge(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfedge*)this) = ((const SHalfedge&)v);
        return *this;
    }
};

template < class SHalfloop>
class SNC_in_place_list_shalfloop
    : public SHalfloop, 
      public In_place_list_base<SNC_in_place_list_shalfloop<SHalfloop> > {
public:
    typedef SNC_in_place_list_shalfloop<SHalfloop> Self;
    //    typedef typename SHalfloop::SHalfloop_handle       SHalfloop_handle;
    //    typedef typename SHalfloop::SHalfloop_const_handle SHalfloop_const_handle;
    SNC_in_place_list_shalfloop() {}
    SNC_in_place_list_shalfloop(const SHalfloop& v)   // down cast
        : SHalfloop(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfloop*)this) = ((const SHalfloop&)v);
        return *this;
    }
};

template < class SFace>
class SNC_in_place_list_sface
    : public SFace, 
      public In_place_list_base<SNC_in_place_list_sface<SFace> > {
public:
    typedef SNC_in_place_list_sface<SFace> Self;
    //    typedef typename SFace::SFace_handle       SFace_handle;
    //    typedef typename SFace::SFace_const_handle SFace_const_handle;
    SNC_in_place_list_sface() {}
    SNC_in_place_list_sface(const SFace& v)   // down cast
        : SFace(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SFace*)this) = ((const SFace&)v);
        return *this;
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_SNC_LIST_H
