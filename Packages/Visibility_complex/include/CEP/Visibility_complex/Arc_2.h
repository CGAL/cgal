#ifndef CGAL_CONVEX_ARC_2
#define CGAL_CONVEX_ARC_2

#include <list>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
// Abstract unsigned convex arc

template < class D_ >
class Arc_base 
{
public:
    // -------------------------------------------------------------------------
    typedef D_                Disk;
    typedef const Disk*   Disk_handle;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    Arc_base() : _object(0) { }
    Arc_base(Disk_handle P) : _object(P) { }
    // -------------------------------------------------------------------------
    Disk_handle   object()               const { return _object; }
    void             set_object(Disk_handle P) { _object = P;    }
    // -------------------------------------------------------------------------
private:
    Disk_handle _object;
};

// -----------------------------------------------------------------------------
// --------------------- Arc_2 general definition -----------------------
// ----- This definition is sufficient if arcs have constant complexity --------
// -----------------------------------------------------------------------------

template < class D_ >
struct Arc_2 : public Arc_base<D_>
{
    // -------------------------------------------------------------------------
    typedef D_                                Disk;
    typedef typename Disk::R                  R;
    typedef typename R::FT                    FT;
    typedef typename R::Point_2               Point_2;
    typedef typename R::Segment_2             Segment_2;
    typedef Arc_base<D_>                      Base;
    typedef typename Base::Disk_handle        Disk_handle;
    // -------------------------------------------------------------------------
    Arc_2() : Base(0) { }
    Arc_2(Disk_handle P) : Base(P) { }
    Arc_2(Disk_handle P,const Point_2& p, const Point_2& q) 
	: Base(P) { }
    // -------------------------------------------------------------------------
  void split (Arc_2& /*tmp*/, const Point_2& /*p*/) {  } 
  void split_cw(Arc_2& /*tmp*/, const Point_2& /*p*/) { }
  void update_begin(const Point_2& /*p*/) { }
  void update_end(const Point_2& /*p*/) { }
  void join (Arc_2& /*y*/) { }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------
// -------------- Arc_2 Specialization for CGAL::Polygon_2 --------------
// -----------------------------------------------------------------------------

#ifdef CGAL_POLYGON_2_H
template < class R_ , class C_ >
class Arc_2 < Polygon_2<R_,C_> >
      : public Arc_base< Polygon_2<R_,C_> >
{
public:
    // -------------------------------------------------------------------------
    typedef R_                                      R;
    typedef typename R_::FT                         FT;
    typedef typename R::Point_2                     Point_2;
    typedef Polygon_2<R_,C_>                        Disk;
    typedef C_                                      Container;
    typedef Arc_base<Disk>                          Base;
    typedef typename Base::Disk_handle              Disk_handle;
    typedef typename Disk::Vertex_const_circulator  Vertex_const_iterator;
    typedef typename Disk::Vertex_const_circulator  Vertex_iterator;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    Arc_2() : Base(0) { }
    Arc_2(Disk_handle P) 
	: Base(P) , first(Vertex_iterator()) , beyond(Vertex_iterator()) { }
    // -------------------------------------------------------------------------
    void split (Arc_2& tmp, const Point_2& p)
    {
	if (first == CGAL_CIRC_NULL) {
	    first = object()->vertices_circulator();
	    while (*first != p) ++first;
	    beyond = first; //++beyond;
	    tmp.first = first; tmp.beyond = beyond;
	}
	else {
	    tmp.beyond = beyond;
	    tmp.first  = first;
	    while (*tmp.first != p) ++tmp.first;
	    beyond = tmp.first; 
	    ++beyond;
	}
    } 
    void split_cw(Arc_2& tmp, const Point_2& p) 
    {
	if (first == CGAL_CIRC_NULL) {
	    first = object()->vertices_circulator();
	    while (*first != p) ++first;
	    beyond = first; //++beyond;
	    tmp.first = first; tmp.beyond = beyond;
	}
	else {
	    tmp.beyond = beyond;
	    tmp.first  = first;
	    while (*tmp.beyond != p) --tmp.beyond;
	    first = tmp.beyond;
	    ++tmp.beyond;
	}
    }
    void update_begin(const Point_2& p) {
	first = beyond; --first;
	while (*first != p) --first;
    }
    void update_end(const Point_2& p) {
	beyond = first;
	while (*beyond != p) ++beyond;
	++beyond;
    }
    void join (Arc_2& y) { beyond  = y.beyond; }
    // -------------------------------------------------------------------------
    // these methods are specific to this specialization
    Vertex_iterator begin()          const { return first;  }
    Vertex_iterator end()            const { return beyond; }
    Vertex_iterator vertices_begin() const { return first;  }
    Vertex_iterator vertices_end()   const { return beyond; }
    // -------------------------------------------------------------------------
private:
    Vertex_iterator first;
    Vertex_iterator beyond;
};
#endif // CGAL_POLYGON_2_H

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
