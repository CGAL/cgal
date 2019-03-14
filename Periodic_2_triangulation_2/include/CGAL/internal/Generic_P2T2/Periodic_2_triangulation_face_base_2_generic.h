#ifndef CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_GENERIC_H
#define CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Dummy_tds_2.h>

#include <CGAL/Periodic_2_offset_2.h>

namespace CGAL
{

template < typename Gt,
           typename Fb = Triangulation_face_base_2<Gt> >
class Periodic_2_triangulation_face_base_2_generic
  : public Fb
{
  typedef Fb                                            Base;
  typedef typename Base::Triangulation_data_structure   Tds;

public:
  typedef Gt                                            Geom_traits;
  typedef Tds                                           Triangulation_data_structure;
  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Face_handle                     Face_handle;
  typedef Periodic_2_offset_2                           Offset;

  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Periodic_2_triangulation_face_base_2_generic<Gt, Fb2>  Other;
  };

public:
  Periodic_2_triangulation_face_base_2_generic() : Fb(), is_canonical(false) { }

  Periodic_2_triangulation_face_base_2_generic(Vertex_handle v0,
                                               Vertex_handle v1,
                                               Vertex_handle v2)
    : Fb(v0, v1, v2), is_canonical(false)
  { }

  Periodic_2_triangulation_face_base_2_generic(Vertex_handle v0,
                                               Vertex_handle v1,
                                               Vertex_handle v2,
                                               Face_handle n0,
                                               Face_handle n1,
                                               Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2), is_canonical(false)
  { }

  /// Periodic functions
  Offset offset(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i < 3 );
    return _off[i];
  }

  bool has_zero_offsets() const
  {
    return (_off[0] == Offset(0,0) &&
            _off[1] == Offset(0,0) &&
            _off[2] == Offset(0,0));
  }

  void set_offsets(const Offset& o0, const Offset& o1, const Offset& o2)
  {
    _off[0] = o0;
    _off[1] = o1;
    _off[2] = o2;
  }

  void set_canonical_flag(const bool b) { is_canonical = b; }
  bool get_canonical_flag() const { return is_canonical; }

private:
  bool is_canonical;
  cpp11::array<Offset, 3> _off;
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_face_base_2_generic<Tds> &)
// non combinatorial information. Default = nothing
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os, const Periodic_2_triangulation_face_base_2_generic<Tds> &)
// non combinatorial information. Default = nothing
{
  return os;
}

} //namespace CGAL

#endif //CGAL_PERIODIC_2_TRIANGULATION_FACE_BASE_2_GENERIC_H
