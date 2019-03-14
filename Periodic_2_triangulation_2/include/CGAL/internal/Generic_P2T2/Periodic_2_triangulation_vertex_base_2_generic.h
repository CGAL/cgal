#ifndef CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_GENERIC_H
#define CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL
{


template < class Gt,
           class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Periodic_2_triangulation_vertex_base_2_generic
  : public Vb
{
  typedef Vb                                            Base;
  typedef typename Vb::Triangulation_data_structure     Tds;

public:
  typedef Gt                                            Geom_traits;
  typedef Tds                                           Triangulation_data_structure;

  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Face_handle                     Face_handle;
  typedef typename Gt::Point_2                          Point;
  typedef Periodic_2_offset_2                           Offset; // @fixme should be given in the traits

  template < typename Tds2 >
  struct Rebind_TDS
  {
    typedef typename Vb::template Rebind_TDS<Tds2>::Other      Vb2;
    typedef Periodic_2_triangulation_vertex_base_2_generic<Gt, Vb2>    Other;
  };

public:
  Periodic_2_triangulation_vertex_base_2_generic() : Base() {}
  Periodic_2_triangulation_vertex_base_2_generic(const Point & p)
    : Base(p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_vertex_base_2_generic(const Point & p, Face_handle f)
    : Base(f, p), _off(), _offset_flag(false) {}
  Periodic_2_triangulation_vertex_base_2_generic(Face_handle f)
    : Base(f), _off(), _offset_flag(false) {}

  const Offset& offset() const
  {
    return _off;
  }

  void set_offset(const Offset& off)
  {
    _off = off;
    _offset_flag = true;
  }

  void clear_offset()
  {
    _offset_flag = false;
    _off = Offset();
  }

  bool get_offset_flag() const
  {
    return _offset_flag;
  }

private:
  /// The offset is needed to be able to copy a triangulation that is
  /// not on the 1-cover.

  /// Normal copying of the vertices would give multiple vertices with
  /// the same location (the periodic copies) and we wouldn't be able
  /// to distinguish them anymore.
  Offset _off;
  /// The flag is used to test whether _off has been set.
  bool _offset_flag;
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_vertex_base_2_generic<Tds> &)
// no combinatorial information.
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os,
           const Periodic_2_triangulation_vertex_base_2_generic<Tds> &)
// no combinatorial information.
{
  return os;
}

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_VERTEX_BASE_2_GENERIC_H
