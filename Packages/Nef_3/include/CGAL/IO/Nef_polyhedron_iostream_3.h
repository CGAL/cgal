#ifndef CGAL_NEF_POLYHEDRON_IOSTREAM_3_H
#define CGAL_NEF_POLYHEDRON_IOSTREAM_3_H

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_io_parser.h>

CGAL_BEGIN_NAMESPACE

template <typename Kernel, typename Items>
std::ostream& operator<<
 (std::ostream& os, Nef_polyhedron_3<Kernel,Items>& NP)
{
  typedef typename Nef_polyhedron_3<Kernel,Items>::SNC_structure SNC_structure;
  CGAL::SNC_io_parser<SNC_structure> O(os, NP.snc(),true);
  O.print();
  return os;
}

template <typename Kernel, typename Items>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_3<Kernel,Items>& NP)
{
  typedef typename Nef_polyhedron_3<Kernel,Items>::SNC_structure SNC_structure;
  CGAL::SNC_io_parser<SNC_structure> I(is, NP.snc());
  I.read();
  NP.pl()->initialize(&NP.snc());
  return is;
}

CGAL_END_NAMESPACE

#endif //CGAL_NEF_POLYHEDRON_IOSTREAM_3_H
