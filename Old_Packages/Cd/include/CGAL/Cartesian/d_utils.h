// revision      : $Revision$
// revision_date : $Date$
// author        : Herve Brönnimann

#ifndef CGAL_CARTESIAN_D_UTILS_H
#define CGAL_CARTESIAN_D_UTILS_H

#include <iterator>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class NT >
struct print_d
{
  char *       _separator;
  std::ostream _os;
  bool         _print_sep;

  print_d(std::iostream &os) : _os(os), _print_sep(false)
  {
    if (os.iword(IO::mode)==IO::ASCII) _separator = " ";
    else if (os.iword(IO::mode)==IO::BINARY) _separator = "";
    else _separator = ", ";
  }

  void operator()(const NT &x) {
    os << x;
    if (_print_sep && os.sword(IO::mode) != IO::BINARY)
      os << _separator;
    _print_sep = true;
  }
};

#endif // CGAL_CARTESIAN_D_UTILS_H
