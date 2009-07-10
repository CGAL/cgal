#ifndef CGAL_PDB_BASIC_H
#define CGAL_PDB_BASIC_H
#include <CGAL/basic.h>
#include <CGAL/auto_link/PDB.h>
#include <CGAL/Tools/utility_macros.h>
#include <CGAL/Tools/Log.h>


#include <boost/version.hpp>
#if BOOST_VERSION < 103400
#include <CGAL/PDB/internal/foreach.h>
#else
#include <boost/foreach.hpp>
#define CGAL_PDB_FOREACH(a,b) BOOST_FOREACH(a,b)
#endif

namespace CGAL {
//! All exposed classes are put in this namespace.
namespace PDB {
}
}


namespace internal {
  template <class T>
  struct EchoType {
    typedef T type;
  };
  template <class T>
  EchoType<T> get_type(T){return EchoType<T>();}
}


#endif
