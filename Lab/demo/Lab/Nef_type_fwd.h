#ifndef NEF_TYPE_FWD_H
#define NEF_TYPE_FWD_H

#ifdef USE_FORWARD_DECL

struct Exact_Kernel;

namespace CGAL {

  class SNC_indexed_items;

  template <typename Kernel_, typename Items_, typename Mark_>
  class Nef_polyhedron_3 ;
}

typedef CGAL::Nef_polyhedron_3<Exact_Kernel,
                               CGAL::SNC_indexed_items,
                               bool> Nef_polyhedron;

#else //  USE_FORWARD_DECL
#include "Nef_type.h"
#endif // USE_FORWARD_DECL

#endif // NEF_TYPE_FWD_H
