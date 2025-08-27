#ifndef HEXMESHING_FUNCTION_ALIAS_H
#define HEXMESHING_FUNCTION_ALIAS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <functional>

namespace CGAL::internal::Hexmeshing {
  using TrimmingFunction = std::function<bool(LCC&, Dart_handle)>;
  // Identifies which 3-cell should be refined
  using MarkingFunction = std::function<bool(LCC&, Dart_handle)>;
  using DetectingFunction = std::function<bool(LCC&, Dart_handle)>;
  using DecideInsideFunction = std::function<bool(Point)>;
}



#endif