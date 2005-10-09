#ifndef CGAL_KDS_INFINITY_OR_MAX_H
#define CGAL_KDS_INFINITY_OR_MAX_H

#include <limits>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
template <class T>
T infinity_or_max() {
  if (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
  else return std::numeric_limits<T>::max();
}

template <class T>
T infinity_or_max(T) {
  if (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
  else return std::numeric_limits<T>::max();
}

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
