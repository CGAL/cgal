#ifndef CGAL_MULTIPASS_DISTANCE_H
#define CGAL_MULTIPASS_DISTANCE_H

#include <iterator>

namespace CGAL {


  template <typename I, typename Cat>
  typename std::iterator_traits<I>::difference_type
  multipass_distance(const I& begin, const I& beyond, const Cat&)
  {
    return std::distance(begin, beyond);
  }

  template <typename I>
  typename std::iterator_traits<I>::difference_type
  multipass_distance(const I& begin, const I& beyond, const std::input_iterator_tag&)
  {
    return 0;
  }


  template <typename I>
  typename std::iterator_traits<I>::difference_type
  multipass_distance(const I& begin, const I& beyond)
  {
    return multipass_distance(begin, beyond, typename std::iterator_traits<I>::iterator_category());
  }


}


#endif // CGAL_MULTIPASS_DISTANCE_H
