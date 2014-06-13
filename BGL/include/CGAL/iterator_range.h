#ifndef CGAL_ITERATOR_RANGE_H
#define CGAL_ITERATOR_RANGE_H

#include <CGAL/tuple.h>
#include <boost/range/irange.hpp>
#include <utility>

namespace CGAL {

  template <typename I>
  class iterator_range
    : public std::pair<I,I>{
    
    typedef std::pair<I,I> Base;

  public:

    typedef I iterator;
    typedef const I const_iterator;

    iterator_range(I b, I e)
      : Base(b,e)
    {}

    iterator_range(const std::pair<I,I>& ip)
      : Base(ip)
    {}

    operator Base() const
    {
      return std::make_pair(begin(),end());
    }

  const I& begin() const
  {
    return first;
  }

  const I& end() const
  {
    return second;
  }
};

} // namespace CGAL

#endif // CGAL_ITERATOR_RANGE_H
