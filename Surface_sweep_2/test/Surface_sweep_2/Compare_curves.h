#ifndef _COMPARE_CURVES_H
#define _COMPARE_CURVES_H

#include <list>
#include <algorithm>

template <typename Traits>
class Equal_pred {
public:
  using Point_2 = typename Traits::Point_2;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;

  Equal_pred(const Traits& traits) : m_traits(traits) {}

  bool operator()(const Point_2& p1, const Point_2& p2)
  { return(m_traits.equal_2_object()(p1, p2)); }

  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  { return(m_traits.equal_2_object()(c1, c2)); }

private:
  const Traits& m_traits;
};


template <typename List, typename Traits>
bool compare_lists(const List& list1, const List& list2, Traits& traits) {
  if(! (list1.size() == list2.size())) {
    std::cerr << "Error: The lists are not of the same lengths ("
              << list1.size() << "," << list2.size() << ")\n";
    return false;
  }

  Equal_pred<Traits> eq(traits);
  auto rc = std::equal(list1.begin(), list1.end(), list2.begin(), eq);
  if (! rc) {
    std::cerr << "Error: The lists do not match\n";
    return false;
  }

  return true;
}


#endif
