#ifndef _COMPARE_CURVE_LIST_H
#define _COMPARE_CURVE_LIST_H

#include <list>
#include <algorithm>


template <class Traits>
class Equal_pred
{
public:
  typedef typename Traits::Point_2                Point_2;
  typedef typename Traits::X_monotone_curve_2     X_monotone_curve_2;
  
  bool operator()(const Point_2& p1, const Point_2& p2)
  {
    return(Traits().equal_2_object()(p1, p2));
  }

  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  {
    return(Traits().equal_2_object()(c1, c2));
  }
};


template <class List, class Traits>
  bool compare_lists(const List& list1, const List& list2, Traits& /*tr*/)
{
  typedef typename List::const_iterator  Iter;
  Iter begin1 = list1.begin();
  Iter end1 = list1.end();

  Iter begin2 = list2.begin();
  Iter end2 = list2.end();
  unsigned int count1 = std::distance(begin1, end1);
  unsigned int count2 = std::distance(begin2, end2);

  if(count1 != count2)
  {
    std::cout << "The lists are not of the same lengths ("
              << count1 << "," << count2 << ")\n";
    return false;
  }

  Equal_pred<Traits> eq;
  return std::equal(begin1, end1, begin2, eq);
}


#endif
