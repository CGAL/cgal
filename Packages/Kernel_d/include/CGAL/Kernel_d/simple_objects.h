#ifndef SIMPLE_OBJECTS_H
#define SIMPLE_OBJECTS_H

CGAL_BEGIN_NAMESPACE

template <class R>
struct Lt_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) == SMALLER; }
};

template <class R>
struct Le_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) != LARGER; }
};

template <class R>
struct Eq_from_method {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { return p1 == p2; }
};

CGAL_END_NAMESPACE
#endif //SIMPLE_OBJECTS_H
