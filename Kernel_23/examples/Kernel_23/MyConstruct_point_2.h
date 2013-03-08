#ifndef MYCONSTRUCT_POINT_2_H
#define MYCONSTRUCT_POINT_2_H

template <typename K, typename OldK>
class MyConstruct_point_2
{
  typedef typename K::RT         RT;
  typedef typename K::Point_2    Point_2;
  typedef typename K::Line_2     Line_2;
  typedef typename Point_2::Rep  Rep;
public:
  typedef Point_2                result_type;

  // Note : the CGAL::Return_base_tag is really internal CGAL stuff.
  // Unfortunately it is needed for optimizing away copy-constructions,
  // due to current lack of delegating constructors in the C++ standard.
  Rep // Point_2
  operator()(CGAL::Return_base_tag, CGAL::Origin o) const
  { return Rep(o); }

  Rep // Point_2
  operator()(CGAL::Return_base_tag, const RT& x, const RT& y) const
  { return Rep(x, y); }

  Rep // Point_2
  operator()(CGAL::Return_base_tag, const RT& x, const RT& y, const RT& w) const
  { return Rep(x, y, w); }

  Point_2
  operator()(const CGAL::Origin&) const
  { return MyPointC2(0, 0, 0); }

  Point_2
  operator()(const RT& x, const RT& y) const
  {
    return MyPointC2(x, y, 0);
  }

  Point_2
  operator()(const Line_2& l) const
  {
    typename OldK::Construct_point_2 base_operator;
    Point_2 p = base_operator(l);
    return p;
  }

  Point_2
  operator()(const Line_2& l, int i) const
  {
    typename OldK::Construct_point_2 base_operator;
    return base_operator(l, i);
  }

  // We need this one, as such a functor is in the Filtered_kernel
  Point_2
  operator()(const RT& x, const RT& y, const RT& w) const
  {
    if(w != 1){
      return MyPointC2(x/w, y/w, 0);
    } else {
      return MyPointC2(x,y, 0);
    }
  }
};

#endif //MYCONSTRUCT_POINT_2_H
