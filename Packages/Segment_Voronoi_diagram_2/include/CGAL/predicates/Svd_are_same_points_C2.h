#ifndef CGAL_SVD_ARE_SAME_POINTS_C2_H
#define CGAL_SVD_ARE_SAME_POINTS_C2_H

CGAL_BEGIN_NAMESPACE

template<class K>
class Svd_are_same_points_C2
{
private:
  typedef typename K::Point_2     Point_2;
  typedef typename K::Segment_2   Segment_2;
  typedef typename K::Site_2      Site_2;

public:
  bool operator()(const Site_2& p, const Site_2& q) const
  {
    return ( p.point().x() == q.point().x() &&
	     p.point().y() == q.point().y() );
  }
};



CGAL_END_NAMESPACE

#endif // CGAL_SVD_ARE_SAME_POINTS_C2_H
