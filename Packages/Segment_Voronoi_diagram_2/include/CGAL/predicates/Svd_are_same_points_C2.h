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
  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;

  Compare_x_2 compare_x_2;
  Compare_y_2 compare_y_2;

  bool are_same(const Point_2& p, const Point_2& q) const
  {
    return
      compare_x_2(p, q) == EQUAL && compare_y_2(p, q) == EQUAL;
  }

public:
  bool operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    if ( !p.is_exact() && !q.is_exact() ) {
      Segment_2 p_supp[2] = { p.supporting_segment(0),
			      p.supporting_segment(1) };
      Segment_2 q_supp[2] = { q.supporting_segment(0),
			      q.supporting_segment(1) };
      bool b[2][2];
      for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	  b[i][j] =
	    (  are_same( p_supp[i].source(), q_supp[j].source() ) &&
	       are_same( p_supp[i].target(), q_supp[j].target() )  ) ||
	    (  are_same( p_supp[i].source(), q_supp[j].target() ) &&
	       are_same( p_supp[i].target(), q_supp[j].source() )  );
	}
      }
      return (b[0][0] && b[1][1]) || (b[0][1] && b[1][0]);
    }

    return are_same(p.point(), q.point());
  }
};



CGAL_END_NAMESPACE

#endif // CGAL_SVD_ARE_SAME_POINTS_C2_H
