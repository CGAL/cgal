namespace CGAL {

template <class NumberType>
struct Cartesian_double_to_Homogeneous
{
  typedef Point_2< Homogeneous< NumberType> >    Point;
  typedef Segment_2< Homogeneous< NumberType> >  Segment;

  Point
  operator()(  const Point_2<Cartesian<double> >& p)
  { return Point( NumberType(p.x()), NumberType(p.y()) ); }

  Segment
  operator()(  const Segment_2<Cartesian<double> >& s)
  { 
    return Segment( Point( NumberType(s.source().x()), 
                           NumberType(s.source().y()) ),
                    Point( NumberType(s.target().x()), 
                           NumberType(s.target().y()) ) );
  }
};

} // namespace CGAL

