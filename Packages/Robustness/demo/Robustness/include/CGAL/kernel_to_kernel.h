#include <LEDA/rat_point.h>

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


template <class NumberType>
struct Cartesian_double_to_Cartesian
{
  typedef Point_2< Cartesian< NumberType> >    Point;
  typedef Segment_2< Cartesian< NumberType> >  Segment;

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

template <class NumberType>
struct Cartesian_float_to_Cartesian
{
  typedef Point_2< Cartesian< NumberType> >    Point;
  typedef Segment_2< Cartesian< NumberType> >  Segment;

  Point
  operator()(  const Point_2<Cartesian<float> >& p)
  { return Point( NumberType(p.x()), NumberType(p.y()) ); }

  Segment
  operator()(  const Segment_2<Cartesian<float> >& s)
  { 
    return Segment( Point( NumberType(s.source().x()), 
                           NumberType(s.source().y()) ),
                    Point( NumberType(s.target().x()), 
                           NumberType(s.target().y()) ) );
  }
};


struct Cartesian_double_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< double> >& s)
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return 
    Segment( Point(to_double(rs.X()), to_double(rs.Y()), to_double(rs.W())), 
             Point(to_double(rt.X()), to_double(rt.Y()), to_double(rt.W())) );
  }
};

struct Cartesian_float_to_H_double_int
{
  typedef Point_2< Homogeneous< double> >    Point;
  typedef Segment_2< Homogeneous< double> >  Segment;

  Segment
  operator()(  const Segment_2< Cartesian< float> >& s)
  {
    leda_rat_point rs =  leda_point(s.source().x(), s.source().y());
    leda_rat_point rt =  leda_point(s.target().x(), s.target().y());

    return 
    Segment( Point(to_double(rs.X()), to_double(rs.Y()), to_double(rs.W())), 
             Point(to_double(rt.X()), to_double(rt.Y()), to_double(rt.W())) );
  }
};

} // namespace CGAL

