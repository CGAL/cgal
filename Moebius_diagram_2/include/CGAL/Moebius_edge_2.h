#ifndef CHR_MOEBIUS_EDGE_2_H
#define CHR_MOEBIUS_EDGE_2_H

#include <CGAL/Moebius_utils.h>

CGAL_BEGIN_NAMESPACE


template <class Gt>
class Moebius_edge_2
{
 public:
  typedef typename Gt::Bare_point_2 Point;
  typedef typename Gt::Line_2 Line;
  typedef typename Gt::Circle_2 Circle;
  typedef typename Gt::Segment_2 Segment;
  typedef typename Gt::Ray_2 Ray;
  typedef typename Gt::Direction_2 Direction;
  
  typedef std::list<Point *> List;
  typedef typename List::iterator Iterator;
  typedef typename List::const_iterator Const_iterator;

  //  typedef Gt::Segment_Circle_2 Arc;

  bool are_ordered (const Point &p, const Point &q) const {
    CGAL_precondition (is_line ());
    return _gt.compare_on_oriented_line_2_object () (line (), p, q);
  }

 public:
  Moebius_edge_2 () { CGAL_assertion (false); }

  Moebius_edge_2 (const Line &bisector) : 
    _is_line (true), _start (), _stop (), _gt () {
    _bisec = (void *) new Line (bisector);
  }

  Moebius_edge_2 (const Circle &bisector) : 
    _is_line (false), _start (), _stop (), _gt () {
    _bisec = (void *) new Circle (bisector);
  }

  void start (const Point &p) {
    TRACE ("        start, ");
    if (_start == NULL) {
      TRACE ("saving (" << p << ")\n");
      _start = new Point (p);
    } else {
      CGAL_assertion (_stop != NULL);
      TRACE ("pushing (" << *_stop);
      _points.push_front (_stop);
      //delete _stop;
      _stop = NULL;
      TRACE (") and (" << p << ")\n");
      _points.push_front (new Point (p));
    }
  }

  void stop (const Point &p) {
    TRACE ("        stop,  ");
    if (_stop == NULL) {
      TRACE ("saving (" << p << ")\n");
      _stop = new Point (p);
    } else {
      CGAL_assertion (_start != NULL);
      TRACE ("pushing (" << p);
      _points.push_front (new Point (p));
      TRACE (") and (" << *_start << ")\n");
      _points.push_front (_start);
      //      delete _start;
      _start = NULL;
    }
  }
  
  bool is_line () const { return _is_line; }
  bool is_circle () const { return ! _is_line; }

  const Circle &circle () const { CGAL_assertion (is_circle ()); return *((Circle *) _bisec); }
  const Line &line () const { CGAL_assertion (is_line ()); return *((Line *) _bisec); }


  List points () const { return _points; }
  Const_iterator begin () const { return _points.begin (); }
  Const_iterator end () const { return _points.end (); }

#ifdef CGAL_USE_QT
  void draw (Qt_widget *widget) const {
    *widget << *this;
  }
#endif

  void finalize () {
    if (_start == NULL && _stop == NULL) return;
    TRACE ("    finalize, ");

    if (_stop == NULL) TRACE ("pushing NULL");
    else TRACE ("pushing ("<<*_stop<<")");
    _points.push_front (_stop);

    if (_start == NULL) TRACE (" and NULL\n");
    else TRACE (" and ("<<*_start<<")\n");
    _points.push_front (_start);

    _start = _stop = NULL;
  }

 private:
  bool _is_line;
  Point *_start;
  Point *_stop;
  Gt _gt;
  void *_bisec;
  List _points;
};

#ifdef CGAL_USE_QT

CGAL_END_NAMESPACE
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
CGAL_BEGIN_NAMESPACE

// workaround : in Qt_widget << Segment, ambiguous call to to_double...
template <class Point>
CGAL::Segment_2<CGAL::Cartesian<double> >
my_segment (const Point &source, const Point &target)
{
  typedef CGAL::Segment_2<CGAL::Cartesian<double> > S;
  typedef CGAL::Point_2<CGAL::Cartesian<double> > P;
  double x1 = CGAL::to_double (source.x());
  double y1 = CGAL::to_double (source.y());
  double x2 = CGAL::to_double (target.x());
  double y2 = CGAL::to_double (target.y());
  return S (P (x1, y1), P (x2, y2));
}


template <class Gt>
Qt_widget& operator<<(Qt_widget& widget, const Moebius_edge_2<Gt>& e)
{
  typedef Moebius_edge_2<Gt> Edge;
  typedef typename Edge::Point Point;
  typedef typename Edge::Ray Ray;
  typedef typename Edge::Direction Direction;

  typedef typename Edge::Const_iterator Iterator;
  Iterator it = e.begin ();
  Iterator end = e.end ();

  TRACE ("Moebius_edge << Qt, ");
  if (e.is_line ()) TRACE ("line: " << e.line ());
  else TRACE ("circle: " << e.circle ());
  TRACE ("\n");
  if (it == end) {
    TRACE ("  whole bisector\n");
    if (e.is_line ()) {
      widget << e.line ();
    } else {
      widget << e.circle ();
    }
  }
  
  while (it != end) {
    Point *start = *(it++);
    CGAL_assertion (it != end);
    Point *stop = *(it++);
    
    if (e.is_line ()) {
      Direction dir (e.line ());
      if (start == NULL) {
	CGAL_assertion (stop != NULL);
	TRACE ("  ray ("<<*stop<<") ("<<-dir<<")\n");
	widget << Ray (*stop, -dir);
      } else if (stop == NULL) {
	TRACE ("  ray ("<<*start<<") ("<<dir<<")\n");
	widget << Ray (*start, dir);
      } else if (e.are_ordered (*start, *stop)) {
	TRACE ("  segment ("<<*start<<") ("<<*stop<<")\n");
	widget << my_segment<Point> (*start, *stop);
      } else {
	TRACE ("  two half lines ("<<*start<<") ("<<*stop<<")\n");
	widget
	  << Ray (*start, dir)
	  << Ray (*stop, -dir);
      }
    } else {
      CGAL_assertion (start != NULL);
      CGAL_assertion (stop != NULL);
      TRACE ("  arc ("<<*start<<") -> ("<<*stop<<")\n");
      qt_draw_arc (widget, e.circle (), *start, *stop);
    }
  }
  return widget;
}
#endif

CGAL_END_NAMESPACE
#endif//CHR_MOEBIUS_EDGE_2_H
