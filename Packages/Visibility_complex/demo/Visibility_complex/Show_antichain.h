#ifndef SHOW_ANTICHAIN_H
#define SHOW_ANTICHAIN_H

#include <list>

#include <CGAL/Polygon_2.h>

#include "Show_lines_base.h"

namespace CGAL {

// Show_antichain is a Qt_widget_styled_layer, because Show_lines_base is
// one.

template <class Antichain>
class Show_antichain : public Show_lines_base {
public:
  typedef typename Antichain::Face_iterator	Face_iterator;
  typedef typename Antichain::Face_handle	Face_handle;
  typedef typename Antichain::Edge_handle	Edge_handle;
  typedef typename Antichain::Vertex Vertex;
  typedef typename Antichain::Vertex_handle Vertex_handle;
  typedef typename Antichain::Gt Gt;
  typedef typename Gt::R Kernel; // WARNING: very dependant from
  // CGAL kernels

  typedef typename Gt::Equal_as_segments Equal_as_segments;

  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef Polygon_2<Kernel> Polygon_2;

  Show_antichain(Antichain* &antichain,
		 Color c=CGAL::YELLOW,
		 int lineWidth = 3,
		 QObject * parent=0, const char * name=0)
    : Show_lines_base(c, lineWidth, parent, name),
      ant(antichain) {};

  void draw()
  {
    widget->lock();

    QColor old_fill_color = widget->fillColor();
    QColor old_color = widget->color();
    int old_line_width = widget->lineWidth();

    QColor c = style()->getColor(color);
    widget->setColor(c);
    widget->setFillColor(c);
    widget->setLineWidth(style()->getInt(width));

    for(Face_iterator it = ant->faces_begin();
	it!=ant->faces_end();
	++it)
      {
	const Vertex_handle& sup = it->sup();
	const Vertex_handle& inf = it->inf();
	const Vertex_handle& top = it->top_edge()->sup();
	const Vertex_handle& bottom = it->bottom_edge()->sup();

	if( sup->type() == Vertex::RL && inf->type() == Vertex::LR )
	  {
	    std::list<Point_2> list;
	    
	    Point_2 point_top;
	    bool point_top_valid = dual(top, point_top);
	    
	    Point_2 point_bottom;
	    bool point_bottom_valid = dual(bottom, point_bottom);

	    Vertex_handle va = top;
	    if( point_top_valid ) list.push_back(point_top);
	    
	    while( va->type() != Vertex::RL )
	      {
		va = va->ccw_target_edge()->sup(); // ccR()

		Point_2 p;
		if( dual(va, p) && 
		    ( !point_top_valid || p.x() >=  point_top.x() ) )
		    list.push_back(p);
	      }
	    CGAL_assertion( va == sup );

	    Point_2 p;
	    if( dual(sup, p) ) list.push_back(p);

	    Polygon_2 poly(list.begin(), list.end());

	    list.clear();

	    va = bottom;
	    if( point_bottom_valid ) list.push_back(point_bottom);
 
	    while( va->type() != Vertex::RL )
	      {
		va = va->ccw_source_edge()->sup(); // ccL()

		Point_2 p;
		if( dual(va, p) && 
		    ( !point_top_valid || p.x() >=  point_bottom.x() ) )
		    list.push_back(p);
	      }
	    CGAL_assertion( va == sup );

	    list.pop_back();
	    // remove the extra "sup" that is counted twice.

	    std::copy(list.rbegin(), list.rend(), std::back_inserter(poly));

	    *widget << poly;
	  }
     }

    widget->setFillColor(old_fill_color);
    widget->setColor(old_color);
    widget->setLineWidth(old_line_width);

    widget->unlock();

   };
private:
  Antichain*	&ant;
};//end class 

} // namespace CGAL

#endif // SHOW_ANTICHAIN_H
