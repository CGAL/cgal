#ifndef SHOW_PSEUDO_TRIANGULATION_H
#define SHOW_PSEUDO_TRIANGULATION_H

#include <CGAL/IO/Qt_widget_layer.h>

using std::endl;

namespace CGAL {

template <class Antichain>
class Show_pseudo_triangulation : public Qt_widget_layer {
public:
  typedef typename Antichain::Vertex_iterator	Vertex_iterator;
  typedef typename Antichain::Vertex Vertex;
  typedef typename Antichain::Gt Gt;
  typedef typename Gt::R Kernel; // WARNING: very dependant from
  // CGAL kernels
  typedef typename Vertex::Bitangent_2 Bitangent;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;

  Show_pseudo_triangulation(Antichain* &antichain,
			    bool in_dual=false,
			    Color c=CGAL::BLUE,
			    int width=1)
    : ant(antichain), _dual(in_dual), color(c), _width(width) {};

  void draw()
  {
    *widget << color;
    if(_dual) *widget << CGAL::PointSize(_width);
    else *widget << CGAL::LineWidth(_width);

    for(Vertex_iterator it = ant->vertices_begin();
	it!=ant->vertices_end();
	++it)
      if(_dual)
	{
	  Line_2 l = it->supporting_line();
	  if(l.b() == 0) return;

	  *widget << dual(l);
	}
      else
	*widget << *it;
	  
  };
private:
  Antichain*	&ant;
  bool _dual;
  Color color;
  int _width;
};//end class 

} // namespace CGAL

#endif // SHOW_PSEUDO_TRIANGULATION_H
