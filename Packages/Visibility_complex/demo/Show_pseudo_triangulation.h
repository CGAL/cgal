#ifndef SHOW_PSEUDO_TRIANGULATION_H
#define SHOW_PSEUDO_TRIANGULATION_H

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {

template <class Antichain>
class Show_pseudo_triangulation : public Qt_widget_layer {
public:
  typedef typename Antichain::Vertex_iterator	Vertex_iterator;
  typedef typename Antichain::Vertex Vertex;
  typedef typename Vertex::Bitangent_2 Bitangent;
  typedef typename Bitangent::R R; // WARNING: very dependant from
  // CGAL kernels
  typedef typename R::Line_2 Line_2;
  typedef typename R::Point_2 Point_2;

  Show_pseudo_triangulation(Antichain* &antichain,
			    bool in_dual=false,
			    Color c=CGAL::BLUE,
			    int linewidth=1)
    : ant(antichain), _dual(in_dual), color(c), width(linewidth) {};

  void draw()
  {
    *widget << color;
    if(_dual) *widget << CGAL::PointSize(width);
    else *widget << CGAL::LineWidth(width);

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
	*widget << (*it);
  };
private:
  Antichain*	&ant;
  bool _dual;
  Color color;
  int width;
};//end class 

} // namespace CGAL

#endif // SHOW_PSEUDO_TRIANGULATION_H
