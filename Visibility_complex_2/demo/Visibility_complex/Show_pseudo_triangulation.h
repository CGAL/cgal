#ifndef SHOW_PSEUDO_TRIANGULATION_H
#define SHOW_PSEUDO_TRIANGULATION_H

#include "Qt_widget_styled_layer.h"

using std::endl;

namespace CGAL {

class Show_pseudo_triangulation_base : public Qt_widget_styled_layer {
  Q_OBJECT
public:
  typedef Qt_widget_styled_layer::Style Style;

  Show_pseudo_triangulation_base(Color c,
				 int size,
				 QObject * parent=0, const char * name=0);
  
  Show_pseudo_triangulation_base(Style* style,
				 QString color_name,
				 QString size_name,
				 QObject * parent=0, const char * name=0);

public slots:
  void setColor(QColor color);
  void setSize(int size);

protected:
  QString color;
  QString size;
}; // end Show_pseudo_triangulation_base  

template <class Antichain>
class Show_pseudo_triangulation : public Show_pseudo_triangulation_base {
public:
  typedef Qt_widget_styled_layer::Style Style;

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
			    int width=1,
			    QObject * parent=0, const char * name=0)
    : Show_pseudo_triangulation_base(c, width, parent, name),
      ant(antichain), _dual(in_dual) {};

  Show_pseudo_triangulation(Antichain* &antichain,
			    bool in_dual,
			    Style* style,
			    QString color_name,
			    QString size_name,
			    QObject * parent=0, const char * name=0)
    : Show_pseudo_triangulation_base(style, color_name, size_name,
				     parent, name),
      ant(antichain), _dual(in_dual) {};

  void draw()
  {
    widget->lock();
    {
      old_color = widget->color();
      if(_dual)
	old_size = widget->pointSize();
      else
	old_size = widget->lineWidth();

      widget->setColor(style()->getColor(color));
      if(_dual)
	widget->setPointSize(style()->getInt(size));
      else
	widget->setLineWidth(style()->getInt(size));
   
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
    }
    widget->unlock();
  };
private:
  Antichain*	&ant;
  bool _dual;

  QColor old_color;
  int old_size;
};//end class 

} // namespace CGAL

#endif // SHOW_PSEUDO_TRIANGULATION_H
