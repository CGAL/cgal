#include "Show_pseudo_triangulation.h"

namespace CGAL {

  Show_pseudo_triangulation_base::
  Show_pseudo_triangulation_base(Color c,
				 int size,
				 QObject * parent, const char * name)
    : Qt_widget_styled_layer(0, parent, name)
  {
    color=tr("Color");
    this->size=tr("Size");
    
    setColor(QColor(c.red(), c.green(), c.blue()));
    setSize(size);
  }

  Show_pseudo_triangulation_base::
  Show_pseudo_triangulation_base(Style* style,
				 QString color_name,
				 QString size_name,
				 QObject * parent,
				 const char * name)
    : Qt_widget_styled_layer(style, parent, name),
      color(color_name),
      size(size_name)
  {}

  void Show_pseudo_triangulation_base::setColor(QColor c) 
  { style()->setColor(color, c); }
  
  void Show_pseudo_triangulation_base::setSize(int size)
  { style()->setInt(this->size, size); }

} // namespace CGAL

// moc_source_file: Show_pseudo_triangulation.h
#include "Show_pseudo_triangulation.moc"
