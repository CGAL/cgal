#include "Show_lines_base.h"

namespace CGAL {

  Show_lines_base::Show_lines_base(Color c,
				   int linewidth,
				   QObject* parent,
				   const char* name)
    : Qt_widget_styled_layer(0, parent, name)
  {
    color=tr("Color");
    width=tr("Point size");

    setColor(QColor(c.red(), c.green(), c.blue()));
    setLineWidth(linewidth);
  }

  Show_lines_base::Show_lines_base(Style* style,
				   QString line_color_name,
				   QString line_width_name,
				   QObject* parent,
				   const char* name)
    : Qt_widget_styled_layer(style, parent, name),
      color(line_color_name),
      width(line_width_name)
  {}

  void Show_lines_base::setColor(QColor c)
  { style()->setColor(color, c); }

  void Show_lines_base::setLineWidth(int line_width)
  { style()->setInt(width, line_width); }

} // namespace CGAL

// moc_source_file: Show_lines_base.h
#include "Show_lines_base.moc"
