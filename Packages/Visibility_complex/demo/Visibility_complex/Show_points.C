#include "Show_points.h"

namespace CGAL {

  void Show_points_base::setColor(QColor c) 
  { 
    color = Color(c.red(),c.green(),c.blue()); 
  }
  
  void Show_points_base::setPointSize(int pointsize)
  { size = pointsize; }

  void Show_points_base::setPointStyle(PointStyle pointstyle)
  { style = pointstyle; }

}

#include "Show_points.moc"
