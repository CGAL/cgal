#ifndef SVD_DRAW_H
#define SVD_DRAW_H



template<class T, class Widget>
void draw_diagram(Widget& widget, const T& svd)
{
  widget << CGAL::BLUE;
#if !defined (__POWERPC__)
  widget << CGAL::PointSize(3);
  widget << CGAL::LineWidth(3);
#endif
  svd.draw_dual(widget);
}




#endif // SVD_DRAW_H
