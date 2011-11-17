#include <QGLViewer/qglviewer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/glu.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;

class Viewer : public QGLViewer
{
  GLUquadricObj *qsphere;
  GLuint dl_nb;
protected :
  virtual void draw();
  virtual void init();
  template <class Kernel>
  void draw_circle_on_unit_sphere(const typename CGAL::Point_3<Kernel>&) const;
  template<class Output_iterator>
  void naive_compute_intersection_points(const std::vector<EPIC::Point_3>&,Output_iterator) const;
};
