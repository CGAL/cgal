#ifndef CGAL_QT_WIDGET_OPENGL_H
#define CGAL_QT_WIDGET_OPENGL_H

// #include <CGAL/Nef_3/OGL_helper.h>
// #include <CGAL/Nef_S2/Sphere_geometry_OGL.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_S2/OGL_base_object.h>

#include <qapplication.h>
#include <qgl.h>
#include <qwidget.h>
#include <qevent.h>
#include <qpopupmenu.h>
#include <qpoint.h>

CGAL_BEGIN_NAMESPACE

class Qt_widget_OpenGL : public QGLWidget {
  Q_OBJECT

    typedef CGAL::OGL::OGL_base_object::Affine_3       Affine_3;
    typedef CGAL::OGL::OGL_base_object::Double_vector  Double_vector;
    //    typedef CGAL::Simple_cartesian<double>::Aff_transformation_3  Affine_3;

protected:
 enum MenuEntriesS2 { ROTATE, SCALE, TRANSLATE, TRANS_Z, RESET_CONTROL, 
		      PERSP, FULLSCREEN, QUIT };

  CGAL::OGL::OGL_base_object* object_;

  int window_width;          // Breite und
  int window_height;         // Hoehe des Fensters
  int window_radius;           // min(width,height) / 2

  bool perspective;
  bool fullscreen;

  int mouse_x, mouse_y;      // Mauskoordinaten linker button
  int interaction;                   // type of interaction in motion fct.
  int motion_mode;           // Bewegen der Maus bei Mouse1 gedrueckt
  int submenu1, submenu2;
  long double dx;                // Translation
  long double dy;                // Translation
  long double dz;                     // Translation in Z
  long double s;                 // Skalierung
  long double init_s;

  Affine_3 rotation;   // Rotation

  long double factor_s;              // Umrechnungsfaktor fuer Skalierung
  long double factor_w;
  long double factor_d;

  QPopupMenu* main;
  QPopupMenu* sub1;
  QPopupMenu* sub2;
  QPopupMenu* sub3;

protected:
  Affine_3 virtual_sphere_transformation( double old_x, double old_y, 
					  double new_x, double new_y);

protected slots:
  void slotControlMenu(int index);
  void slotRenderMenu(int index);
  void slotOptionsMenu(int index);
  void slotFullscreen();
  void slotPerspective();

protected:
  virtual void mouseMoveEvent( QMouseEvent* event);
  virtual void mousePressEvent( QMouseEvent* event);
  virtual void mouseReleaseEvent( QMouseEvent* event);

public:
  Qt_widget_OpenGL(int width, int height, double scale);
  virtual void paintGL();
  virtual void initializeGL();
  virtual void resizeGL(int w, int h);
};

CGAL_END_NAMESPACE
#endif // CGAL_QT_WIDGET_OPENGL_H
