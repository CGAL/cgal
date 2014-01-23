#ifndef SCENE_C3T3_ITEM_H
#define SCENE_C3T3_ITEM_H

#include <CGAL/config.h>
#ifndef CGAL_CPP11_OVERRIDE
#  define CGAL_CPP11_OVERRIDE
#endif

#include "Scene_item_with_display_list.h"
#include "Scene_c3t3_item_config.h"
#include "C3t3_type.h"
#include <iostream>

#include <QtCore/qglobal.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <CGAL/glu.h>

struct Scene_c3t3_item_priv;

class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item : public Scene_item_with_display_list
{
  friend struct Scene_c3t3_item_priv;

  Q_OBJECT
  Q_PROPERTY(int numberOfVertices READ number_of_vertices)
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;
  typedef C3t3::Triangulation Tr;

  Scene_c3t3_item();
  Scene_c3t3_item(const C3t3& c3t3, Scene_item* parent = 0);
  ~Scene_c3t3_item() CGAL_CPP11_OVERRIDE;

  bool load_binary(std::istream&);
  bool save_binary(std::ostream&) const;
  bool export_boundary_to_OFF(std::ostream&) const;

  inline const Scene_item* parent() const;
  inline Scene_item* parent();

  inline const C3t3& c3t3() const;
  inline C3t3& c3t3();

  inline bool manipulatable() const;
  inline ManipulatedFrame* manipulatedFrame();

  inline void setPosition(float x, float y, float z);
  inline void setNormal(float x, float y, float z);

  Kernel::Plane_3 plane() const;

  inline bool isFinite() const CGAL_CPP11_OVERRIDE;
  inline bool isEmpty() const CGAL_CPP11_OVERRIDE;

  Bbox bbox() const CGAL_CPP11_OVERRIDE;

  inline Scene_c3t3_item* clone() const;

  QString toolTip() const CGAL_CPP11_OVERRIDE;
  QPixmap graphicalToolTip() const CGAL_CPP11_OVERRIDE;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const;

  QMenu* contextMenu();

  void draw() const CGAL_CPP11_OVERRIDE;
  void draw_edges() const CGAL_CPP11_OVERRIDE;
  void draw_points() const CGAL_CPP11_OVERRIDE;
  // override Scene_item_with_display_list::draw() to deal with the cutting
  // plane
  void changed() CGAL_CPP11_OVERRIDE;

  void setColor(QColor c) CGAL_CPP11_OVERRIDE;

  // rebuild histogram
  inline void update_histogram();

  void direct_draw() const; // draw nothing
  void direct_draw_edges() const; // draw all edges

  void draw_triangles(bool with_cut_plane = false,
                      bool in_draw_edges = false) const;
  void draw_cut_tetrahedra() const;
  void draw_sphere(const Tr::Point& p) const;

  int number_of_vertices() const;

public slots:
  void frameIsModified();
  void show_spheres(bool b);
  void show_cut_plane(bool b);
  void show_tetrahedra(bool b);
  void show_triangles_duals(bool b);

  void new_spheres_drawn_radius();
  void new_spheres_drawn_radius(double);
  void restore_spheres_drawn_radii();

  void recenter_cut_plane();
  void reset_cut_plane();

  void selectedPoint(double,
                     double,
                     double);
private:
  void build_histogram();
  QColor get_histogram_color(const double v) const;

private:
  static void draw_triangle(const Kernel::Point_3& pa,
                            const Kernel::Point_3& pb,
                            const Kernel::Point_3& pc);

  double complex_diag() const;

private:
  QPixmap histogram_;
  Scene_c3t3_item_priv* d; // d-pointer

  Scene_item* parent_item;
  C3t3 c3t3_;
  qglviewer::ManipulatedFrame* frame;
  mutable GLuint sphere_display_list;
  mutable GLUquadric* quadric;
  bool draw_spheres;
  double spheres_drawn_radius;
  bool draw_tets;
};


inline
const Scene_item*
Scene_c3t3_item::parent() const
{
  return parent_item;
}

inline
Scene_item*
Scene_c3t3_item::parent()
{
  return parent_item;
}

inline
const C3t3&
Scene_c3t3_item::c3t3() const
{
  return c3t3_;
}

inline
C3t3&
Scene_c3t3_item::c3t3()
{
  return c3t3_;
}

inline
bool
Scene_c3t3_item::manipulatable() const
{
  return true;
}

inline
Scene_c3t3_item::ManipulatedFrame*
Scene_c3t3_item::manipulatedFrame()
{
  return frame;
}

inline
void
Scene_c3t3_item::setPosition(float x, float y, float z)
{
  frame->setPosition(x, y, z);
}

inline
void
Scene_c3t3_item::setNormal(float x, float y, float z)
{
  frame->setOrientation(x, y, z, 0.f);
}

inline
bool
Scene_c3t3_item::isFinite() const
{
  return true;
}

inline
bool
Scene_c3t3_item::isEmpty() const
{
  return c3t3().triangulation().number_of_vertices() == 0;
}

inline
Scene_c3t3_item*
Scene_c3t3_item::clone() const
{
  return 0;
}

inline
bool
Scene_c3t3_item::supportsRenderingMode(RenderingMode m) const
{
  return (m != Gouraud); // CHECK THIS!
}



#endif // SCENE_C3T3_ITEM_H
