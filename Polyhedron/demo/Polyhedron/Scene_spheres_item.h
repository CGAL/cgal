#ifndef SCENE_SPHERES_ITEM_H
#define SCENE_SPHERES_ITEM_H
#include "Scene_basic_objects_config.h"
#include "create_sphere.h"

#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <QOpenGLShaderProgram>

#include <QList>
#include <vector>
struct Scene_spheres_item_priv;

/* This item contains spheres and associated colors. They are kept in a Spheres_container,
 * sorted by the value of their "index". This item also has an internal picking mechanism that
 * colorizes all the spheres that has the same index as the one that has been picked.
 * The picking is only usable if several indices exist.
 * If all the spheres have the index 0, they can have independant colors (generally used by the items that
 * have a Scene_spheres_item child).
*/
class SCENE_BASIC_OBJECTS_EXPORT Scene_spheres_item
    : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Sphere_3<Kernel> Sphere;
  typedef std::pair<Sphere, CGAL::Color> Sphere_pair;
  typedef std::vector<std::vector<Sphere_pair> > Spheres_container;

  Scene_spheres_item(Scene_group_item* parent, std::size_t max_index = 0, bool planed = false, bool pickable = true);

  ~Scene_spheres_item();

  bool isFinite() const Q_DECL_OVERRIDE{ return true; }
  bool isEmpty() const Q_DECL_OVERRIDE{ return false; }
  Scene_item* clone() const Q_DECL_OVERRIDE{return 0;}
  QString toolTip() const Q_DECL_OVERRIDE;
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE{
    return (m == Gouraud || m == Wireframe);
  }
  void compute_bbox() const Q_DECL_OVERRIDE;
  void add_sphere(const Sphere &sphere, std::size_t index = 0, CGAL::Color = CGAL::Color(120,120,120));
  void clear_spheres();
  void setPrecision(int prec);
  void gl_initialization(CGAL::Three::Viewer_interface* viewer);
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;
  void setPlane(Kernel::Plane_3 p_plane);
  void setToolTip(QString s);
  void setColor(QColor c) Q_DECL_OVERRIDE;
  bool save(const std::string &file_name) const;


  void initializeBuffers(Viewer_interface *) const Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
  bool eventFilter(QObject *watched, QEvent *event)Q_DECL_OVERRIDE;
Q_SIGNALS:
  void on_color_changed();
  void picked(std::size_t id) const;
  void destroyMe();
protected:
  friend struct Scene_spheres_item_priv;
  Scene_spheres_item_priv* d;
};

#endif // SCENE_SPHERES_ITEM_H
