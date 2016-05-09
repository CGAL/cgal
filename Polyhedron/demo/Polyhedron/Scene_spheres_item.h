#ifndef SCENE_SPHERES_ITEM_H
#define SCENE_SPHERES_ITEM_H
#include "Scene_basic_objects_config.h"
#include "create_sphere.h"

#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <QOpenGLShaderProgram>

#include <QList>
#include <vector>
struct Scene_spheres_item_priv;
class SCENE_BASIC_OBJECTS_EXPORT Scene_spheres_item
    : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef std::pair<CGAL::Sphere_3<Kernel>*, CGAL::Color> Sphere_pair ;

  Scene_spheres_item(Scene_group_item* parent, bool planed = false);

  ~Scene_spheres_item();

  bool isFinite() const { return false; }
  bool isEmpty() const { return false; }
  Scene_item* clone() const {return 0;}
  QString toolTip() const {return QString();}
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Gouraud || m == Wireframe);
  }
  void compute_bbox() const { _bbox = Bbox(); }
  void add_sphere(CGAL::Sphere_3<Kernel>* sphere, CGAL::Color = CGAL::Color(120,120,120));
  void remove_sphere(CGAL::Sphere_3<Kernel>* sphere);
  void clear_spheres();
  void setPrecision(int prec);

  void draw(CGAL::Three::Viewer_interface* viewer) const;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
  void invalidateOpenGLBuffers();
  void computeElements() const;
  void setPlane(Kernel::Plane_3 p_plane);

protected:
  friend struct Scene_spheres_item_priv;
  Scene_spheres_item_priv* d;
};

#endif // SCENE_SPHERES_ITEM_H
