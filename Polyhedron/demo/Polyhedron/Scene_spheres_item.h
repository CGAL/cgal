#ifndef SCENE_SPHERES_ITEM_H
#define SCENE_SPHERES_ITEM_H
#include "Scene_basic_objects_config.h"
#include "create_sphere.h"
#include "CGAL/Three/Scene_group_item.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <QOpenGLShaderProgram>

#include <QList>
#include<vector>

class SCENE_BASIC_OBJECTS_EXPORT Scene_spheres_item
    : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef std::pair<CGAL::Sphere_3<Kernel>*, CGAL::Color> sphere_pair ;

  Scene_spheres_item(Scene_group_item* parent, bool planed = false)
    :CGAL::Three::Scene_item(NbOfVbos,NbOfVaos)
    ,precision(36)
    ,has_plane(planed)

  {
    setParent(parent);
    create_flat_and_wire_sphere(1.0f,vertices,normals, edges);
  }

  ~Scene_spheres_item() {
    Q_FOREACH(sphere_pair sphere, spheres)
      delete sphere.first;
  }

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
  void setPrecision(int prec) { precision = prec; }

  void draw(CGAL::Three::Viewer_interface* viewer) const;
  void draw_edges(CGAL::Three::Viewer_interface* viewer) const;
  void invalidateOpenGLBuffers(){are_buffers_filled = false;}
  void computeElements() const;
  void setPlane(Kernel::Plane_3 p_plane) { plane = p_plane; }

private:
  enum Vbos
  {
    Vertices = 0,
    Edge_vertices,
    Normals,
    Center,
    Radius,
    Color,
    Edge_color,
    NbOfVbos
  };
  enum Vaos
  {
    Facets = 0,
    Edges,
    NbOfVaos
  };


  int precision;
  mutable CGAL::Plane_3<Kernel> plane;
  bool has_plane;

  QList<sphere_pair> spheres;
  mutable std::vector<float> vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> edges;
  mutable std::vector<float> colors;
  mutable std::vector<float> edges_colors;
  mutable std::vector<float> centers;
  mutable std::vector<float> radius;
  mutable QOpenGLShaderProgram *program;
  mutable int nb_centers;
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;

};

#endif // SCENE_SPHERES_ITEM_H
