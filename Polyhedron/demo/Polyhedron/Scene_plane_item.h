#ifndef SCENE_PLANE_ITEM_H
#define SCENE_PLANE_ITEM_H


#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include "Scene_basic_objects_config.h"

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <CGAL/Three/Viewer_interface.h>

#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_epic;
typedef Kernel_epic::Plane_3 Plane_3;

class SCENE_BASIC_OBJECTS_EXPORT Scene_plane_item 
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_plane_item(const CGAL::Three::Scene_interface* scene_interface);
  ~Scene_plane_item();

  double scene_diag() const {
    const Scene_item::Bbox& bbox = scene->bbox();
    const double& xdelta = bbox.xmax()-bbox.xmin();
    const double& ydelta = bbox.ymax()-bbox.ymin();
    const double& zdelta = bbox.zmax()-bbox.zmin();
    const double diag = std::sqrt(xdelta*xdelta +
                            ydelta*ydelta +
                            zdelta*zdelta);
    return diag * 0.7;
  }
  bool isFinite() const { return false; }
  bool isEmpty() const { return false; }
  void compute_bbox() const { _bbox = Bbox(); }
  bool manipulatable() const;
  ManipulatedFrame* manipulatedFrame();

  Scene_plane_item* clone() const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe || m == Flat || m == FlatPlusEdges);
  }
  virtual void draw(CGAL::Three::Viewer_interface*) const;
 virtual void drawEdges(CGAL::Three::Viewer_interface* viewer)const;
  Plane_3 plane() const;

public Q_SLOTS:
  virtual void invalidateOpenGLBuffers();

  void setPosition(float x, float y, float z);
  
  void setPosition(double x, double y, double z);
  
  void setNormal(float x, float y, float z);

  void setNormal(double x, double y, double z);
  void flipPlane();
  void setClonable(bool b = true);

  void setManipulatable(bool b = true);
protected:

  const CGAL::Three::Scene_interface* scene;
  bool manipulable;
  bool can_clone;
  qglviewer::ManipulatedFrame* frame;

  enum VAOs {
      Facets = 0,
      Edges,
      NbOfVaos
  };
  enum VBOs {
      Facets_vertices = 0,
      Edges_vertices,
      NbOfVbos
  };

  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_quad;
  mutable GLint sampler_location;
  mutable bool smooth_shading;
  mutable QOpenGLShaderProgram *program;

  void initializeBuffers(CGAL::Three::Viewer_interface*)const;
  void compute_normals_and_vertices(void) const;
};

#endif // SCENE_PLANE_ITEM_H
