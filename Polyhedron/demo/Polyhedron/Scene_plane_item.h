#ifndef SCENE_PLANE_ITEM_H
#define SCENE_PLANE_ITEM_H


#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Scene_interface.h>

#include "Scene_basic_objects_config.h"

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>
#include <CGAL/Three/Viewer_interface.h>

#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_epic;
typedef Kernel_epic::Plane_3 Plane_3;

class SCENE_BASIC_OBJECTS_EXPORT Scene_plane_item
  : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public:
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_plane_item(const CGAL::Three::Scene_interface* scene_interface);
  ~Scene_plane_item();

  double scene_diag() const {
    /*If no item is visible, scene->bbox is 0,0,0,0,0,0 and the texture is empty.
    To avoid that, we need to compute the scene's bbox if the items were visible.
    {
*/
    Scene_item::Bbox bbox = scene->bbox();
    if(bbox == Scene_item::Bbox(std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()))
      bbox = Scene_item::Bbox(0,0,0,0,0,0);
    if(bbox == Scene_item::Bbox(0,0,0,0,0,0))
    {
      for(int id = 0; id< scene->numberOfEntries(); ++id)
      {
        if(scene->item(id)->isFinite() && !scene->item(id)->isEmpty())
          bbox = bbox + scene->item(id)->bbox();
      }
    }
    else
    {
      bbox = scene->bbox();
    }
    //}
    const double& xdelta = bbox.xmax()-bbox.xmin();
    const double& ydelta = bbox.ymax()-bbox.ymin();
    const double& zdelta = bbox.zmax()-bbox.zmin();
    const double diag = std::sqrt(xdelta*xdelta +
                            ydelta*ydelta +
                            zdelta*zdelta);
    return diag * 0.7;
  }
  bool isFinite() const Q_DECL_OVERRIDE { return false; }
  bool isEmpty() const Q_DECL_OVERRIDE { return false; }
  void compute_bbox() const Q_DECL_OVERRIDE { _bbox = Bbox(); }
  bool manipulatable() const Q_DECL_OVERRIDE;
  ManipulatedFrame* manipulatedFrame() Q_DECL_OVERRIDE;
  QMenu* contextMenu() Q_DECL_OVERRIDE;

  Scene_plane_item* clone() const Q_DECL_OVERRIDE ;

  QString toolTip() const Q_DECL_OVERRIDE ;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE{
    return (m == Wireframe || m == Flat || m == FlatPlusEdges);
  }
  virtual void draw(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
  virtual void drawEdges(CGAL::Three::Viewer_interface* viewer)const Q_DECL_OVERRIDE;
  Plane_3 plane(CGAL::qglviewer::Vec offset = CGAL::qglviewer::Vec(0,0,0)) const;

public Q_SLOTS:
  virtual void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;

  void setPosition(float x, float y, float z);

  void setPosition(double x, double y, double z);

  void setNormal(float x, float y, float z);

  void setNormal(double x, double y, double z);
  void flipPlane();
  void setClonable(bool b = true);

  void setManipulatable(bool b = true);
  void setPlaneOrientation();
protected:

  const CGAL::Three::Scene_interface* scene;
  bool manipulable;
  bool can_clone;
  CGAL::qglviewer::ManipulatedFrame* frame;

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
  mutable std::size_t nb_quads;
  mutable std::size_t nb_lines;
  mutable GLint sampler_location;
  mutable bool smooth_shading;
  mutable QOpenGLShaderProgram *program;

  void initializeBuffers(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
};

#endif // SCENE_PLANE_ITEM_H
