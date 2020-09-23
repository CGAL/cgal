#ifndef POINT_SET_ITEM_H
#define POINT_SET_ITEM_H
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Scene_item_with_properties.h>
#include <CGAL/Three/Scene_zoomable_item_interface.h>
#include "Scene_points_with_normal_item_config.h"
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include "Kernel_type.h"
#include "Point_set_3.h"
#include <iostream>
#include <CGAL/config.h>

struct Scene_points_with_normal_item_priv;
// point set
typedef Point_set_3<Kernel> Point_set;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;

class QMenu;
class QAction;

// This class represents a point set in the OpenGL scene
class SCENE_POINTS_WITH_NORMAL_ITEM_EXPORT Scene_points_with_normal_item
  : public CGAL::Three::Scene_item_rendering_helper,
    public CGAL::Three::Scene_item_with_properties,
    public CGAL::Three::Scene_zoomable_item_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Scene_zoomable_item_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.ZoomInterface/1.0")

public:
  Scene_points_with_normal_item();
  Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy);

  Scene_points_with_normal_item(const SMesh& input_mesh);

  ~Scene_points_with_normal_item();
  Scene_points_with_normal_item* clone() const Q_DECL_OVERRIDE;

  // Is selection empty?
  virtual bool isSelectionEmpty() const;

  // Function to override the context menu
  QMenu* contextMenu() Q_DECL_OVERRIDE;

  // IO
#ifdef CGAL_LINKED_WITH_LASLIB
  bool read_las_point_set(std::istream& in);
  bool write_las_point_set(std::ostream& out) const;
#endif // LAS
  bool read_ply_point_set(std::istream& in);
  bool write_ply_point_set(std::ostream& out, bool binary) const;
  bool read_off_point_set(std::istream& in);
  bool write_off_point_set(std::ostream& out) const;
  bool read_xyz_point_set(std::istream& in);
  bool write_xyz_point_set(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const Q_DECL_OVERRIDE;

  virtual void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE;

  virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  virtual void drawPoints(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;

  // Gets wrapped point set
  Point_set*       point_set();
  const Point_set* point_set() const;

  // Gets PLY comments (empty if point set not originated from PLY input)
  std::string& comments();
  const std::string& comments() const;

  // Gets dimensions
  virtual bool isFinite() const Q_DECL_OVERRIDE { return true; }
  virtual bool isEmpty() const Q_DECL_OVERRIDE;
  virtual void compute_bbox() const Q_DECL_OVERRIDE;

  // computes the local point spacing (aka radius) of each point
  void computes_local_spacing(int k);

  bool has_normals() const;
  void copyProperties(Scene_item *) Q_DECL_OVERRIDE;
  int getNormalSliderValue();
  int getPointSliderValue();
  void computeElements() const Q_DECL_OVERRIDE;
  void initializeBuffers(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;

public Q_SLOTS:
  // Delete selection
  virtual void deleteSelection();
  // Invert selection
  void invertSelection();
  // Select all points
  void selectAll();
  // Reset selection mark
  void resetSelection();
  //Select duplicated points
  void selectDuplicates();
  //Set the status of the slider as `pressed`
  void pointSliderPressed();
  //Set the status of the slider as `released`
  void pointSliderReleased();
  void itemAboutToBeDestroyed(Scene_item *) Q_DECL_OVERRIDE;
  void setPointSize(int size);
  void setNormalSize(int size);
  void resetColors();
// Data
protected:
  friend struct Scene_points_with_normal_item_priv;
  Scene_points_with_normal_item_priv* d;

public:
 void zoomToPosition(const QPoint &, CGAL::Three::Viewer_interface *)const Q_DECL_OVERRIDE;
}; // end class Scene_points_with_normal_item


#endif // POINT_SET_ITEM_H
