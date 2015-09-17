#ifndef POINT_SET_ITEM_H
#define POINT_SET_ITEM_H
#include "Scene_item.h"
#include "Scene_points_with_normal_item_config.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>


// point set
typedef Point_set_3<Kernel> Point_set;
typedef Point_set::UI_point UI_point; // type of points in Point_set_3

class QMenu;
class QAction;

// This class represents a point set in the OpenGL scene
class SCENE_POINTS_WITH_NORMAL_ITEM_EXPORT Scene_points_with_normal_item
  : public Scene_item
{
  Q_OBJECT

public:
  Scene_points_with_normal_item();
  Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy);
  Scene_points_with_normal_item(const Polyhedron& p);
  ~Scene_points_with_normal_item();
  Scene_points_with_normal_item* clone() const;

  // Is selection empty?
  virtual bool isSelectionEmpty() const;

  // Function to override the context menu
  QMenu* contextMenu();

  // IO
  bool read_off_point_set(std::istream& in);
  bool write_off_point_set(std::ostream& out) const;
  bool read_xyz_point_set(std::istream& in);
  bool write_xyz_point_set(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  virtual void invalidate_buffers();

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;

  virtual void draw_edges(Viewer_interface* viewer) const;
  virtual void draw_points(Viewer_interface*) const;

  virtual void draw_splats(Viewer_interface*) const;
  
  // Gets wrapped point set
  Point_set*       point_set();
  const Point_set* point_set() const;

  // Gets dimensions
  virtual bool isFinite() const { return true; }
  virtual bool isEmpty() const;
  virtual Bbox bbox() const;

  virtual void setRenderingMode(RenderingMode m);

  // computes the local point spacing (aka radius) of each point
  void computes_local_spacing(int k);

  bool has_normals() const;
  void set_has_normals(bool b);

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

// Data
private:
  Point_set* m_points;
  bool m_has_normals;
  QAction* actionDeleteSelection;
  QAction* actionResetSelection;
  QAction* actionSelectDuplicatedPoints;


  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> positions_selected_points;
  mutable std::vector<double> normals;
  mutable std::size_t nb_points;
  mutable std::size_t nb_selected_points;
  mutable std::size_t nb_lines;

  mutable QOpenGLShaderProgram *program;

  using Scene_item::initialize_buffers;
  void initialize_buffers(Viewer_interface *viewer) const;

  void compute_normals_and_vertices() const;


}; // end class Scene_points_with_normal_item


#endif // POINT_SET_ITEM_H
