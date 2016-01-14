#ifndef POINT_SET_ITEM_H
#define POINT_SET_ITEM_H
#include  <CGAL/Three/Scene_item.h>

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
  : public CGAL::Three::Scene_item
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
  bool read_ply_point_set(std::istream& in);
  bool write_ply_point_set(std::ostream& out) const;
  bool read_off_point_set(std::istream& in);
  bool write_off_point_set(std::ostream& out) const;
  bool read_xyz_point_set(std::istream& in);
  bool write_xyz_point_set(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  virtual void invalidate_buffers();

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;

  virtual void draw_edges(CGAL::Three::Viewer_interface* viewer) const;
  virtual void draw_points(CGAL::Three::Viewer_interface*) const;

  virtual void draw_splats(CGAL::Three::Viewer_interface*) const;
  
  // Gets wrapped point set
  Point_set*       point_set();
  const Point_set* point_set() const;

  // Gets dimensions
  virtual bool isFinite() const { return true; }
  virtual bool isEmpty() const;
  virtual void compute_bbox() const;

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
  std::vector<QColor> m_colors;
  bool m_has_normals;
  QAction* actionDeleteSelection;
  QAction* actionResetSelection;
  QAction* actionSelectDuplicatedPoints;

  enum VAOs {
      Edges=0,
      ThePoints,
      Selected_points,
      NbOfVaos = Selected_points+1
  };
  enum VBOs {
      Edges_vertices = 0,
      Points_vertices,
      Selected_points_vertices,
      NbOfVbos = Selected_points_vertices+1
  };

  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> positions_selected_points;
  mutable std::vector<double> normals;
  mutable std::size_t nb_points;
  mutable std::size_t nb_selected_points;
  mutable std::size_t nb_lines;

  mutable QOpenGLShaderProgram *program;

  using CGAL::Three::Scene_item::initialize_buffers;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;

  void compute_normals_and_vertices() const;


}; // end class Scene_points_with_normal_item


template <typename PlyReadNumber>
class Custom_ply_interpreter
{
  Point_set* m_points;
  std::vector<QColor>& m_colors;
    
public:
  Custom_ply_interpreter (Point_set* points,
                          std::vector<QColor>& colors)
    : m_points (points), m_colors (colors)
  { }

  bool is_applicable (const std::vector<PlyReadNumber*>& readers)
  {
    bool x_found = false, y_found = false, z_found = false;
    for (std::size_t i = 0; i < readers.size (); ++ i)
      if (readers[i]->name () == "x")
        x_found = true;
      else if (readers[i]->name () == "y")
        y_found = true;
      else if (readers[i]->name () == "z")
        z_found = true;

    return x_found && y_found && z_found;
  }
      
  void operator() (const std::vector<PlyReadNumber*>& readers)
  {
    double x, y, z, nx, ny, nz;
    int r = -1, g = -1, b = -1;
    for (std::size_t i = 0; i < readers.size (); ++ i)
      if (readers[i]->name () == "x")
        readers[i]->assign (x);
      else if (readers[i]->name () == "y")
        readers[i]->assign (y);
      else if (readers[i]->name () == "z")
        readers[i]->assign (z);
      else if (readers[i]->name () == "nx")
        readers[i]->assign (nx);
      else if (readers[i]->name () == "ny")
        readers[i]->assign (ny);
      else if (readers[i]->name () == "nz")
        readers[i]->assign (nz);
      else if (readers[i]->name () == "red"
               || readers[i]->name () == "r")
        readers[i]->assign (r);
      else if (readers[i]->name () == "green"
               || readers[i]->name () == "g")
        readers[i]->assign (g);
      else if (readers[i]->name () == "blue"
               || readers[i]->name () == "b")
        readers[i]->assign (b);

    m_points->push_back (UI_point (typename Kernel::Point_3 (x, y, z),
                                   typename Kernel::Vector_3 (nx, ny, nz)));
    if (r != -1 && g != -1 && b != -1)
      m_colors.push_back (QColor (r, g, b));
  }

};



#endif // POINT_SET_ITEM_H
