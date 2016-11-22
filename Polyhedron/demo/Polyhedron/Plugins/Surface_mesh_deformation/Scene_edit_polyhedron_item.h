#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H
//#define CGAL_PROFILE 
#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"


#include <CGAL/Three/Scene_group_item.h>

#include "Scene_polyhedron_item_k_ring_selection.h"
#include "Travel_isolated_components.h"

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <iostream>
#include <fstream>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/camera.h>

#include "ui_Deform_mesh.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <boost/function_output_iterator.hpp>
#include <QGLBuffer>
#include <QGLShader>
#include <QGLShaderProgram>


typedef Polyhedron::Vertex_handle                             Vertex_handle;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::face_descriptor      face_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;
class Scene_spheres_item;
namespace PMP = CGAL::Polygon_mesh_processing;
struct Array_based_vertex_point_map
{
public:
  typedef vertex_descriptor            key_type;
  typedef Polyhedron::Traits::Point_3  value_type;
  typedef const value_type&                  reference;
  typedef boost::read_write_property_map_tag category;
  Array_based_vertex_point_map(std::vector<double>* positions) : positions(positions) {}
  std::vector<double>* positions;
};


inline
Array_based_vertex_point_map::reference
get(Array_based_vertex_point_map,
  Array_based_vertex_point_map::key_type key) {
    return key->point();
}

inline
void
put(Array_based_vertex_point_map pmap,
  Array_based_vertex_point_map::key_type key,
  Array_based_vertex_point_map::value_type val)
{
  key->point() = val; // to make things easy (ray selection after deformation, save to polyhedron after close etc),
  // I also change point() of vertex together with positions list
  // So that we do not need to pmap everywhere other than draw
  if (key->id() == std::size_t(-1))
  {
    key->id() = pmap.positions->size() / 3;
    pmap.positions->push_back(val.x());
    pmap.positions->push_back(val.y());
    pmap.positions->push_back(val.z());
  }
  else
  {
  std::size_t pos = key->id() * 3;
  (*pmap.positions)[pos] = val.x();
  (*pmap.positions)[pos+1] = val.y();
  (*pmap.positions)[pos+2] = val.z();
  }
}

typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
  ,CGAL::Default, CGAL::Default, CGAL::Default, 
  Array_based_vertex_point_map> Deform_mesh;


typedef Deform_mesh::Point  Point;

/// For storing associated data with a group of control vertices
class Control_vertices_data
{
public:
  std::vector<vertex_descriptor> ctrl_vertices_group;
  qglviewer::ManipulatedFrame* frame;  // manframe assoc with a group of control vertices
  qglviewer::Vec frame_initial_center; // initial center of frame
  CGAL::Three::Scene_interface::Bbox bbox;          // bbox of control vertices inside group
  qglviewer::Vec rot_direction;        // vector for constraint rotation
private:
  std::vector<qglviewer::Vec> initial_positions;
  Deform_mesh* deform_mesh;

public:
  Control_vertices_data(Deform_mesh* deform_mesh, qglviewer::ManipulatedFrame* frame = 0)
    : frame(frame), bbox(0,0,0,0,0,0), rot_direction(0.,0.,1.), deform_mesh(deform_mesh)
  { }
  void refresh()
  {
    for(std::vector<vertex_descriptor>::iterator it = ctrl_vertices_group.begin(); it != ctrl_vertices_group.end(); ) {
      if(!deform_mesh->is_control_vertex(*it)) {
        it = ctrl_vertices_group.erase(it);
      }
      else { ++it; }
    }

    reset_initial_positions();
    frame_initial_center = calculate_initial_center();
    bbox = calculate_initial_bbox();

    bool oldState = frame->blockSignals(true); // do not let it Q_EMIT modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    // frame->setOrientation(qglviewer::Quaternion());
    frame->setPosition(frame_initial_center);
    frame->blockSignals(oldState);
  }
  void set_target_positions()
  {
    std::vector<vertex_descriptor>::iterator hb = ctrl_vertices_group.begin();
    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it, ++hb)
    {
      qglviewer::Vec dif_from_initial_center = (*it) - frame_initial_center;
      qglviewer::Vec rotated = frame->orientation() * dif_from_initial_center;
      qglviewer::Vec rotated_and_translated = rotated + frame->position();

      deform_mesh->set_target_position(*hb, Point(rotated_and_translated.x, rotated_and_translated.y, rotated_and_translated.z) );
    }
  }
  qglviewer::Vec calculate_initial_center() const
  {
    qglviewer::Vec center_acc(0, 0, 0);
    if (initial_positions.empty()) { return center_acc; }

    for (std::vector<qglviewer::Vec>::const_iterator it = initial_positions.begin();
         it != initial_positions.end(); ++it)
    {
      center_acc += (*it);
    }
    return center_acc / initial_positions.size();
  }

private:
  void reset_initial_positions()
  {
    initial_positions.clear();
    
    for(std::vector<vertex_descriptor>::iterator hb = ctrl_vertices_group.begin(); hb != ctrl_vertices_group.end(); ++hb)
    {
      qglviewer::Vec point((*hb)->point().x(), (*hb)->point().y(), (*hb)->point().z() );
      initial_positions.push_back(point);
    }
  }
  CGAL::Three::Scene_interface::Bbox calculate_initial_bbox()
  {    
    if(initial_positions.empty()) {return CGAL::Three::Scene_interface::Bbox(0,0,0,0,0,0); }

    const qglviewer::Vec& p_i = *(initial_positions.begin());
    CGAL::Three::Scene_interface::Bbox bbox(p_i.x, p_i.y, p_i.z, p_i.x, p_i.y, p_i.z);

    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it)
    {
      const qglviewer::Vec& p_i = (*it);
      CGAL::Three::Scene_interface::Bbox bbox_it(p_i.x, p_i.y, p_i.z, p_i.x, p_i.y, p_i.z);
      bbox = bbox + bbox_it;
    }
    return bbox;
  }
};

// To hold pressing states together
struct Mouse_keyboard_state_deformation
{
  bool ctrl_pressing;
  bool shift_pressing;
  bool left_button_pressing;
  bool right_button_pressing;

  Mouse_keyboard_state_deformation() 
    : ctrl_pressing(false), shift_pressing(false), left_button_pressing(false), right_button_pressing(false)
  { }
};

typedef std::list<Control_vertices_data> Ctrl_vertices_group_data_list;
struct Scene_edit_polyhedron_item_priv;
// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_EXPORT Scene_edit_polyhedron_item 
  : public CGAL::Three::Scene_group_item {
  Q_OBJECT
public:  
  /// Create an Scene_edit_polyhedron_item from a Scene_polyhedron_item.
  /// The ownership of the polyhedron is moved to the new edit_polyhedron
  /// item.
  Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item, Ui::DeformMesh* ui_widget, QMainWindow* mw);
  ~Scene_edit_polyhedron_item();

  /// Returns 0, so that one cannot clone an "edit polyhedron" item.
  Scene_edit_polyhedron_item* clone() const;

  // Function for displaying meta-data of the item
  QString toolTip() const;

  void setColor(QColor c);
  void setName(QString n);
  void setVisible(bool b);
  void setRenderingMode(RenderingMode m);
  
  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const { 
    return m == Gouraud; 
  }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  void draw() const{}
  void draw(CGAL::Three::Viewer_interface*) const;
  void drawEdges(CGAL::Three::Viewer_interface*) const;
  void draw_bbox(const CGAL::Three::Scene_interface::Bbox&) const;
  void draw_ROI_and_control_vertices(CGAL::Three::Viewer_interface *viewer) const;
  void draw_frame_plane(QGLViewer *) const;

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  /// Returns a Scene_polyhedron_item from the edit polyhedron item, and
  /// transfer the ownership of the polyhedron to it.
  /// The item 'this' must be destroy just after a call to this function.
  Scene_polyhedron_item* to_polyhedron_item();

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  void compute_bbox() const;
  Bbox bbox() const{return Scene_item::bbox();}

  int get_k_ring();
  void set_k_ring(int v);

  // take mouse events from viewer, main-window does not work
  // take keyboard events from main-window, which is more stable
  bool eventFilter(QObject *target, QEvent *event);
  void update_frame_plane();
  void ShowAsSphere(bool b);

public Q_SLOTS:
  void reset_spheres();
  void updateDeform();
  void change();

  void invalidateOpenGLBuffers();
  void selected(const std::set<Polyhedron::Vertex_handle>& m);

  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

  void deform(); // deform the mesh
  void remesh();

protected:
  friend struct Scene_edit_polyhedron_item_priv;
  Scene_edit_polyhedron_item_priv* d;


public:
  // Deformation related functions //
  bool insert_control_vertex(vertex_descriptor v);

  bool insert_roi_vertex(vertex_descriptor v);
  
  bool erase_control_vertex(vertex_descriptor v);

  bool erase_roi_vertex(vertex_descriptor v);

  void set_all_vertices_as_roi()
  {
    vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = vertices(*polyhedron()); vb != ve; ++vb)
    {
      insert_roi_vertex(*vb);
    }   
  }

  void clear_roi();

  void create_ctrl_vertices_group();

  void delete_ctrl_vertices_group(bool create_new = true);

  void prev_ctrl_vertices_group();
  void next_ctrl_vertices_group();
  void pivoting_end();

  void pivoting_begin();

  void save_roi(const char* file_name) const;
  void read_roi(const char* file_name);
  void overwrite_deform_object();

  void reset_deform_object();
  struct Is_selected {
    Deform_mesh* dm;
    Is_selected(Deform_mesh* dm) : dm(dm) {}
    bool count(Vertex_handle vh) const {
      return dm->is_roi_vertex(vh);
    }
  };

  boost::optional<std::size_t> get_minimum_isolated_component();
  struct Select_roi_output {
    Select_roi_output(Deform_mesh* dm) : dm(dm) { }
    void operator()(Vertex_handle vh) {
      dm->insert_roi_vertex(vh);
    }
    Deform_mesh* dm;
  };

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) ;

protected:
  // Deformation related functions //
  void print_message(const QString& /*message*/)
  {
    // std::cout << message.toStdString() << std::endl;
  }

  bool is_there_any_ctrl_vertices_group(Ctrl_vertices_group_data_list::iterator& hgb, Ctrl_vertices_group_data_list::iterator& hge);

  bool is_there_any_ctrl_vertices_group();

  bool is_there_any_ctrl_vertices();
  void refresh_all_group_centers();
  bool activate_closest_manipulated_frame(int x, int y);

  bool keyPressEvent(QKeyEvent* e);

  void update_normals();

  double scene_diag() const {
    const double& xdelta = bbox().xmax() - bbox().xmin();
    const double& ydelta = bbox().ymax() - bbox().ymin();
    const double& zdelta = bbox().zmax() - bbox().zmin();
    const double diag = std::sqrt(xdelta*xdelta +
      ydelta*ydelta +
      zdelta*zdelta);
    return diag * 0.5;
  }

}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
