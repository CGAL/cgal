#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H
//#define CGAL_PROFILE 
#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"


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


typedef boost::graph_traits<SMesh>::vertex_descriptor sm_vertex_descriptor;
typedef boost::graph_traits<SMesh>::edge_descriptor sm_edge_descriptor;
typedef boost::graph_traits<SMesh>::face_descriptor sm_face_descriptor;
typedef boost::graph_traits<SMesh>::halfedge_descriptor sm_halfedge_descriptor;
//sets/gets the ID of a Mesh vertex descriptor in a different manner if Mesh is a Polyhedron or a SMesh
struct Id_setter{
  typedef boost::graph_traits<SMesh>::vertex_descriptor sm_vd;
  typedef boost::graph_traits<Polyhedron>::vertex_descriptor p_vd;
  typedef boost::graph_traits<SMesh>::face_descriptor sm_fd;
  typedef boost::graph_traits<Polyhedron>::face_descriptor p_fd;

  Polyhedron* poly;

  SMesh* sm;
  boost::property_map< SMesh, boost::vertex_index_t >::type im;

  Id_setter(Polyhedron* poly)
    :poly(poly), sm(NULL){}

  Id_setter(SMesh* sm)
    :poly(NULL), sm(sm)
  {
    im = get(boost::vertex_index, *sm);
  }

  std::size_t get_id(p_vd vd)
  {
    return vd->id();
  }
  std::size_t get_id(sm_vd vd)
  {
    return static_cast<std::size_t>(im[vd]);
  }

  void set_id(p_vd vd, std::size_t id)
  {
    vd->id() = id;
  }
  //update the visu map
  void set_id(sm_vd, std::size_t)
  {
  }

  void set_id(p_fd fd, std::size_t id)
  {
    boost::property_map<Polyhedron, CGAL::face_index_t>::type fim
        = get(CGAL::face_index, *poly);
    put(fim, fd, id);
  }
  //cannot set a surface_mesh id, but it is only use din case the id = -1 which cannot happen so it's ok
  void set_id(sm_fd, std::size_t)
  {
      return;
  }
};




typedef boost::graph_traits<Polyhedron>::vertex_descriptor    poly_vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      poly_vertex_iterator;

class Scene_spheres_item;
namespace PMP = CGAL::Polygon_mesh_processing;
template<typename Mesh>
struct Array_based_vertex_point_map
{
public:
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor     key_type;
  typedef Kernel::Point_3                                           value_type;
  typedef const value_type&                                         reference;
  typedef boost::read_write_property_map_tag                        category;
  std::vector<double>* positions;
  Mesh* mesh;
  Id_setter* id_setter;
  Array_based_vertex_point_map(std::vector<double>* positions, Mesh* mesh, Id_setter* id_setter) : positions(positions), mesh(mesh), id_setter(id_setter) {}

};


template<typename Mesh> inline
typename Array_based_vertex_point_map<Mesh>::reference
get(Array_based_vertex_point_map<Mesh> map,
  typename Array_based_vertex_point_map<Mesh>::key_type key) {
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VertexPointMap;
  VertexPointMap pmap = get(boost::vertex_point, *map.mesh);
  return get(pmap, key);

}

template<typename Mesh> inline
void
put(Array_based_vertex_point_map<Mesh> map,
  typename Array_based_vertex_point_map<Mesh>::key_type key,
  typename Array_based_vertex_point_map<Mesh>::value_type val)
{
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(boost::vertex_point, *map.mesh);
  put(vpmap, key, val);// to make things easy (ray selection after deformation, save to polyhedron after close etc),
  // I also change point() of vertex together with positions list
  // So that we do not need to pmap everywhere other than draw
  if (map.id_setter->get_id(key) == std::size_t(-1))
  {
    map.id_setter->set_id(key, static_cast<std::size_t>(map.positions->size() / 3));
    map.positions->push_back(val.x());
    map.positions->push_back(val.y());
    map.positions->push_back(val.z());
  }
  else
  {
    std::size_t pos = map.id_setter->get_id(key) * 3;
    if(pos < map.positions->size()-1)
    {
      (*map.positions)[pos] = val.x();
      (*map.positions)[pos+1] = val.y();
      (*map.positions)[pos+2] = val.z();
    }
  }
}

typedef CGAL::Surface_mesh_deformation<Polyhedron, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
  ,CGAL::Default, CGAL::Default, CGAL::Default, 
  Array_based_vertex_point_map<Polyhedron> > Deform_mesh;

typedef CGAL::Surface_mesh_deformation<SMesh, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
  ,CGAL::Default, CGAL::Default, CGAL::Default,
  Array_based_vertex_point_map<SMesh> > Deform_sm_mesh;

typedef Deform_mesh::Point  Point;
typedef Deform_sm_mesh::Point  SM_Point;

/// For storing associated data with a group of control vertices
template<typename Mesh>
class Control_vertices_data
{
public:
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor mesh_vd;
  typedef typename CGAL::Surface_mesh_deformation<Mesh, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
    ,CGAL::Default, CGAL::Default, CGAL::Default,
    Array_based_vertex_point_map<Mesh> > M_Deform_mesh;

  std::vector<mesh_vd> ctrl_vertices_group;
  qglviewer::ManipulatedFrame* frame;  // manframe assoc with a group of control vertices
  qglviewer::Vec frame_initial_center; // initial center of frame
  CGAL::Three::Scene_interface::Bbox bbox;          // bbox of control vertices inside group
  qglviewer::Vec rot_direction;        // vector for constraint rotation
private:
  std::vector<qglviewer::Vec> initial_positions;
  M_Deform_mesh* deform_mesh;

public:
  Control_vertices_data(M_Deform_mesh* deform_mesh, qglviewer::ManipulatedFrame* frame = 0)
    : frame(frame), bbox(0,0,0,0,0,0), rot_direction(0.,0.,1.), deform_mesh(deform_mesh)
  { }

  void refresh(Mesh *mesh)
  {
    for(typename std::vector<mesh_vd>::iterator it = ctrl_vertices_group.begin(); it != ctrl_vertices_group.end(); ) {
      if(!deform_mesh->is_control_vertex(*it)) {
        it = ctrl_vertices_group.erase(it);
      }
      else { ++it; }
    }

    reset_initial_positions(mesh);
    frame_initial_center = calculate_initial_center();
    bbox = calculate_initial_bbox();

    bool oldState = frame->blockSignals(true); // do not let it Q_EMIT modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    // frame->setOrientation(qglviewer::Quaternion());
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    frame->setPosition(frame_initial_center+offset);
    frame->blockSignals(oldState);
  }
  void set_target_positions()
  {
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    typename std::vector<mesh_vd>::iterator hb = ctrl_vertices_group.begin();
    for( std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it, ++hb)
    {
      qglviewer::Vec dif_from_initial_center = (*it) - frame_initial_center;
      qglviewer::Vec rotated = frame->orientation() * dif_from_initial_center;
      qglviewer::Vec rotated_and_translated = rotated + frame->position();

      deform_mesh->set_target_position(*hb, typename M_Deform_mesh::Point(rotated_and_translated.x-offset.x, rotated_and_translated.y-offset.y, rotated_and_translated.z-offset.z) );
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
  void reset_initial_positions(Mesh* mesh)
  {
    initial_positions.clear();
    typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VertexPointMap;
    VertexPointMap pmap = get(boost::vertex_point, *mesh);
    for(typename std::vector<mesh_vd>::iterator hb = ctrl_vertices_group.begin(); hb != ctrl_vertices_group.end(); ++hb)
    {

      typename M_Deform_mesh::Point p = get(pmap, (*hb));
      qglviewer::Vec point(p.x(), p.y(), p.z() );
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

typedef std::list<Control_vertices_data<Polyhedron> > Ctrl_vertices_poly_group_data_list;
typedef std::list<Control_vertices_data<SMesh> > Ctrl_vertices_sm_group_data_list;
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
  Scene_edit_polyhedron_item(Scene_surface_mesh_item* sm_item, Ui::DeformMesh* ui_widget, QMainWindow* mw);
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
  template<typename Mesh>
  void draw_frame_plane(QGLViewer *, Mesh *mesh) const;

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;
  Scene_polyhedron_item* poly_item() const;
  // Get wrapped Surface_mesh
  SMesh*       surface_mesh();
  const SMesh* surface_mesh() const;
  Scene_surface_mesh_item* sm_item() const;

  /// Returns a Scene_polyhedron_item from the edit polyhedron item, and
  /// transfer the ownership of the polyhedron to it.
  /// The item 'this' must be destroy just after a call to this function.
  Scene_polyhedron_item* to_polyhedron_item();
  Scene_surface_mesh_item* to_sm_item();

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
  bool hasPolyhedronItem()const;

public Q_SLOTS:
  void reset_spheres();
  void updateDeform();
  void change();

  void invalidateOpenGLBuffers();
  void selected(const std::set<boost::graph_traits<Scene_polyhedron_item::Face_graph>::vertex_descriptor>&);
  void selected(const std::set<boost::graph_traits<Scene_surface_mesh_item::Face_graph>::vertex_descriptor>& );

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
  template<typename Mesh>
  bool insert_control_vertex(typename boost::graph_traits<Mesh>::vertex_descriptor v, Mesh* mesh);

  template<typename Mesh>
  bool insert_roi_vertex(typename boost::graph_traits<Mesh>::vertex_descriptor v, Mesh* mesh);

  //for calls from the plugin
  bool insert_roi_vertex(poly_vertex_descriptor vh);
  bool insert_roi_vertex(sm_vertex_descriptor vh);

  template<typename Mesh>
  bool erase_control_vertex(typename boost::graph_traits<Mesh>::vertex_descriptor v, Mesh* mesh);

  template<typename Mesh>
  bool erase_roi_vertex(typename boost::graph_traits<Mesh>::vertex_descriptor v, Mesh* mesh);

  void set_all_vertices_as_roi();

  void clear_roi();

  void create_ctrl_vertices_group();

  void delete_ctrl_vertices_group( bool create_new = true);

  void pivoting_end();
  void pivoting_begin();

  void prev_ctrl_vertices_group();
  void next_ctrl_vertices_group();

  void save_roi(const char* file_name) const;
  void read_roi(const char* file_name);
  void overwrite_deform_object();

  void reset_deform_object();
  template<typename Mesh>
  struct Is_selected {
    typedef typename CGAL::Surface_mesh_deformation<Mesh, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
      ,CGAL::Default, CGAL::Default, CGAL::Default,
      Array_based_vertex_point_map<Mesh> > M_Deform_mesh;

    M_Deform_mesh* dm;
    Is_selected(M_Deform_mesh* dm) : dm(dm) {}
    bool count(poly_vertex_descriptor vh) const {
      return dm->is_roi_vertex(vh);
    }

    bool count(sm_vertex_descriptor vh) const {
      return dm->is_roi_vertex(vh);
    }
  };

  boost::optional<std::size_t> get_minimum_isolated_component();
  template<typename Mesh>
  struct Select_roi_output{
    typedef typename CGAL::Surface_mesh_deformation<Mesh, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP
    ,CGAL::Default, CGAL::Default, CGAL::Default,
    Array_based_vertex_point_map<Mesh> > M_Deform_mesh;

    Select_roi_output(M_Deform_mesh* dm) : dm(dm) { }
    void operator()(poly_vertex_descriptor vh) {
      dm->insert_roi_vertex(vh);
    }
    void operator()(sm_vertex_descriptor vd) {
      dm->insert_roi_vertex(vd);
    }
    M_Deform_mesh* dm;
  };

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) ;

protected:
  // Deformation related functions //
  void print_message(const QString& /*message*/)
  {
    // std::cout << message.toStdString() << std::endl;
  }

  bool is_there_any_ctrl_vertices_group(Ctrl_vertices_poly_group_data_list::iterator& hgb,
                                        Ctrl_vertices_poly_group_data_list::iterator& hge);

  bool is_there_any_ctrl_vertices_group(Ctrl_vertices_sm_group_data_list::iterator& hgb,
                                        Ctrl_vertices_sm_group_data_list::iterator& hge);

  bool is_there_any_ctrl_vertices_group();

  bool is_there_any_ctrl_vertices();
  void refresh_all_group_centers();
  template<typename Mesh>
  bool activate_closest_manipulated_frame(int x, int y, Mesh* mesh);

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
  Id_setter * id_setter; //needed as a class member because of the deform_meshes
}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
