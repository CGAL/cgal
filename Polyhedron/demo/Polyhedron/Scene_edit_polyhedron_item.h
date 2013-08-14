#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H
//#define CGAL_PROFILE 
#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_item_k_ring_selection.h"
#include "Travel_isolated_components.h"

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/compute_normal.h>

#include <iostream>
#include <fstream>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/camera.h>

#include "ui_Deform_mesh.h"
#include <CGAL/Deform_mesh.h> 
#include <boost/function_output_iterator.hpp>

typedef Polyhedron::Vertex_handle Vertex_handle;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor		vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor		  edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
typedef boost::graph_traits<Polyhedron>::in_edge_iterator		  in_edge_iterator;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;
        
    reference operator[](key_type key) const { return key->id(); }
};

struct Array_based_vertex_point_map
{
public:
  typedef vertex_descriptor            key_type;
  typedef Polyhedron::Traits::Point_3  value_type;
  typedef value_type&  reference;
  typedef boost::read_write_property_map_tag category;
  Array_based_vertex_point_map(std::vector<double>* positions) : positions(positions) {}
  std::vector<double>* positions;
};


Array_based_vertex_point_map::value_type
get(Array_based_vertex_point_map,
  Array_based_vertex_point_map::key_type key) {
    return key->point();
}

void
put(Array_based_vertex_point_map pmap,
  Array_based_vertex_point_map::key_type key,
  Array_based_vertex_point_map::value_type val) {
  key->point() = val; // to make things easy (ray selection after deformation, save to polyhedron after close etc),
  // I also change point() of vertex together with positions list
  // So that we do not need to pmap everywhere other than draw
  std::size_t pos = key->id() * 3;
  (*pmap.positions)[pos] = val.x();
  (*pmap.positions)[pos+1] = val.y();
  (*pmap.positions)[pos+2] = val.z();
}

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; 
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map; 

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP
  ,CGAL::Default, CGAL::Default, CGAL::Default, 
  Array_based_vertex_point_map> Deform_mesh;


typedef Deform_mesh::Point  Point;

/// For storing associated data with handle group
class Handle_group_data
{
public:
  std::vector<vertex_descriptor> handle_group;
  qglviewer::ManipulatedFrame* frame;  // manframe assoc with handle group
  qglviewer::Vec frame_initial_center; // initial center of frame
  Scene_interface::Bbox bbox;          // bbox of handles inside group  
  qglviewer::Vec rot_direction;        // vector for constraint rotation
private:
  std::vector<qglviewer::Vec> initial_positions;
  Deform_mesh* deform_mesh;

public:
  Handle_group_data(Deform_mesh* deform_mesh, qglviewer::ManipulatedFrame* frame = 0)
    : frame(frame), bbox(0,0,0,0,0,0), rot_direction(0.,0.,1.), deform_mesh(deform_mesh)
  { }
  void refresh()
  {
    for(std::vector<vertex_descriptor>::iterator it = handle_group.begin(); it != handle_group.end(); ) {
      if(!deform_mesh->is_control(*it)) {
        it = handle_group.erase(it);
      }
      else { ++it; }
    }

    reset_initial_positions();
    frame_initial_center = calculate_initial_center();
    bbox = calculate_initial_bbox();

    bool oldState = frame->blockSignals(true); // do not let it emit modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    frame->setOrientation(qglviewer::Quaternion());
    frame->setPosition(frame_initial_center);
    frame->blockSignals(oldState);
  }
  void assign()
  {
    std::vector<vertex_descriptor>::iterator hb = handle_group.begin();
    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it, ++hb)
    {
      qglviewer::Vec dif_from_initial_center = (*it) - frame_initial_center;
      qglviewer::Vec rotated = frame->orientation() * dif_from_initial_center;
      qglviewer::Vec rotated_and_translated = rotated + frame->position();

      deform_mesh->assign(*hb, Point(rotated_and_translated.x, rotated_and_translated.y, rotated_and_translated.z) );
    }
  }

private:
  void reset_initial_positions()
  {
    initial_positions.clear();
    
    for(std::vector<vertex_descriptor>::iterator hb = handle_group.begin(); hb != handle_group.end(); ++hb)
    {
      qglviewer::Vec point((*hb)->point().x(), (*hb)->point().y(), (*hb)->point().z() );
      initial_positions.push_back(point);
    }
  }
  qglviewer::Vec calculate_initial_center()
  {
    qglviewer::Vec center_acc(0, 0, 0);
    if(initial_positions.empty()) {return center_acc; }

    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it)
    {
      center_acc += (*it);
    }
    return center_acc / initial_positions.size();
  }
  Scene_interface::Bbox calculate_initial_bbox()
  {    
    if(initial_positions.empty()) {return Scene_interface::Bbox(0,0,0,0,0,0); }

    const qglviewer::Vec& p_i = *(initial_positions.begin());
    Scene_interface::Bbox bbox(p_i.x, p_i.y, p_i.z, p_i.x, p_i.y, p_i.z);

    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it)
    {
      const qglviewer::Vec& p_i = (*it);
      Scene_interface::Bbox bbox_it(p_i.x, p_i.y, p_i.z, p_i.x, p_i.y, p_i.z);
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

// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_EXPORT Scene_edit_polyhedron_item 
  : public Scene_item {
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
  void draw() const;
  void draw_edges() const;
  void draw_bbox(const Scene_interface::Bbox& bb ) const;
  void gl_draw_edge(double px, double py, double pz,
                          double qx, double qy, double qz) const;
  void gl_draw_point(const Point& p) const;

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
  Bbox bbox() const;

  int get_k_ring()       { return k_ring_selector.k_ring; }
  void set_k_ring(int v) { k_ring_selector.k_ring = v; }

  // take mouse events from viewer, main-window does not work
  // take keyboard events from main-window, which is more stable
  bool eventFilter(QObject *target, QEvent *event);
  
protected:
  void timerEvent(QTimerEvent *event);
  void draw_ROI_and_handles() const;

public slots:
  void changed();
  void selected(const std::map<Polyhedron::Vertex_handle, int>& m)
  {
    bool any_changes = false;
    for(std::map<vertex_descriptor, int>::const_iterator it = m.begin(); it != m.end(); ++it)
    {
      vertex_descriptor vh = it->first;
      bool changed = false;
      if(ui_widget->ROIRadioButton->isChecked()) {
        if(ui_widget->InsertRadioButton->isChecked()) { changed = insert_roi(vh); }
        else          { changed = erase_roi(vh);  }
      }
      else {
        if(ui_widget->InsertRadioButton->isChecked()) { changed = insert_handle(vh); }
        else          { changed = erase_handle(vh);  }
      }
      any_changes |= changed;
    }
    if(any_changes) { emit itemChanged(); }
  }

  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);

  void deform(); // deform the mesh
// members
private:
  Ui::DeformMesh* ui_widget;
  Scene_polyhedron_item* poly_item;
  // For drawing
  std::vector<double> positions;
  std::vector<unsigned int> tris;
  std::vector<unsigned int> edges;
  std::vector<double> normals;

  Deform_mesh deform_mesh;
typedef std::list<Handle_group_data> Handle_group_data_list;
  Handle_group_data_list::iterator active_group;
  Handle_group_data_list handle_frame_map; // keep list of handle_groups with assoc data

  double length_of_axis; // for drawing axis at a handle group

  // by interleaving 'viewer's events (check constructor), keep followings:
  Mouse_keyboard_state_deformation state;

  //For constraint rotation
  qglviewer::LocalConstraint rot_constraint;
  bool is_rot_free;

  bool own_poly_item; //indicates if the poly_item should be deleted by the destructor
  Scene_polyhedron_item_k_ring_selection k_ring_selector;

public:
  // Deformation related functions //
  bool insert_handle(vertex_descriptor v)
  {
    if(!is_there_any_handle_group()) {
      print_message("There is no handle group, create one!");
      return false; 
    } // no handle group to insert

    bool inserted = deform_mesh.insert_control(v);
    if(inserted) {
      active_group->handle_group.push_back(v);
      active_group->refresh();
    }
    return inserted;
  }

  bool insert_roi(vertex_descriptor v)
  {
    return deform_mesh.insert_roi(v);
  }
  
  bool erase_handle(vertex_descriptor v)
  {
    if(deform_mesh.erase_control(v)) // API should be safe enough to do that (without checking empty handle group etc.)
    {
      refresh_all_handle_centers(); // since we don't know which handle group that v erased of, refresh all
      return true;
    }

    print_message("Selected vertex is not handle!");
    return false;     
  }

  bool erase_roi(vertex_descriptor v)
  {
    erase_handle(v); // erase from handles
    return deform_mesh.erase_roi(v);
  }

  void set_all_vertices_as_roi()
  {
    vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = boost::vertices(*polyhedron()); vb != ve; ++vb)
    {
      insert_roi(*vb);
    }   
  }

  void clear_roi()
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      delete it->frame;
    }
    handle_frame_map.clear();
    deform_mesh.reset();

    create_handle_group(); // create one new handle group
  } 

  void create_handle_group()
  {
    qglviewer::ManipulatedFrame* new_frame = new qglviewer::ManipulatedFrame();
    new_frame->setRotationSensitivity(2.0f);

    Handle_group_data hgd(&deform_mesh, new_frame);
    handle_frame_map.push_back(hgd);
    hgd.refresh();

    active_group = --handle_frame_map.end();

    connect(new_frame, SIGNAL(modified()), this, SLOT(deform()));  // OK we are deforming via timer,
    // but it makes demo more responsive if we also add this signal
    emit itemChanged();

    print_message("A new empty handle group is created.");
  }

  void delete_handle_group(bool create_new = true)
  {
    if(!is_there_any_handle_group()) { 
      print_message("There is no handle group to be deleted!");
      return; 
    } // no handle group

    // delete group representative    
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it == active_group)
      {
        delete it->frame;        
        for(std::vector<vertex_descriptor>::iterator v_it = it->handle_group.begin(); v_it != it->handle_group.end(); ++v_it) {
          deform_mesh.erase_control(*v_it);
        }
        handle_frame_map.erase(it);
        break;
      }
    }

    // assign another handle_group to active_group
    Handle_group_data_list::iterator hgb, hge;
    if( is_there_any_handle_group(hgb, hge) )
    { 
      active_group = hgb; 
    } // no handle group 
    else if(create_new)
    { 
      create_handle_group(); 
    }
  }

  void prev_handle_group()
  {
    Handle_group_data_list::iterator hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) ) {
      print_message("There is no handle group to iterate on!");
      return; 
    }
    // shift
    if(hgb == active_group) { active_group = --hge; }
    else                    {--active_group; }    
  }

  void next_handle_group()
  {
    Handle_group_data_list::iterator hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) ) {
      print_message("There is no handle group to iterate on!");
      return; 
    }
    // shift
    if(--hge == active_group) { active_group = hgb; }
    else                      {++active_group; }    
  }

  void pivoting_end()
  {       
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      //update constraint rotation vector, set only for the last group
      it->rot_direction = it->frame->rotation().rotate( qglviewer::Vec(0.,0.,1.) );
      //translate center of the frame
      qglviewer::Vec vec= it->frame->position();
      it->refresh();
      it->frame_initial_center = vec;
      it->frame->setPosition(vec);
    }
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame->blockSignals(false);
    }
  }

  void pivoting_begin()
  {
    is_rot_free=true;
    rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);

    // just block signals to prevent deformation
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame->blockSignals(true);
    }
  }

  void save_roi(const char* file_name) const
  { 
    std::ofstream out(file_name);
    // save roi
    int hc = 0;
    Deform_mesh::Roi_const_iterator rb, re;
    for(boost::tie(rb, re) = deform_mesh.roi_vertices(); rb != re; ++rb) { ++hc; }
    out << hc << std::endl;
    for(boost::tie(rb, re) = deform_mesh.roi_vertices(); rb != re; ++rb)
    {
      out << (*rb )->id() << " ";
    }
    out << std::endl;
    // save handles
    
    out << handle_frame_map.size() << std::endl; // handle count
    for(Handle_group_data_list::const_iterator hgb = handle_frame_map.begin(); hgb != handle_frame_map.end(); ++hgb) {

      out << hgb->handle_group.size() << std::endl;
      for(std::vector<vertex_descriptor>::const_iterator hb = hgb->handle_group.begin(); hb != hgb->handle_group.end(); ++hb) 
      {
        out << (*hb)->id() << " ";
      }
      out << std::endl;
    }
  }

  void read_roi(const char* file_name)
  { 
    clear_roi();
    delete_handle_group(false);

    // put vertices to vector
    std::vector<vertex_descriptor> all_vertices;
    all_vertices.reserve(boost::num_vertices(deform_mesh.halfedge_graph()));
    vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = boost::vertices(deform_mesh.halfedge_graph()); vb != ve; ++vb) {
      all_vertices.push_back(*vb);
    }
    // read roi
    std::ifstream in(file_name);
    int roi_size;
    in >> roi_size;
    while(roi_size-- > 0)
    {
      std::size_t v_id;
      in >> v_id;
      insert_roi(all_vertices[v_id]);
    }
    // read handles
    int handle_group_size;
    in >> handle_group_size;
    while(handle_group_size-- > 0)
    {
      create_handle_group();
      int handle_size;
      in >> handle_size;      
      while(handle_size-- > 0) 
      {                    
        std::size_t v_id;
        in >> v_id;
        insert_handle(all_vertices[v_id]);
      }
    }
  }

  void overwrite_deform_object()
  {
    deform_mesh.overwrite_original_positions();

    refresh_all_handle_centers();
  }

  struct Is_selected {
    Deform_mesh& dm;
    Is_selected(Deform_mesh& dm) : dm(dm) {}
    bool is_selected(Vertex_handle vh) const {
      return dm.is_roi(vh);
    }
  };

  boost::optional<std::size_t> get_minimum_isolated_component() {
    Travel_isolated_components::Minimum_visitor visitor;
    Travel_isolated_components().travel<Vertex_handle>
      (polyhedron()->vertices_begin(), polyhedron()->vertices_end(), 
       polyhedron()->size_of_vertices(), Is_selected(deform_mesh), visitor);
    return visitor.minimum;
  }

  struct Select_roi_output {
    Select_roi_output(Deform_mesh* dm) : dm(dm) { }
    void operator()(Vertex_handle vh) {
      dm->insert_roi(vh);
    }
    Deform_mesh* dm;
  };

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    typedef boost::function_output_iterator<Select_roi_output> Output_iterator;
    Output_iterator out(&deform_mesh);

    Travel_isolated_components::Selection_visitor<Output_iterator> visitor(threshold, out);
    Travel_isolated_components().travel<Vertex_handle>
      (polyhedron()->vertices_begin(), polyhedron()->vertices_end(), 
      polyhedron()->size_of_vertices(), Is_selected(deform_mesh), visitor);

    if(visitor.any_inserted) { emit itemChanged(); }
    return visitor.minimum_visitor.minimum;
  }
protected:
  // Deformation related functions //
  void print_message(const QString& /*message*/)
  {
    // std::cout << message.toStdString() << std::endl;
  }

  bool is_there_any_handle_group(Handle_group_data_list::iterator& hgb, Handle_group_data_list::iterator& hge)
  {
    hgb = handle_frame_map.begin(); hge = handle_frame_map.end();
    return hgb != hge;
  }

  bool is_there_any_handle_group()
  {
    Handle_group_data_list::iterator hgb, hge;
    return is_there_any_handle_group(hgb, hge);
  }

  bool is_there_any_handle()
  {
    Handle_group_data_list::iterator hgb, hge;
    if(!is_there_any_handle_group(hgb, hge)) { return false; } // there isn't any handle group

    for(; hgb != hge; ++hgb) // check inside handle groups
    {
      if(!hgb->handle_group.empty()) { return true; }
    }
    return false;
  }

  void refresh_all_handle_centers()
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    { it->refresh(); }
  }

  bool activate_closest_manipulated_frame(int x, int y)
  {
    if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing) ) 
    { // user is deforming currently don't change the state 
      return false;  
    }
    if(handle_frame_map.empty()) { return false; }

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    if(!state.ctrl_pressing) 
    {   
      if(viewer->manipulatedFrame() == NULL) 
      { return false;}
      viewer->setManipulatedFrame(NULL);    
      return true;
    }
    
    // now find closest frame and make it active manipulated frame
    Handle_group_data_list::iterator min_it = handle_frame_map.begin();    
    const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(min_it->frame->position());
    float min_dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);

    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->frame->position());
      float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
      if(dist < min_dist) {
        min_dist = dist;
        min_it = it;
      }
    }

    //set rotation constraint for the manipulated frame
    if (!is_rot_free){
      rot_constraint.setRotationConstraintDirection(min_it->rot_direction);
      rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::AXIS);
      min_it->frame->setConstraint(&rot_constraint);
    }
    else
      rot_constraint.setRotationConstraintType(qglviewer::AxisPlaneConstraint::FREE);

    if(viewer->manipulatedFrame() == min_it->frame)
    { return false; }
    viewer->setManipulatedFrame(min_it->frame);

    return true;
  }

  bool keyPressEvent(QKeyEvent* e);

  void update_normals() {
    Deform_mesh::Roi_const_iterator rb, re;
    for(boost::tie(rb, re) = deform_mesh.roi_vertices(); rb != re; ++rb) {
      std::size_t id = (*rb)->id();
      const Polyhedron::Traits::Vector_3& n = 
        compute_vertex_normal<Polyhedron::Vertex, Polyhedron::Traits>(**rb);
      normals[id*3] = n.x();
      normals[id*3+1] = n.y(); 
      normals[id*3+2] = n.z(); 
    }
  }
protected:
  GLUquadric* quadric; // for drawing spheres
}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
