#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H
//#define CGAL_PROFILE 
#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <queue>

#include <QColor>
#include <QList>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/camera.h>

#include <CGAL/glu.h>

#include "ui_Deform_mesh.h"
#include <CGAL/Deform_mesh.h> 

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

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; 
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map; 

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP> Deform_mesh;


typedef Deform_mesh::Point  Point;
typedef Deform_mesh::Vector Vector;

/// For storing associated data with handle group
struct Handle_group_data
{
  Deform_mesh::Handle_group handle_group;
  qglviewer::ManipulatedFrame* frame;  // manframe assoc with handle group
  qglviewer::Vec frame_initial_center; // initial center of frame
  Scene_interface::Bbox bbox;          // bbox of handles inside group  


  std::vector<qglviewer::Vec> initial_positions;
  Deform_mesh* deform_mesh;

  Handle_group_data(Deform_mesh::Handle_group handle_group, Deform_mesh* deform_mesh, qglviewer::ManipulatedFrame* frame = 0)
    : handle_group(handle_group), frame(frame), bbox(0,0,0,0,0,0), deform_mesh(deform_mesh)
  { }

  void refresh()
  {
    reset_initial_positions();
    frame_initial_center = calculate_initial_center();
    bbox = calculate_initial_bbox();

    bool oldState = frame->blockSignals(true); // do not let it emit modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    frame->setOrientation(qglviewer::Quaternion());
    frame->setPosition(frame_initial_center);
    frame->blockSignals(oldState);
  }

  void reset_initial_positions()
  {
    initial_positions.clear();
    Deform_mesh::Handle_iterator hb, he;
    for(boost::tie(hb, he) = deform_mesh->handles(handle_group); hb != he; ++hb)
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

  void assign()
  {
    Deform_mesh::Handle_iterator hb, he;
    boost::tie(hb, he) = deform_mesh->handles(handle_group);
    for(std::vector<qglviewer::Vec>::iterator it = initial_positions.begin(); it != initial_positions.end(); ++it, ++hb)
    {
      qglviewer::Vec dif_from_initial_center = (*it) - frame_initial_center;
      qglviewer::Vec rotated = frame->orientation() * dif_from_initial_center;
      qglviewer::Vec rotated_and_translated = rotated + frame->position();

      deform_mesh->assign(*hb, Point(rotated_and_translated.x, rotated_and_translated.y, rotated_and_translated.z) );
    }
  }
};

// To hold pressing states together
struct Mouse_keyboard_state
{
  bool ctrl_pressing;
  bool shift_pressing;
  bool left_button_pressing;
  bool right_button_pressing;

  Mouse_keyboard_state() 
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
  Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item, Ui::DeformMesh* ui_widget);
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
  bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals); }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  void draw() const;
  void draw_bbox(const Scene_interface::Bbox& bb ) const;
  void gl_draw_edge(double px, double py, double pz,
                          double qx, double qy, double qz) const;
  void gl_draw_point(const Point& p) const;

  bool manipulatable() const { return true; }
  qglviewer::ManipulatedFrame* manipulatedFrame();

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  /// Returns a Scene_polyhedron_item from the edit polyhedron item, and
  /// transfer the ownership of the polyhedron to it.
  /// The item 'this' must be destroy just after a call to this function.
  Scene_polyhedron_item* to_polyhedron_item() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  // take mouse events from viewer, main-window does not work
  // take keyboard events from main-window, which is more stable
  bool eventFilter(QObject *target, QEvent *event);
  
protected:
  void timerEvent(QTimerEvent *event);

public slots:
  void changed();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
  
  void vertex_has_been_selected(void* void_ptr); // a vertex is selected by shift + left_click
  void deform();                                 // deform the mesh
  
signals:
  void mesh_deformed(Scene_edit_polyhedron_item* edit_item);        // emits when deformation is completed (currently handled by plugin)
  void mesh_repaint_needed(Scene_edit_polyhedron_item* edit_item);  // emits when there is a need for repaint (currently handled by plugin)
public:

  Ui::DeformMesh* ui_widget;

  Scene_polyhedron_item* poly_item;
  qglviewer::ManipulatedFrame* frame;

  Deform_mesh deform_mesh;
  Deform_mesh::Handle_group active_group;
typedef std::list<Handle_group_data> Handle_group_data_list;
  Handle_group_data_list handle_frame_map; // keep list of handle_groups with assoc data

  double length_of_axis; // for drawing axis at a handle group

  // by interleaving 'viewer's events (check constructor), keep followings:
  Mouse_keyboard_state state;

public:
  // Deformation related functions //
  void insert_handle(vertex_descriptor v)
  {
    if(!is_there_any_handle_group()) {
      print_message("There is no handle group, create one!");
      return; 
    } // no handle group to insert

    deform_mesh.insert_handle(active_group, v);

    get_data(active_group).refresh();
  }

  void insert_roi(vertex_descriptor v)
  {
    deform_mesh.insert_roi(v);
  }
  
  void erase_handle(vertex_descriptor v)
  {
    if(deform_mesh.erase_handle(v)) // API should be safe enough to do that (without checking empty handle group etc.)
    {
      refresh_all_handle_centers(); // since we don't know which handle group that v erased of, refresh all
    }

    if(!is_there_any_handle_group()) { 
      print_message("There is no handle group!");
      return; 
    }        
  }

  void erase_roi(vertex_descriptor v)
  {
    erase_handle(v); // erase from handles
    deform_mesh.erase_roi(v);
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
    active_group = deform_mesh.create_handle_group();

    qglviewer::ManipulatedFrame* new_frame = new qglviewer::ManipulatedFrame();
    new_frame->setRotationSensitivity(2.0f);

    Handle_group_data hgd(active_group, &deform_mesh, new_frame);
    handle_frame_map.push_back(hgd);
    hgd.refresh();

    connect(new_frame, SIGNAL(modified()), this, SLOT(deform()));  // OK we are deforming via timer,
    // but it makes demo more responsive if we also add this signal
    emit mesh_repaint_needed(this);

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
      if(it->handle_group == active_group)
      {
        delete it->frame;
        handle_frame_map.erase(it); 
        break;
      }
    }
    deform_mesh.erase_handle(active_group);

    // assign another handle_group to active_group
    Deform_mesh::Handle_group hgb, hge;
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
    Deform_mesh::Handle_group hgb, hge;
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
    Deform_mesh::Handle_group hgb, hge;
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
    hc = 0;
    Deform_mesh::Const_handle_group hgb, hge;
    for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb) 
    { ++hc; }
    out << hc << std::endl; // handle count
    for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb) {
      hc = 0;
      Deform_mesh::Handle_const_iterator hb, he;
      for(boost::tie(hb, he) = deform_mesh.handles(hgb); hb != he; ++hb) { ++hc; }
      out << hc << std::endl;
      for(boost::tie(hb, he) = deform_mesh.handles(hgb); hb != he; ++hb) 
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

  bool show_roi() const
  { return ui_widget->ShowROICheckBox->isChecked(); }

  bool show_as_sphere() const
  { return ui_widget->ShowAsSphereCheckBox->isChecked(); }

  void overwrite_deform_object()
  {
    deform_mesh.overwrite_original_positions();

    refresh_all_handle_centers();
  }


protected:
  // Deformation related functions //
  void print_message(const QString& /*message*/)
  {
    // std::cout << message.toStdString() << std::endl;
  }

  bool is_there_any_handle_group(Deform_mesh::Handle_group& hgb, Deform_mesh::Handle_group& hge)
  {
    boost::tie(hgb, hge) = deform_mesh.handle_groups();
    return hgb != hge;
  }

  bool is_there_any_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    return is_there_any_handle_group(hgb, hge);
  }

  bool is_there_any_handle()
  {
    Deform_mesh::Handle_group hgb, hge;
    if(!is_there_any_handle_group(hgb, hge)) { return false; } // there isn't any handle group

    for(; hgb != hge; ++hgb) // check inside handle groups
    {
      Deform_mesh::Handle_iterator hb, he;
      boost::tie(hb, he) = deform_mesh.handles(hgb);
      if(hb != he) { return true; }      
    }
    return false;
  }

  void process_selection(vertex_descriptor v, int k_ring, bool is_roi, bool is_insert, bool use_euclidean)
  {
    // std::cout << "Process k-ring: " << k_ring << " roi: " << is_roi << " insert: " << is_insert << std::endl;

    std::map<vertex_descriptor, int> neighs = use_euclidean ?
      extract_k_ring_with_distance(*polyhedron(), v, k_ring) :
      extract_k_ring(*polyhedron(), v, k_ring);

    for(std::map<vertex_descriptor, int>::iterator it = neighs.begin(); it != neighs.end(); ++it)
    {
      vertex_descriptor vh = it->first;
      if(is_roi) {
        if(is_insert) { insert_roi(vh); }
        else          { erase_roi(vh);  }
      }
      else {
        if(is_insert) { insert_handle(vh); }
        else          { erase_handle(vh);  }
      }
    }
  }

  void refresh_all_handle_centers()
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->refresh();
    }
  }

  Handle_group_data& get_data(Deform_mesh::Handle_group hg)
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it->handle_group == hg)
      {
        return *it;
      }
    }
    return *handle_frame_map.end(); // crash
    // this can't happen since every created handle group is inserted in handle_frame_map
  }

  Handle_group_data get_data(Deform_mesh::Const_handle_group hg) const
  {
    for(Handle_group_data_list::const_iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it->handle_group == hg)
      {
        return *it;
      }
    }
    return *handle_frame_map.end(); // crash 
    // this can't happen since every created handle group is inserted in handle_frame_map
  }

  std::map<vertex_descriptor, int> extract_k_ring(const Polyhedron &P, vertex_descriptor v, int k)
  {
    std::map<vertex_descriptor, int>  D;
    std::queue<vertex_descriptor>     Q;
    Q.push(v); D[v] = 0;

    int dist_v;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      }   
    }
    return D;
  }

  std::map<vertex_descriptor, int> extract_k_ring_with_distance(const Polyhedron &P, vertex_descriptor root, int k)
  {
    // tuple: <vertex, geodistance, euclidian distane>
    std::map<vertex_descriptor, int>  D;
    std::queue<vertex_descriptor>     Q;
    Q.push(root); D[root] = 0;

    int dist_v;
    double max_distance = 0.0;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          max_distance = (std::max)( (new_v->point() - root->point()).squared_length(), max_distance);
          Q.push(new_v);
        }
      }   
    }
    // now Q contains all nonprocessed 
    while( !Q.empty() ) {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        double distance = (new_v->point() - root->point()).squared_length();
        if(distance < max_distance)
        {
          if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
            Q.push(new_v);
          }
        }
      }   
    }

    return D;
  }

  void activate_closest_manipulated_frame(int x, int y)
  {
    if(state.ctrl_pressing && (state.left_button_pressing || state.right_button_pressing) ) 
    { // user is deforming currently don't change the state 
      return;  
    }

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    if(! state.ctrl_pressing) 
    {     
      viewer->setManipulatedFrame(frame);    
      return;
    }
    
    if(handle_frame_map.empty()) { return; }

    // now find closest frame and make it active manipulated frame
    std::list<Handle_group_data>::iterator min_it = handle_frame_map.begin();    
    const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(min_it->frame->position());
    float min_dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);

    for(std::list<Handle_group_data>::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->frame->position());
      float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
      if(dist < min_dist) {
        min_dist = dist;
        min_it = it;
      }
    }

    viewer->setManipulatedFrame(min_it->frame);  
  }

protected:
  GLUquadric* quadric; // for drawing spheres
}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
