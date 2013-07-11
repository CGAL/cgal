#ifndef SCENE_POLYHEDRON_SELECTION_ITEM_H
#define SCENE_POLYHEDRON_SELECTION_ITEM_H

#include "Scene_polyhedron_selection_item_config.h"
#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_type.h"
#include "ui_Selection_widget.h"
#include "opengl_tools.h"
#include <CGAL/gl_render.h>

#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QMouseEvent>

#include <queue>
#include <vector>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
// Bunch of structs used in Scene_polyhedron_selection_item //

// For iterating on one ring neighbor
template<class T>
struct One_ring_iterator {};

template<>
struct One_ring_iterator<Polyhedron::Vertex_handle> {
  One_ring_iterator(Polyhedron::Vertex_handle v) 
    : circ(v->vertex_begin()), end(circ), first(true) { }

  operator bool() { return first || circ != end; }
  operator Polyhedron::Vertex_handle() { return circ->opposite()->vertex(); }
  One_ring_iterator& operator++() { 
    first = false;
    ++circ; 
    return *this;
  }

  Polyhedron::Halfedge_around_vertex_circulator circ;
  Polyhedron::Halfedge_around_vertex_circulator end;
  bool first;
  // to be used in One_ring_iterator<Halfedge_handle>
  operator Polyhedron::Halfedge_handle() { 
    return &*circ < &*circ->opposite() ? circ : circ->opposite(); 
  }
};

template<>
struct One_ring_iterator<Polyhedron::Facet_handle> {
  One_ring_iterator(Polyhedron::Facet_handle f)
    : circ(f->facet_begin()), end(circ), first(true)
  {
    iterate_to_non_border(); // move it to valid location
  }

  operator bool() { return first || circ != end; }
  operator Polyhedron::Facet_handle() {
    CGAL_assertion(!circ->opposite()->is_border());
    return circ->opposite()->facet(); 
  }
  One_ring_iterator& operator++() {
    first = false;
    ++circ;
    if(circ != end) { iterate_to_non_border(); }
    return *this;
  }
  
  void iterate_to_non_border() {
    while(circ->opposite()->is_border()) {
      first = false;
      ++circ;
      if(circ == end) { break; }
    }
  }
  Polyhedron::Halfedge_around_facet_circulator circ;
  Polyhedron::Halfedge_around_facet_circulator end;
  bool first;
};

template<>
struct One_ring_iterator<Polyhedron::Halfedge_handle> {
  One_ring_iterator(Polyhedron::Halfedge_handle h)
    : it_1(h->vertex()), it_2(h->opposite()->vertex())
  { }

  operator bool() { return it_1 || it_2; }
  operator Polyhedron::Halfedge_handle() {
     if(it_1) { return it_1; }
     return it_2;
  }
  One_ring_iterator& operator++() {
    it_1 ? ++it_1 : ++it_2;
    return *this;
  }

  One_ring_iterator<Polyhedron::Vertex_handle> it_1;
  One_ring_iterator<Polyhedron::Vertex_handle> it_2;
};

// to be used in get_minimum_isolated_component function
struct Minimum_visitor
{
  template<class HandleType>
  void operator()(const std::vector<HandleType>& C) {
    if(!minimum) { minimum = C.size(); }
    else         { minimum = (std::min)(*minimum, C.size()); }
  }

  boost::optional<std::size_t> minimum;
};

// to be used in select_isolated_components function
template<class OutputIterator>
struct Selection_visitor
{
  Selection_visitor(std::size_t threshold_size, OutputIterator out) 
    : threshold_size(threshold_size), out(out), any_inserted(false) { }

  template<class HandleType>
  void operator()(const std::vector<HandleType>& C) {
    if(C.size() <= threshold_size) {
      any_inserted = true;
      out = std::copy(C.begin(), C.end(), out);
    }
    else {
      minimum_visitor(C);
    }
  }

  std::size_t     threshold_size;
  OutputIterator  out;
  bool            any_inserted;
  Minimum_visitor minimum_visitor; // hold minimum of NOT inserted components
};

// Hold mouse keyboard state together
struct Mouse_keyboard_state
{
  Mouse_keyboard_state() : shift_pressing(false), left_button_pressing(false) { }
  bool shift_pressing, left_button_pressing;
};

// Wrapper for holding selected entities
template <class Entity, class Listener>
class Selection_set : public std::set<Entity>
{
private:
  typedef std::set<Entity>                Base;
  typedef Selection_set<Entity, Listener> Self;

public:
  Selection_set(Listener* listener = NULL) : listener(listener)
  { }

  struct Selection_set_inserter {
    Selection_set_inserter(Self* container) : container(container) { }
    void operator()(const Entity & e) const {
      container->insert(e);
    }
    Self* container;
  };
// types from base
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
// functions from base
  using Base::begin;
  using Base::end;
  using Base::empty;
  using Base::clear;

  bool insert(const Entity& entity) {
    if(listener != NULL) { listener->inserted(entity); } // no matter entity is inserted or not, call notifier
    return Base::insert(entity).second;
  }

  bool erase(const Entity& entity) {
    if(listener != NULL) { listener->erased(entity); } // no matter entity is erased or not, call notifier
    return Base::erase(entity) != 0;
  }

  bool is_selected(const Entity& entity) {
    return Base::find(entity) != end();
  }

// members
  Listener* listener;
};

//
struct Minimum_address_halfedge_iterator {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  Minimum_address_halfedge_iterator(Halfedge_iterator hb, Halfedge_iterator he)
    : hb(hb), he(he), current(hb) { }

  Minimum_address_halfedge_iterator& operator++() {
    ++current;
    while(current != he && (&*current > &*current->opposite())) {
      ++current;
    }
    return *this;
  }
  operator Polyhedron::Halfedge_handle() { return current; }
  bool operator!=(const Minimum_address_halfedge_iterator& other) {
    return current != other.he;
  }
  Halfedge_iterator operator->() {return current;}

  Halfedge_iterator hb, he, current;
};
//////////////////////////////////////////////////////////////////////////

class SCENE_POLYHEDRON_SELECTION_ITEM_EXPORT Scene_polyhedron_selection_item 
  : public Scene_polyhedron_item_decorator
{
  Q_OBJECT

public:
  typedef Polyhedron::Vertex_handle Vertex_handle;
  typedef Polyhedron::Facet_handle  Facet_handle;
  typedef Polyhedron::Halfedge_handle Halfedge_handle;

  Scene_polyhedron_selection_item(Scene_polyhedron_item* poly_item, Ui::Selection* ui_widget) 
    : Scene_polyhedron_item_decorator(poly_item, false), 
      ui_widget(ui_widget),
      selected_vertices(this),
      selected_facets(this),
      selected_edges(this)
  { 
    connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_facet(void*)), this, SLOT(facet_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_edge(void*)), this, SLOT(edge_has_been_selected(void*)));
    poly_item->enable_facets_picking(true);
    poly_item->set_color_vector_read_only(true);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    
    // id field is used in operations
    poly_item->update_vertex_indices();
    poly_item->update_facet_indices();
    poly_item->update_halfedge_indices();
  }
// drawing
  void draw() const {
    draw_selected_vertices();
    draw_selected_facets();
    draw_selected_edges();
  }

  void draw_selected_vertices() const {
    GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);

    CGAL::GL::Color color;
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);
    color.set_rgb_color(0, 1.f, 0);

    ::glBegin(GL_POINTS);
    for(Selection_set_vertex::iterator 
      it = selected_vertices.begin(),
      end = selected_vertices.end();
    it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();

    if(enable_back_lighting) { glEnable(GL_LIGHTING); }
  }

  void draw_selected_facets() const {
    CGAL::GL::Color color;
    color.set_rgb_color(0.f,1.f,0.f);

    GLfloat offset_factor;
    GLfloat offset_units;
    ::glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
    ::glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);

    ::glPolygonOffset(-1.f, 1.f);
    ::glBegin(GL_TRIANGLES);
    for(Selection_set_facet::iterator
          it = selected_facets.begin(),
          end = selected_facets.end();
        it != end; ++it)
    {
      const Kernel::Vector_3 n =
        compute_facet_normal<Polyhedron::Facet,Kernel>(**it);
      ::glNormal3d(n.x(),n.y(),n.z());

      Polyhedron::Halfedge_around_facet_circulator
        he = (*it)->facet_begin(),
        cend = he;

      CGAL_For_all(he,cend)
      {
        const Kernel::Point_3& p = he->vertex()->point();
        ::glVertex3d(p.x(),p.y(),p.z());
      }
    }
    ::glEnd();
    ::glPolygonOffset(offset_factor, offset_units);
  }

  void draw_selected_edges() const {
    GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);

    CGAL::GL::Color color;
    color.set_rgb_color(0.f,1.f,0.f);
    ::glLineWidth(3.f);
    ::glBegin(GL_LINES);
    for(Selection_set_edge::iterator it = selected_edges.begin(); it != selected_edges.end(); ++it) {
      const Kernel::Point_3& a = (*it)->vertex()->point();
      const Kernel::Point_3& b = (*it)->opposite()->vertex()->point();
      ::glVertex3d(a.x(),a.y(),a.z());
      ::glVertex3d(b.x(),b.y(),b.z());
    }
    ::glEnd();

    if(enable_back_lighting) { glEnable(GL_LIGHTING); }
  }

  bool supportsRenderingMode(RenderingMode m) const { return (m==Flat); }
  //void save_roi(const char* file_name) const { 
  //  std::ofstream out(file_name);
  //  // save roi
  //  for(std::set<Vertex_handle>::iterator it = selected_vertices.begin();
  //    it != selected_vertices.end(); ++it) {
  //      out << (*it)->id() << " ";
  //  }
  //  out.close();
  //}
  //void load_roi(const char* file_name) {
  //  // put vertices to vector
  //  std::vector<Polyhedron::Vertex_handle> all_vertices;
  //  all_vertices.reserve(polyhedron()->size_of_vertices());
  //  Polyhedron::Vertex_iterator vb(polyhedron()->vertices_begin()), ve(polyhedron()->vertices_end());
  //  for( ;vb != ve; ++vb) {
  //    all_vertices.push_back(vb);
  //  }
  //  // read roi
  //  std::ifstream in(file_name);
  //  std::size_t idx;
  //  while(in >> idx) {
  //    selected_vertices.insert(all_vertices[idx]);
  //  }
  //  in.close();
  //}

  void select_all() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      Polyhedron::Vertex_iterator vb(polyhedron()->vertices_begin()), ve(polyhedron()->vertices_end());
      for( ;vb != ve; ++vb) {
        selected_vertices.insert(vb);
      }
    }
    else if(ui_widget->Selection_type_combo_box->currentIndex() == 1){
      Polyhedron::Facet_iterator fb(polyhedron()->facets_begin()), fe(polyhedron()->facets_end());
      for( ;fb != fe; ++fb) {
        selected_facets.insert(fb);
      }
    }
    else {
      Minimum_address_halfedge_iterator hb(polyhedron()->halfedges_begin(), polyhedron()->halfedges_end()), he(hb);
      for( ;hb != he; ++hb) {
        selected_edges.insert(hb);
      }
    }
    emit itemChanged();
  }

  void clear() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      selected_vertices.clear();
    }
    else if(ui_widget->Selection_type_combo_box->currentIndex() == 1) {
      selected_facets.clear();
    }
    else {
      selected_edges.clear();
    }
    emit itemChanged();
  }

  void get_minimum_isolated_component() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      get_minimum_isolated_vertex_component();
    }
    else if(ui_widget->Selection_type_combo_box->currentIndex() == 1) {
      get_minimum_isolated_facet_component();
    }
    else {
      get_minimum_isolated_edge_component();
    }
  }

  void select_marked_edges() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      selected_vertices.clear();
    } else {
      selected_facets.clear();
    }
    for(Polyhedron::Edge_iterator
          eit = polyhedron()->edges_begin(),
          end = polyhedron()->edges_end(); eit != end; ++eit)
    {
      if(!eit->is_feature_edge()) continue;
      if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
        selected_vertices.insert(eit->vertex());
        selected_vertices.insert(eit->opposite()->vertex());
      } else {
        selected_facets.insert(eit->face());
        selected_facets.insert(eit->opposite()->face());
      }
    }
    for(int i = 0; i < ui_widget->neighb_size->value(); ++i) {
      if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
        std::set<Vertex_handle> new_set;//(selected_vertices);
        BOOST_FOREACH(Vertex_handle v, selected_vertices)
        {
          Polyhedron::Halfedge_around_vertex_circulator
            he_circ = v->vertex_begin(), end = he_circ;
          if(he_circ != NULL) do {
              new_set.insert(he_circ->opposite()->vertex());
            } while (++he_circ != end);
        }
        BOOST_FOREACH(Vertex_handle v, new_set) {
          selected_vertices.insert(v);
        }
      } else {
        std::set<Facet_handle> new_set;//(selected_facets);
        BOOST_FOREACH(Facet_handle f, selected_facets)
        {
          Polyhedron::Halfedge_around_facet_circulator
            he_circ = f->facet_begin(), end = he_circ;
          if(he_circ != NULL) do {
              new_set.insert(he_circ->opposite()->facet());
            } while (++he_circ != end);
        }
        BOOST_FOREACH(Facet_handle f, new_set) {
          selected_facets.insert(f);
        }
      }
    }
  }

  void get_minimum_isolated_vertex_component() {
    Minimum_visitor visitor;
    travel_not_selected_connected_components<Vertex_handle>
      (visitor, polyhedron()->vertices_begin(), polyhedron()->vertices_end(), polyhedron()->size_of_vertices());
    if(visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum);
    }
  }
  void get_minimum_isolated_facet_component() {
    Minimum_visitor visitor;
    travel_not_selected_connected_components<Facet_handle>
      (visitor, polyhedron()->facets_begin(), polyhedron()->facets_end(), polyhedron()->size_of_facets());
    if(visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum);
    }
  }
  void get_minimum_isolated_edge_component() {
    Minimum_visitor visitor;
    Minimum_address_halfedge_iterator hb(polyhedron()->halfedges_begin(), polyhedron()->halfedges_end()), he(hb);
    travel_not_selected_connected_components<Halfedge_handle>
      (visitor, hb, he, polyhedron()->size_of_halfedges()); // do not divide by 2
    if(visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum);
    }
  }

  void select_isolated_components() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      select_isolated_vertex_components();
    }
    else if(ui_widget->Selection_type_combo_box->currentIndex() == 1){
      select_isolated_facet_components();
    }
    else {
      select_isolated_edge_components();
    }
  }

  void select_isolated_vertex_components() {
    typedef boost::function_output_iterator<
      Selection_set_vertex::Selection_set_inserter> Output_iterator;
    Output_iterator out(&selected_vertices);

    Selection_visitor<Output_iterator> visitor(ui_widget->Threshold_size_spin_box->value() , out);
    travel_not_selected_connected_components<Vertex_handle>
      (visitor, polyhedron()->vertices_begin(), polyhedron()->vertices_end(), polyhedron()->size_of_vertices());

    if(visitor.any_inserted) { emit itemChanged(); }
    if(visitor.minimum_visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum_visitor.minimum);
    }
  }
  void select_isolated_facet_components() {
    typedef boost::function_output_iterator<
      Selection_set_facet::Selection_set_inserter> Output_iterator;
    Output_iterator out(&selected_facets);

    Selection_visitor<Output_iterator> visitor(ui_widget->Threshold_size_spin_box->value() , out);
    travel_not_selected_connected_components<Facet_handle>
      (visitor, polyhedron()->facets_begin(), polyhedron()->facets_end(), polyhedron()->size_of_facets());

    if(visitor.any_inserted) { emit itemChanged(); }
    if(visitor.minimum_visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum_visitor.minimum);
    }
  }

  void select_isolated_edge_components() {
    typedef boost::function_output_iterator<
      Selection_set_edge::Selection_set_inserter> Output_iterator;
    Output_iterator out(&selected_edges);
    Minimum_address_halfedge_iterator hb(polyhedron()->halfedges_begin(), polyhedron()->halfedges_end()), he(hb);
    Selection_visitor<Output_iterator> visitor(ui_widget->Threshold_size_spin_box->value() , out);
    travel_not_selected_connected_components<Halfedge_handle>
      (visitor, hb, he, polyhedron()->size_of_halfedges()); // do not divide by 2

    if(visitor.any_inserted) { emit itemChanged(); }
    if(visitor.minimum_visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum_visitor.minimum);
    }
  }

  void changed_with_poly_item() {
    poly_item->changed();
    emit itemChanged();
  }

  bool is_selected(Vertex_handle v) { return selected_vertices.is_selected(v); }
  bool is_selected(Facet_handle f) { return selected_facets.is_selected(f); }
  bool is_selected(Halfedge_handle h) { return selected_edges.is_selected(h); }

  bool insert(Vertex_handle v) { return selected_vertices.insert(v); }
  bool insert(Facet_handle f) { return selected_facets.insert(f); }
  bool insert(Halfedge_handle h) { return selected_edges.insert(h); }

  bool erase(Vertex_handle v) { return selected_vertices.erase(v); }
  bool erase(Facet_handle f) { return selected_facets.erase(f); }
  bool erase(Halfedge_handle h) { return selected_edges.erase(h); }

public slots:
  void changed() {
    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
  }
  void vertex_has_been_selected(void* void_ptr) 
  {
    if(!visible() || ui_widget->Selection_type_combo_box->currentIndex() != 0) { return; }
    has_been_selected( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
  }
  void facet_has_been_selected(void* void_ptr)
  {
    if(!visible() || ui_widget->Selection_type_combo_box->currentIndex() != 1) { return; }
    has_been_selected( static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet() );
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    if(!visible() || ui_widget->Selection_type_combo_box->currentIndex() != 2) { return; }
    has_been_selected( static_cast<Polyhedron::Halfedge*>(void_ptr)->opposite()->opposite() );
  }

signals:
  void inserted(Vertex_handle v);
  void inserted(Facet_handle f);
  void inserted(Halfedge_handle f);

  void erased(Vertex_handle v);
  void erased(Facet_handle f);
  void erased(Halfedge_handle f);

protected:
  template<class HandleType>
  void has_been_selected(HandleType clicked) {
    bool is_insert = ui_widget->Insertion_radio_button->isChecked();
    int k_ring = ui_widget->Brush_size_spin_box->value();
    std::map<HandleType, int> selection = extract_k_ring(*polyhedron(), clicked, k_ring);

    bool any_change = false;
    if(is_insert) {
      for(typename std::map<HandleType, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= insert(it->first);
      }
    }else {
      for(typename std::map<HandleType, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= erase(it->first);
      }
    }
    if(any_change) { emit itemChanged(); }
  }

  bool eventFilter(QObject* /*target*/, QEvent *event)
  {
    // This filter is both filtering events from 'viewer' and 'main window'
    // key events
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }
    // mouse events
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      if(mouse_event->button() == Qt::LeftButton) {
        state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
      }   
    }

    // use mouse move event for paint-like selection
    if(event->type() == QEvent::MouseMove &&
      (state.shift_pressing && state.left_button_pressing) )    
    { // paint with mouse move event 
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(mouse_event->pos(), found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
      }
    }//end MouseMove
    return false;
  }

  template<class HandleType>
  std::map<HandleType, int> extract_k_ring(const Polyhedron &P, HandleType v, int k)
  {
    std::map<HandleType, int>  D;
    std::queue<HandleType>     Q;
    Q.push(v); D[v] = 0;

    int dist_v;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      v = Q.front();
      Q.pop();

      for(One_ring_iterator<HandleType> circ(v); circ; ++circ)
      {
        HandleType new_v = circ;
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      }
    }
    return D;
  }

  template<class HandleType, class Visitor, class InputIterator>
  void travel_not_selected_connected_components(
    Visitor& visitor, InputIterator begin, InputIterator end, std::size_t size) 
  {
    std::vector<bool> mark(size, false);
    for( ;begin != end; ++begin) 
    {
      if(mark[begin->id()] || is_selected(begin)) { continue; }

      std::vector<HandleType> C;
      C.push_back(begin);
      mark[begin->id()] = true;
      std::size_t current_index = 0;

      while(current_index < C.size()) {
        HandleType current = C[current_index++];

        for(One_ring_iterator<HandleType> circ(current); circ; ++circ)
        {
          HandleType nv = circ;
          if(!mark[nv->id()] && !is_selected(nv)) {
            mark[nv->id()] = true; 
            C.push_back(nv);
          }
        }
      } // while(!Q.empty())

      visitor(C);
    }
  }

  typedef Selection_set<Vertex_handle, Scene_polyhedron_selection_item> Selection_set_vertex;
  typedef Selection_set<Facet_handle, Scene_polyhedron_selection_item> Selection_set_facet;
  typedef Selection_set<Halfedge_handle, Scene_polyhedron_selection_item> Selection_set_edge;
  friend class Selection_set_vertex;
  friend class Selection_set_facet;
  friend class Selection_set_edge;
// members
  Ui::Selection* ui_widget;
  Mouse_keyboard_state state;
  //qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())
public:
  Selection_set_vertex selected_vertices;
  Selection_set_facet  selected_facets;
  Selection_set_edge   selected_edges;
};

#endif 
