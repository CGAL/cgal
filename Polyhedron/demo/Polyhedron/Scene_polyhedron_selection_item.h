#ifndef SCENE_POLYHEDRON_SELECTION_ITEM_H
#define SCENE_POLYHEDRON_SELECTION_ITEM_H

#include "Scene_polyhedron_selection_item_config.h"
#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_type.h"
#include "ui_Selection_widget.h"
#include "opengl_tools.h"

#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QMouseEvent>

#include <queue>
#include <vector>
#include <fstream>

// Bunch of structs used in Scene_polyhedron_selection_item //

// For iterating on one ring neighbor
template<class T>
struct One_ring_iterator {};

template<>
struct One_ring_iterator<Polyhedron::Vertex_handle> {
  One_ring_iterator(Polyhedron::Vertex_handle v) : circ(v->vertex_begin()), end(circ) { }

  operator bool() { return circ != end; }
  operator Polyhedron::Vertex_handle() { return circ->opposite()->vertex(); }
  One_ring_iterator& operator++() { ++circ; return *this;}

  Polyhedron::Halfedge_around_vertex_circulator circ;
  Polyhedron::Halfedge_around_vertex_circulator end;
};

template<>
struct One_ring_iterator<Polyhedron::Facet_handle> {
  One_ring_iterator(Polyhedron::Facet_handle f) : circ(f->facet_begin()), end(circ) { }

  operator bool() { return circ != end; }
  operator Polyhedron::Facet_handle() { return circ->opposite()->facet(); }
  One_ring_iterator& operator++() { ++circ; return *this;}

  Polyhedron::Halfedge_around_facet_circulator circ;
  Polyhedron::Halfedge_around_facet_circulator end;
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
//////////////////////////////////////////////////////////////////////////

class SCENE_POLYHEDRON_SELECTION_ITEM_EXPORT Scene_polyhedron_selection_item 
  : public Scene_polyhedron_item_decorator
{
  Q_OBJECT

public:
  typedef Polyhedron::Vertex_handle Vertex_handle;
  typedef Polyhedron::Facet_handle Facet_handle;

  Scene_polyhedron_selection_item(Scene_polyhedron_item* poly_item, Ui::Selection* ui_widget) 
    : Scene_polyhedron_item_decorator(poly_item), ui_widget(ui_widget)
  { 
    connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_facet(void*)), this, SLOT(facet_has_been_selected(void*)));
    poly_item->enable_facets_picking(true);
    poly_item->set_color_vector_read_only(true);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    
    // id field is used in operations
    poly_item->update_vertex_indices();
    poly_item->update_facet_indices();
    // for rendering
    std::size_t id = 0;
    for (Polyhedron::Facet_iterator fit = polyhedron()->facets_begin(), 
      fit_end = polyhedron()->facets_end(); fit != fit_end; ++fit, ++id)
    { fit->set_patch_id(id); }
    poly_item->color_vector().resize(id, poly_item->color());
  }

  void draw_edges() const {
    poly_item->draw_edges();
    if(rendering_mode == Wireframe) {
      draw_selected_vertices();
    }
  }  
  void draw() const {
    poly_item->direct_draw();
    draw_selected_vertices();
  }

  void draw_selected_vertices() const {
    if(selected_vertices.empty() || !ui_widget->Show_selection_check_box->isChecked())
    { return; }

    GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);

    CGAL::GL::Color color;
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);
    color.set_rgb_color(0, 1.f, 0);

    ::glBegin(GL_POINTS);
    for(std::set<Vertex_handle>::iterator 
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

  void save_roi(const char* file_name) const { 
    std::ofstream out(file_name);
    // save roi
    for(std::set<Vertex_handle>::iterator it = selected_vertices.begin();
      it != selected_vertices.end(); ++it) {
        out << (*it)->id() << " ";
    }
    out.close();
  }
  void load_roi(const char* file_name) {
    // put vertices to vector
    std::vector<Polyhedron::Vertex_handle> all_vertices;
    all_vertices.reserve(polyhedron()->size_of_vertices());
    Polyhedron::Vertex_iterator vb(polyhedron()->vertices_begin()), ve(polyhedron()->vertices_end());
    for( ;vb != ve; ++vb) {
      all_vertices.push_back(vb);
    }
    // read roi
    std::ifstream in(file_name);
    std::size_t idx;
    while(in >> idx) {
      selected_vertices.insert(all_vertices[idx]);
    }
    in.close();
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

        One_ring_iterator<HandleType> circ(current);
        do {
          HandleType nv = circ;
          if(!mark[nv->id()] && !is_selected(nv)) {
            mark[nv->id()] = true; 
            C.push_back(nv);
          }
        } while(++circ);
      } // while(!Q.empty())

      visitor(C);
    }
  }

  void select_all() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      Polyhedron::Vertex_iterator vb(polyhedron()->vertices_begin()), ve(polyhedron()->vertices_end());
      for( ;vb != ve; ++vb) {
        selected_vertices.insert(vb);
      }
    }
    else {
      Polyhedron::Facet_iterator fb(polyhedron()->facets_begin()), fe(polyhedron()->facets_end());
      for( ;fb != fe; ++fb) {
        selected_facets.insert(fb);
        poly_item->color_vector()[fb->patch_id()] = QColor(0,255,0);
      }
    }
    emit itemChanged();
  }

  void clear() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      selected_vertices.clear();
    }
    else {
      selected_facets.clear();
      std::fill(poly_item->color_vector().begin(), poly_item->color_vector().end(), poly_item->color());
    }
    emit itemChanged();
  }

  void get_minimum_isolated_component() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      get_minimum_isolated_vertex_component();
    }
    else {
      get_minimum_isolated_facet_component();
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

  void select_isolated_components() {
    if(ui_widget->Selection_type_combo_box->currentIndex() == 0) {
      select_isolated_vertex_components();
    }
    else {
      select_isolated_facet_components();
    }
  }

  void select_isolated_vertex_components() {
    typedef std::insert_iterator<std::set<Vertex_handle> > Output_iterator;
    Output_iterator out(selected_vertices, selected_vertices.begin());

    Selection_visitor<Output_iterator> visitor(ui_widget->Threshold_size_spin_box->value() , out);
    travel_not_selected_connected_components<Vertex_handle>
      (visitor, polyhedron()->vertices_begin(), polyhedron()->vertices_end(), polyhedron()->size_of_vertices());

    if(visitor.any_inserted) { emit itemChanged(); }
    if(visitor.minimum_visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum_visitor.minimum);
    }
  }

  void select_isolated_facet_components() {
    typedef std::insert_iterator<std::set<Facet_handle> > Output_iterator;
    Output_iterator out(selected_facets, selected_facets.begin());

    Selection_visitor<Output_iterator> visitor(ui_widget->Threshold_size_spin_box->value() , out);
    travel_not_selected_connected_components<Facet_handle>
      (visitor, polyhedron()->facets_begin(), polyhedron()->facets_end(), polyhedron()->size_of_facets());

    if(visitor.any_inserted) { emit itemChanged(); }
    if(visitor.minimum_visitor.minimum) {
      ui_widget->Threshold_size_spin_box->setValue(*visitor.minimum_visitor.minimum);
    }
  }

  void changed_with_poly_item() {
    poly_item->changed();
    emit itemChanged();
  }

  bool is_selected(Vertex_handle v) { return selected_vertices.find(v) != selected_vertices.end(); }
  bool is_selected(Facet_handle f) { return selected_facets.find(f) != selected_facets.end(); }
public slots:
  void changed() {
    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
  }
  void vertex_has_been_selected(void* void_ptr) 
  {
    if(ui_widget->Selection_type_combo_box->currentIndex() != 0) { return; }
    Vertex_handle clicked_vertex = static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex();
    bool is_insert = ui_widget->Insertion_radio_button->isChecked();
    int k_ring = ui_widget->Brush_size_spin_box->value();
    std::map<Vertex_handle, int> selection = extract_k_ring(*polyhedron(), clicked_vertex, k_ring);

    bool any_change = false;
    if(is_insert) {
      for(std::map<Vertex_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= selected_vertices.insert(it->first).second;
      }
    }else {
      for(std::map<Vertex_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= selected_vertices.erase(it->first) != 0;
      }
    }
    if(any_change) { emit itemChanged(); }
  }
  void facet_has_been_selected(void* void_ptr)
  {
    if(ui_widget->Selection_type_combo_box->currentIndex() != 1) { return; }
    Facet_handle clicked_facet = static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet();
    bool is_insert = ui_widget->Insertion_radio_button->isChecked();
    int k_ring = ui_widget->Brush_size_spin_box->value();
    std::map<Facet_handle, int> selection = extract_k_ring(*polyhedron(), clicked_facet, k_ring);

    std::size_t size_before = selected_facets.size();
    if(is_insert) {
      for(std::map<Facet_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        if(selected_facets.insert(it->first).second) {
          std::cout << "inserted" << std::endl;
          poly_item->color_vector()[it->first->patch_id()] = QColor(0,255,0);
        }
      }
    }else {
      for(std::map<Facet_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        if(selected_facets.erase(it->first) != 0) {
          poly_item->color_vector()[it->first->patch_id()] = poly_item->color();
        }
      }
    }
    
    if(size_before != selected_facets.size() /* any change */) {
      emit itemChanged(); 
    }
  }

protected:
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

    if(!poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

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

      One_ring_iterator<HandleType> circ(v);
      do {
        HandleType new_v = circ;
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      } while(++circ);
    }
    return D;
  }

// members
  Ui::Selection* ui_widget;
  Mouse_keyboard_state state;
public:
  std::set<Vertex_handle> selected_vertices;
  std::set<Facet_handle>  selected_facets;
};

#endif 
