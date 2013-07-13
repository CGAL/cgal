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

// For iterating on one ring neighbor, sample usage:
// Halfedge_handle h; //or Vertex_handle or Facet_handle
// for(One_ring_iterator<Halfedge_handle> circ(h); circ; ++circ) {
//   Halfedge_handle h_neighbor = circ;
// }
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

  // for back_insert_iterator
  void push_back(const Entity& entity) { insert(entity); }
// members
  Listener* listener;
};

// To iterate on each minimum address halfedge as edge
struct Minimum_address_halfedge_iterator {
  typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
  Minimum_address_halfedge_iterator() {}
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
  typedef Polyhedron::Vertex_handle   Vertex_handle;
  typedef Polyhedron::Facet_handle    Facet_handle;
  typedef Polyhedron::Halfedge_handle Halfedge_handle;
  typedef Polyhedron::Vertex_iterator Vertex_iterator;
  typedef Polyhedron::Facet_iterator  Facet_iterator;

  enum ACTIVE_HANDLE_TYPE { VERTEX = 0, FACET = 1, EDGE = 2 };
  // To be used inside loader
  Scene_polyhedron_selection_item() : Scene_polyhedron_item_decorator(NULL, false) 
  { }

  Scene_polyhedron_selection_item(
    Scene_polyhedron_item* poly_item, 
    ACTIVE_HANDLE_TYPE aht,
    bool is_insert,
    bool k_ring) 
    : Scene_polyhedron_item_decorator(poly_item, false),
      active_handle_type(aht),
      is_insert(is_insert),
      k_ring(k_ring),
      selected_vertices(this),
      selected_facets(this),
      selected_edges(this)
  { init(); }

protected: 
  void init() 
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
  
  typedef Selection_set<Vertex_handle, Scene_polyhedron_selection_item> Selection_set_vertex;
  typedef Selection_set<Facet_handle, Scene_polyhedron_selection_item> Selection_set_facet;
  typedef Selection_set<Halfedge_handle, Scene_polyhedron_selection_item> Selection_set_edge;

  friend class Selection_set<Vertex_handle, Scene_polyhedron_selection_item>;
  friend class Selection_set<Facet_handle, Scene_polyhedron_selection_item>;
  friend class Selection_set<Halfedge_handle, Scene_polyhedron_selection_item>;

// get various informations to corresponding type
// should admit that this part is a little bit overkill (at least there is no code duplication)
  template<class HandleType>
  struct Selection_traits {};

  template<>
  struct Selection_traits<Vertex_handle> 
  {
    typedef Selection_set_vertex Container;
    typedef Vertex_iterator Iterator;
    Selection_traits(Scene_polyhedron_selection_item* item) : item(item) { }

    Container& container() { return item->selected_vertices; }
    Vertex_iterator iterator_begin() { return item->polyhedron()->vertices_begin(); }
    Vertex_iterator iterator_end() { return item->polyhedron()->vertices_end(); }
    std::size_t size() { return item->polyhedron()->size_of_vertices(); }

    Scene_polyhedron_selection_item* item;
  };

  template<>
  struct Selection_traits<Facet_handle> 
  {
    typedef Selection_set_facet Container;
    typedef Facet_iterator Iterator;
    Selection_traits(Scene_polyhedron_selection_item* item) : item(item) { }

    Container& container() { return item->selected_facets; }
    Facet_iterator iterator_begin() { return item->polyhedron()->facets_begin(); }
    Facet_iterator iterator_end() { return item->polyhedron()->facets_end(); }
    std::size_t size() { return item->polyhedron()->size_of_facets(); }

    Scene_polyhedron_selection_item* item;
  };

  template<>
  struct Selection_traits<Halfedge_handle> 
  {
    typedef Selection_set_edge Container;
    typedef Minimum_address_halfedge_iterator Iterator;
    Selection_traits(Scene_polyhedron_selection_item* item) : item(item) { }

    Container& container() { return item->selected_edges; }
    Minimum_address_halfedge_iterator iterator_begin() 
    { return Minimum_address_halfedge_iterator(item->polyhedron()->halfedges_begin(),
              item->polyhedron()->halfedges_end()); }
    Minimum_address_halfedge_iterator iterator_end() { return iterator_begin(); }
    std::size_t size() { return item->polyhedron()->size_of_halfedges(); }

    Scene_polyhedron_selection_item* item;
  };

public:
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

  bool save(const std::string& file_name) const {
    std::ofstream out(file_name);
    if(!out) { return false; }

    for(Selection_set_vertex::const_iterator it = selected_vertices.begin(); it != selected_vertices.end(); ++it) 
    { out << (*it)->id() << " "; }
    out << std::endl;

    for(Selection_set_facet::const_iterator it = selected_facets.begin(); it != selected_facets.end(); ++it) 
    { out << (*it)->id() << " "; }
    out << std::endl;

    for(Selection_set_edge::const_iterator it = selected_edges.begin(); it != selected_edges.end(); ++it) 
    { out << (*it)->id() << " "; }
    out << std::endl;
    return true;
  }
  bool load(const std::string& file_name) {
    file_name_holder = file_name;
    return true;
  }
  bool actual_load() {
    init();

    std::vector<Vertex_handle> all_vertices;
    all_vertices.reserve(polyhedron()->size_of_vertices());
    Polyhedron::Vertex_iterator vb(polyhedron()->vertices_begin()), ve(polyhedron()->vertices_end());
    for(;vb != ve; ++vb) { all_vertices.push_back(vb); }
    
    std::vector<Facet_handle> all_facets;
    all_facets.reserve(polyhedron()->size_of_facets());
    Polyhedron::Facet_iterator fb(polyhedron()->facets_begin()), fe(polyhedron()->facets_end());
    for(;fb != fe; ++fb) { all_facets.push_back(fb); }

    std::vector<Halfedge_handle> all_halfedges;
    all_facets.reserve(polyhedron()->size_of_halfedges());
    Polyhedron::Halfedge_iterator hb(polyhedron()->halfedges_begin()), he(polyhedron()->halfedges_end());
    for(;hb != he; ++hb) { all_halfedges.push_back(hb); }

    std::ifstream in(file_name_holder);
    if(!in) { return false; }

    std::string line;
    std::size_t id;

    if(!std::getline(in, line)) { return false; }
    std::istringstream vertex_line(line);
    while(vertex_line >> id) {
      if(id >= all_vertices.size()) { return false; }
      selected_vertices.insert(all_vertices[id]);
    }

    if(!std::getline(in, line)) { return false; }
    std::istringstream facet_line(line);
    while(facet_line >> id) {
      if(id >= all_facets.size()) { return false; }
      selected_facets.insert(all_facets[id]);
    }

    if(!std::getline(in, line)) { return false; }
    std::istringstream edge_line(line);
    while(edge_line >> id) {
      if(id >= all_halfedges.size()) { return false; }
      Halfedge_handle h = all_halfedges[id];
      h = &*h < &*h->opposite() ? h : h->opposite();
      selected_edges.insert(h);
    }
    return true;
  }

  void select_all() {
    switch(active_handle_type) {
    case VERTEX:
      select_all<Vertex_handle>(); break;
    case FACET:
      select_all<Facet_handle>(); break;
    case EDGE:
      select_all<Halfedge_handle>(); break;
    }
  }
  template<class HandleType>
  void select_all() {
    typedef Selection_traits<HandleType> Tr;
    Tr tr(this);
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      tr.container().insert(it);
    }
    emit itemChanged();
  }
  
  void clear() {
    switch(active_handle_type) {
    case VERTEX:
      clear<Vertex_handle>(); break;
    case FACET:
      clear<Facet_handle>(); break;
    case EDGE:
      clear<Halfedge_handle>(); break;
    }
  }
  template<class HandleType>
  void clear() {
    Selection_traits<HandleType> tr(this);
    tr.container().clear();
    emit itemChanged();
  }

  void select_marked_edges(int neighb_size) {
    clear();
    for(Polyhedron::Edge_iterator
          eit = polyhedron()->edges_begin(),
          end = polyhedron()->edges_end(); eit != end; ++eit)
    {
      if(!eit->is_feature_edge()) continue;
      if(active_handle_type == VERTEX) {
        selected_vertices.insert(eit->vertex());
        selected_vertices.insert(eit->opposite()->vertex());
      } else if(active_handle_type == FACET) {
        selected_facets.insert(eit->face());
        selected_facets.insert(eit->opposite()->face());
      }
      else {
        selected_edges.insert(&*eit < &*eit->opposite() ? eit : eit->opposite());
      }
    }
    for(int i = 0; i < neighb_size; ++i) {
      if(active_handle_type == VERTEX) {
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
      } else if(active_handle_type == FACET){
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
      else {
        std::set<Halfedge_handle> new_set;//(selected_facets);
        BOOST_FOREACH(Halfedge_handle h, selected_edges)
        {
          for(One_ring_iterator<Halfedge_handle> circ(h); circ; ++circ) {
            new_set.insert(circ);
          }
        }
        BOOST_FOREACH(Halfedge_handle h, new_set) {
          selected_edges.insert(h);
        }
      }
    }
    emit itemChanged();
  }

  boost::optional<std::size_t> get_minimum_isolated_component() {
    switch(active_handle_type) {
    case VERTEX:
      return get_minimum_isolated_component<Vertex_handle>();
    case FACET:
      return get_minimum_isolated_component<Facet_handle>();
    default:
      return get_minimum_isolated_component<Halfedge_handle>();
    }
  }
  template<class HandleType>
  boost::optional<std::size_t> get_minimum_isolated_component() {
    Selection_traits<HandleType> tr(this);
    Minimum_visitor visitor;
    travel_not_selected_connected_components<HandleType>
      (visitor, tr.iterator_begin(), tr.iterator_end(), tr.size());
    return visitor.minimum;
  }

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    switch(active_handle_type) {
    case VERTEX:
      return select_isolated_components<Vertex_handle>(threshold);
    case FACET:
      return select_isolated_components<Facet_handle>(threshold);
    default:
      return select_isolated_components<Halfedge_handle>(threshold);
    }
  }
  template<class HandleType>
  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    typedef Selection_traits<HandleType> Tr;
    Tr tr(this);
    typedef std::back_insert_iterator<typename Tr::Container> Output_iterator;
    Output_iterator out(tr.container());

    Selection_visitor<Output_iterator> visitor(threshold , out);
    travel_not_selected_connected_components<HandleType>
      (visitor, tr.iterator_begin(), tr.iterator_end(), tr.size());

    if(visitor.any_inserted) { emit itemChanged(); }
    return visitor.minimum_visitor.minimum;
  }

  void changed_with_poly_item() {
    poly_item->changed();
    emit itemChanged();
  }

  virtual bool isEmpty() const { 
    if(poly_item == NULL) { return true; }
    return Scene_polyhedron_item_decorator::isEmpty();
  }

public slots:
  void changed() {
    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
  }
  void vertex_has_been_selected(void* void_ptr) 
  {
    if(!visible() || active_handle_type != VERTEX) { return; }
    has_been_selected( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
  }
  void facet_has_been_selected(void* void_ptr)
  {
    if(!visible() || active_handle_type != FACET) { return; }
    has_been_selected( static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet() );
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    if(!visible() || active_handle_type != EDGE) { return; }
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
    Selection_traits<HandleType> tr(this);
    std::map<HandleType, int> selection = extract_k_ring(*polyhedron(), clicked, k_ring);

    bool any_change = false;
    if(is_insert) {
      for(typename std::map<HandleType, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= tr.container().insert(it->first);
      }
    }else {
      for(typename std::map<HandleType, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= tr.container().erase(it->first);
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
    Selection_traits<HandleType> tr(this);

    std::vector<bool> mark(size, false);
    for( ;begin != end; ++begin) 
    {
      if(mark[begin->id()] || tr.container().is_selected(begin)) { continue; }

      std::vector<HandleType> C;
      C.push_back(begin);
      mark[begin->id()] = true;
      std::size_t current_index = 0;

      while(current_index < C.size()) {
        HandleType current = C[current_index++];

        for(One_ring_iterator<HandleType> circ(current); circ; ++circ)
        {
          HandleType nv = circ;
          if(!mark[nv->id()] && !tr.container().is_selected(nv)) {
            mark[nv->id()] = true; 
            C.push_back(nv);
          }
        }
      } // while(!Q.empty())

      visitor(C);
    }
  }

// members
  Mouse_keyboard_state state;
  std::string file_name_holder;
public:
// action state
  ACTIVE_HANDLE_TYPE active_handle_type;
  bool is_insert;
  int  k_ring;
// selection
  Selection_set_vertex selected_vertices;
  Selection_set_facet  selected_facets;
  Selection_set_edge   selected_edges;
};

#endif 
