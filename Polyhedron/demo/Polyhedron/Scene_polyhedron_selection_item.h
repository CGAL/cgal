#ifndef SCENE_POLYHEDRON_SELECTION_ITEM_H
#define SCENE_POLYHEDRON_SELECTION_ITEM_H
#include "Scene_polyhedron_selection_item_config.h"
#include "Scene_polyhedron_item_k_ring_selection.h"
#include "Travel_isolated_components.h"

#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_type.h"
#include "opengl_tools.h"
#include <CGAL/gl_render.h>
#include <CGAL/orient_polygon_soup.h>

#include <fstream>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
// Wrapper for holding selected entities
template <class Entity, class Base>
class Selection_set : Base
{
public:
// types from base
  typedef typename Base::iterator         iterator;
  typedef typename Base::const_iterator   const_iterator;
  typedef typename Base::const_reference  const_reference;
  typedef typename Base::value_type       value_type;
// functions from base
  using Base::begin;
  using Base::end;
  using Base::size;
  using Base::clear;
  using Base::empty;

  bool insert(const Entity& entity) {
    Entity e = edge_filter(entity);
    return Base::insert(e).second;
  }
  bool erase(const Entity& entity) {
    Entity e = edge_filter(entity);
    return Base::erase(e) != 0;
  }
  bool is_selected(const Entity& entity) const {
    Entity e = edge_filter(entity);
    return Base::find(e) != end();
  }
  // for back_insert_iterator
  void push_back(const Entity& entity) {
    Entity e = edge_filter(entity);
    insert(e); 
  }

  Polyhedron::Halfedge_handle edge_filter(Polyhedron::Halfedge_handle h) const {
    return &*h < &*h->opposite() ? h : h->opposite();
  }
  template<class E> E edge_filter(E e) const { return e; }
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

template<class HandleType, class SelectionItem>
struct Selection_traits {};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::Vertex_handle, SelectionItem> 
{
  typedef typename SelectionItem::Selection_set_vertex Container;
  typedef typename SelectionItem::Vertex_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_vertices; }
  Iterator iterator_begin() { return item->polyhedron()->vertices_begin(); }
  Iterator iterator_end() { return item->polyhedron()->vertices_end(); }
  std::size_t size() { return item->polyhedron()->size_of_vertices(); }
  void update_indices() { item->polyhedron_item()->update_vertex_indices(); }

  SelectionItem* item;
};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::Facet_handle, SelectionItem> 
{
  typedef typename SelectionItem::Selection_set_facet Container;
  typedef typename SelectionItem::Facet_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_facets; }
  Iterator iterator_begin() { return item->polyhedron()->facets_begin(); }
  Iterator iterator_end() { return item->polyhedron()->facets_end(); }
  std::size_t size() { return item->polyhedron()->size_of_facets(); }
  void update_indices() { item->polyhedron_item()->update_facet_indices(); }

  SelectionItem* item;
};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::Halfedge_handle, SelectionItem> 
{
  typedef typename SelectionItem::Selection_set_edge Container;
  typedef Minimum_address_halfedge_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_edges; }
  Iterator iterator_begin() 
  { return Minimum_address_halfedge_iterator(item->polyhedron()->halfedges_begin(),
  item->polyhedron()->halfedges_end()); }
  Iterator iterator_end() { return iterator_begin(); }
  std::size_t size() { return item->polyhedron()->size_of_halfedges(); }
  void update_indices() { item->polyhedron_item()->update_halfedge_indices(); }
  SelectionItem* item;
};

//////////////////////////////////////////////////////////////////////////

class SCENE_POLYHEDRON_SELECTION_ITEM_EXPORT Scene_polyhedron_selection_item 
  : public Scene_polyhedron_item_decorator
{
  Q_OBJECT

friend class Polyhedron_demo_selection_plugin;

public:
  typedef Polyhedron::Vertex_handle   Vertex_handle;
  typedef Polyhedron::Facet_handle    Facet_handle;
  typedef Polyhedron::Halfedge_handle Halfedge_handle;
  typedef Polyhedron::Vertex_iterator Vertex_iterator;
  typedef Polyhedron::Facet_iterator  Facet_iterator;
  typedef Scene_polyhedron_item_k_ring_selection::Active_handle Active_handle;
  // To be used inside loader
  Scene_polyhedron_selection_item() 
    : Scene_polyhedron_item_decorator(NULL, false)
  { }

  Scene_polyhedron_selection_item(Scene_polyhedron_item* poly_item, QMainWindow* mw) 
    : Scene_polyhedron_item_decorator(NULL, false)
  { init(poly_item, mw); }

protected: 
  void init(Scene_polyhedron_item* poly_item, QMainWindow* mw)
  {
    this->poly_item = poly_item;
    connect(poly_item, SIGNAL(item_is_about_to_be_changed()), this, SLOT(poly_item_changed())); 
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Vertex_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Vertex_handle, int>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Facet_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Facet_handle, int>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Halfedge_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Halfedge_handle, int>&)));

    k_ring_selector.init(poly_item, mw, Active_handle::VERTEX, -1);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    mw->installEventFilter(this);

    facet_color = QColor(87,87,87);
    edge_color = QColor(173,35,35);
    vertex_color = QColor(255,205,243);
  }

  Active_handle::Type get_active_handle_type() 
  { return k_ring_selector.active_handle_type; }
  void set_active_handle_type(Active_handle::Type aht) 
  { k_ring_selector.active_handle_type = aht; }

  int get_k_ring() { return k_ring_selector.k_ring; }
  void set_k_ring(int k) { k_ring_selector.k_ring = k; }

  bool get_is_insert() { return is_insert; }
  void set_is_insert(bool i) { is_insert = i; }
  
public:
  //typedef Selection_set<Vertex_handle, std::set<Vertex_handle> >      Selection_set_vertex;
  //typedef Selection_set<Facet_handle, std::set<Facet_handle> >        Selection_set_facet;
  //typedef Selection_set<Halfedge_handle, std::set<Halfedge_handle> >  Selection_set_edge;
  typedef Selection_set<Vertex_handle, 
    boost::unordered_set<Vertex_handle, CGAL::Handle_hash_function> >    Selection_set_vertex;
  typedef Selection_set<Facet_handle,
    boost::unordered_set<Facet_handle, CGAL::Handle_hash_function> >     Selection_set_facet;
  typedef Selection_set<Halfedge_handle, 
    boost::unordered_set<Halfedge_handle, CGAL::Handle_hash_function> >  Selection_set_edge;

// drawing
  void draw() const {
    draw_selected_vertices();
    draw_selected_facets();
    draw_selected_edges();
  }
  void draw_edges() const { }

  void draw_selected_vertices() const {
    GLboolean enable_back_lighting = glIsEnabled(GL_LIGHTING);
    glDisable(GL_LIGHTING);

    
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);
    CGALglcolor(vertex_color);

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
    CGALglcolor(facet_color);

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

    CGALglcolor(edge_color);
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

  bool isEmpty() const {
    return selected_vertices.empty() && selected_edges.empty() && selected_facets.empty();
  }
  Bbox bbox() const 
  {
    boost::optional<CGAL::Bbox_3> item_bbox;

    for(Selection_set_vertex::const_iterator v_it = selected_vertices.begin(); 
        v_it != selected_vertices.end(); ++v_it) {
      *item_bbox = item_bbox ? *item_bbox + (*v_it)->point().bbox() : (*v_it)->point().bbox();
    }

    for(Selection_set_edge::const_iterator e_it = selected_edges.begin(); 
        e_it != selected_edges.end(); ++e_it) {
        CGAL::Bbox_3 e_bbox = (*e_it)->vertex()->point().bbox();
        e_bbox = e_bbox + (*e_it)->opposite()->vertex()->point().bbox();
        *item_bbox = item_bbox ? *item_bbox + e_bbox : e_bbox;
    }

    for(Selection_set_facet::const_iterator f_it = selected_facets.begin(); 
        f_it != selected_facets.end(); ++f_it) {

        Polyhedron::Halfedge_around_facet_circulator he = (*f_it)->facet_begin(), cend = he;
        CGAL_For_all(he,cend) {
          *item_bbox = item_bbox ? *item_bbox + he->vertex()->point().bbox() : he->vertex()->point().bbox();
        }
    }

    if(!item_bbox) { return Bbox(); }
    return Bbox(item_bbox->xmin(),item_bbox->ymin(),item_bbox->zmin(),
                item_bbox->xmax(),item_bbox->ymax(),item_bbox->zmax());
  }

  bool save(const std::string& file_name) const {
    // update id fields before using
    if(selected_vertices.size() > 0) { poly_item->update_vertex_indices();   }
    if(selected_facets.size() > 0)   { poly_item->update_facet_indices();    }
    if(selected_edges.size() > 0)    { poly_item->update_halfedge_indices(); }

    std::ofstream out(file_name.c_str());
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
  // this function is called by selection_plugin, since at the time of the call of load(...) 
  // we do not have access to selected polyhedron item
  bool actual_load(Scene_polyhedron_item* poly_item, QMainWindow* mw) 
  {
    init(poly_item, mw);

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

    std::ifstream in(file_name_holder.c_str());
    if(!in) { return false; }

    std::string line;
    std::size_t id;

    if(!std::getline(in, line)) { return true; }
    std::istringstream vertex_line(line);
    while(vertex_line >> id) {
      if(id >= all_vertices.size()) { return false; }
      selected_vertices.insert(all_vertices[id]);
    }

    if(!std::getline(in, line)) { return true; }
    std::istringstream facet_line(line);
    while(facet_line >> id) {
      if(id >= all_facets.size()) { return false; }
      selected_facets.insert(all_facets[id]);
    }

    if(!std::getline(in, line)) { return true; }
    std::istringstream edge_line(line);
    while(edge_line >> id) {
      if(id >= all_halfedges.size()) { return false; }
      Halfedge_handle h = all_halfedges[id];
      selected_edges.insert(h);
    }
    return true;
  }

  // select all of `active_handle_type`(vertex, facet or edge)
  void select_all() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      select_all<Vertex_handle>(); break;
    case Active_handle::FACET:
      select_all<Facet_handle>(); break;
    case Active_handle::EDGE:
      select_all<Halfedge_handle>(); break;
    }
  }
  // select all of vertex, facet or edge (use Vertex_handle, Facet_handle, Halfedge_handle as template argument)
  template<class HandleType>
  void select_all() {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      tr.container().insert(it);
    }
    emit itemChanged();
  }

  // clear all of `active_handle_type`(vertex, facet or edge)
  void clear() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      clear<Vertex_handle>(); break;
    case Active_handle::FACET:
      clear<Facet_handle>(); break;
    case Active_handle::EDGE:
      clear<Halfedge_handle>(); break;
    }
  }
  // select all of vertex, facet or edge (use Vertex_handle, Facet_handle, Halfedge_handle as template argument)
  template<class HandleType>
  void clear() {

    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
    tr.container().clear();
    emit itemChanged();
  }

  boost::optional<std::size_t> get_minimum_isolated_component() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      return get_minimum_isolated_component<Vertex_handle>();
    case Active_handle::FACET:
      return get_minimum_isolated_component<Facet_handle>();
    default:
      return get_minimum_isolated_component<Halfedge_handle>();
    }
  }
  template<class HandleType> // use Vertex_handle, Facet_handle, Halfedge_handle
  boost::optional<std::size_t> get_minimum_isolated_component() {
    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
    tr.update_indices();
    Travel_isolated_components::Minimum_visitor visitor;
    Travel_isolated_components().travel<HandleType>
      (tr.iterator_begin(), tr.iterator_end(), tr.size(), tr.container(), visitor);
    return visitor.minimum;
  }

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      return select_isolated_components<Vertex_handle>(threshold);
    case Active_handle::FACET:
      return select_isolated_components<Facet_handle>(threshold);
    default:
      return select_isolated_components<Halfedge_handle>(threshold);
    }
  }
  template<class HandleType> // use Vertex_handle, Facet_handle, Halfedge_handle
  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    tr.update_indices();
    typedef std::back_insert_iterator<typename Tr::Container> Output_iterator;
    Output_iterator out(tr.container());

    Travel_isolated_components::Selection_visitor<Output_iterator> visitor(threshold , out);
    Travel_isolated_components().travel<HandleType>
      (tr.iterator_begin(), tr.iterator_end(), tr.size(), tr.container(), visitor);

    if(visitor.any_inserted) { emit itemChanged(); }
    return visitor.minimum_visitor.minimum;
  }

  void expand_or_shrink(int steps) {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      expand_or_shrink<Vertex_handle>(steps); break;
    case Active_handle::FACET:
      expand_or_shrink<Facet_handle>(steps); break;
    default:
      expand_or_shrink<Halfedge_handle>(steps);
    }
  }

  template<class HandleType>
  void expand_or_shrink(int steps) {
   // It is good for large values of `steps`
    if(steps == 0) { return; }
    bool expand_req = steps > 0;
    steps = std::abs(steps);
    expand_req ? expand<HandleType>(steps) : shrink<HandleType>(steps);
  }

  template<class HandleType>
  void shrink(unsigned int steps) {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    tr.update_indices();
    std::vector<bool> mark(tr.size());
    std::vector<HandleType> to_be_shrink;
    std::vector<HandleType> next;
    to_be_shrink.reserve(tr.container().size());
    for(typename Tr::Container::iterator it = tr.container().begin(); it != tr.container().end(); ++it) {
      for(One_ring_iterator<HandleType> circ(*it); circ; ++circ) {
        if(!tr.container().is_selected(circ) && !mark[HandleType(circ)->id()]) {
          to_be_shrink.push_back(circ);
          mark[HandleType(circ)->id()] = true;
        }
      }
    }

    while(steps-- > 0) {
      for(typename std::vector<HandleType>::iterator it = to_be_shrink.begin(); 
        it != to_be_shrink.end(); ++it)
      {
        for(One_ring_iterator<HandleType> circ(*it); circ; ++circ) {
          HandleType ht = circ;
          if(tr.container().is_selected(ht) && !mark[ht->id()]) {
            next.push_back(ht);
            mark[ht->id()] = true;
          }
        }
      }

      to_be_shrink.swap(next);
      next.clear();
    }

    bool any_change = false;
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      if(mark[it->id()]) {
        any_change |= tr.container().erase(it);
      }
    }
    if(any_change) { emit itemChanged(); }
  }

  template<class HandleType>
  void expand(unsigned int steps) {

    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    tr.update_indices();
    std::vector<bool> mark(tr.size());

    std::vector<HandleType> to_be_expand;
    std::vector<HandleType> next;
    to_be_expand.reserve(tr.container().size());
    for(typename Tr::Container::iterator it = tr.container().begin(); it != tr.container().end(); ++it) {
      to_be_expand.push_back(*it);
    }

    while(steps-- > 0) {
      for(typename std::vector<HandleType>::iterator it = to_be_expand.begin(); 
        it != to_be_expand.end(); ++it)
      {
        for(One_ring_iterator<HandleType> circ(*it); circ; ++circ) {
          HandleType ht = circ;

          if(!tr.container().is_selected(ht) && !mark[ht->id()]) {
            next.push_back(ht);
            mark[ht->id()] = true;
          }
        }
      }

      to_be_expand.swap(next);
      next.clear();
    }

    bool any_change = false;
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      if(mark[it->id()]) {
        any_change |= tr.container().insert(it);
      }
    }
    if(any_change) { emit itemChanged(); }
  }

  void erase_selected_facets() {
    if(selected_facets.empty()) {return;}
    // no-longer-valid vertices and edges will be handled when item_about_to_be_changed() 

    // erase facets from poly
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      polyhedron()->erase_facet((*fb)->halfedge());
    }
    selected_facets.clear();
    changed_with_poly_item();
  }

  bool export_selected_facets_as_polyhedron(Polyhedron* out) {
    // Note: might be a more performance wise solution
    // assign sequential id to vertices neighbor to selected facets
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      Polyhedron::Halfedge_around_facet_circulator hb((*fb)->facet_begin()), hend(hb);
      do {
        hb->vertex()->id() = 0;
      } while(++hb != hend);
    }
    // construct point vector
    std::vector<Polyhedron::Point_3> points;
    points.reserve(selected_facets.size());
    std::size_t counter = 1;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      Polyhedron::Halfedge_around_facet_circulator hb((*fb)->facet_begin()), hend(hb);
      do {
        if(hb->vertex()->id() == 0) {
          hb->vertex()->id() = counter++; 
          points.push_back(hb->vertex()->point());
        }
      } while(++hb != hend);
    }
    // construct polygon vector
    std::vector<std::vector<std::size_t> > polygons(selected_facets.size());
    counter = 0;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb, ++counter) {
      Polyhedron::Halfedge_around_facet_circulator hb((*fb)->facet_begin()), hend(hb);
      do {
        polygons[counter].push_back(hb->vertex()->id() -1);
      } while(++hb != hend);
    }
    CGAL::Polygon_soup_to_polyhedron_3<Polyhedron::HalfedgeDS, Polyhedron::Point_3> builder(points, polygons);
    out->delegate(builder);
    return out->size_of_vertices() > 0;
  }

  void changed_with_poly_item() {
    // no need to update indices
    poly_item->changed();
    emit itemChanged();
  }

public slots:
  void changed() {
    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
  }
  // slots are called by signals of polyhedron_k_ring_selector
  void selected(const std::map<Polyhedron::Vertex_handle, int>& m)
  { has_been_selected(m); }
  void selected(const std::map<Polyhedron::Facet_handle, int>& m)
  { has_been_selected(m); }
  void selected(const std::map<Polyhedron::Halfedge_handle, int>& m)
  { has_been_selected(m); }
  void poly_item_changed() {
    remove_erased_handles<Vertex_handle>();
    remove_erased_handles<Halfedge_handle>();
    remove_erased_handles<Facet_handle>();
  }

protected:
  bool eventFilter(QObject* /*target*/, QEvent * gen_event)
  {
    if(!visible() || !k_ring_selector.state.shift_pressing) { return false; }
    if(gen_event->type() == QEvent::Wheel)
    {
      QWheelEvent *event = static_cast<QWheelEvent*>(gen_event);
      int steps = event->delta() / 120;
      expand_or_shrink(steps);
      return true;
    }
    return false;
  }

  template<class HandleType>
  void remove_erased_handles() {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    if(tr.container().empty()) { return;}

    std::vector<HandleType> exists;
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      if(tr.container().is_selected(it)) {
        exists.push_back(it);
      }
    }
    tr.container().clear();
    for(typename std::vector<HandleType>::iterator it = exists.begin(); it != exists.end(); ++it) {
      tr.container().insert(*it);
    }
  }

  template<class HandleType>
  void has_been_selected(const std::map<HandleType, int>& selection) 
  {
    if(!visible()) { return; }
    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);

    bool any_change = false;
    if(is_insert) {
      for(typename std::map<HandleType, int>::const_iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= tr.container().insert(it->first);
      }
    }else {
      for(typename std::map<HandleType, int>::const_iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= tr.container().erase(it->first);
      }
    }
    if(any_change) { emit itemChanged(); }
  }

// members
  std::string file_name_holder;
  Scene_polyhedron_item_k_ring_selection k_ring_selector;
  // action state
  bool is_insert;

public:
// selection
  Selection_set_vertex selected_vertices;
  Selection_set_facet  selected_facets;
  Selection_set_edge   selected_edges; // stores one halfedge for each pair (halfedge with minimum address)
// 
  QColor vertex_color, facet_color, edge_color;
};

#endif 
