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

// Wrapper for holding selected entities
template <class Entity>
class Selection_set : std::set<Entity>
{
private:
  typedef std::set<Entity>                Base;

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
    return Base::insert(entity).second;
  }
  bool erase(const Entity& entity) {
    return Base::erase(entity) != 0;
  }
  bool is_selected(const Entity& entity) const {
    return Base::find(entity) != end();
  }
  // for back_insert_iterator
  void push_back(const Entity& entity) { insert(entity); }
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

  Scene_polyhedron_selection_item(
    Scene_polyhedron_item* poly_item, 
    Active_handle::Type aht,
    bool is_insert,
    int k_ring,
    QMainWindow* mw) 
    : Scene_polyhedron_item_decorator(NULL, false)
  { init(poly_item, aht, is_insert, k_ring, mw); }

protected: 
  void init(Scene_polyhedron_item* poly_item, Active_handle::Type aht, bool is_insert, int k_ring, QMainWindow* mw)
  {
    this->poly_item = poly_item;
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Vertex_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Vertex_handle, int>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Facet_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Facet_handle, int>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::map<Polyhedron::Halfedge_handle, int>&)), this, 
      SLOT(selected(const std::map<Polyhedron::Halfedge_handle, int>&)));

    k_ring_selector.init(poly_item, mw, aht, k_ring);
    set_is_insert(is_insert);

    facet_color = QColor(87,87,87);
    edge_color = QColor(173,35,35);
    vertex_color = QColor(255,205,243);
  }

public:
  typedef Selection_set<Vertex_handle> Selection_set_vertex;
  typedef Selection_set<Facet_handle> Selection_set_facet;
  typedef Selection_set<Halfedge_handle> Selection_set_edge;

  Active_handle::Type get_active_handle_type() 
  { return k_ring_selector.active_handle_type; }
  void set_active_handle_type(Active_handle::Type aht) 
  { k_ring_selector.active_handle_type = aht; }

  int get_k_ring() { return k_ring_selector.k_ring; }
  void set_k_ring(int k) { k_ring_selector.k_ring = k; }

  bool get_is_insert() { return is_insert; }
  void set_is_insert(bool i) { is_insert = i; }

// drawing
  void draw() const {
    draw_selected_vertices();
    draw_selected_facets();
    draw_selected_edges();
  }
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
  bool actual_load(Scene_polyhedron_item* poly_item, Active_handle::Type aht, bool is_insert, int k_ring, QMainWindow* mw)
  {
    init(poly_item, aht, is_insert, k_ring, mw);

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
      h = &*h < &*h->opposite() ? h : h->opposite();
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

  void select_marked_edges(int neighb_size) {
    clear();

    for(Polyhedron::Edge_iterator
          eit = polyhedron()->edges_begin(),
          end = polyhedron()->edges_end(); eit != end; ++eit)
    {
      if(!eit->is_feature_edge()) continue;
      if(get_active_handle_type() == Active_handle::VERTEX) {
        selected_vertices.insert(eit->vertex());
        selected_vertices.insert(eit->opposite()->vertex());
      } else if(get_active_handle_type() == Active_handle::FACET) {
        selected_facets.insert(eit->face());
        selected_facets.insert(eit->opposite()->face());
      }
      else {
        selected_edges.insert(&*eit < &*eit->opposite() ? eit : eit->opposite());
      }
    }
    for(int i = 0; i < neighb_size; ++i) {
      if(get_active_handle_type() == Active_handle::VERTEX) {
        std::set<Vertex_handle> new_set;
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
      } else if(get_active_handle_type() == Active_handle::FACET){
        std::set<Facet_handle> new_set;
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
        std::set<Halfedge_handle> new_set;
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

  void erase_selected_facets() {
    if(selected_facets.empty()) {return;}
    // erase will-be-erased vertices and edges from selection
    for(Selection_set_vertex::iterator vb = selected_vertices.begin(); vb != selected_vertices.end(); ) {
      Polyhedron::Halfedge_around_vertex_circulator hvb((*vb)->vertex_begin()), hvbend(hvb);
      bool erase = true;
      do {
        if(!hvb->is_border() && !selected_facets.is_selected(hvb->facet()) ) {
          erase = false;
          break; 
        }
      } while(++hvb != hvbend);

      if(erase) {
        Vertex_handle v_erase = *vb;
        ++vb;
        selected_vertices.erase(v_erase);
      }
      else { ++vb; }
    }

    for(Selection_set_edge::iterator eb = selected_edges.begin(); eb != selected_edges.end(); ) {
      bool first_selected = (*eb)->is_border() || selected_facets.is_selected((*eb)->facet());
      bool second_selected = (*eb)->opposite()->is_border() || selected_facets.is_selected((*eb)->opposite()->facet());

      if(first_selected && second_selected) {
        Halfedge_handle h_erase = *eb;
        ++eb;
        selected_edges.erase(h_erase);
      }
      else { ++eb; }
    }

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
    std::map<Vertex_handle, std::size_t> index_map;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      Polyhedron::Halfedge_around_facet_circulator hb((*fb)->facet_begin()), hend(hb);
      do {
        index_map.insert(std::make_pair(hb->vertex(), index_map.size()));
      } while(++hb != hend);
    }
    // construct point vector
    std::vector<Polyhedron::Point_3> points(index_map.size());
    for(std::map<Vertex_handle, std::size_t>::iterator it = index_map.begin(); it != index_map.end(); ++it) {
      points[it->second] = it->first->point();
    }
    // construct polygon vector
    std::vector<std::vector<std::size_t> > polygons(selected_facets.size());
    std::size_t counter = 0;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb, ++counter) {
      Polyhedron::Halfedge_around_facet_circulator hb((*fb)->facet_begin()), hend(hb);
      do {
        polygons[counter].push_back(index_map[hb->vertex()]);
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

  virtual bool isEmpty() const { 
    if(poly_item == NULL) { return true; }
    return Scene_polyhedron_item_decorator::isEmpty();
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

protected:
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
