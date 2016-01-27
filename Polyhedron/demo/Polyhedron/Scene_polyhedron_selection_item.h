#ifndef SCENE_POLYHEDRON_SELECTION_ITEM_H
#define SCENE_POLYHEDRON_SELECTION_ITEM_H
#include "opengl_tools.h"
#include "Scene_polyhedron_selection_item_config.h"
#include "Scene_polyhedron_item_k_ring_selection.h"
#include "Travel_isolated_components.h"

#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_type.h"

#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include "Polyhedron_demo_detect_sharp_edges.h"

#include <fstream>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace PMP = CGAL::Polygon_mesh_processing;


template<class HandleType, class SelectionItem>
struct Selection_traits {};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::Vertex_handle, SelectionItem> 
{
  typedef typename SelectionItem::Selection_set_vertex Container;
  typedef boost::graph_traits<Polyhedron>::vertex_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_vertices; }
  Iterator iterator_begin() { return vertices(*item->polyhedron()).first; }
  Iterator iterator_end() { return vertices(*item->polyhedron()).second; }
  std::size_t size() { return item->polyhedron()->size_of_vertices(); }
  void update_indices() { item->polyhedron_item()->update_vertex_indices(); }
  std::size_t id(typename SelectionItem::Vertex_handle  vh) {return vh->id();}

  template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
  static
  OutputIterator
  reduce_selection(
    const VertexRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsVertexSelectedPMap is_selected,
    OutputIterator out)
  {
    return reduce_vertex_selection(selection, graph, k, is_selected, out);
  }
  template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
  static
  OutputIterator
  expand_selection(
    const VertexRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsVertexSelectedPMap is_selected,
    OutputIterator out)
  {
    return expand_vertex_selection(selection, graph, k, is_selected, out);
  }

  SelectionItem* item;
};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::Facet_handle, SelectionItem>
{
  typedef typename SelectionItem::Selection_set_facet Container;
  typedef boost::graph_traits<Polyhedron>::face_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_facets; }
  Iterator iterator_begin() { return faces(*item->polyhedron()).first; }
  Iterator iterator_end() { return faces(*item->polyhedron()).second; }
  std::size_t size() { return item->polyhedron()->size_of_facets(); }
  void update_indices() { item->polyhedron_item()->update_facet_indices(); }
  std::size_t id(typename SelectionItem::Facet_handle fh) {return fh->id();}

  template <class FaceRange, class HalfedgeGraph, class IsFaceSelectedPMap, class OutputIterator>
  static
  OutputIterator
  reduce_selection(
    const FaceRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsFaceSelectedPMap is_selected,
    OutputIterator out)
  {
    return reduce_face_selection(selection, graph, k, is_selected, out);
  }
  template <class FaceRange, class HalfedgeGraph, class IsFaceSelectedPMap, class OutputIterator>
  static
  OutputIterator
  expand_selection(
    const FaceRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsFaceSelectedPMap is_selected,
    OutputIterator out)
  {
    return expand_face_selection(selection, graph, k, is_selected, out);
  }

  SelectionItem* item;
};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::edge_descriptor, SelectionItem> 
{
  typedef typename SelectionItem::Selection_set_edge Container;
  typedef boost::graph_traits<Polyhedron>::edge_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_edges; }
  Iterator iterator_begin() { return edges(*item->polyhedron()).first; }
  Iterator iterator_end() { return edges(*item->polyhedron()).second; }
  std::size_t size() { return item->polyhedron()->size_of_halfedges()/2; }
  void update_indices() { item->polyhedron_item()->update_halfedge_indices(); }
  std::size_t id(boost::graph_traits<Polyhedron>::edge_descriptor ed) {return ed.halfedge()->id()/2;}

  template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
  static
  OutputIterator
  reduce_selection(
    const EdgeRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsEdgeSelectedPMap is_selected,
    OutputIterator out)
  {
    return reduce_edge_selection(selection, graph, k, is_selected, out);
  }
  template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
  static
  OutputIterator
  expand_selection(
    const EdgeRange& selection,
    HalfedgeGraph& graph,
    unsigned int k,
    IsEdgeSelectedPMap is_selected,
    OutputIterator out)
  {
    return expand_edge_selection(selection, graph, k, is_selected, out);
  }

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
  typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
  typedef Polyhedron::Vertex_iterator Vertex_iterator;
  typedef Polyhedron::Facet_iterator  Facet_iterator;
  typedef Scene_polyhedron_item_k_ring_selection::Active_handle Active_handle;
  // To be used inside loader
  Scene_polyhedron_selection_item() 
    : Scene_polyhedron_item_decorator(NULL, false)
    {
        for(int i=0; i<3; i++)
        {
            addVaos(i);
            vaos[i]->create();
        }

        for(int i=0; i<7; i++)
        {
            buffers[i].create();
        }
        nb_facets = 0;
        nb_points = 0;
        nb_lines = 0;
    }

  Scene_polyhedron_selection_item(Scene_polyhedron_item* poly_item, QMainWindow* mw) 
    : Scene_polyhedron_item_decorator(NULL, false)
    {
        nb_facets = 0;
        nb_points = 0;
        nb_lines = 0;

        for(int i=0; i<3; i++)
        {
            addVaos(i);
            vaos[i]->create();
        }

        for(int i=0; i<7; i++)
        {
            buffers[i].create();
        }
        init(poly_item, mw);
        invalidateOpenGLBuffers();
    }

   ~Scene_polyhedron_selection_item()
    {
    }

protected: 
  void init(Scene_polyhedron_item* poly_item, QMainWindow* mw)
  {
    this->poly_item = poly_item;
    connect(poly_item, SIGNAL(item_is_about_to_be_changed()), this, SLOT(poly_item_changed())); 
    connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Vertex_handle>&)), this,
      SLOT(selected(const std::set<Polyhedron::Vertex_handle>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::set<Polyhedron::Facet_handle>&)), this,
      SLOT(selected(const std::set<Polyhedron::Facet_handle>&)));
    connect(&k_ring_selector, SIGNAL(selected(const std::set<edge_descriptor>&)), this,
      SLOT(selected(const std::set<edge_descriptor>&)));
    connect(&k_ring_selector, SIGNAL(endSelection()), this,SLOT(endSelection()));
    connect(&k_ring_selector, SIGNAL(toogle_insert(bool)), this,SLOT(toggle_insert(bool)));
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
  typedef boost::unordered_set<Vertex_handle, CGAL::Handle_hash_function>    Selection_set_vertex;
  typedef boost::unordered_set<Facet_handle, CGAL::Handle_hash_function>      Selection_set_facet;
  typedef boost::unordered_set<edge_descriptor, CGAL::Handle_hash_function>    Selection_set_edge;

  Polyhedron* polyhedron()
  {
    return this->poly_item->polyhedron();
  }

    using Scene_polyhedron_item_decorator::draw;
    virtual void draw(CGAL::Three::Viewer_interface*) const;
    virtual void draw_edges() const { }
    virtual void draw_edges(CGAL::Three::Viewer_interface*) const;
    virtual void draw_points(CGAL::Three::Viewer_interface*) const;

  bool supportsRenderingMode(RenderingMode m) const { return (m==Flat); }

  bool isEmpty() const {
    return selected_vertices.empty() && selected_edges.empty() && selected_facets.empty();
  }
  void selection_changed(bool b)
  {
      QGLViewer* v = *QGLViewer::QGLViewerPool().begin();
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      if(!viewer)
          return;
      if(!b)
        viewer->setBindingSelect();
      else
        viewer->setNoBinding();
  }
  void compute_bbox() const
  {
    // Workaround a bug in g++-4.8.3:
    //   http://stackoverflow.com/a/21755207/1728537
    // Using boost::make_optional to copy-initialize 'item_bbox' hides the
    //   warning about '*item_bbox' not being initialized.
    // -- Laurent Rineau, 2014/10/30
    boost::optional<CGAL::Bbox_3> item_bbox
      = boost::make_optional(false, CGAL::Bbox_3());

    for(Selection_set_vertex::const_iterator v_it = selected_vertices.begin(); 
        v_it != selected_vertices.end(); ++v_it) {

      if(item_bbox) { *item_bbox = *item_bbox + (*v_it)->point().bbox(); }
      else          {  item_bbox = (*v_it)->point().bbox(); }
    }

    for(Selection_set_edge::const_iterator e_it = selected_edges.begin(); 
        e_it != selected_edges.end(); ++e_it) {
        CGAL::Bbox_3 e_bbox = e_it->halfedge()->vertex()->point().bbox();
        e_bbox = e_bbox + e_it->halfedge()->opposite()->vertex()->point().bbox();
        if(item_bbox) { *item_bbox = *item_bbox + e_bbox; }
        else          {  item_bbox = e_bbox; }
    }

    for(Selection_set_facet::const_iterator f_it = selected_facets.begin(); 
        f_it != selected_facets.end(); ++f_it) {

        Polyhedron::Halfedge_around_facet_circulator he = (*f_it)->facet_begin(), cend = he;
        CGAL_For_all(he,cend) {
          if(item_bbox) { *item_bbox = *item_bbox + he->vertex()->point().bbox(); }
          else          {  item_bbox = he->vertex()->point().bbox(); }
        }
    }

    if(!item_bbox) { _bbox = this->poly_item->bbox(); return;}
    _bbox = Bbox(item_bbox->xmin(),item_bbox->ymin(),item_bbox->zmin(),
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
    { out << it->id() << " "; }
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

    std::vector<edge_descriptor> all_edges(edges(*polyhedron()).first, edges(*polyhedron()).second);

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
      if(id >= all_edges.size()) { return false; }
      selected_edges.insert(all_edges[id]);
    }
    return true;
  }

  // select all of `active_handle_type`(vertex, facet or edge)
  void select_all() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      select_all<Vertex_handle>(); break;
    case Active_handle::FACET:
    case Active_handle::CONNECTED_COMPONENT:
      select_all<Facet_handle>(); break;
    case Active_handle::EDGE:
      selected_edges.insert(edges(*polyhedron()).first, edges(*polyhedron()).second);
      invalidateOpenGLBuffers();
      QGLViewer* v = *QGLViewer::QGLViewerPool().begin();
      v->update();

    }
  }
  // select all of vertex, facet or edge (use Vertex_handle, Facet_handle, edge_descriptor as template argument)
  template<class HandleType>
  void select_all() {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      tr.container().insert(*it);
    }
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  // clear all of `active_handle_type`(vertex, facet or edge)
  void clear() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      clear<Vertex_handle>(); break;
    case Active_handle::FACET:
    case Active_handle::CONNECTED_COMPONENT:
      clear<Facet_handle>(); break;
    case Active_handle::EDGE:
      clear<edge_descriptor>(); break;
    }
  }
  // select all of vertex, facet or edge (use Vertex_handle, Facet_handle, edge_descriptor as template argument)
  template<class HandleType>
  void clear() {

    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
    tr.container().clear();
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void clear_all(){
    clear<Vertex_handle>();
    clear<Facet_handle>();
    clear<edge_descriptor>();
  }

  boost::optional<std::size_t> get_minimum_isolated_component() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      return get_minimum_isolated_component<Vertex_handle>();
    case Active_handle::FACET:
      return get_minimum_isolated_component<Facet_handle>();
    default:
      return get_minimum_isolated_component<edge_descriptor>();
    }
  }
  template<class HandleType> // use Vertex_handle, Facet_handle, edge_descriptor
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
      return select_isolated_components<edge_descriptor>(threshold);
    }
  }
  template<class HandleType> // use Vertex_handle, Facet_handle, edge_descriptor
  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    tr.update_indices();
    typedef std::insert_iterator<typename Tr::Container> Output_iterator;
    Output_iterator out(tr.container(), tr.container().begin());

    Travel_isolated_components::Selection_visitor<Output_iterator> visitor(threshold , out);
    Travel_isolated_components().travel<HandleType>
      (tr.iterator_begin(), tr.iterator_end(), tr.size(), tr.container(), visitor);

    if(visitor.any_inserted) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
    return visitor.minimum_visitor.minimum;
  }

  void expand_or_reduce(int steps) {
    if (steps>0)
    {
      switch(get_active_handle_type()) {
        case Active_handle::VERTEX:
          expand_selection<Vertex_handle>(steps);
        break;
        case Active_handle::FACET:
        case Active_handle::CONNECTED_COMPONENT:
          expand_selection<Facet_handle>(steps);
        break;
        case Active_handle::EDGE:
          expand_selection<edge_descriptor>(steps);
      }
    }
    else
    {
      switch(get_active_handle_type()) {
        case Active_handle::VERTEX:
          reduce_selection<Vertex_handle>(-steps);
        break;
        case Active_handle::FACET:
        case Active_handle::CONNECTED_COMPONENT:
          reduce_selection<Facet_handle>(-steps);
        break;
        case Active_handle::EDGE:
          reduce_selection<edge_descriptor>(-steps);
      }
    }
  }

  template <class Handle>
  struct Is_selected_property_map{
    std::vector<bool>* is_selected_ptr;
    Is_selected_property_map()
      : is_selected_ptr(NULL) {}
    Is_selected_property_map(std::vector<bool>& is_selected)
      : is_selected_ptr( &is_selected) {}

    template<class H>
    std::size_t id(H h){ return h->id(); }
    std::size_t id(edge_descriptor ed) { return ed.halfedge()->id()/2; }

    friend bool get(Is_selected_property_map map, Handle h)
    {
      CGAL_assertion(map.is_selected_ptr!=NULL);
      return (*map.is_selected_ptr)[map.id(h)];
    }

    friend void put(Is_selected_property_map map, Handle h, bool b)
    {
      CGAL_assertion(map.is_selected_ptr!=NULL);
      (*map.is_selected_ptr)[map.id(h)]=b;
    }
  };

  template <class Handle>
  struct Index_map
  {
    typedef Handle  key_type;
    typedef std::size_t     value_type;
    typedef value_type&    reference;
    typedef boost::read_write_property_map_tag category;

    friend value_type get(Index_map, Handle h)
    {
      return h->id();
    }
    friend void put(Index_map, Handle h, value_type i)
    {
      h->id() = i;
    }
  };

  template <class Handle>
  void expand_selection(unsigned int steps)
  {
    typedef Selection_traits<Handle, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    tr.update_indices();
    std::vector<bool> mark(tr.size(),false);

    BOOST_FOREACH(Handle h,tr.container())
      mark[tr.id(h)]=true;

    Tr::expand_selection(
      tr.container(),
      *this->poly_item->polyhedron(),
      steps,
      Is_selected_property_map<Handle>(mark),
      CGAL::Emptyset_iterator()
    );

    bool any_change = false;
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      if(mark[tr.id(*it)]) {
        any_change |= tr.container().insert(*it).second;
      }
    }
    if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
  }

  template <class Handle>
  void reduce_selection(unsigned int steps) {

    typedef Selection_traits<Handle, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    tr.update_indices();
    std::vector<bool> mark(tr.size(),false);

    BOOST_FOREACH(Handle h,tr.container())
      mark[tr.id(h)]=true;

    Tr::reduce_selection(
      tr.container(),
      *this->poly_item->polyhedron(),
      steps,
      Is_selected_property_map<Handle>(mark),
      CGAL::Emptyset_iterator()
    );

    bool any_change = false;
    for(typename Tr::Iterator it = tr.iterator_begin() ; it != tr.iterator_end(); ++it) {
      if(!mark[tr.id(*it)]) {
        any_change |= (tr.container().erase(*it)!=0);
      }
    }
    if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
  }

  void erase_selected_facets() {
    if(selected_facets.empty()) {return;}
    // no-longer-valid vertices and edges will be handled when item_about_to_be_changed() 

    // erase facets from poly
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      polyhedron()->erase_facet((*fb)->halfedge());
    }
    selected_facets.clear();
    invalidateOpenGLBuffers();
    changed_with_poly_item();
  }

  void keep_connected_components() {
    if (selected_facets.empty()) { return; }

    Selection_traits<Polyhedron::Face_handle, Scene_polyhedron_selection_item> trf(this);
    trf.update_indices();
    Selection_traits<Polyhedron::Vertex_handle, Scene_polyhedron_selection_item> trv(this);
    trv.update_indices();

    //Selection_traits<edge_descriptor, Scene_polyhedron_selection_item> tre(this);
    //tre.update_indices();
    //std::vector<bool> mark(polyhedron()->size_of_halfedges() / 2, false);
    //BOOST_FOREACH(edge_descriptor e, selected_edges)
    //  mark[tre.id(e)] = true;

    PMP::keep_connected_components(*polyhedron()
      , trf.container()
      , PMP::parameters::face_index_map(Index_map<Polyhedron::Face_handle>())
      .vertex_index_map(Index_map<Polyhedron::Vertex_handle>()));
//      .edge_is_constrained_map(Is_selected_property_map<edge_descriptor>(mark)));

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
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<Polyhedron>(
      points, polygons, *out);

    return out->size_of_vertices() > 0;
  }

  void select_sharp_edges(const double angle)
  {
    CGAL::detect_sharp_edges(polyhedron(), angle);

    BOOST_FOREACH(edge_descriptor e, edges(*polyhedron()))
    {
      Polyhedron::Halfedge_handle h = halfedge(e, *polyhedron());
      if (h->is_feature_edge())
        selected_edges.insert(e);
    }
    invalidateOpenGLBuffers();
  }

  void changed_with_poly_item() {
    // no need to update indices
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
    Q_EMIT itemChanged();
  }

Q_SIGNALS:
  void simplicesSelected(CGAL::Three::Scene_item*);

public Q_SLOTS:
  void invalidateOpenGLBuffers() {

    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
      //  poly_item->invalidateOpenGLBuffers();
        are_buffers_filled = false;
        compute_bbox();
  }
  // slots are called by signals of polyhedron_k_ring_selector
  void selected(const std::set<Polyhedron::Vertex_handle>& m)
  { has_been_selected(m); }
  void selected(const std::set<Polyhedron::Facet_handle>& m)
  { has_been_selected(m); }
  void selected(const std::set<edge_descriptor>& m)
  { has_been_selected(m); }
  void poly_item_changed() {
    remove_erased_handles<Vertex_handle>();
    remove_erased_handles<edge_descriptor>();
    remove_erased_handles<Facet_handle>();
  }
  void endSelection(){
    Q_EMIT simplicesSelected(this);
  }
  void toggle_insert(bool b)
  {
   is_insert = b;
  }

protected:
  bool eventFilter(QObject* /*target*/, QEvent * gen_event)
  {
    if(!visible() || !k_ring_selector.state.shift_pressing) { return false; }
    if(gen_event->type() == QEvent::Wheel)
    {
      QWheelEvent *event = static_cast<QWheelEvent*>(gen_event);
      int steps = event->delta() / 120;
      expand_or_reduce(steps);
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
      if(tr.container().count(*it)) {
        exists.push_back(*it);
      }
    }
    tr.container().clear();
    for(typename std::vector<HandleType>::iterator it = exists.begin(); it != exists.end(); ++it) {
      tr.container().insert(*it);
    }
  }

  template<typename HandleRange>
  bool treat_selection(const HandleRange& selection)
  {
    typedef typename HandleRange::value_type HandleType;
    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);

    bool any_change = false;
    if(is_insert) {
      BOOST_FOREACH(HandleType h, selection)
        any_change |= tr.container().insert(h).second;
    }
    else{
      BOOST_FOREACH(HandleType h, selection)
        any_change |= (tr.container().erase(h)!=0);
    }
    if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
    return any_change;
  }

  Facet_handle face(Facet_handle fh)
  { return fh; }
  Facet_handle face(Vertex_handle)
  { return boost::graph_traits<Polyhedron>::null_face(); }
  Facet_handle face(edge_descriptor)
  { return boost::graph_traits<Polyhedron>::null_face(); }

  template<class HandleType>
  void has_been_selected(const std::set<HandleType>& selection)
  {
    if(!visible()) { return; }


    if (get_active_handle_type() == Active_handle::CONNECTED_COMPONENT)
    {
      Selection_traits<edge_descriptor,
                       Scene_polyhedron_selection_item> tr(this);
      tr.update_indices();
      std::vector<bool> mark(tr.size(), false);
      BOOST_FOREACH(edge_descriptor e, selected_edges)
        mark[tr.id(e)] = true;
      std::vector<Facet_handle> selected_cc;
      CGAL::Polygon_mesh_processing::connected_component(
        face(*selection.begin()),
        *polyhedron(),
        std::back_inserter(selected_cc),
        CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(
          Is_selected_property_map<edge_descriptor>(mark)));
       treat_selection(selected_cc);
    }
    else
    {
      treat_selection(selection);
    }
  }

public:
  Is_selected_property_map<edge_descriptor>
    selected_edges_pmap(std::vector<bool>& mark)
  {
    Selection_traits<edge_descriptor,
      Scene_polyhedron_selection_item> tr(this);
    tr.update_indices();

    for (unsigned int i = 0; i < mark.size(); ++i)
      mark[i] = false;

    BOOST_FOREACH(edge_descriptor e, selected_edges)
      mark[tr.id(e)] = true;

    return Is_selected_property_map<edge_descriptor>(mark);
  }

protected:
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

private:

  mutable std::vector<float> positions_facets;
  mutable std::vector<float> normals;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_points;
  mutable std::size_t nb_facets;
  mutable std::size_t nb_points;
  mutable std::size_t nb_lines;
  mutable QOpenGLShaderProgram *program;
  using CGAL::Three::Scene_item::initialize_buffers;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_elements() const;

};

#endif 
