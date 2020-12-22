#ifndef SCENE_POLYHEDRON_SELECTION_ITEM_H
#define SCENE_POLYHEDRON_SELECTION_ITEM_H

#include "Scene_polyhedron_selection_item_config.h"
#include "Plugins/PMP/Scene_facegraph_item_k_ring_selection.h"
#include "Travel_isolated_components.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_item_decorator.h"
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Three/Scene_print_item_interface.h>

#include "Polyhedron_demo_detect_sharp_edges.h"

// Laurent Rineau, 2016/04/07: that header should not be included here, but
// only in the .cpp file. But that header file does contain the body of a
// few member functions.
#include <CGAL/Three/Viewer_interface.h>

#include <fstream>
#include <boost/unordered_set.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/Qt/manipulatedFrame.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_face_graph_item::Face_graph Face_graph;
typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::type VPmap;
typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::const_type constVPmap;
typedef Scene_face_graph_item::Vertex_selection_map Vertex_selection_map;
typedef Scene_face_graph_item::Face_selection_map Face_selection_map;

typedef boost::graph_traits<Face_graph>::vertex_descriptor fg_vertex_descriptor;
typedef boost::graph_traits<Face_graph>::halfedge_descriptor fg_halfedge_descriptor;
typedef boost::graph_traits<Face_graph>::edge_descriptor fg_edge_descriptor;
typedef boost::graph_traits<Face_graph>::face_descriptor fg_face_descriptor;

typedef boost::graph_traits<Face_graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Face_graph>::halfedge_iterator halfedge_iterator;
typedef boost::graph_traits<Face_graph>::edge_iterator edge_iterator;
typedef boost::graph_traits<Face_graph>::face_iterator face_iterator;

template<class HandleType, class SelectionItem>
struct Selection_traits {};

template<class SelectionItem>
struct Selection_traits<typename SelectionItem::fg_vertex_descriptor, SelectionItem>
{
  typedef typename SelectionItem::Selection_set_vertex Container;
  typedef vertex_iterator Iterator;

  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_vertices; }
  Iterator iterator_begin() { return vertices(*item->polyhedron()).first; }
  Iterator iterator_end() { return vertices(*item->polyhedron()).second; }
  std::size_t size() { return num_vertices(*item->polyhedron()); }
  std::size_t id(typename SelectionItem::fg_vertex_descriptor  vh)
  {
    return get(get(boost::vertex_index, *item->polyhedron()), vh);
  }

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
struct Selection_traits<typename SelectionItem::fg_face_descriptor, SelectionItem>
{
  typedef typename SelectionItem::Selection_set_facet Container;
  typedef face_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_facets; }
  Iterator iterator_begin() { return faces(*item->polyhedron()).first; }
  Iterator iterator_end() { return faces(*item->polyhedron()).second; }
  std::size_t size() { return num_faces(*item->polyhedron()); }
  std::size_t id(typename SelectionItem::fg_face_descriptor fh)
{
  return get(get(boost::face_index, *item->polyhedron()), fh);
}

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
struct Selection_traits<typename SelectionItem::fg_edge_descriptor, SelectionItem>
{
  typedef typename SelectionItem::Selection_set_edge Container;
  typedef edge_iterator Iterator;
  Selection_traits(SelectionItem* item) : item(item) { }

  Container& container() { return item->selected_edges; }
  Iterator iterator_begin() { return edges(*item->polyhedron()).first; }
  Iterator iterator_end() { return edges(*item->polyhedron()).second; }
  std::size_t size() { return num_edges(*item->polyhedron()); }
  std::size_t id(boost::graph_traits<Face_graph>::edge_descriptor ed)
  {
    return get(boost::halfedge_index, *item->polyhedron(), halfedge(ed,*item->polyhedron()))/2;
  }

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
struct Scene_polyhedron_selection_item_priv;
class SCENE_POLYHEDRON_SELECTION_ITEM_EXPORT Scene_polyhedron_selection_item
  : public Scene_polyhedron_item_decorator,
    public CGAL::Three::Scene_print_item_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Scene_print_item_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PrintInterface/1.0")

friend class Polyhedron_demo_selection_plugin;

public:
  typedef boost::graph_traits<Face_graph>::face_descriptor fg_face_descriptor;
  typedef boost::graph_traits<Face_graph>::edge_descriptor fg_edge_descriptor;
  typedef boost::graph_traits<Face_graph>::halfedge_descriptor fg_halfedge_descriptor;
  typedef boost::graph_traits<Face_graph>::vertex_descriptor fg_vertex_descriptor;

  typedef Scene_facegraph_item_k_ring_selection::Active_handle Active_handle;

  enum SelectionType{
    Vertex=0x1,
    Edge=0x2,
    Facet=0x4,
    None=0x8
  };
  Q_DECLARE_FLAGS(SelectionTypes, SelectionType)
  void common_constructor();
  Scene_polyhedron_selection_item() ;
  Scene_polyhedron_selection_item(Scene_face_graph_item* poly_item, QMainWindow* mw);
  ~Scene_polyhedron_selection_item();
  void inverse_selection();
  void setPathSelection(bool b);
  //For ID printing
  void printPrimitiveId(QPoint, CGAL::Three::Viewer_interface*);
  bool printVertexIds() const;
  bool printEdgeIds() const;
  bool printFaceIds() const;
  void printAllIds();
  bool testDisplayId(double, double, double, CGAL::Three::Viewer_interface*)const;
  bool shouldDisplayIds(CGAL::Three::Scene_item *current_item) const;
  QString defaultSaveName() const
  {
    QString res = name();
    res.remove(" (selection)");
    return res;
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  void computeElements() const;
  void setKeepSelectionValid(SelectionTypes type);

protected:
  void init(Scene_face_graph_item* poly_item, QMainWindow* mw);

  Active_handle::Type get_active_handle_type()
  { return k_ring_selector.active_handle_type; }
  void set_active_handle_type(Active_handle::Type aht)
  { k_ring_selector.active_handle_type = aht; }
  void set_lasso_mode(bool b)
  { k_ring_selector.set_lasso_mode(b); }
  int get_k_ring() { return k_ring_selector.k_ring; }
  void set_k_ring(int k) { k_ring_selector.k_ring = k; }

  bool get_is_insert() { return is_insert; }
  void set_is_insert(bool i) { is_insert = i; }

public:
  typedef boost::unordered_set<fg_vertex_descriptor, CGAL::Handle_hash_function>    Selection_set_vertex;
  typedef boost::unordered_set<fg_face_descriptor, CGAL::Handle_hash_function>      Selection_set_facet;
  typedef boost::unordered_set<fg_edge_descriptor, CGAL::Handle_hash_function>    Selection_set_edge;

  Vertex_selection_map vertex_selection_map()
  {
    return this->poly_item->vertex_selection_map();
  }

  Face_graph* polyhedron()
  {
    return this->poly_item->polyhedron();
  }

  const Face_graph* polyhedron() const
  {
    return this->poly_item->polyhedron();
  }

    using Scene_polyhedron_item_decorator::draw;
    virtual void draw(CGAL::Three::Viewer_interface*) const;
    virtual void drawEdges() const { }
    virtual void drawEdges(CGAL::Three::Viewer_interface*) const;
    virtual void drawPoints(CGAL::Three::Viewer_interface*) const;

  bool supportsRenderingMode(RenderingMode m) const { return (m==Flat); }

  bool isEmpty() const {
    return selected_vertices.empty() && selected_edges.empty() && selected_facets.empty();
  }

  void reset_numbers();

  void compute_bbox() const
  {
    // Workaround a bug in g++-4.8.3:
    //   http://stackoverflow.com/a/21755207/1728537
    // Using boost::make_optional to copy-initialize 'item_bbox' hides the
    //   warning about '*item_bbox' not being initialized.
    // -- Laurent Rineau, 2014/10/30

    constVPmap vpm = get(CGAL::vertex_point, *polyhedron());

    boost::optional<CGAL::Bbox_3> item_bbox
      = boost::make_optional(false, CGAL::Bbox_3());


    for(Selection_set_vertex::const_iterator v_it = selected_vertices.begin();
        v_it != selected_vertices.end(); ++v_it) {

      if(item_bbox) { *item_bbox = *item_bbox + get(vpm,*v_it).bbox(); }
      else          {  item_bbox = get(vpm,*v_it).bbox(); }
    }

    for(Selection_set_edge::const_iterator e_it = selected_edges.begin();
        e_it != selected_edges.end(); ++e_it) {
      CGAL::Bbox_3 e_bbox = get(vpm,target(halfedge(*e_it,*polyhedron()),*polyhedron())).bbox();
      e_bbox = e_bbox + get(vpm,target(opposite(halfedge(*e_it,*polyhedron()),*polyhedron()),*polyhedron())).bbox();
        if(item_bbox) { *item_bbox = *item_bbox + e_bbox; }
        else          {  item_bbox = e_bbox; }
    }

    for(Selection_set_facet::const_iterator f_it = selected_facets.begin();
        f_it != selected_facets.end(); ++f_it) {
      fg_face_descriptor fd = *f_it;
      for(fg_halfedge_descriptor he : halfedges_around_face(halfedge(fd,*polyhedron()),*polyhedron())){
        if(item_bbox) { *item_bbox = *item_bbox + get(vpm,target(he,*polyhedron())).bbox(); }
        else          {  item_bbox = get(vpm,target(he,*polyhedron())).bbox(); }

      }
    }

    if(!item_bbox)
    {
      setBbox(this->poly_item->bbox());
      return;
    }
    setBbox(Bbox(item_bbox->xmin(),item_bbox->ymin(),item_bbox->zmin(),
                 item_bbox->xmax(),item_bbox->ymax(),item_bbox->zmax()));
  }

  bool save(const std::string& file_name) const {
    // update id fields before using
    if(selected_vertices.size() > 0
    ||selected_facets.size() > 0
    || (selected_edges.size() > 0 &&
        selected_vertices.empty() ))   { poly_item->face_graph()->collect_garbage(); }

    std::ofstream out(file_name.c_str());
    if(!out) { return false; }

    for(Selection_set_vertex::const_iterator it = selected_vertices.begin(); it != selected_vertices.end(); ++it)
      { out << get(boost::vertex_index, *polyhedron(), *it) << " "; }
    out << "\n";

    for(Selection_set_facet::const_iterator it = selected_facets.begin(); it != selected_facets.end(); ++it)
      { out << get(boost::face_index, *polyhedron(), *it) << " "; }
    out << "\n";

    for(Selection_set_edge::const_iterator it = selected_edges.begin(); it != selected_edges.end(); ++it)
    {
      fg_edge_descriptor ed = *it;
      out << get(boost::vertex_index, *polyhedron(), source(ed,*polyhedron())) << " "
          << get(boost::vertex_index, *polyhedron(),target(ed,*polyhedron())) << " ";
    }

    out << "\n";
    return true;
  }


  bool load(const std::string& file_name) {
    file_name_holder = file_name;
    return true;
  }


  // this function is called by selection_plugin, since at the time of the call of load(...)
  // we do not have access to selected polyhedron item
  bool actual_load(Scene_face_graph_item* poly_item, QMainWindow* mw)
  {

    init(poly_item, mw);

    std::vector<fg_vertex_descriptor> all_vertices;
    all_vertices.reserve(num_vertices(*polyhedron()));

    for(fg_vertex_descriptor vb : vertices(*polyhedron()))
      { all_vertices.push_back(vb); }

    std::vector<fg_face_descriptor> all_facets;
    all_facets.reserve(num_faces(*polyhedron()));
    for(fg_face_descriptor fb : faces(*polyhedron()))
      { all_facets.push_back(fb); }

    std::vector<fg_edge_descriptor> all_edges(edges(*polyhedron()).first, edges(*polyhedron()).second);

    std::ifstream in(file_name_holder.c_str());
    if(!in) { return false; }

    std::string line;
    std::size_t id, id2;

    if(!std::getline(in, line)) { compute_normal_maps(); return true; }
    std::istringstream vertex_line(line);
    while(vertex_line >> id) {
      if(id >= all_vertices.size()) { return false; }
      selected_vertices.insert(all_vertices[id]);
    }

    if(!std::getline(in, line)) { compute_normal_maps(); return true; }
    std::istringstream facet_line(line);
    while(facet_line >> id) {
      if(id >= all_facets.size()) { return false; }
      selected_facets.insert(all_facets[id]);
    }

    if(!std::getline(in, line)) { compute_normal_maps(); return true; }
    std::istringstream edge_line(line);
    while(edge_line >> id >> id2) {
      if(id >= all_edges.size() || id2 >= all_edges.size()) { return false; }
      fg_vertex_descriptor s = all_vertices[id];
      fg_vertex_descriptor t = all_vertices[id2];
      fg_halfedge_descriptor hd;
      bool exists;
      boost::tie(hd,exists) = halfedge(s,t,*polyhedron());
      if(! exists) { return false; }
      selected_edges.insert(edge(hd,*polyhedron()));
    }
    compute_normal_maps();
    return true;
  }

  //adds the content of temp_selection to the current selection
  void add_to_selection();
  // select all of `active_handle_type`(vertex, facet or edge)
  void select_all() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      select_all<fg_vertex_descriptor>(); break;
    case Active_handle::FACET:
    case Active_handle::CONNECTED_COMPONENT:
      select_all<fg_face_descriptor>(); break;
    case Active_handle::EDGE:
    case Active_handle::PATH:
      selected_edges.insert(edges(*polyhedron()).first, edges(*polyhedron()).second);
      invalidateOpenGLBuffers();
      CGAL::QGLViewer* v = *CGAL::QGLViewer::QGLViewerPool().begin();
      v->update();
      break;
    }
  }

  void select_boundary();
  void select_all_NT();
  // select all of vertex, facet or edge (use fg_vertex_descriptor, fg_face_descriptor, fg_edge_descriptor as template argument)
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
      clear<fg_vertex_descriptor>(); break;
    case Active_handle::FACET:
    case Active_handle::CONNECTED_COMPONENT:
      clear<fg_face_descriptor>(); break;
    case Active_handle::EDGE:
    case Active_handle::PATH:
      clear<fg_edge_descriptor>(); break;
    }
  }
  // select all of vertex, facet or edge (use fg_vertex_descriptor, fg_face_descriptor, fg_edge_descriptor as template argument)
  template<class HandleType>
  void clear() {

    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
    tr.container().clear();
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void clear_all(){
    clear<fg_vertex_descriptor>();
    clear<fg_face_descriptor>();
    clear<fg_edge_descriptor>();
  }

  boost::optional<std::size_t> get_minimum_isolated_component() {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      return get_minimum_isolated_component<fg_vertex_descriptor>();
    case Active_handle::FACET:
      return get_minimum_isolated_component<fg_face_descriptor>();
    default:
      return get_minimum_isolated_component<fg_edge_descriptor>();
    }
  }
  template<class HandleType> // use fg_vertex_descriptor, fg_face_descriptor, fg_edge_descriptor
  boost::optional<std::size_t> get_minimum_isolated_component() {
    Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
    Travel_isolated_components<Face_graph>::Minimum_visitor visitor;
    Travel_isolated_components<Face_graph>(*polyhedron()).travel<HandleType>
      (tr.iterator_begin(), tr.iterator_end(), tr.size(), tr.container(), visitor);
    return visitor.minimum;
  }

  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    switch(get_active_handle_type()) {
    case Active_handle::VERTEX:
      return select_isolated_components<fg_vertex_descriptor>(threshold);
    case Active_handle::FACET:
      return select_isolated_components<fg_face_descriptor>(threshold);
    default:
      return select_isolated_components<fg_edge_descriptor>(threshold);
    }
  }
  template<class HandleType> // use fg_vertex_descriptor, fg_face_descriptor, fg_edge_descriptor
  boost::optional<std::size_t> select_isolated_components(std::size_t threshold) {
    typedef Selection_traits<HandleType, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);
    typedef std::insert_iterator<typename Tr::Container> Output_iterator;
    Output_iterator out(tr.container(), tr.container().begin());

    Travel_isolated_components<Face_graph>::Selection_visitor<Output_iterator> visitor(threshold , out);
    Travel_isolated_components<Face_graph>(*polyhedron()).travel<HandleType>
      (tr.iterator_begin(), tr.iterator_end(), tr.size(), tr.container(), visitor);

    if(visitor.any_inserted) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
    return visitor.minimum_visitor.minimum;
  }

  void expand_or_reduce(int steps) {
    if (steps>0)
    {
      switch(get_active_handle_type()) {
        case Active_handle::VERTEX:
          expand_selection<fg_vertex_descriptor, boost::vertex_index_t>(steps);
        break;
        case Active_handle::FACET:
        case Active_handle::CONNECTED_COMPONENT:
          expand_selection<fg_face_descriptor, boost::face_index_t>(steps);
        break;
        case Active_handle::EDGE:
          expand_selection<fg_edge_descriptor, boost::edge_index_t>(steps);
        break;
        case Active_handle::PATH:
        break;
      }
    }
    else
    {
      switch(get_active_handle_type()) {
        case Active_handle::VERTEX:
          reduce_selection<fg_vertex_descriptor, boost::vertex_index_t>(-steps);
        break;
        case Active_handle::FACET:
        case Active_handle::CONNECTED_COMPONENT:
          reduce_selection<fg_face_descriptor, boost::face_index_t>(-steps);
        break;
        case Active_handle::EDGE:
          reduce_selection<fg_edge_descriptor, boost::edge_index_t>(-steps);
        break;
        case Active_handle::PATH:
        break;
      }
    }
  }

  template <class Handle, class IDmap>
  struct Is_selected_property_map{
    std::vector<bool>* is_selected_ptr;
    IDmap idmap;
    Is_selected_property_map()
      : is_selected_ptr(NULL) {}
    Is_selected_property_map(std::vector<bool>& is_selected, IDmap idmap)
      : is_selected_ptr( &is_selected), idmap(idmap) {}

    template<class H>
    std::size_t id(H h){ return get(idmap,h); }

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

  template <typename SelectionSet>
  struct Is_constrained_map
  {
    SelectionSet* m_set_ptr;

    typedef typename SelectionSet::key_type    key_type;
    typedef bool                               value_type;
    typedef bool                               reference;
    typedef boost::read_write_property_map_tag category;

    Is_constrained_map()
      : m_set_ptr(NULL)
    {}
    Is_constrained_map(SelectionSet* set_)
      : m_set_ptr(set_)
    {}
    friend bool get(const Is_constrained_map& map, const key_type& k)
    {
      CGAL_assertion(map.m_set_ptr != NULL);
      return map.m_set_ptr->count(k);
    }
    friend void put(Is_constrained_map& map, const key_type& k, const value_type b)
    {
      CGAL_assertion(map.m_set_ptr != NULL);
      if (b)  map.m_set_ptr->insert(k);
      else    map.m_set_ptr->erase(k);
    }
  };

  template <class Handle>
  struct Index_map
  {
    typedef Handle                             key_type;
    typedef std::size_t                        value_type;
    typedef value_type&                        reference;
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

  template <class Handle, class Tag>
  void expand_selection(unsigned int steps)
  {
    typedef Selection_traits<Handle, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    std::vector<bool> mark(tr.size(),false);

    for(Handle h :tr.container())
    {
      std::size_t id = tr.id(h);
      mark[id]=true;
    }

    typedef typename boost::property_map<Face_graph,Tag>::type PM;
    Tr::expand_selection(
      tr.container(),
      *this->poly_item->polyhedron(),
      steps,
      Is_selected_property_map<Handle,PM>(mark, get(Tag(),*this->poly_item->polyhedron())),
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

  template <class Handle, class Tag>
  void reduce_selection(unsigned int steps) {

    typedef Selection_traits<Handle, Scene_polyhedron_selection_item> Tr;
    Tr tr(this);

    std::vector<bool> mark(tr.size(),false);

    for(Handle h :tr.container())
      mark[tr.id(h)]=true;

    typedef typename boost::property_map<Face_graph,Tag>::type PM;
    Tr::reduce_selection(
      tr.container(),
      *this->poly_item->polyhedron(),
      steps,
      Is_selected_property_map<Handle,PM>(mark, get(Tag(),*this->poly_item->polyhedron())),
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

    for (Selection_set_edge::iterator eit = selected_edges.begin(); eit != selected_edges.end();)
    {
      if(//both incident faces will be erased
         (selected_facets.find(face(halfedge(*eit,*polyhedron()),*polyhedron())) != selected_facets.end()
          && selected_facets.find(face(opposite(halfedge(*eit,*polyhedron()),*polyhedron()),*polyhedron())) != selected_facets.end())
         //OR eit is a boundary edge and its incident face will be erased
         || (is_border(halfedge(*eit,*polyhedron()),*polyhedron())
             && (selected_facets.find(face(halfedge(*eit,*polyhedron()),*polyhedron())) != selected_facets.end()
                 || selected_facets.find(face(opposite(halfedge(*eit,*polyhedron()),*polyhedron()),*polyhedron())) != selected_facets.end())))
      {
        fg_edge_descriptor tmp = *eit;
        ++eit;
        selected_edges.erase(tmp);
      }
      else
        ++eit;
    }

    // erase facets from poly
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      CGAL::Euler::remove_face(halfedge(*fb,*polyhedron()), *polyhedron());
    }
    selected_facets.clear();
    invalidateOpenGLBuffers();
    changed_with_poly_item();
  }


  // This function writes into the id() of a Polyhedron
  void keep_connected_components() {
    if (selected_facets.empty()) { return; }

    Selection_traits<fg_face_descriptor, Scene_polyhedron_selection_item> trf(this);
    Selection_traits<fg_vertex_descriptor, Scene_polyhedron_selection_item> trv(this);

    PMP::keep_connected_components(*polyhedron(), trf.container());
    changed_with_poly_item();
  }

  bool export_selected_facets_as_polyhedron(Face_graph* out) {
    // Note: might be a more performance wise solution
    // assign sequential id to vertices neighbor to selected facets
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      for(fg_halfedge_descriptor hb : halfedges_around_face(halfedge(*fb,*polyhedron()),*polyhedron())){
        put(vertex_selection_map(),target(hb,*polyhedron()), 0);
      }
    }
    // construct point vector
    std::vector<EPICK::Point_3> points;
    points.reserve(selected_facets.size());
    VPmap vpm = get(CGAL::vertex_point, *polyhedron());
    unsigned int counter = 1;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb) {
      for(fg_halfedge_descriptor hb : halfedges_around_face(halfedge(*fb,*polyhedron()),*polyhedron())){
        if(get(vertex_selection_map(), target(hb,*polyhedron())) == 0) {
          put(vertex_selection_map(),target(hb,*polyhedron()), counter++);
          points.push_back(get(vpm,target(hb,*polyhedron())));
        }
      }
    }
    // construct polygon vector
    std::vector<std::vector<std::size_t> > polygons(selected_facets.size());
    counter = 0;
    for(Selection_set_facet::iterator fb = selected_facets.begin(); fb != selected_facets.end(); ++fb, ++counter) {
      for(fg_halfedge_descriptor hb : halfedges_around_face(halfedge(*fb,*polyhedron()),*polyhedron()))
      {
        polygons[counter].push_back(get(vertex_selection_map(),target(hb,*polyhedron())) -1);
      }
    }
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<Face_graph>(
      points, polygons, *out);

    return num_vertices(*out) > 0;
  }

  void select_sharp_edges(const double angle)
  {
    CGAL::detect_sharp_edges(polyhedron(), angle);

    boost::property_map<Face_graph,CGAL::edge_is_feature_t>::type is_feature = get(CGAL::edge_is_feature,*polyhedron());
    for(fg_edge_descriptor e : edges(*polyhedron()))
    {
      if (get(is_feature,e))
        selected_edges.insert(e);
    }
    invalidateOpenGLBuffers();
  }

  void changed_with_poly_item() {
    // no need to update indices
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
    compute_normal_maps();
    Q_EMIT itemChanged();
  }

  void setItemIsMulticolor(bool b) {
    poly_item->setItemIsMulticolor(b);
    poly_item->computeItemColorVectorAutomatically(b);
  }

  void selection_changed(bool);
  void updateDisplayedIds(QEvent *e);

Q_SIGNALS:
  void updateInstructions(QString);
  void simplicesSelected(CGAL::Three::Scene_item*);
  void isCurrentlySelected(Scene_facegraph_item_k_ring_selection*);
  void printMessage(QString);

public Q_SLOTS:
  void connectNewViewer(QObject* o)
  {
      o->installEventFilter(this);
  }
  void update_poly();
  void on_Ctrlz_pressed();
  void on_Ctrlu_pressed();
  void emitTempInstruct();
  void resetIsTreated();
  void save_handleType();
  void set_operation_mode(int mode);
  void set_highlighting(bool b);
  void invalidateOpenGLBuffers();
  void validateMoveVertex();
  void compute_normal_maps();
  void clearHL();
  QString toolTip() const;

  // slots are called by signals of polyhedron_k_ring_selector
  void selected(const std::set<fg_vertex_descriptor>& m)
  { has_been_selected(m); }
  void selected(const std::set<fg_face_descriptor>& m)
  { has_been_selected(m); }
  void selected(const std::set<fg_edge_descriptor>& m)
  { has_been_selected(m); }

  void selected_HL(const std::set<fg_vertex_descriptor>& m);
  void selected_HL(const std::set<fg_face_descriptor>& m);
  void selected_HL(const std::set<fg_edge_descriptor>& m);
  void poly_item_changed();
  void endSelection(){
    Q_EMIT simplicesSelected(this);
  }
  void toggle_insert(bool b)
  {
   is_insert = b;
  }

  void updateTick();
  void moveVertex();
protected:
  bool eventFilter(QObject* /*target*/, QEvent * gen_event)
  {
    if(gen_event->type() == QEvent::KeyPress
            && static_cast<QKeyEvent*>(gen_event)->key()==Qt::Key_Z)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(gen_event);
      if(keyEvent->modifiers().testFlag(Qt::ControlModifier)){
        on_Ctrlz_pressed();
        return true;
      }
    }
    else if(gen_event->type() == QEvent::KeyPress
            && static_cast<QKeyEvent*>(gen_event)->key()==Qt::Key_U)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(gen_event);
      if(keyEvent->modifiers().testFlag(Qt::ControlModifier))
      {
        on_Ctrlu_pressed();
        return true;
      }
    }

    updateDisplayedIds(gen_event);

    if(!visible() || !k_ring_selector.state.shift_pressing) { return false; }
    if(gen_event->type() == QEvent::Wheel)
    {
      QWheelEvent *event = static_cast<QWheelEvent*>(gen_event);
      int steps = event->angleDelta().y()/120;
      if(do_process)
      {
        expand_or_reduce(steps);
        do_process = false;
        QTimer::singleShot(0,this, [this](){do_process = true;});
      }
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
  void join_vertex(Scene_polyhedron_selection_item::fg_edge_descriptor ed)
  {
    CGAL::Euler::join_vertex(halfedge(ed, *polyhedron()),*polyhedron());
  }


  void selectPath(fg_vertex_descriptor vh);

//Generic class
  template<typename HandleRange>
  bool treat_selection(const HandleRange&)
  {
    qDebug()<<"ERROR : unknown HandleRange";
    return false;
  }

template<typename HandleRange>
  bool treat_classic_selection(const HandleRange& selection);
//Specialization for set<fg_vertex_descriptor>
  bool treat_selection(const std::set<fg_vertex_descriptor>& selection);
  bool treat_selection(const std::set<fg_edge_descriptor>& selection);
  bool treat_selection(const std::set<fg_face_descriptor>& selection);
  bool treat_selection(const std::vector<fg_face_descriptor>& selection);

  fg_face_descriptor get_face(fg_face_descriptor fh)
  { return fh; }
  fg_face_descriptor get_face(fg_vertex_descriptor)
  { return boost::graph_traits<Face_graph>::null_face(); }
  fg_face_descriptor get_face(fg_edge_descriptor)
  { return boost::graph_traits<Face_graph>::null_face(); }

  template<class HandleType>
  void has_been_selected(const std::set<HandleType>& selection)
  {
    if(!visible()) { return; }


    if (get_active_handle_type() == Active_handle::CONNECTED_COMPONENT)
    {
      Selection_traits<fg_edge_descriptor,
                       Scene_polyhedron_selection_item> tr(this);
      std::vector<bool> mark(tr.size(), false);
      for(fg_edge_descriptor e : selected_edges)
        mark[tr.id(e)] = true;
      std::vector<fg_face_descriptor> selected_cc;
      typedef typename boost::property_map<Face_graph,boost::edge_index_t>::type PM;
      CGAL::Polygon_mesh_processing::connected_component(
        get_face(*selection.begin()),
        *polyhedron(),
        std::back_inserter(selected_cc),
        CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(
                                                                           Is_selected_property_map<fg_edge_descriptor,PM>(mark, get(boost::edge_index,*polyhedron()))));
       treat_selection(selected_cc);
    }
    else
    {
      treat_selection(selection);
    }
  }

  typedef boost::property_map<Face_graph,boost::edge_index_t>::type Face_graph_edge_index_pm;

public:
  Is_selected_property_map<fg_edge_descriptor, Face_graph_edge_index_pm>
    selected_edges_pmap(std::vector<bool>& mark)
  {
    Selection_traits<fg_edge_descriptor,
      Scene_polyhedron_selection_item> tr(this);

    for (unsigned int i = 0; i < mark.size(); ++i)
      mark[i] = false;

    for(fg_edge_descriptor e : selected_edges)
      mark[tr.id(e)] = true;

    return Is_selected_property_map<fg_edge_descriptor, Face_graph_edge_index_pm>(mark, get(boost::edge_index,*polyhedron()));
  }

  Is_constrained_map<Selection_set_edge> constrained_edges_pmap()
  {
    return Is_constrained_map<Selection_set_edge>(&selected_edges);
  }

  Is_constrained_map<Selection_set_vertex> constrained_vertices_pmap()
  {
    return Is_constrained_map<Selection_set_vertex>(&selected_vertices);
  }

protected:
  // members
  std::string file_name_holder;
  Scene_facegraph_item_k_ring_selection k_ring_selector;
  // action state
  bool is_insert;
  bool do_process;

public:
// selection
  Selection_set_vertex selected_vertices;
  Selection_set_facet  selected_facets;
  Selection_set_edge   selected_edges; // stores one halfedge for each pair (halfedge with minimum address)

  Selection_set_vertex fixed_vertices;
  Selection_set_vertex temp_selected_vertices;
  Selection_set_facet  temp_selected_facets;
  Selection_set_edge   temp_selected_edges; // stores one halfedge for each pair (halfedge with minimum address)
  Selection_set_vertex HL_selected_vertices;
  Selection_set_facet  HL_selected_facets;
  Selection_set_edge   HL_selected_edges; // stores one halfedge for each pair (halfedge with minimum address)

protected :
  friend struct Scene_polyhedron_selection_item_priv;
  Scene_polyhedron_selection_item_priv *d;

  struct Update_indices_visitor
  {
    Selection_set_vertex& m_vertices;
    Selection_set_edge&   m_edges;
    Selection_set_facet&  m_facets;
    FaceGraph& m_mesh;

    Update_indices_visitor(Selection_set_vertex& vertices,
                           Selection_set_edge&   edges,
                           Selection_set_facet&  facets,
                           FaceGraph& mesh)
      : m_vertices(vertices), m_edges(edges), m_facets(facets), m_mesh(mesh)
    {}

    template<typename V2V, typename E2E, typename F2F>
    void operator()(const V2V& v2v, const E2E& e2e, const F2F& f2f)
    {
      //in *2* maps,
      //left is old simplex, right is new simplex

      Selection_set_vertex new_vertices;
      Selection_set_edge   new_edges;
      Selection_set_facet  new_facets;

      for(vertex_descriptor v : m_vertices)
      {
        if(v2v[v] != boost::graph_traits<SMesh>::null_vertex()
           && int(v2v[v])  < static_cast<int>(m_mesh.number_of_vertices()))
          new_vertices.insert(v2v[v]);
      }
      m_vertices.clear();
      m_vertices.insert(new_vertices.begin(), new_vertices.end());

      for (fg_edge_descriptor e : m_edges)
      {
        halfedge_descriptor h = halfedge(e, m_mesh);
        if(e2e[h] != boost::graph_traits<SMesh>::null_halfedge()
           && int(e2e[h])  < static_cast<int>(m_mesh.number_of_halfedges()))
          new_edges.insert(edge(e2e[h], m_mesh));
      }
      m_edges.clear();
      m_edges.insert(new_edges.begin(), new_edges.end());

      for (face_descriptor f : m_facets)
        if(f2f[f] != boost::graph_traits<SMesh>::null_face()
           && int(f2f[f])  < static_cast<int>(m_mesh.number_of_faces()))
          new_facets.insert(f2f[f]);
      m_facets.clear();
      m_facets.insert(new_facets.begin(), new_facets.end());
    }
  };


public:
  //statistics
  enum STATS {
    NB_VERTICES = 0,
    NB_CONNECTED_COMPOS,
    NB_BORDER_EDGES,
    IS_PURE_TRIANGLE,
    IS_PURE_QUAD,
    NB_DEGENERATED_FACES,
    HOLES,
    AREA,
    VOLUME,
    SELFINTER,
    NB_FACETS,
    MIN_AREA,
    MAX_AREA,
    MED_AREA,
    MEAN_AREA,
    MIN_ALTITUDE,
    MIN_ASPECT_RATIO,
    MAX_ASPECT_RATIO,
    MEAN_ASPECT_RATIO,
    GENUS,
    NB_EDGES,
    MIN_LENGTH,
    MAX_LENGTH,
    MID_LENGTH,
    MEAN_LENGTH,
    NB_NULL_LENGTH,
    MIN_ANGLE,
    MAX_ANGLE,
    MEAN_ANGLE
  };

  bool has_stats()const {return true;}
  QString computeStats(int type);
  CGAL::Three::Scene_item::Header_data header() const ;

  void set_num_faces(const std::size_t n);
};

#endif
