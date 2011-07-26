#include "Scene_edit_polyhedron_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <boost/foreach.hpp>
#include <algorithm>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <QVariant>
#include <set>

#include <QObject>
#include <QMenu>
#include <QAction>
#include <CGAL/gl_render.h>

#include <QGLViewer/manipulatedFrame.h>

typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef std::set<Vertex_handle> Selected_vertices;
typedef Selected_vertices::iterator Selected_vertices_it;
typedef Polyhedron_vertex_deformation_index_map<Polyhedron> Vertex_index_map;
typedef Polyhedron_edge_deformation_length_map<Polyhedron> Edge_length_map;
typedef boost::iterator_property_map<std::vector<double>::iterator, Vertex_index_map> Dist_pmap;

struct Scene_edit_polyhedron_item_priv {
  Scene_polyhedron_item* poly_item;
  int handlesRegionSize;
  int interestRegionSize;
  qglviewer::ManipulatedFrame* frame;
  Selected_vertices selected_handles;        
  Selected_vertices handles_vertices;
  Selected_vertices roi_vertices;
  Selected_vertices roi_minus_handles_vertices;
  Vertex_handle selected_vertex;
  Kernel::Point_3 orig_pos;
  Kernel::Point_3 last_pos;
  Vertex_index_map* vertex_index_map;
  Edge_length_map* edge_length_map;
  std::vector<double> geodesic_distance;
  Dist_pmap* dist_pmap;
}; // end struct Scene_edit_polyhedron_item_priv

Scene_edit_polyhedron_item::Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item)
  : d(new Scene_edit_polyhedron_item_priv)
{
  d->poly_item = poly_item;
  d->handlesRegionSize = 0;
  d->frame = new ManipulatedFrame();
  d->frame->setProperty("item", QVariant::fromValue<QObject*>(this));
  if(!connect(poly_item, SIGNAL(selected_vertex(void*)),
              this, SLOT(vertex_has_been_selected(void*))))
    std::cerr << __FILE__ << ": connection failed!\n";
  poly_item->enable_facets_picking(true);

  d->vertex_index_map = new Vertex_index_map(*poly_item->polyhedron());
  int idx = 0;
  for ( Vertex_handle vh = poly_item->polyhedron()->vertices_begin(); vh != poly_item->polyhedron()->vertices_end(); vh++ )
  {
    boost::put(*d->vertex_index_map, vh, idx++);
  }
  d->edge_length_map = new Edge_length_map(*poly_item->polyhedron());
  for ( Halfedge_handle eh = poly_item->polyhedron()->edges_begin(); eh != poly_item->polyhedron()->edges_end(); eh++ )
  {
    Kernel::Vector_3 edge = eh->vertex()->point() - eh->opposite()->vertex()->point();
    double edge_length = std::sqrt(edge.squared_length());
    boost::put(*d->edge_length_map, eh, edge_length);
  }
  d->geodesic_distance.resize(boost::num_vertices(*poly_item->polyhedron()), 0);
  d->dist_pmap = new Dist_pmap(d->geodesic_distance.begin(), *d->vertex_index_map);
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
  delete d->frame;
  delete d->vertex_index_map;
  delete d->edge_length_map;
  delete d->dist_pmap;
  delete d;
}

Scene_edit_polyhedron_item* 
Scene_edit_polyhedron_item::clone() const {
  return 0;
}

QString 
Scene_edit_polyhedron_item::toolTip() const
{
  if(!d->poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(d->poly_item->polyhedron()->size_of_vertices())
    .arg(d->poly_item->polyhedron()->size_of_halfedges()/2)
    .arg(d->poly_item->polyhedron()->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

#include "opengl_tools.h"

void Scene_edit_polyhedron_item::draw() const {
  d->poly_item->direct_draw();
  if(!d->handles_vertices.empty() || !d->roi_vertices.empty() || !d->selected_handles.empty()) {
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);
    CGAL::GL::Color color;
    color.set_rgb_color(1.f, 0, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
      it = d->selected_handles.begin(),
      end = d->selected_handles.end();
    it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();
    color.set_rgb_color(1.f, 0.5f, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
          it = d->handles_vertices.begin(),
          end = d->handles_vertices.end();
        it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();
    color.set_rgb_color(0, 1.f, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
          it = d->roi_vertices.begin(),
          end = d->roi_vertices.end();
        it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();
  }
}

Polyhedron* 
Scene_edit_polyhedron_item::polyhedron()       { return d->poly_item->polyhedron(); }
const Polyhedron* 
Scene_edit_polyhedron_item::polyhedron() const { return d->poly_item->polyhedron(); }

bool
Scene_edit_polyhedron_item::isEmpty() const {
  return d->poly_item->isEmpty();
}

Scene_edit_polyhedron_item::Bbox
Scene_edit_polyhedron_item::bbox() const {
  return d->poly_item->bbox();
}


void
Scene_edit_polyhedron_item::
changed()
{
  d->poly_item->changed();
  Scene_item::changed();
  d->last_pos = current_position();
}

void 
Scene_edit_polyhedron_item::select(double orig_x,
                                   double orig_y,
                                   double orig_z,
                                   double dir_x,
                                   double dir_y,
                                   double dir_z)
{
  Scene_item::select(orig_x,
                     orig_y,
                     orig_z,
                     dir_x,
                     dir_y,
                     dir_z);
  d->poly_item->select(orig_x,
                       orig_y,
                       orig_z,
                       dir_x,
                       dir_y,
                       dir_z);
}

Scene_polyhedron_item* Scene_edit_polyhedron_item::to_polyhedron_item() const {
  return d->poly_item;
}

void 
Scene_edit_polyhedron_item::setHandlesRegionSize(int i) {
  if(i >= 0) {
    d->handlesRegionSize = i;
    if(d->selected_vertex != Vertex_handle()) {
      vertex_has_been_selected(&*d->selected_vertex);
    }
  }
}

void 
Scene_edit_polyhedron_item::setInterestRegionSize(int i) {
  if(i >= 0) {
    d->interestRegionSize = i;
    if(d->selected_vertex != Vertex_handle()) {
      vertex_has_been_selected(&*d->selected_vertex);
    }
  }
}

qglviewer::ManipulatedFrame* 
Scene_edit_polyhedron_item::manipulatedFrame() {
  return d->frame;
}

struct Get_vertex_handle : public CGAL::Modifier_base<Polyhedron::HDS>
{
  Polyhedron::Vertex* vertex_ptr;
  Vertex_handle vh;
  void operator()(Polyhedron::HDS& hds) {
    vh = hds.vertex_handle(vertex_ptr);
  }
};

Selected_vertices extend_once(Selected_vertices selected_vertices)
{
  std::set<Vertex_handle> original_set = selected_vertices;
  BOOST_FOREACH(Vertex_handle v, original_set) {
    Polyhedron::Halfedge_around_vertex_circulator 
      he_it = v->vertex_begin(), he_it_end(he_it);
    if(he_it != 0) {
      do {
        const Vertex_handle other_v = he_it->opposite()->vertex();
        if( selected_vertices.find(other_v) == selected_vertices.end() )
        {
          selected_vertices.insert(other_v);
        }
      } while(++he_it != he_it_end);
    }
  }
  return selected_vertices;
}

// extend k-neighboring vertices 
Selected_vertices extend_k_ring(Selected_vertices selected_vertices, int k)
{
  std::vector<Vertex_handle> selected_vertices_vector;
  selected_vertices_vector.insert(selected_vertices_vector.begin(), selected_vertices.begin(), selected_vertices.end());

  int idx_lv = 0;    // pointing the neighboring vertices on current level
  int idx_lv_end;
  for ( int lv = 0; lv < k; lv++ )
  {
    idx_lv_end = selected_vertices_vector.size();
    for (; idx_lv < idx_lv_end; idx_lv++)
    {
      Vertex_handle v = selected_vertices_vector[idx_lv];
      Polyhedron::Halfedge_around_vertex_circulator 
        he_it = v->vertex_begin(), he_it_end(he_it);
      if(he_it != 0) 
      {
        do {
          const Vertex_handle other_v = he_it->opposite()->vertex();
          std::vector<Vertex_handle>::iterator 
            it = std::find(selected_vertices_vector.begin(), selected_vertices_vector.end(), other_v);
          if (it == selected_vertices_vector.end())
          {
            selected_vertices_vector.push_back(other_v);
            selected_vertices.insert(other_v);
          }        
        }  while(++he_it != he_it_end);    
      }     
    }
  }
  
  return selected_vertices;
}


// extend vertices inside specific radius, so that the selected region is close circle
Selected_vertices extend_circle(Selected_vertices selected_vertices, double radius, Dist_pmap dist_pmap)
{
  std::vector<Vertex_handle> selected_vertices_vector;
  selected_vertices_vector.insert(selected_vertices_vector.begin(), selected_vertices.begin(), selected_vertices.end());
  bool new_vertex_selected = true; 
  int idx_lv = 0;    // pointing the neighboring vertices on current level
  int idx_lv_end;
  while (new_vertex_selected)
  {
    new_vertex_selected = false;
    idx_lv_end = selected_vertices_vector.size();
    for (; idx_lv < idx_lv_end; idx_lv++)
    {
      Vertex_handle v = selected_vertices_vector[idx_lv];
      Polyhedron::Halfedge_around_vertex_circulator 
        he_it = v->vertex_begin(), he_it_end(he_it);
      if(he_it != 0) {
        do {
          const Vertex_handle other_v = he_it->opposite()->vertex();
          if(  boost::get(dist_pmap, other_v) <= radius )
          {
            std::vector<Vertex_handle>::iterator it = std::find(selected_vertices_vector.begin(), selected_vertices_vector.end(), other_v);
            if (it == selected_vertices_vector.end())
            {
              selected_vertices_vector.push_back(other_v);
              selected_vertices.insert(other_v);
              new_vertex_selected = true;
            }        
          }
        } while(++he_it != he_it_end);
      }
    }
  }

  return selected_vertices;
}

void Scene_edit_polyhedron_item::vertex_has_been_selected(void* void_ptr) {
  Polyhedron* poly = d->poly_item->polyhedron();

  // Need a modifier to get access to the HDS, to get the vertex handle
  // from the vertex pointer.

  Get_vertex_handle get_vertex_handle;  
  get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);

  poly->delegate(get_vertex_handle);
  Vertex_handle vh = get_vertex_handle.vh;
  if (d->selected_vertex != vh)
  {
    d->selected_vertex = vh;
    // compute geodesic distances relative to selected_vertex
    boost::dijkstra_shortest_paths( *poly, vh, 
      boost::vertex_index_map (*d->vertex_index_map).
      weight_map (*d->edge_length_map).
      distance_map (*d->dist_pmap));
  }

  // std::cerr << "Selected vertex: " << void_ptr << " = " << vh->point()
  //           << std::endl;

  d->selected_handles.clear();
  d->selected_handles.insert(vh);

  // compute the k-neighborhood of vh, with k==d->handlesRegionSize.
  d->selected_handles = extend_k_ring(d->selected_handles, d->handlesRegionSize);

  d->roi_vertices = d->selected_handles;
  std::cerr << d->handlesRegionSize << " " << d->interestRegionSize << std::endl;
  d->roi_vertices = extend_k_ring(d->roi_vertices, d->interestRegionSize-d->handlesRegionSize);

  // add geodesic distance constraints into handles and ROI vertices
  double radius = 0;
  BOOST_FOREACH(Vertex_handle v, d->selected_handles)
  {
    double dist = boost::get(*d->dist_pmap, v);
    if ( dist> radius ) radius = dist;
  }
  d->selected_handles = extend_circle(d->selected_handles, radius, *d->dist_pmap);

  radius = 0;
  BOOST_FOREACH(Vertex_handle v, d->roi_vertices)
  {
    double dist = boost::get(*d->dist_pmap, v);
    if ( dist> radius ) radius = dist;
  }
  d->roi_vertices = extend_circle(d->roi_vertices, radius, *d->dist_pmap);

  // compute the difference
  // YX: not really need this
  d->roi_minus_handles_vertices.clear();
  std::set_difference(d->roi_vertices.begin(),
                      d->roi_vertices.end(),
                      d->selected_handles.begin(),
                      d->selected_handles.end(),
                      std::inserter(d->roi_minus_handles_vertices,
                                    d->roi_minus_handles_vertices.begin()));

  const Kernel::Point_3& p = vh->point();
  d->orig_pos = p;
  d->last_pos = p;
  d->frame->setPosition(qglviewer::Vec(p.x(), p.y(), p.z()));
  connect(d->frame, SIGNAL(modified()),
          this, SIGNAL(modified()));       
  emit begin_edit();
}

Vertex_handle 
Scene_edit_polyhedron_item::selected_vertex() const {
  return d->selected_vertex;
}

QList<Vertex_handle>
Scene_edit_polyhedron_item::selected_handles() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->selected_handles) {
    result << vh;
  }
  return result;
}

QList<Vertex_handle>
Scene_edit_polyhedron_item::handles_vertices() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->handles_vertices) {
    result << vh;
  }
  return result;
}

QList<Vertex_handle>
Scene_edit_polyhedron_item::vertices_in_region_of_interest() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->roi_vertices) {
    result << vh;
  }
  return result;
}

void Scene_edit_polyhedron_item::clear_roi_vertices() {
  d->roi_vertices.clear();
  d->roi_minus_handles_vertices.clear();
}

void Scene_edit_polyhedron_item::insert_roi(Polyhedron::Vertex_handle vh) {
  d->roi_vertices.insert(vh);
}

void Scene_edit_polyhedron_item::clear_handles_vertices() {
  d->handles_vertices.clear();
}

void Scene_edit_polyhedron_item::clear_selected_handles() {
  d->selected_handles.clear();
  d->selected_vertex = Vertex_handle();
}

void Scene_edit_polyhedron_item::insert_handle(Polyhedron::Vertex_handle vh) {
  d->handles_vertices.insert(vh);
}


Kernel::Point_3 Scene_edit_polyhedron_item::current_position() const {
  const qglviewer::Vec vec = d->frame->position();
  return Kernel::Point_3(vec.x, vec.y, vec.z);
}

Kernel::Point_3 Scene_edit_polyhedron_item::last_position() const {
  return d->last_pos;
}

Kernel::Point_3 Scene_edit_polyhedron_item::original_position() const {
  return d->orig_pos;
}

void Scene_edit_polyhedron_item::setVisible(bool b) {
  d->poly_item->setVisible(b);
  Scene_item::setVisible(b);
}
void Scene_edit_polyhedron_item::setColor(QColor c) {
  d->poly_item->setColor(c);
  Scene_item::setColor(c);
}
void Scene_edit_polyhedron_item::setName(QString n) {
  Scene_item::setName(n);
  n.replace(" (edit)", "");
  d->poly_item->setName(n);
}
void Scene_edit_polyhedron_item::setRenderingMode(RenderingMode m) {
  d->poly_item->setRenderingMode(m);
  Scene_item::setRenderingMode(m);
}

#include "Scene_edit_polyhedron_item.moc"
