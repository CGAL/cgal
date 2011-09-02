#include "Scene_edit_polyhedron_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <boost/foreach.hpp>
#include <algorithm>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

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
typedef std::vector<std::pair<Kernel::Point_3, Kernel::Point_3> > Transform_vectors;
typedef Selected_vertices::iterator Selected_vertices_it;
typedef Polyhedron_vertex_deformation_index_map<Polyhedron> Vertex_index_map;
typedef Polyhedron_edge_deformation_length_map<Polyhedron> Edge_length_map;
typedef boost::iterator_property_map<std::vector<double>::iterator, Vertex_index_map> Dist_pmap;

#define PI 3.14159265359


struct Scene_edit_polyhedron_item_priv {
  Scene_polyhedron_item* poly_item;
  int handlesRegionSize;
  int interestRegionSize;
  bool geodesicCircle;
  bool sharpFeature;
  int usageScenario;
  bool selected_vertex_changed;
  bool selected_handles_moved;
  qglviewer::ManipulatedFrame* frame;
  Selected_vertices selected_handles;        
  Selected_vertices non_selected_handles;
  Selected_vertices selected_roi;
  Selected_vertices non_selected_roi;
  Vertex_handle selected_vertex;
  Transform_vectors selected_vectors;
  Transform_vectors non_selected_vectors;
  Kernel::Point_3 orig_pos;
  Kernel::Point_3 last_pos;
  Vertex_index_map* vertex_index_map;
  Edge_length_map* edge_length_map;
  std::vector<double> geodesic_distance;
  Dist_pmap* dist_pmap;
  std::vector<bool> is_sharp_vertices;
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
    put(*d->vertex_index_map, vh, idx++);
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
  find_sharp_vertices_1();
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
  if(!d->non_selected_handles.empty() || !d->non_selected_roi.empty() 
    || !d->selected_handles.empty() || !d->selected_roi.empty() ) {
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
          it = d->non_selected_handles.begin(),
          end = d->non_selected_handles.end();
        it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();

    color.set_rgb_color(0, 1.f, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
      it = d->selected_roi.begin(),
      end = d->selected_roi.end();
    it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();
    color.set_rgb_color(0, 1.f, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
          it = d->non_selected_roi.begin(),
          end = d->non_selected_roi.end();
        it != end; ++it)
    {
      const Kernel::Point_3& p = (*it)->point();
      ::glVertex3d(p.x(), p.y(), p.z());
    }
    ::glEnd();

    color.set_rgb_color(1.f, 0, 0);
    ::glBegin(GL_LINES);
    for (std::size_t i = 0; i < d->selected_vectors.size(); i++)
    {
      const Kernel::Point_3& p0 = d->selected_vectors[i].first;
      const Kernel::Point_3& p1 = d->selected_vectors[i].second;

      ::glVertex3d(p0.x(), p0.y(), p0.z());
      ::glVertex3d(p1.x(), p1.y(), p1.z());
    }
    ::glEnd();
    color.set_rgb_color(1.f, 0.5f, 0);
    ::glBegin(GL_LINES);
    for (std::size_t i = 0; i < d->non_selected_vectors.size(); i++)
    {
      const Kernel::Point_3& p0 = d->non_selected_vectors[i].first;
      const Kernel::Point_3& p1 = d->non_selected_vectors[i].second;

      ::glVertex3d(p0.x(), p0.y(), p0.z());
      ::glVertex3d(p1.x(), p1.y(), p1.z());
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

void
Scene_edit_polyhedron_item::setGeodesicCircle(bool status) {
  d->geodesicCircle = status;
  if(d->selected_vertex != Vertex_handle()) {
    vertex_has_been_selected(&*d->selected_vertex);
  }
}

void
Scene_edit_polyhedron_item::setSharpFeature(bool status) {
  d->sharpFeature = status;
  if(d->selected_vertex != Vertex_handle()) {
    vertex_has_been_selected(&*d->selected_vertex);
  }
}

void
Scene_edit_polyhedron_item::setUsageScenario(int i) {
  d->usageScenario = i;
}

void
Scene_edit_polyhedron_item::setSelectedVertexChanged(bool status) {
  d->selected_vertex_changed = status;
}

void
Scene_edit_polyhedron_item::setSelectedHandlesMoved(bool status) {
  d->selected_handles_moved = status;
}

qglviewer::ManipulatedFrame* 
Scene_edit_polyhedron_item::manipulatedFrame() {
  return d->frame;
}

void
Scene_edit_polyhedron_item::setSelectedVector(Kernel::Vector_3 translation_last)
{
  for (std::size_t i = 0; i < d->selected_vectors.size(); i++)
  {
    Kernel::Point_3 old_position = d->selected_vectors[i].second;
    Kernel::Point_3 new_position = old_position + translation_last;
    d->selected_vectors[i].second = new_position;
  }
}

struct Get_vertex_handle : public CGAL::Modifier_base<Polyhedron::HDS>
{
  Polyhedron::Vertex* vertex_ptr;
  Vertex_handle vh;
  void operator()(Polyhedron::HDS& hds) {
    vh = hds.vertex_handle(vertex_ptr);
  }
};

double Scene_edit_polyhedron_item::dihedral_angle(edge_descriptor e)
{
  Polyhedron* poly = d->poly_item->polyhedron();
  vertex_descriptor v0 = boost::target(e, *poly);
  vertex_descriptor v1 = boost::source(e, *poly);
  // Only one triangle for border edges
  if (boost::get(CGAL::edge_is_border, *poly, e)||boost::get(CGAL::edge_is_border, *poly, CGAL::opposite_edge(e, *poly)))
  {
    return -1;
  }
  else
  {
    edge_descriptor e_cw = CGAL::next_edge_cw(e, *poly);
    vertex_descriptor v2 = boost::source(e_cw, *poly);     
    edge_descriptor e_ccw = CGAL::next_edge_ccw(e, *poly);
    vertex_descriptor v3 = boost::source(e_ccw, *poly);
    
    Kernel::Vector_3 e01 = v1->point() - v0->point();
    Kernel::Vector_3 e02 = v2->point() - v0->point();
    Kernel::Vector_3 e03 = v3->point() - v0->point();
    // compute dihedral angle between v0_v1_v2 and v0_v1_v3
    // cross(e01, e02)
    Kernel::Vector_3 n012( e01[1]*e02[2]-e01[2]*e02[1], e01[2]*e02[0]-e01[0]*e02[2], e01[0]*e02[1]-e01[1]*e02[0] );
    // cross(e01, e03)
    Kernel::Vector_3 n013( e01[1]*e03[2]-e01[2]*e03[1], e01[2]*e03[0]-e01[0]*e03[2], e01[0]*e03[1]-e01[1]*e03[0] );
    // n012*n013 / (|n012|*|n013|)
    double cos_angle = (n012[0]*n013[0]+n012[1]*n013[1]+n012[2]*n013[2]) / std::sqrt(n012.squared_length()*n013.squared_length());
    return acos(cos_angle);
  }
}

void
Scene_edit_polyhedron_item::find_sharp_vertices()
{
  Polyhedron* poly = d->poly_item->polyhedron();
  d->is_sharp_vertices.clear();
  d->is_sharp_vertices.resize(poly->size_of_vertices(), false);
  std::vector<double> dihedral_angles;
  vertex_iterator vb, ve;
  for ( boost::tie(vb,ve) = boost::vertices(*poly); vb != ve; vb++ )
  {
    // compute dihedral angles
    in_edge_iterator eb, ee;
    dihedral_angles.clear();
    double sum = 0;
    for ( boost::tie(eb,ee) = boost::in_edges(*vb, *poly); eb != ee; eb++ )
    {
      double angle = dihedral_angle(*eb);
      if (angle != -1)
      {
        dihedral_angles.push_back(angle);
        sum += angle;
      }
    }
    // compute variance of angles
    double ave_angle = sum/dihedral_angles.size();
    double var = 0;
    for ( std::size_t i = 0; i < dihedral_angles.size(); i++ )
    {
      var += (dihedral_angles[i] - ave_angle)*(dihedral_angles[i] - ave_angle)/dihedral_angles.size();
    }
    if (var > 0.2)
    {
      d->is_sharp_vertices[ get(*d->vertex_index_map, *vb) ] = true ;
    }
    
  }
}

void
Scene_edit_polyhedron_item::find_sharp_vertices_1()
{
  Polyhedron* poly = d->poly_item->polyhedron();
  d->is_sharp_vertices.clear();
  d->is_sharp_vertices.resize(poly->size_of_vertices(), false);
  vertex_iterator vb, ve;
  for ( boost::tie(vb,ve) = boost::vertices(*poly); vb != ve; vb++ )
  {
    // compute dihedral angles
    in_edge_iterator eb, ee;
    for ( boost::tie(eb,ee) = boost::in_edges(*vb, *poly); eb != ee; eb++ )
    {
      double angle = dihedral_angle(*eb);
      if (angle < PI*3.0/4.0)
      {
        d->is_sharp_vertices[ get(*d->vertex_index_map, *vb) ] = true;
        break;
      }
    }
    
  }
}

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


// extend vertices until reaching the sharp edges
Selected_vertices extend_sharp_edge(Selected_vertices selected_vertices, std::vector<bool> is_sharp_vertices, Vertex_index_map vertex_index_map)
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
          int idx = get(vertex_index_map, other_v);
          if( !is_sharp_vertices[idx] )
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
    // re-compute geodesic distances relative to selected_vertex
    d->edge_length_map = new Edge_length_map(*poly);
    for ( Halfedge_handle eh = poly->edges_begin(); eh != poly->edges_end(); eh++ )
    {
      Kernel::Vector_3 edge = eh->vertex()->point() - eh->opposite()->vertex()->point();
      double edge_length = std::sqrt(edge.squared_length());
      boost::put(*d->edge_length_map, eh, edge_length);
    }
    boost::dijkstra_shortest_paths( *poly, vh, 
      boost::vertex_index_map (*d->vertex_index_map).
      weight_map (*d->edge_length_map).
      distance_map (*d->dist_pmap));
  }

  // set new handles and ROI regions
  Selected_vertices new_handles;
  new_handles.insert(vh);
  new_handles = extend_k_ring(new_handles, d->handlesRegionSize);
  Selected_vertices new_roi;
  if (d->sharpFeature)
  {
    new_roi.insert(vh);
    new_roi = extend_sharp_edge(new_roi, d->is_sharp_vertices, *d->vertex_index_map);
  }
  else
  {
    new_roi = new_handles;
    new_roi = extend_k_ring( new_roi, d->interestRegionSize - d->handlesRegionSize );
    if (d->geodesicCircle)
    {
      double radius = 0;
      BOOST_FOREACH(Vertex_handle v, new_handles)
      {
        double dist = boost::get(*d->dist_pmap, v);
        if ( dist> radius ) radius = dist;
      }
      new_handles = extend_circle(new_handles, radius, *d->dist_pmap);

      radius = 0;
      BOOST_FOREACH(Vertex_handle v, new_roi)
      {
        double dist = boost::get(*d->dist_pmap, v);
        if ( dist> radius ) radius = dist;
      }
      new_roi = extend_circle(new_roi, radius, *d->dist_pmap);
    }
  }
  
  // multiple handles scenario
  if (d->usageScenario == 1)
  {
    if (d->selected_vertex != vh)        // selected_vertex changed
    {
      BOOST_FOREACH(Vertex_handle v, new_handles)
      {
        if ( d->non_selected_handles.find(v) != d->non_selected_handles.end() ||
             d->selected_handles.find(v) != d->selected_handles.end() )
        {
          std::cerr << "New handles rejected: overlapped!\n";
          return;
        }
      }
      BOOST_FOREACH(Vertex_handle v, d->selected_handles)
      {
        // add the latest handles into non_selected_handles
        if ( d->non_selected_handles.find(v) == d->non_selected_handles.end() )
        {
          d->non_selected_handles.insert(v);
        }  
      }
      BOOST_FOREACH(Vertex_handle v, d->selected_roi)
      {
        // add the latest ROI into non_selected_roi
        if ( d->non_selected_roi.find(v) == d->non_selected_roi.end() )
        {
          d->non_selected_roi.insert(v);
        }  
      }
      for (std::size_t i = 0; i < d->selected_vectors.size(); i++)
      {
        std::pair<Kernel::Point_3, Kernel::Point_3> v = d->selected_vectors[i];
        Transform_vectors::iterator it = find(d->non_selected_vectors.begin(), d->non_selected_vectors.end(), v );
        // add the latest transform vectors into non_selected_vectors
        if ( it == d->non_selected_vectors.end() )
        {
          d->non_selected_vectors.push_back(v);
        }
      }
      
    }
    else
    {
      BOOST_FOREACH(Vertex_handle v, new_handles)
      {
        if ( d->non_selected_handles.find(v) != d->non_selected_handles.end() )
        {
          std::cerr << "New handles rejected: overlapped!\n";
          return;
        }
      }
      if (d->selected_handles_moved)
      {
        std::cerr << "Update of old handles rejected: already moved!\n";
        return;
      }
    }
    d->selected_vectors.clear();
    std::pair<Kernel::Point_3, Kernel::Point_3> v(vh->point(), vh->point());
    d->selected_vectors.push_back(v);  
  }

  d->selected_vertex = vh;
  d->selected_handles = new_handles;
  d->selected_handles_moved = false;

  d->selected_roi = new_roi;
  std::cerr << d->handlesRegionSize << " " << d->interestRegionSize << std::endl;

  
  const Kernel::Point_3& p = vh->point();
  d->orig_pos = p;
  d->last_pos = p;
  d->frame->setPosition(qglviewer::Vec(p.x(), p.y(), p.z()));
  connect(d->frame, SIGNAL(modified()),
          this, SIGNAL(modified()));       
  emit begin_edit();
}

void Scene_edit_polyhedron_item::vertex_has_been_selected_2(void* void_ptr) {
  Polyhedron* poly = d->poly_item->polyhedron();

  // Need a modifier to get access to the HDS, to get the vertex handle
  // from the vertex pointer.

  Get_vertex_handle get_vertex_handle;  
  get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);

  poly->delegate(get_vertex_handle);
  Vertex_handle vh = get_vertex_handle.vh;
  if (d->selected_vertex != vh)
  {
    d->selected_vertex_changed = true;
    d->selected_vertex = vh;
  }
  else
  {
    d->selected_vertex_changed = false;
  }
  if (d->selected_handles_moved)      // re-compute variance of dihedral angles for each vertex
  {
    find_sharp_vertices_1();
  }
  if (d->selected_vertex_changed || d->selected_handles_moved)        // selected_vertex is changed or moved
  {  
    // compute geodesic distances relative to selected_vertex
    d->edge_length_map = new Edge_length_map(*poly);
    for ( Halfedge_handle eh = poly->edges_begin(); eh != poly->edges_end(); eh++ )
    {
      Kernel::Vector_3 edge = eh->vertex()->point() - eh->opposite()->vertex()->point();
      double edge_length = std::sqrt(edge.squared_length());
      boost::put(*d->edge_length_map, eh, edge_length);
    }
    boost::dijkstra_shortest_paths( *poly, vh, 
      boost::vertex_index_map (*d->vertex_index_map).
      weight_map (*d->edge_length_map).
      distance_map (*d->dist_pmap));

    if (d->usageScenario == 1)
    {
      BOOST_FOREACH(Vertex_handle v, d->selected_handles)
      {
        // add the latest handles into non_selected_handles
        if ( d->non_selected_handles.find(v) == d->non_selected_handles.end() )
        {
          d->non_selected_handles.insert(v);
        }  
      }
      BOOST_FOREACH(Vertex_handle v, d->selected_roi)
      {
        // add the latest ROI into non_selected_roi
        if ( d->non_selected_roi.find(v) == d->non_selected_roi.end() )
        {
          d->non_selected_roi.insert(v);
        }  
      }
     
    }
  }

  d->selected_handles.clear();
  d->selected_handles.insert(vh);
  // compute the k-neighborhood of vh, with k==d->handlesRegionSize.
  d->selected_handles = extend_k_ring(d->selected_handles, d->handlesRegionSize);
  d->selected_handles_moved = false;

  d->selected_roi = d->selected_handles;
  std::cerr << d->handlesRegionSize << " " << d->interestRegionSize << std::endl;
  d->selected_roi = extend_k_ring( d->selected_roi, d->interestRegionSize - d->handlesRegionSize );

  //d->selected_roi = extend_sharp_edge(d->selected_roi, 0.2, d->dihedral_angle_variance, *d->vertex_index_map);
  // add geodesic distance constraints into handles and ROI vertices
  if (d->geodesicCircle)
  {
    double radius = 0;
    BOOST_FOREACH(Vertex_handle v, d->selected_handles)
    {
      double dist = boost::get(*d->dist_pmap, v);
      if ( dist> radius ) radius = dist;
    }
    d->selected_handles = extend_circle(d->selected_handles, radius, *d->dist_pmap);

    radius = 0;
    BOOST_FOREACH(Vertex_handle v, d->selected_roi)
    {
      double dist = boost::get(*d->dist_pmap, v);
      if ( dist> radius ) radius = dist;
    }
    d->selected_roi = extend_circle(d->selected_roi, radius, *d->dist_pmap);
  }


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
Scene_edit_polyhedron_item::non_selected_handles() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->non_selected_handles) {
    result << vh;
  }
  return result;
}

QList<Vertex_handle>
Scene_edit_polyhedron_item::selected_roi() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->selected_roi) {
    result << vh;
  }
  return result;
}

std::pair<Kernel::Point_3, Kernel::Point_3>
Scene_edit_polyhedron_item::selected_vector() const {
  std::pair<Kernel::Point_3, Kernel::Point_3> result;
  if (!d->selected_vectors.empty())
  {
    result = d->selected_vectors[0];
  }
  return result;
}

QList<Vertex_handle>
Scene_edit_polyhedron_item::non_selected_roi() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->non_selected_roi) {
    result << vh;
  }
  return result;
}

void Scene_edit_polyhedron_item::clear_non_selected_roi() {
  d->non_selected_roi.clear();
}

void Scene_edit_polyhedron_item::clear_selected_roi() {
  d->selected_roi.clear();
}


void Scene_edit_polyhedron_item::clear_non_selected_handles() {
  d->non_selected_handles.clear();
}

void Scene_edit_polyhedron_item::clear_selected_handles() {
  d->selected_handles.clear();
  d->selected_vertex = Vertex_handle();
}

void Scene_edit_polyhedron_item::clear_selected_vectors() {
  d->selected_vectors.clear();
}

void Scene_edit_polyhedron_item::clear_non_selected_vectors() {
  d->non_selected_vectors.clear();
}

int Scene_edit_polyhedron_item::usage_scenario() {
  return d->usageScenario;
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
