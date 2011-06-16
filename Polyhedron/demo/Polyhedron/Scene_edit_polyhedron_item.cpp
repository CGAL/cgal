#include "Scene_edit_polyhedron_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include <boost/foreach.hpp>

#include <QVariant>
#include <set>

#include <QObject>
#include <QMenu>
#include <QAction>
#include <CGAL/gl_render.h>

#include <QGLViewer/manipulatedFrame.h>

typedef Polyhedron::Vertex_handle Vertex_handle;
typedef std::set<Vertex_handle> Selected_vertices;
typedef Selected_vertices::iterator Selected_vertices_it;

struct Scene_edit_polyhedron_item_priv {
  Scene_polyhedron_item* poly_item;
  int zone_size;
  qglviewer::ManipulatedFrame* frame;
  Selected_vertices selected_vertices;
  Vertex_handle selected_vertex;
  Kernel::Point_3 orig_pos;
}; // end struct Scene_edit_polyhedron_item_priv

Scene_edit_polyhedron_item::Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item)
  : d(new Scene_edit_polyhedron_item_priv)
{
  d->poly_item = poly_item;
  d->zone_size = 0;
  d->frame = new ManipulatedFrame();
  d->frame->setProperty("item", QVariant::fromValue<QObject*>(this));
  if(!connect(poly_item, SIGNAL(selected_vertex(void*)),
              this, SLOT(vertex_has_been_selected(void*))))
    std::cerr << __FILE__ << ": connection failed!\n";
  poly_item->enable_facets_picking(true);
}

Scene_edit_polyhedron_item::~Scene_edit_polyhedron_item()
{
  delete d->frame;
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
  if(!d->selected_vertices.empty()) {
    CGAL::GL::Point_size point_size; point_size.set_point_size(5);
    CGAL::GL::Color color; color.set_rgb_color(0, 0, 0);
    ::glBegin(GL_POINTS);
    for(Selected_vertices_it 
          it = d->selected_vertices.begin(),
          end = d->selected_vertices.end();
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
  d->orig_pos = current_position();
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
Scene_edit_polyhedron_item::setZoneSize(int i) {
  if(i >= 0) {
    std::cerr << "item \"" << qPrintable(name()) 
              << "\".setZoneSize(" << i << ")\n";
    d->zone_size = i;
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

void Scene_edit_polyhedron_item::vertex_has_been_selected(void* void_ptr) {
  Polyhedron* poly = d->poly_item->polyhedron();

  // Need a modifier to get access to the HDS, to get the vertex handle
  // from the vertex pointer.

  Get_vertex_handle get_vertex_handle;  
  get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);

  poly->delegate(get_vertex_handle);
  Vertex_handle vh = get_vertex_handle.vh;

  std::cerr << "Selected vertex: " << void_ptr << " = " << vh->point()
            << std::endl;
  d->selected_vertices.clear();

  d->selected_vertices.insert(vh);

  std::cerr << "d->zone_size = " << d->zone_size << std::endl;
  // Naive way to compute the k-neighborhood of vh, with k==d->zone_size.
  for(int i = 0; i < d->zone_size; ++i) {
    std::set<Vertex_handle> selected_vertices;
    for(Selected_vertices_it 
          it = d->selected_vertices.begin(),
          end = d->selected_vertices.end();
        it != end; ++it) {
      selected_vertices.insert(*it);
    }
    BOOST_FOREACH(Vertex_handle v, selected_vertices) {
      Polyhedron::Halfedge_around_vertex_circulator 
        he_it = v->vertex_begin(), he_it_end(he_it);
      if(he_it != 0) {
        do {
          const Vertex_handle other_v = he_it->opposite()->vertex();
          if( d->selected_vertices.find(other_v) == d->selected_vertices.end() )
          {
            d->selected_vertices.insert(other_v);
          }
        } while(++he_it != he_it_end);
      }
    }
  }
  const Kernel::Point_3& p = vh->point();
  d->orig_pos = p;
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
Scene_edit_polyhedron_item::selected_vertices() const {
  QList<Vertex_handle> result;
  BOOST_FOREACH(Vertex_handle vh, d->selected_vertices) {
    result << vh;
  }
  return result;
}

Kernel::Point_3 Scene_edit_polyhedron_item::current_position() const {
  const qglviewer::Vec vec = d->frame->position();
  return Kernel::Point_3(vec.x, vec.y, vec.z);
}

Kernel::Point_3 Scene_edit_polyhedron_item::original_position() const {
  return d->orig_pos;
}

#include "Scene_edit_polyhedron_item.moc"
