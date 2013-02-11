#include <CGAL/AABB_intersections.h>
#include "Scene_polyhedron_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <QVariant>
#include <list>

typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel, Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;

const char* aabb_property_name = "Scene_polyhedron_item aabb tree";

Input_facets_AABB_tree* get_aabb_tree(Scene_polyhedron_item* item)
{
  QVariant aabb_tree_property = item->property(aabb_property_name);
  if(aabb_tree_property.isValid()) {
    void* ptr = aabb_tree_property.value<void*>();
    return static_cast<Input_facets_AABB_tree*>(ptr);
  }
  else {
    Polyhedron* poly = item->polyhedron();
    if(poly) {
      Input_facets_AABB_tree* tree = 
        new Input_facets_AABB_tree(poly->facets_begin(),
                                   poly->facets_end());
      item->setProperty(aabb_property_name, 
                        QVariant::fromValue<void*>(tree));
      return tree;
    }
    else return 0;
  }
}

void delete_aabb_tree(Scene_polyhedron_item* item)
{
  QVariant aabb_tree_property = item->property(aabb_property_name);
  if(aabb_tree_property.isValid()) {
    void* ptr = aabb_tree_property.value<void*>();
    Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
    if(tree) {
      delete tree;
      tree = 0;
    }
    item->setProperty(aabb_property_name, QVariant());
  }
}

#include <QObject>
#include <QMenu>
#include <QAction>
#include <CGAL/gl_render.h>

Scene_polyhedron_item::Scene_polyhedron_item()
  : Scene_item_with_display_list(),
    poly(new Polyhedron),
    show_only_feature_edges_m(false),
    facet_picking_m(false),
    erase_next_picked_facet_m(false),
    plugin_has_set_color_vector_m(false)
{
  //init();
}

Scene_polyhedron_item::Scene_polyhedron_item(Polyhedron* const p)
  : Scene_item_with_display_list(),
    poly(p),
    show_only_feature_edges_m(false),
    facet_picking_m(false),
    erase_next_picked_facet_m(false),
    plugin_has_set_color_vector_m(false)
{
  init();
}

Scene_polyhedron_item::Scene_polyhedron_item(const Polyhedron& p)
  : Scene_item_with_display_list(),
    poly(new Polyhedron(p)),
    show_only_feature_edges_m(false),
    facet_picking_m(false),
    erase_next_picked_facet_m(false),
    plugin_has_set_color_vector_m(false)
{
  init();
}

// Scene_polyhedron_item::Scene_polyhedron_item(const Scene_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     poly(new Polyhedron(*item.poly)),
//     show_only_feature_edges_m(false)
// {
// }

Scene_polyhedron_item::~Scene_polyhedron_item()
{
  delete_aabb_tree(this);
  delete poly;
}

#include "Color_map.h"

void
Scene_polyhedron_item::
init()
{
  typedef Polyhedron::Facet_iterator Facet_iterator;
  
  if ( !plugin_has_set_color_vector_m )
  {
    // Fill indices map and get max subdomain value
    int max = 0;
    for(Facet_iterator fit = poly->facets_begin(), end = poly->facets_end() ;
        fit != end; ++fit)
    {
      max = (std::max)(max, fit->patch_id());
    }
    
    colors_.clear();
    compute_color_map(this->color(), max + 1, 
                      std::back_inserter(colors_));
  }
}


Scene_polyhedron_item* 
Scene_polyhedron_item::clone() const {
  return new Scene_polyhedron_item(*poly);
}

// Load polyhedron from .OFF file
bool
Scene_polyhedron_item::load(std::istream& in)
{
  in >> *poly;
  
  if ( in && !isEmpty() )
  {
    changed();
    return true;
  }
  return false;
}

// Write polyhedron to .OFF file
bool 
Scene_polyhedron_item::save(std::ostream& out) const
{
  out.precision(13);
  out << *poly;
  return out;
}

QString 
Scene_polyhedron_item::toolTip() const
{
  if(!poly)
    return QString();

  return QObject::tr("<p>Polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4</p>")
    .arg(this->name())
    .arg(poly->size_of_vertices())
    .arg(poly->size_of_halfedges()/2)
    .arg(poly->size_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

QMenu* Scene_polyhedron_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_polyhedron_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.trolltech.com/lastest/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {

    QAction* actionShowOnlyFeatureEdges = 
      menu->addAction(tr("Show only &feature edges"));
    actionShowOnlyFeatureEdges->setCheckable(true);
    actionShowOnlyFeatureEdges->setObjectName("actionShowOnlyFeatureEdges");
    connect(actionShowOnlyFeatureEdges, SIGNAL(toggled(bool)),
            this, SLOT(show_only_feature_edges(bool)));

    QAction* actionPickFacets = 
      menu->addAction(tr("Facets picking"));
    actionPickFacets->setCheckable(true);
    actionPickFacets->setObjectName("actionPickFacets");
    connect(actionPickFacets, SIGNAL(toggled(bool)),
            this, SLOT(enable_facets_picking(bool)));

    QAction* actionEraseNextFacet = 
      menu->addAction(tr("Erase next picked facet"));
    actionEraseNextFacet->setCheckable(true);
    actionEraseNextFacet->setObjectName("actionEraseNextFacet");
    connect(actionEraseNextFacet, SIGNAL(toggled(bool)),
            this, SLOT(set_erase_next_picked_facet(bool)));

    menu->setProperty(prop_name, true);
  }
  QAction* action = menu->findChild<QAction*>("actionPickFacets");
  if(action) action->setChecked(facet_picking_m);
  action = menu->findChild<QAction*>("actionEraseNextFacet");
  if(action) action->setChecked(erase_next_picked_facet_m);
  return menu;
}

void Scene_polyhedron_item::show_only_feature_edges(bool b)
{
  show_only_feature_edges_m = b;
  emit itemChanged();
}

void Scene_polyhedron_item::enable_facets_picking(bool b)
{
  facet_picking_m = b;
}

void Scene_polyhedron_item::set_erase_next_picked_facet(bool b)
{
  erase_next_picked_facet_m = b;
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::direct_draw() const {
  gl_render_facets(*poly,colors_);
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_polyhedron_item::direct_draw_edges() const {
  typedef Kernel::Point_3		Point;
  typedef Polyhedron::Edge_iterator	Edge_iterator;

  ::glBegin(GL_LINES);
  Edge_iterator he;
  if(!show_only_feature_edges_m) {
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
      if(he->is_feature_edge()) continue;
      const Point& a = he->vertex()->point();
      const Point& b = he->opposite()->vertex()->point();
      ::glVertex3d(a.x(),a.y(),a.z());
      ::glVertex3d(b.x(),b.y(),b.z());
    }
  }
  ::glColor3d(1.0, 0.0, 0.0);
  for(he = poly->edges_begin();
    he != poly->edges_end();
    he++)
  {
    if(!he->is_feature_edge()) continue;
    const Point& a = he->vertex()->point();
    const Point& b = he->opposite()->vertex()->point();
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
  }
  ::glEnd();
}

Polyhedron* 
Scene_polyhedron_item::polyhedron()       { return poly; }
const Polyhedron* 
Scene_polyhedron_item::polyhedron() const { return poly; }

bool
Scene_polyhedron_item::isEmpty() const {
  return (poly == 0) || poly->empty();
}

Scene_polyhedron_item::Bbox
Scene_polyhedron_item::bbox() const {
  const Kernel::Point_3& p = *(poly->points_begin());
  CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  for(Polyhedron::Point_iterator it = poly->points_begin();
      it != poly->points_end();
      ++it) {
    bbox = bbox + it->bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void
Scene_polyhedron_item::
changed()
{
  delete_aabb_tree(this);
  init();
  Base::changed();
}

void 
Scene_polyhedron_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
  if(facet_picking_m) {
    typedef Input_facets_AABB_tree Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    Tree* aabb_tree = get_aabb_tree(this);
    if(aabb_tree) {
      const Kernel::Point_3 ray_origin(orig_x, orig_y, orig_z);
      const Kernel::Vector_3 ray_dir(dir_x, dir_y, dir_z);
      const Kernel::Ray_3 ray(ray_origin, ray_dir);

      typedef std::list<Object_and_primitive_id> Intersections;
      Intersections intersections;

      aabb_tree->all_intersections(ray, std::back_inserter(intersections));

      Intersections::iterator closest = intersections.begin();
      if(closest != intersections.end()) {
        const Kernel::Point_3* closest_point = 
          CGAL::object_cast<Kernel::Point_3>(&closest->first);

        for(Intersections::iterator 
              it = boost::next(intersections.begin()),
              end = intersections.end();
            it != end; ++it)
        {
          if(! closest_point) {
            closest = it;
          }
          else {
            const Kernel::Point_3* it_point = 
              CGAL::object_cast<Kernel::Point_3>(&it->first);
            if(it_point && 
               (ray_dir * (*it_point - *closest_point)) < 0)
            {
              closest = it;
              closest_point = it_point;
            }
          }
        }
        if(closest_point) {
          Polyhedron::Facet_handle selected_fh = closest->second;

          // The computation of the nearest vertex may be costly.  Only
          // do it if some objects are connected to the signal
          // 'selected_vertex'.
          if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
          {
            Polyhedron::Halfedge_around_facet_circulator 
              he_it = selected_fh->facet_begin(),
              around_end = he_it;

            Polyhedron::Vertex_handle v = he_it->vertex(), nearest_v = v;

            Kernel::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                        v->point());

            while(++he_it != around_end) {
              v = he_it->vertex();
              Kernel::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                              v->point());
              if(new_sq_dist < sq_dist) {
                sq_dist = new_sq_dist;
                nearest_v = v;
              }
            }
            std::cerr << "Selected vertex: " << v->point() << std::endl;
            emit selected_vertex((void*)(&*nearest_v));
          }
          
          emit selected_facet((void*)(&*selected_fh));
          if(erase_next_picked_facet_m) {
            polyhedron()->erase_facet(selected_fh->halfedge());
            polyhedron()->normalize_border();
            set_erase_next_picked_facet(false);
            changed();
            emit itemChanged();
          }
          std::cerr << "Facet selected. patch_id="
                    << selected_fh->patch_id() << std::endl;
        }
      }
    }
  }
  Base::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
}

void Scene_polyhedron_item::update_vertex_indices()
{
  std::size_t id=0;
  for (Polyhedron::Vertex_iterator vit = polyhedron()->vertices_begin(), 
                                   vit_end = polyhedron()->vertices_end(); vit != vit_end; ++vit)
  {
    vit->id()=id++;
  }
}
void Scene_polyhedron_item::update_facet_indices()
{
  std::size_t id=0;
  for (Polyhedron::Facet_iterator  fit = polyhedron()->facets_begin(), 
                                   fit_end = polyhedron()->facets_end(); fit != fit_end; ++fit)
  {
    fit->id()=id++;
  }  
}
void Scene_polyhedron_item::update_halfedge_indices()
{
  std::size_t id=0;
  for (Polyhedron::Halfedge_iterator hit = polyhedron()->halfedges_begin(), 
                                     hit_end = polyhedron()->halfedges_end(); hit != hit_end; ++hit)
  {
    hit->id()=id++;
  }
}

#include "Scene_polyhedron_item.moc"
