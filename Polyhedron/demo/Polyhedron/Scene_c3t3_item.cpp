#include <CGAL/Mesh_3/io_signature.h>

#include "Show_point_dialog.h"

#include <boost/foreach.hpp>

#include "Scene_c3t3_item.h"
#include <string>
#include <iostream>
#include <CGAL/IO/io.h>
#include <QInputDialog>
#include <QMenu>
#include <QVariant>
#include <QPainter>

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include "DisplayList.h"
#include <CGAL/gl.h>

#include <CGAL/IO/File_binary_mesh_3.h>

namespace {
  void CGALglcolor(QColor c)
  {
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

struct Scene_c3t3_item_priv {
  Scene_c3t3_item_priv(Scene_c3t3_item* item)
    : item(item),
      draw_triangles_duals(false),
      draw_clipping_plane(true),
      draw_triangles_display_list(4) {}

  typedef std::map<C3t3::Surface_patch_index, QColor> Color_map;
  typedef C3t3::Triangulation Tr;

  qglviewer::ManipulatedFrame* newManipulatedFrame() const {
    qglviewer::ManipulatedFrame* frame =
      new qglviewer::ManipulatedFrame();
    QObject::connect(frame, SIGNAL(modified()),
                     item, SLOT(frameIsModified()));
    return frame;
  }

  void begin_OpenGL_clip_plane() {
    Kernel::Plane_3 plane = item->plane();
    if(item->draw_tets) {
      // If one draws the tetrahedra of the cut plane, then offset the
      // OpenGL cut plane a bit, to make surfaces go further than the cut
      // plane.
      Kernel::Point_3 point = plane.point();
      Kernel::Vector_3 normal = plane.orthogonal_vector();
      // const double len = std::sqrt(normal*normal);
      // if(len != 0) normal = normal * ( item->len_diagonal() * 0.5 / len );
      // point = point + 0.02 * normal;
      plane = Kernel::Plane_3(point, normal);
    }
    GLdouble clip_plane[4];
    clip_plane[0] = -plane.a();
    clip_plane[1] = -plane.b();
    clip_plane[2] = -plane.c();
    clip_plane[3] = -plane.d();

    ::glClipPlane(GL_CLIP_PLANE0, clip_plane);
    ::glEnable(GL_CLIP_PLANE0);
  }

  void end_OpenGL_clip_plane() {
    ::glDisable(GL_CLIP_PLANE0);
  }

  void draw_cut_plane() {
    if(draw_clipping_plane) {
      ::glPushMatrix();
      ::glMultMatrixd(item->frame->matrix());
      QGLViewer::drawGrid((float)item->complex_diag());
      ::glPopMatrix();
    }
  }

  void draw_spheres() {
    // create the display list for a single sphere
    if(item->sphere_display_list == 0) {
      item->sphere_display_list = glGenLists(1);
      if(item->sphere_display_list == 0)
        std::cerr << "ERROR: Cannot create display list!\n";
      if(item->quadric == 0)
        item->quadric = gluNewQuadric();
      if(item->quadric == 0)
        std::cerr << "ERROR: Cannot create GLU quadric!\n";
      glNewList(item->sphere_display_list, GL_COMPILE);
      gluSphere(item->quadric, 1., 10, 10);
      glEndList();
      if(glGetError() != GL_NO_ERROR)
        std::cerr << gluErrorString(glGetError());
    }
    if(this->spheres_display_list.begin_draw()) {
      // force wireframe for protecting spheres
      GLint polygon_mode[2];
      ::glGetIntegerv(GL_POLYGON_MODE, &polygon_mode[0]);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      for(Tr::Finite_vertices_iterator
            vit = item->c3t3().triangulation().finite_vertices_begin(),
            end =  item->c3t3().triangulation().finite_vertices_end();
          vit != end; ++vit)
      {
        typedef Tr::Vertex_handle Vertex_handle;
        std::vector<Vertex_handle> incident_vertices;
        item->c3t3().triangulation().incident_vertices(vit, std::back_inserter(incident_vertices));
        bool red = vit->is_special();
        for(std::vector<Vertex_handle>::const_iterator
              vvit = incident_vertices.begin(), end = incident_vertices.end();
            vvit != end; ++vvit)
        {
          if(Kernel::Sphere_3(vit->point().point(),
                              vit->point().weight()).bounded_side((*vvit)->point().point())
             == CGAL::ON_BOUNDED_SIDE)
            red = true;
        }
        if(red)
          CGALglcolor(Qt::red);
        else
          CGALglcolor(item->color().darker(250));
        item->draw_sphere(vit->point());
      }
      ::glPolygonMode(GL_FRONT_AND_BACK, polygon_mode[0]);
    }
    this->spheres_display_list.end_draw();
  }

public:
  Scene_c3t3_item* item;
  Color_map color_map;

  bool draw_triangles_duals;
  bool draw_clipping_plane;

  OpenGL_display_lists duals_display_list;
  OpenGL_display_lists spheres_display_list;
  OpenGL_display_lists cut_display_list;
  OpenGL_display_lists draw_triangles_display_list;
}; //end struct Scene_c3t3_item_priv

void Scene_c3t3_item::frameIsModified() {
  d->cut_display_list.invalidate();
  d->draw_triangles_display_list.invalidate(2);
  d->draw_triangles_display_list.invalidate(3);
}

Scene_c3t3_item::Scene_c3t3_item()
  : d(new Scene_c3t3_item_priv(this))
  , parent_item(0)
  , c3t3_()
  , frame(d->newManipulatedFrame())
  , sphere_display_list(0)
  , quadric(0)
{
  draw_spheres = false;
  spheres_drawn_radius = 0.;
  draw_tets = false;
}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3, Scene_item* parent)
  : d(new Scene_c3t3_item_priv(this))
  , parent_item(parent)
  , c3t3_(c3t3)
  , frame(d->newManipulatedFrame())
  , sphere_display_list(0)
  , quadric(0)
{
  draw_spheres = false;
  spheres_drawn_radius = 0.;
  draw_tets = false;
  reset_cut_plane();
  changed();
}

QMenu* Scene_c3t3_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_c3t3_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.trolltech.com/lastest/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {

    menu->addSeparator();

    QAction* newRadiiAction =
      menu->addAction(tr("Display spheres with new radius..."));
    QAction* restoreRadiiAction =
      menu->addAction(tr("Display spheres with original radii"));

    newRadiiAction->setProperty("restoreRadiiAction",
                                qVariantFromValue((QObject*)restoreRadiiAction));

    restoreRadiiAction->setEnabled(false);
    connect(newRadiiAction, SIGNAL(triggered()),
            this, SLOT(new_spheres_drawn_radius()));
    connect(restoreRadiiAction, SIGNAL(triggered()),
            this, SLOT(restore_spheres_drawn_radii()));

    menu->addSeparator();

    QAction* actionShowCutPlane =
      menu->addAction(tr("Show cut &plane"));
    actionShowCutPlane->setCheckable(true);
    actionShowCutPlane->setObjectName("actionShowCutPlane");
    connect(actionShowCutPlane, SIGNAL(toggled(bool)),
            this, SLOT(show_cut_plane(bool)));

    QAction* recenterCutPlaneAction =
      menu->addAction(tr("Change the cut plane center"));
    connect(recenterCutPlaneAction, SIGNAL(triggered()),
            this, SLOT(recenter_cut_plane()));

    QAction* resetCutPlaneAction =
      menu->addAction(tr("Reset the cut plane position"));
    connect(resetCutPlaneAction, SIGNAL(triggered()),
            this, SLOT(reset_cut_plane()));

    menu->addSeparator();

    QAction* actionShowTets= menu->addAction(tr("Show &tetrahedra"));
    actionShowTets->setCheckable(true);
    actionShowTets->setObjectName("actionShowTets");
    connect(actionShowTets, SIGNAL(toggled(bool)),
            this, SLOT(show_tetrahedra(bool)));

    QAction* actionShowDuals= menu->addAction(tr("Show triangles &duals"));
    actionShowDuals->setCheckable(true);
    actionShowDuals->setObjectName("actionShowDuals");
    connect(actionShowDuals, SIGNAL(toggled(bool)),
            this, SLOT(show_triangles_duals(bool)));
    connect(actionShowTets, SIGNAL(toggled(bool)),
            actionShowDuals, SLOT(setEnabled(bool)));

    QAction* actionShowSpheres =
      menu->addAction(tr("Show protecting &spheres"));
    actionShowSpheres->setCheckable(true);
    actionShowSpheres->setObjectName("actionShowSpheres");
    connect(actionShowSpheres, SIGNAL(toggled(bool)),
            this, SLOT(show_spheres(bool)));

    menu->setProperty(prop_name, true);
  }

  QAction* action = menu->findChild<QAction*>("actionShowTets");
  if(action) action->setChecked(draw_tets);
  action = menu->findChild<QAction*>("actionShowSpheres");
  if(action) action->setChecked(draw_spheres);
  action = menu->findChild<QAction*>("actionShowCutPlane");
  if(action) action->setChecked(d->draw_clipping_plane);

  return menu;
}

void Scene_c3t3_item::recenter_cut_plane()
{
  Show_point_dialog dialog(0);
  int i = dialog.exec();
  if( i == QDialog::Accepted &&
      dialog.has_correct_coordinates() )
  {
    frame->setPosition((float)dialog.get_x(),
                       (float)dialog.get_y(),
                       (float)dialog.get_z());
    emit changed();
  }
}

void Scene_c3t3_item::restore_spheres_drawn_radii()
{
  QAction* restoreRadiiAction = qobject_cast<QAction*>(sender());
  if(restoreRadiiAction) restoreRadiiAction->setEnabled(false);

  spheres_drawn_radius = 0.;
  emit itemChanged();
}

void Scene_c3t3_item::new_spheres_drawn_radius(double d)
{
  if(d > 0) {
    spheres_drawn_radius = d;
    emit itemChanged();
  }
}

void Scene_c3t3_item::new_spheres_drawn_radius()
{
  bool ok = true;
  double d = QInputDialog::getDouble(NULL,
                                     tr("Display spheres with new radius..."),
                                     tr("Radius:"),
                                     spheres_drawn_radius, // value
                                     0.,          // min
                                     2147483647., // max
                                     10,          // decimals
                                     &ok);

  if(ok && d >= 0) {
    QAction* newRadiiAction = qobject_cast<QAction*>(sender());
    if(newRadiiAction) {
      QVariant v = newRadiiAction->property("restoreRadiiAction");
      QAction* restoreRadiiAction = qobject_cast<QAction*>(v.value<QObject*>());
      if(restoreRadiiAction)
        restoreRadiiAction->setEnabled(d>0);
    }
    new_spheres_drawn_radius(d);
  }
}

void
Scene_c3t3_item::reset_cut_plane() {
  const Bbox& bbox = this->bbox();
  const float xcenter = static_cast<float>((bbox.xmax+bbox.xmin)/2.);
  const float ycenter = static_cast<float>((bbox.ymax+bbox.ymin)/2.);
  const float zcenter = static_cast<float>((bbox.zmax+bbox.zmin)/2.);

  frame->setPosition(qglviewer::Vec(xcenter, ycenter, zcenter));
}


Scene_c3t3_item::~Scene_c3t3_item()
{
  if(frame != 0) delete frame;
  if(quadric != 0)
    gluDeleteQuadric(quadric);
  if(sphere_display_list  != 0)
    glDeleteLists(sphere_display_list, 1);
  delete d;
}

void
Scene_c3t3_item::setColor(QColor c) {
  Scene_item::setColor(c);
  changed();
}

bool
Scene_c3t3_item::save_binary(std::ostream& os) const
{
  return CGAL::Mesh_3::save_binary_file(os, c3t3());
}

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// bool
// Scene_c3t3_item::export_boundary_to_OFF(std::ostream& os) const {
//   // CGAL::output_surface_facets_to_off(os, this->c3t3());
//   std::cerr << "export to OFF" << std::endl;
//   CGAL::output_oriented_surface_facets_to_off(os, this->c3t3().triangulation());
//   return os.good();
// }

bool
Scene_c3t3_item::load_binary(std::istream& is)
{
  if(!CGAL::Mesh_3::load_binary_file(is, c3t3())) return false;
  // if(!c3t3().triangulation().is_valid()) std::cerr << "INVALID\n";
  if(is && frame == 0) {
    frame = d->newManipulatedFrame();
  }
  reset_cut_plane();
  if(is.good()) {
    changed();
    return true;
  }
  else
    return false;
}

Kernel::Plane_3
Scene_c3t3_item::plane() const
{
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n =
  frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
}


Scene_c3t3_item::Bbox
Scene_c3t3_item::bbox() const
{
  if(isEmpty())
    return Bbox();
  else {
    CGAL::Bbox_3 result = c3t3().triangulation().finite_vertices_begin()->point().bbox();
    for(Tr::Finite_vertices_iterator
        vit = ++c3t3().triangulation().finite_vertices_begin(),
        end = c3t3().triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      result = result + vit->point().bbox();
    }
    return Bbox(result.xmin(), result.ymin(), result.zmin(),
                result.xmax(), result.ymax(), result.zmax());
  }
}


QString
Scene_c3t3_item::toolTip() const
{
  int number_of_tets = 0;
  for(Tr::Finite_cells_iterator
      cit = c3t3().triangulation().finite_cells_begin(),
      end = c3t3().triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if( c3t3().is_in_complex(cit) )
      ++number_of_tets;
  }
  return tr("<p><b>%5</b><br /><i>3D complex in a 3D triangulation</i></p>"
            "<p>Number of vertices: %1<br />"
            "Number of surface facets: %2<br />"
            "Number of volume tetrahedra: %3<br />"
            "Number of different colors: %4</p>")
    .arg(c3t3().triangulation().number_of_vertices())
    .arg(c3t3().number_of_facets())
    .arg(number_of_tets)
    .arg(d->color_map.size())
    .arg(name());
}

#include "Color_map.h"

void Scene_c3t3_item::changed() {
  Scene_item_with_display_list::changed();
  d->spheres_display_list.invalidate();
  d->cut_display_list.invalidate();
  d->duals_display_list.invalidate();
  d->draw_triangles_display_list.invalidate(0);
  d->draw_triangles_display_list.invalidate(1);
  d->draw_triangles_display_list.invalidate(2);
  d->draw_triangles_display_list.invalidate(3);
  // std::set<int> indices;
  d->color_map.clear();
  for(C3t3::Facet_iterator
      fit = c3t3().facets_begin(),
      end = c3t3().facets_end();
      fit != end; ++fit)
  {
    const int index = fit->first->surface_index(fit->second);
    d->color_map[index]; // create the map entry
    // indices.insert(index);
  }
  // .. and same for cells
  for (C3t3::Cells_in_complex_iterator
         cit = c3t3().cells_in_complex_begin(),
         end = c3t3().cells_in_complex_end();
       cit !=  end; ++cit)
  {
    const int index = cit->subdomain_index();
    d->color_map[index]; // create the map entry
  }
  // std::cerr << "C3t3 indices:\n";
  // for(std::set<int>::const_iterator it = indices.begin();
  //     it != indices.end(); ++it) {
  //   std::cerr << *it << "\n";
  // }
  std::vector<QColor> colors;
  colors.reserve(d->color_map.size());
  compute_color_map(this->color(), d->color_map.size(),
                    std::back_inserter(colors));
  std::size_t i = 0;
  for(Scene_c3t3_item_priv::Color_map::iterator
        it = d->color_map.begin(),
        end = d->color_map.end();
      it != end; ++it)
  {
    it->second = colors[i++]; // start with colors[1], ignore the first color
  }

  // Rebuild histogram
  build_histogram();
}

void
Scene_c3t3_item::draw_points() const
{
  const Tr& tr = c3t3().triangulation();

  ::glBegin(GL_POINTS);
  for(Tr::Finite_vertices_iterator
        vit = tr.finite_vertices_begin(),
        end = tr.finite_vertices_end();
      vit != end; ++vit)
  {
    const Kernel::Point_3 p = vit->point().point();
    ::glVertex3d(p.x(),p.y(),p.z());
  }
  ::glEnd();

  if(this->draw_spheres) {
    d->draw_spheres();
  } // end if(draw_sphere)
}

void
Scene_c3t3_item::draw_triangles(bool with_cut_plane,
                                bool in_draw_edges) const
{
  const Kernel::Plane_3& plane = this->plane();

  unsigned int i = 2 * (with_cut_plane) + in_draw_edges;

  if(in_draw_edges || d->draw_triangles_display_list.begin_draw(i)) {
    ::glBegin(GL_TRIANGLES);
    for(C3t3::Facet_iterator
          fit = c3t3().facets_begin(),
          end = c3t3().facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const int color_index = fit->first->surface_index(fit->second);
      if(in_draw_edges)
        CGALglcolor(d->color_map[color_index].lighter(50));
      else
        CGALglcolor(d->color_map[color_index]);

      const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
      const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
      const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
      if(with_cut_plane) {
        typedef Kernel::Oriented_side Side;
        using CGAL::ON_ORIENTED_BOUNDARY;
        const Side sa = plane.oriented_side(pa);
        const Side sb = plane.oriented_side(pb);
        const Side sc = plane.oriented_side(pc);
        if(sa != ON_ORIENTED_BOUNDARY &&
           sb != ON_ORIENTED_BOUNDARY &&
           sc != ON_ORIENTED_BOUNDARY &&
           sb == sa && sc == sa )
        {
          draw_triangle(pa, pb, pc);
        }
      } else {
        draw_triangle(pa, pb, pc);
      }
    }
    ::glEnd();
  }
  if(!in_draw_edges) d->draw_triangles_display_list.end_draw(i);
}

void
Scene_c3t3_item::direct_draw_edges() const {
  draw_triangles(false, true);
}

void
Scene_c3t3_item::direct_draw() const {
}

void
Scene_c3t3_item::draw_edges() const
{
  double current_color[4];
  ::glGetDoublev(GL_CURRENT_COLOR, current_color);

  d->begin_OpenGL_clip_plane();
  Scene_item_with_display_list::draw_edges(); // draw edges with display list
  d->end_OpenGL_clip_plane();
  if(draw_tets) {
    draw_cut_tetrahedra();
  }

  ::glColor4dv(current_color);
  if(draw_spheres) {
    d->draw_spheres();
  } // end if(draw_sphere)
}

void
Scene_c3t3_item::draw_cut_tetrahedra() const
{
  if(frame->isManipulated()) return;
  double current_color[4];
  ::glGetDoublev(GL_CURRENT_COLOR, current_color);

  std::set<Tr::Cell_handle> cells;

  if(d->cut_display_list.begin_draw()) {
    const Kernel::Plane_3& plane = this->plane();

    ::glBegin(GL_TRIANGLES);
    for(Tr::Finite_cells_iterator
          cit = c3t3().triangulation().finite_cells_begin(),
          end = c3t3().triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if(! c3t3().is_in_complex(cit) )
        continue;

      const Kernel::Point_3& pa = cit->vertex(0)->point();
      const Kernel::Point_3& pb = cit->vertex(1)->point();
      const Kernel::Point_3& pc = cit->vertex(2)->point();
      const Kernel::Point_3& pd = cit->vertex(3)->point();
      typedef Kernel::Oriented_side Side;
      using CGAL::ON_ORIENTED_BOUNDARY;
      const Side sa = plane.oriented_side(pa);
      const Side sb = plane.oriented_side(pb);
      const Side sc = plane.oriented_side(pc);
      const Side sd = plane.oriented_side(pd);

      if( sa == ON_ORIENTED_BOUNDARY ||
          sb == ON_ORIENTED_BOUNDARY ||
          sc == ON_ORIENTED_BOUNDARY ||
          sd == ON_ORIENTED_BOUNDARY ||
          sb != sa || sc != sa || sd != sa)
      {
        if(d->draw_triangles_duals) {
          cells.insert(cit);
        }
        CGALglcolor(d->color_map[cit->subdomain_index()]);
        draw_triangle(pa, pb, pc);
        draw_triangle(pa, pb, pd);
        draw_triangle(pa, pc, pd);
        draw_triangle(pb, pc, pd);
      }
    }
    ::glEnd();
  }
  d->cut_display_list.end_draw();
  if(d->draw_triangles_duals) {
    const Tr& tr = c3t3().triangulation();
    if(d->duals_display_list.begin_draw())
    {
      ::glBegin(GL_LINES);
      BOOST_FOREACH(Tr::Cell_handle ch, cells) {
        for(int i = 0; i < 4; ++i) {
          Tr::Cell_handle n = ch->neighbor(i);
          const int j = n->index(ch);
          const Tr::Vertex_handle v = n->vertex(j);
          if(tr.is_infinite(v)) continue;
          if(cells.find(n) != cells.end()) continue;
          const Tr::Point p1 = tr.dual(ch);
          const Tr::Point p2 = tr.dual(n);
          ::glVertex3d(p1.x(), p1.y(), p1.z());
          ::glVertex3d(p2.x(), p2.y(), p2.z());
        }
      }
      ::glEnd();
    }
    d->duals_display_list.end_draw();
  }
  ::glColor4dv(current_color);
}

// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif

void
Scene_c3t3_item::draw() const
{
  if(isEmpty()) {
    return;
  }
  double current_color[4];
  ::glGetDoublev(GL_CURRENT_COLOR, current_color);

  GLboolean two_side;
  ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  d->begin_OpenGL_clip_plane();
  draw_triangles(true);
  d->end_OpenGL_clip_plane();

  if(draw_tets)
  {
    draw_cut_tetrahedra();
  }
  if(!two_side)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  if(draw_spheres) {
    d->draw_spheres();
  } // end if(draw_sphere)

  CGALglcolor(this->color());
  d->draw_cut_plane();

  ::glColor4dv(current_color);
}

#undef darker

void
Scene_c3t3_item::draw_sphere(const Tr::Point& p) const
{
  if(p.weight() > 0) {
    glPushMatrix();
    glTranslated(CGAL::to_double(p.point().x()),
                 CGAL::to_double(p.point().y()),
                 CGAL::to_double(p.point().z()));

    const GLdouble r =
      (spheres_drawn_radius != 0) ?
      spheres_drawn_radius :
      CGAL::to_double(CGAL_NTS sqrt(p.weight()));

    glScaled(r, r, r);
    glCallList(sphere_display_list);
    glPopMatrix();
  }
}




QPixmap
Scene_c3t3_item::graphicalToolTip() const
{
  if ( ! histogram_.isNull() )
  {
    return histogram_;
  }
  else
  {
    const_cast<Scene_c3t3_item&>(*this).build_histogram();
    return histogram_;
  }
}

template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value);

void
Scene_c3t3_item::build_histogram()
{
  // Create an histogram_ and display it
  const int height = 140;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height-top_margin*2;
  const int width = 402;
  const int cell_width = 2;
  const int text_margin = 3;
  const int text_height = 20;

  histogram_ = QPixmap(width,height+text_height);
  histogram_.fill(QColor(192,192,192));

  QPainter painter(&histogram_);
  painter.setPen(Qt::black);
  painter.setBrush(QColor(128,128,128));
  //painter.setFont(QFont("Arial", 30));

  // Build histogram_ data
  double min_value, max_value;
  std::vector<int> histo_data = create_histogram(c3t3(),min_value,max_value);

  // Get maximum value (to normalize)
  int max_size = 0;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it )
  {
    max_size = (std::max)(max_size,*it);
  }

  // colored histogram
  int j = 0;

  // draw
  int i=left_margin;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it, i+=cell_width )
  {
    int line_height = static_cast<int>( std::ceil(static_cast<double>(drawing_height) *
      static_cast<double>(*it)/static_cast<double>(max_size)) + .5);

    painter.fillRect(i,
                     drawing_height+top_margin-line_height,
                     cell_width,
                     line_height,
                     get_histogram_color(j++));
  }

  // draw bottom horizontal line
  painter.setPen(Qt::blue);

  painter.drawLine(QPoint(left_margin, drawing_height + top_margin),
                   QPoint(left_margin + static_cast<int>(histo_data.size())*cell_width,
                          drawing_height + top_margin));


  // draw min value and max value
  const int min_tr_width = static_cast<int>( 2*(std::floor(min_value)*cell_width + left_margin) );
  const int max_tr_width = static_cast<int>(
    2*((histo_data.size()-std::floor(max_value))*cell_width + left_margin) );
  const int tr_y = drawing_height + top_margin + text_margin;

  painter.setPen(get_histogram_color(min_value));
  QRect min_text_rect (0, tr_y, min_tr_width, text_height);
  painter.drawText(min_text_rect, Qt::AlignCenter, tr("%1").arg(min_value,0,'f',1));

  painter.setPen(get_histogram_color(max_value));
  QRect max_text_rect (width - max_tr_width, tr_y, max_tr_width, text_height);
  painter.drawText(max_text_rect, Qt::AlignCenter, tr("%1").arg(max_value,0,'f',1));
}


template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value)
{
  typedef typename C3t3::Triangulation::Point Point_3;

  std::vector<int> histo(181,0);

  min_value = 180.;
  max_value = 0.;

  for (typename C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit)
  {
    if( !c3t3.is_in_complex(cit))
      continue;

    const Point_3& p0 = cit->vertex(0)->point();
    const Point_3& p1 = cit->vertex(1)->point();
    const Point_3& p2 = cit->vertex(2)->point();
    const Point_3& p3 = cit->vertex(3)->point();

    double a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0,p1,p2,p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p2, p1, p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p3, p1, p2)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p2, p0, p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p3, p0, p2)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p2, p3, p0, p1)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

  }

  return histo;
}


QColor
Scene_c3t3_item::get_histogram_color(const double v) const
{
  if ( v < 5 )            { return Qt::red; }
  else if ( v < 10 )      { return QColor(215,108,0); }
  else if ( v < 15 )      { return QColor(138,139,0); }
  else if ( v < 165 )     { return Qt::darkGreen; }
  else if ( v < 170 )     { return QColor(138,139,1); }
  else if ( v < 175 )     { return QColor(215,108,0); }
  else /* 175<v<=180 */   { return Qt::red; }
}

void
Scene_c3t3_item::show_spheres(bool b)
{
  draw_spheres = b;
  emit itemChanged();
}


void
Scene_c3t3_item::show_tetrahedra(bool b)
{
  draw_tets = b;
  emit itemChanged();
}

void
Scene_c3t3_item::show_cut_plane(bool b)
{
  d->draw_clipping_plane = b;
  emit itemChanged();
}

void
Scene_c3t3_item::show_triangles_duals(bool b)
{
  d->draw_triangles_duals = b;
  emit itemChanged();
}


void
Scene_c3t3_item::draw_triangle(const Kernel::Point_3& pa,
                               const Kernel::Point_3& pb,
                               const Kernel::Point_3& pc)
{
  Kernel::Vector_3 n = cross_product(pb - pa, pc -pa);
  n = n / CGAL::sqrt(n*n);

  ::glNormal3d(n.x(),n.y(),n.z());

  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glVertex3d(pc.x(),pc.y(),pc.z());
}

double
Scene_c3t3_item::complex_diag() const
{
  const Bbox& bbox = this->bbox();
  const double xdelta = bbox.xmax-bbox.xmin;
  const double ydelta = bbox.ymax-bbox.ymin;
  const double zdelta = bbox.zmax-bbox.zmin;
  const double diag = std::sqrt(xdelta*xdelta +
                                ydelta*ydelta +
                                zdelta*zdelta);
  return diag * 0.7;
}

int
Scene_c3t3_item::number_of_vertices() const {
  return c3t3().triangulation().number_of_vertices();
}

void
Scene_c3t3_item::selectedPoint(double x,
                               double y,
                               double z)
{
  const Tr::Vertex_handle v =
    c3t3().triangulation().nearest_power_vertex(Kernel::Point_3(x, y, z));
  if(Tr::Vertex_handle() != v) {
    std::cerr << "c3t3_item: nearest vertex: " << v->point()
              << "   dim=" << v->in_dimension()
              << " index=" << v->index() << std::endl;
  }
}

#include "Scene_c3t3_item.moc"
