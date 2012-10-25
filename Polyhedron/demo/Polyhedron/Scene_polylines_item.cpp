#include "Scene_polylines_item.h"

#include <CGAL/bounding_box.h>
#include <CGAL/gl.h>
#include <CGAL/glu.h>
#include <QMenu>
#include <QAction>

#include <QInputDialog>

namespace {
  void CGALglcolor(QColor c, int dv = 0)
  {
    if ( 0 != dv )
    {
// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
      c = c.darker(dv);
#undef darker
    }
    ::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

class Scene_polylines_item_private {
public:
  typedef Scene_polylines_item::K K;
  typedef K::Point_3 Point_3;

  Scene_polylines_item_private() :
    draw_extremities(false),
    spheres_drawn_radius(0),
    sphere_display_list(0),
    quadric(0)
  {}

  ~Scene_polylines_item_private()
  {
    if(quadric != 0)
      gluDeleteQuadric(quadric);
    if(sphere_display_list  != 0)
      glDeleteLists(sphere_display_list, 1);
  }

  void draw_sphere(const K::Point_3&, double) const;
  void draw_spheres(const Scene_polylines_item*) const;

  bool draw_extremities;
  double spheres_drawn_radius;
private:
  mutable GLuint sphere_display_list;
  mutable GLUquadric* quadric;
};

Scene_polylines_item::Scene_polylines_item() 
  : d(new Scene_polylines_item_private())
{

}

Scene_polylines_item::~Scene_polylines_item()
{
  delete d;
}

bool
Scene_polylines_item::isEmpty() const {
  return polylines.empty();
}

Scene_interface::Bbox 
Scene_polylines_item::bbox() const {
  if(isEmpty())
    return Bbox();
  std::list<Point_3> boxes;
  for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
      it != polylines.end();
      ++it){
    if(it->begin() != it->end()) {
      Iso_cuboid_3 cub = CGAL::bounding_box(it->begin(), it->end());
      boxes.push_back((cub.min)());
      boxes.push_back((cub.max)());
    }
  }
  Iso_cuboid_3 bbox = 
    boxes.begin() != boxes.end() ?
    CGAL::bounding_box(boxes.begin(), boxes.end()) :
    Iso_cuboid_3();

  return Bbox(bbox.xmin(),
              bbox.ymin(),
              bbox.zmin(),
              bbox.xmax(),
              bbox.ymax(),
              bbox.zmax());
}

Scene_polylines_item* 
Scene_polylines_item::clone() const {
  Scene_polylines_item* item = new Scene_polylines_item;
  item->polylines = polylines;
  QVariant metadata_variant = property("polylines metadata");
  if(metadata_variant.type() == QVariant::StringList)
  {
    item->setProperty("polylines metadata", metadata_variant);
  }
  return item;
}

QString
Scene_polylines_item::toolTip() const {
  QString s = 
    tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
       "<i>Polylines</i></p>"
       "<p>Number of polylines: %4</p>")
    .arg(this->name())
    .arg(this->renderingModeName())
    .arg(this->color().name())
    .arg(polylines.size());
  if(d->draw_extremities) {
    s += tr("<p>Legende of endpoints colors: <ul>"
            "<li>black: one incident polyline</li>"
            "<li>green: two incident polylines</li>"
            "<li>blue: three incident polylines</li>"
            "<li>red: four incident polylines</li>"
            "<li>fuchsia: five or more incident polylines</li>"
            "</ul></p>");
  }
  return s;
}

bool
Scene_polylines_item::supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe || 
            m == FlatPlusEdges ||
            m == Points);
}

// Shaded OpenGL drawing: only draw spheres
void
Scene_polylines_item::draw() const {
  if(d->draw_extremities)
    d->draw_spheres(this);
}

// Wireframe OpenGL drawing
void 
Scene_polylines_item::draw_edges() const {
  CGALglcolor(this->color());
  ::glBegin(GL_LINES);
  for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
      it != polylines.end();
      ++it){
    if(it->empty()) continue;
    for(size_t i = 0, end = it->size()-1;
        i < end; ++i)
    {
      const Point_3& a = (*it)[i];
      const Point_3& b = (*it)[i+1];
      ::glVertex3d(a.x(), a.y(), a.z());
      ::glVertex3d(b.x(), b.y(), b.z());
    }
  }
  ::glEnd();
  if(d->draw_extremities)
  {
    d->draw_spheres(this);
  }
}

void 
Scene_polylines_item::draw_points() const {
  ::glBegin(GL_POINTS);
  // draw all points but endpoints
  for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
      it != polylines.end();
      ++it)
  {
    if(it->empty()) continue;
    for(size_t i = 1, end = it->size()-1;
        i < end; ++i)
    {
      const Point_3& a = (*it)[i];
      ::glVertex3d(a.x(), a.y(), a.z());
    }
  }
  ::glEnd();

  ::glColor3d(1., 0., 0.); //red
  // draw endpoints
  ::glBegin(GL_POINTS);
  for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
      it != polylines.end();
      ++it){
    if(it->empty()) continue;
    const Point_3& a = (*it)[0];
    const Point_3& b = (*it)[it->size()-1];
    ::glVertex3d(a.x(), a.y(), a.z());
    ::glVertex3d(b.x(), b.y(), b.z());
  }
  ::glEnd();
}

void
Scene_polylines_item_private::
draw_spheres(const Scene_polylines_item* item) const {
  // FIRST, count the number of incident cycles and polylines
  // for all extremities.
  typedef std::map<Point_3, int> Point_to_int_map;
  typedef Point_to_int_map::iterator iterator;
  Point_to_int_map corner_polyline_nb;

  { // scope to fill corner_polyline_nb'
    Point_to_int_map corner_cycles_nb;

    for(std::list<std::vector<Point_3> >::const_iterator 
          it = item->polylines.begin(),
          end = item->polylines.end();
        it != end; ++it)
    {
      const K::Point_3& a = *it->begin();
      const K::Point_3& b = *it->rbegin();
      if(a == b) {
        if ( it->size()>1 )
          ++corner_cycles_nb[a];
        else
          ++corner_polyline_nb[a];  
      }
      else {
        ++corner_polyline_nb[a];
        ++corner_polyline_nb[b];
      }
    }
    // THEN, ignore points that are incident to one cycle only.
    for(iterator 
          c_it = corner_cycles_nb.begin(),
          end = corner_cycles_nb.end();
        c_it != end; ++c_it)
    {
      const Point_3& a = c_it->first;

      iterator p_it = corner_polyline_nb.find(a);

      // If the point 'a'=c_it->first has only incident cycles...
      if(p_it == corner_polyline_nb.end()) {
        // ...then count it as a corner only if it has two incident cycles
        // or more.
        if(c_it->second > 1) {
          corner_polyline_nb[a] = c_it->second;
        }
      } else {
        // else add the number of cycles.
        p_it->second += c_it->second;
      }
    }
  }
  // At this point, 'corner_polyline_nb' gives the multiplicity of all
  // corners.
  for(iterator 
        p_it = corner_polyline_nb.begin(),
        end = corner_polyline_nb.end();
      p_it != end; ++p_it)
  {
    switch(p_it->second) {
    case 1: 
      ::glColor3d(0.0, 0.0, 0.0); // black
      break;
    case 2:
      ::glColor3d(0.0, 0.8, 0.0); // green
      break;
    case 3:
      ::glColor3d(0.0, 0.0, 0.8); // blue
      break;
    case 4:
      ::glColor3d(0.8, 0.0, 0.0); //red
      break;
    default:
      ::glColor3d(0.8, 0.0, 0.8); //fuschia
    }
    this->draw_sphere(p_it->first, this->spheres_drawn_radius);
  }
}

void 
Scene_polylines_item_private::draw_sphere(const K::Point_3& p,
                                          double r) const 
{
  if(sphere_display_list == 0) {
    sphere_display_list = glGenLists(1);
    if(sphere_display_list == 0)
      std::cerr << "ERROR: Cannot create display list!\n";
    if(quadric == 0)
      quadric = gluNewQuadric();
    if(quadric == 0)
      std::cerr << "ERROR: Cannot create GLU quadric!\n";
    glNewList(sphere_display_list, GL_COMPILE);
    gluSphere(quadric, 1., 10, 10);
    glEndList();
    if(glGetError() != GL_NO_ERROR)
      std::cerr << gluErrorString(glGetError());
  }
  glPushMatrix();
  glTranslated(CGAL::to_double(p.x()),
               CGAL::to_double(p.y()),
               CGAL::to_double(p.z()));

  glScaled(r, r, r);
  glCallList(sphere_display_list);
  glPopMatrix();
}

QMenu* Scene_polylines_item::contextMenu() 
{
  const char* prop_name = "Menu modified by Scene_polylines_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.trolltech.com/lastest/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {
    menu->addSeparator();
    // TODO: add actions to display corners
    QAction* action = menu->addAction(tr("Display corners with radius..."));
    connect(action, SIGNAL(triggered()),
            this, SLOT(change_corner_radii()));

    QAction* actionSmoothPolylines = 
      menu->addAction(tr("Smooth polylines"));
    actionSmoothPolylines->setObjectName("actionSmoothPolylines");
    connect(actionSmoothPolylines, SIGNAL(triggered()),this, SLOT(smooth()));
    menu->setProperty(prop_name, true);
  }
  return menu;
}

void Scene_polylines_item::change_corner_radii() {
  bool ok = true;
  double proposed_radius = d->spheres_drawn_radius;
  if(proposed_radius == 0) {
    Scene_interface::Bbox b = bbox();
    proposed_radius = (std::max)(b.xmax - b.xmin,
                                 proposed_radius);
    proposed_radius = (std::max)(b.ymax - b.ymin,
                                 proposed_radius);
    proposed_radius = (std::max)(b.zmax - b.zmin,
                                 proposed_radius);
    proposed_radius /= 100;
  }
  double r = QInputDialog::getDouble(NULL, 
                                     tr("Display corners with new radius..."),
                                     tr("Radius:"),
                                     proposed_radius, // value
                                     0.,          // min 
                                     2147483647., // max
                                     10,          // decimals
                                     &ok);
  if(ok) {
    change_corner_radii(r);
  }
}

void Scene_polylines_item::change_corner_radii(double r) {
  if(r >= 0) {
    d->spheres_drawn_radius = r;
    d->draw_extremities = (r > 0);
    this->changed();
    emit itemChanged();
  }
}

void Scene_polylines_item::split_at_sharp_angles()
{
  typedef Polylines_container Bare_polyline_container;
  typedef Polyline Bare_polyline;
  Polylines_container& bare_polylines = polylines;

  int counter = 0;
  for(Bare_polyline_container::iterator
        bare_polyline_it = bare_polylines.begin();
      bare_polyline_it != bare_polylines.end(); // the end changes
      // during the loop
      ++counter /* bare_polyline_it is incremented in the loop */)
  {
    Bare_polyline_container::iterator current_polyline_it = 
        bare_polyline_it;
    Bare_polyline& bare_polyline = *bare_polyline_it;
    Bare_polyline::iterator it = boost::next(bare_polyline.begin());

    if(boost::next(bare_polyline.begin()) == bare_polyline.end())
    {
      std::cerr << "WARNING: Isolated point in polylines\n";
      bare_polyline_it = bare_polylines.erase(bare_polyline_it);
      continue;
    }
    else 
      ++bare_polyline_it;
    if(it != bare_polyline.end()) {
      for(; it != boost::prior(bare_polyline.end()); ++it) {
        const Point_3 pv = *it;
        const Point_3 pa = *boost::prior(it);
        const Point_3 pb = *boost::next(it);
        const K::Vector_3 av = pv - pa;
        const K::Vector_3 bv = pv - pb;
        const K::FT sc_prod = av * bv;
        if( sc_prod >= 0 ||
            (sc_prod < 0 && 
             CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) ) 
        {
#ifdef PROTECTION_DEBUG
          std::cerr << "Split polyline (small angle) "
                    <<  std::acos(sqrt(CGAL::square(sc_prod) /
                                       ((av*av) * (bv*bv)))) * 180 /CGAL_PI
                    << " degres\n";
#endif
          Bare_polyline new_polyline;
          std::copy(it, bare_polyline.end(), 
                    std::back_inserter(new_polyline));
          
          if(*bare_polyline.begin() == *bare_polyline.rbegin()) {
            // if the polyline is a cycle, test if its beginning is a sharp
            // angle...
            const Point_3 pv = *bare_polyline.begin();
            const Point_3 pa = *boost::prior(boost::prior(bare_polyline.end()));
            const Point_3 pb = *boost::next(bare_polyline.begin());
            const K::Vector_3 av = pv - pa;
            const K::Vector_3 bv = pv - pb;
            const K::FT sc_prod = av * bv;
            if( sc_prod >= 0 ||
                (sc_prod < 0 && 
                 CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) ) 
            {
              // if its beginning is a sharp angle, then split
              bare_polyline.erase(boost::next(it), bare_polyline.end());
            }
            else {
              // ...if not, modifies its beginning
              std::copy(boost::next(bare_polyline.begin()),
                        boost::next(it),
                        std::back_inserter(new_polyline));
              bare_polylines.erase(current_polyline_it);
            }
          }
          else {
            bare_polyline.erase(boost::next(it), bare_polyline.end());
          }
          bare_polylines.push_back(new_polyline);
          break;
        }
      }
    }
  }
  emit itemChanged();
}

void
Scene_polylines_item::merge(Scene_polylines_item* other_item) {
  if(other_item == 0) return;
  std::copy(other_item->polylines.begin(),
            other_item->polylines.end(),
            std::back_inserter(polylines));
  QVariant other_metadata_variant = other_item->property("polylines metadata");
  if(other_metadata_variant.type() == QVariant::StringList)
  {
    QStringList metadata = property("polylines metadata").toStringList();
    metadata.append(other_metadata_variant.toStringList());
    setProperty("polylines metadata", metadata);
  }
  changed();
}

#include "Scene_polylines_item.moc"
