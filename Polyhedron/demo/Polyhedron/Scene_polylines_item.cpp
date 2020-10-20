#include "Scene_polylines_item.h"
#include "Scene_spheres_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/bounding_box.h>
#include <CGAL/Three/Three.h>
#include <QMenu>
#include <QSlider>
#include <QWidgetAction>
#include <QAction>
#include <QInputDialog>
#include <QApplication>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

using namespace CGAL::Three;
typedef Edge_container Ec;
typedef Point_container Pc;
typedef Viewer_interface Vi;

struct Scene_polylines_item_private {
    typedef Scene_polylines_item::K K;
    typedef K::Point_3 Point_3;

    Scene_polylines_item_private(Scene_polylines_item *parent) :
        draw_extremities(false),
        spheres_drawn_square_radius(0)
    {
      line_Slider = new QSlider(Qt::Horizontal);
      line_Slider->setValue(CGAL::Three::Three::getDefaultLinesWidth());
      line_Slider->setMaximum(2000);
      line_Slider->setMinimum(1);
      item = parent;
      invalidate_stats();
    }

    void invalidate_stats()
    {
      nb_vertices = 0;
      nb_edges = 0;
      min_length = (std::numeric_limits<double>::max)();
      max_length = 0;
      mean_length = 0;
      computed_stats = false;
    }

    mutable Scene_spheres_item *spheres;
    mutable std::vector<float> positions_lines;
    mutable std::size_t nb_lines;
    typedef std::map<Point_3, int> Point_to_int_map;
    typedef Point_to_int_map::iterator iterator;
    void computeSpheres();
    void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
    void computeElements() const;
    bool draw_extremities;
    double spheres_drawn_square_radius;
    Scene_polylines_item *item;
    mutable std::size_t nb_vertices;
    mutable std::size_t nb_edges;
    mutable double min_length;
    mutable double max_length;
    mutable double mean_length;
    mutable bool computed_stats;
    QSlider* line_Slider;
};


void
Scene_polylines_item_private::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  float lineWidth[2];
  if(!viewer->isOpenGL_4_3())
    viewer->glGetFloatv(GL_LINE_WIDTH_RANGE, lineWidth);
  else
  {
    lineWidth[0] = 0;
    lineWidth[1] = 10;
  }
  line_Slider->setMaximum(lineWidth[1]);

  item->getEdgeContainer(0)->initializeBuffers(viewer);
  item->getEdgeContainer(0)->setFlatDataSize(nb_lines);
  item->getPointContainer(0)->initializeBuffers(viewer);
  item->getPointContainer(0)->setFlatDataSize(nb_lines);

  positions_lines.clear();
  positions_lines.shrink_to_fit();
}

void
Scene_polylines_item_private::computeElements() const
{
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first())->offset();
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_lines.resize(0);
    double mean = 0;
    bool all_equal=true;
    //Fills the VBO with the lines
    for(std::list<std::vector<Point_3> >::const_iterator it = item->polylines.begin();
        it != item->polylines.end();
        ++it)
    {
      if(it->empty()) continue;
      if(it->front() == it->back())
        nb_vertices += it->size() - 1;
      else
        nb_vertices += it->size();

      for(size_t i = 0, end = it->size()-1;
          i < end; ++i)
      {
            const Point_3& a = (*it)[i];
            const Point_3& b = (*it)[i+1];
            if(a!=b)
              all_equal = false;
            if(!computed_stats)
            {
              ++nb_edges;
                double length = CGAL::sqrt(
                      (a.x()-b.x()) * (a.x()-b.x()) +
                      (a.y()-b.y()) * (a.y()-b.y()) +
                      (a.z()-b.z()) * (a.z()-b.z()) );
                if(max_length < length)
                  max_length = length;
                if(min_length > length)
                  min_length = length;
                mean += length;
            }

            positions_lines.push_back(a.x()+offset.x);
            positions_lines.push_back(a.y()+offset.y);
            positions_lines.push_back(a.z()+offset.z);

            positions_lines.push_back(b.x()+offset.x);
            positions_lines.push_back(b.y()+offset.y);
            positions_lines.push_back(b.z()+offset.z);
        }

    }
    if(all_equal)
      item->setPointsMode();
    if(!computed_stats)
      mean_length = mean/double(nb_edges);
    computed_stats = true;

    item->getEdgeContainer(0)->allocate(
          Ec::Vertices, positions_lines.data(),
          static_cast<int>(positions_lines.size()*sizeof(float)));
    item->getPointContainer(0)->allocate(
          Pc::Vertices, positions_lines.data(),
          static_cast<int>(positions_lines.size()*sizeof(float)));

    item->setBuffersFilled(true);
    nb_lines = positions_lines.size();
    QApplication::restoreOverrideCursor();
}

void
Scene_polylines_item_private::computeSpheres()
{
  const CGAL::qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  K::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);

      spheres->clear_spheres();
      QApplication::setOverrideCursor(Qt::WaitCursor);
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
      //Finds the centers of the spheres and their color
      int s_id = 0;
      for(iterator
          p_it = corner_polyline_nb.begin(),
          end = corner_polyline_nb.end();
          p_it != end; ++p_it)
      {
          const K::Point_3& center = p_it->first;
          unsigned char colors[3];
          switch(p_it->second) {
          case 1:
              colors[0] = 0; // black
              colors[1] = 0;
              colors[2] = 0;
              break;
          case 2:
              colors[0] = 0; // green
              colors[1] = 200;
              colors[2] = 0;
              break;
          case 3:
              colors[0] = 0; // blue
              colors[1] = 0;
              colors[2] = 200;
              break;
          case 4:
              colors[0] = 200; //red
              colors[1] = 0;
              colors[2] = 0;
              break;
          default:
              colors[0] = 200; //fuschia
              colors[1] = 0;
              colors[2] = 200;
          }

          CGAL::Color c(colors[0], colors[1], colors[2]);
          spheres->add_sphere(K::Sphere_3(center+offset, spheres_drawn_square_radius),s_id++, c);
      }
      spheres->setToolTip(
            QString("<p>Legend of endpoints colors: <ul>"
                    "<li>black: one incident polyline</li>"
                    "<li>green: two incident polylines</li>"
                    "<li>blue: three incident polylines</li>"
                    "<li>red: four incident polylines</li>"
                    "<li>fuchsia: five or more incident polylines</li>"
                    "</ul></p>"
                    "<p>Tip: To erase this item, set its radius to 0 or less, or select it in the itme list and press the Delete key. </p>")
                            );
      spheres->computeElements();
      QApplication::restoreOverrideCursor();
}

Scene_polylines_item::Scene_polylines_item()
    :CGAL::Three::Scene_group_item("unnamed")
    ,d(new Scene_polylines_item_private(this))
{
    setRenderingMode(FlatPlusEdges);
    d->nb_lines = 0;
    d->spheres = NULL;
    setEdgeContainer(0,
                     new Ec(Three::mainViewer()->isOpenGL_4_3() ? Vi::PROGRAM_SOLID_WIREFRAME
                                                                : Vi::PROGRAM_NO_SELECTION
                                                                  , false));
    setPointContainer(0,
                      new Pc(Vi::PROGRAM_NO_SELECTION, false));

}

Scene_polylines_item::~Scene_polylines_item()
{
  delete d->line_Slider;
  delete d;

}

bool
Scene_polylines_item::isEmpty() const {
    return polylines.empty();
}

void
Scene_polylines_item::compute_bbox() const {
    typedef K::Iso_cuboid_3 Iso_cuboid_3;

    if(isEmpty())
    {
        setBbox(Bbox());
        return;
    }
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

    setBbox(Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax()));
}

Scene_item::Bbox Scene_polylines_item::bbox() const
{
  return Scene_item_rendering_helper::bbox();
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
    if(polylines.size() == 1 )
    {
      double length = 0;
      for(std::size_t i=1; i<polylines.front().size(); ++i)
      {
        K::Vector_3 vec(polylines.front()[i-1], polylines.front()[i]);
        length += CGAL::sqrt(vec.squared_length());
      }
      QString vertices_info = tr("<p>Number of vertices: %1<br />"
                                 "Polyline's length: %2</p>")
          .arg(polylines.front().size())
          .arg(length);
     s.append(vertices_info);

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
Scene_polylines_item::draw(CGAL::Three::Viewer_interface* viewer) const {
  if(!visible())
    return;
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  if(d->draw_extremities)
  {
    Scene_group_item::draw(viewer);
  }
}

// Wireframe OpenGL drawing
void
Scene_polylines_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
  if(!visible())
    return;

  if(renderingMode() == Wireframe
     || renderingMode() == FlatPlusEdges)
  {
    if(!isInit(viewer))
      initGL(viewer);
    if ( getBuffersFilled() &&
         ! getBuffersInit(viewer))
    {
      initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }
    if(!getBuffersFilled())
    {
      computeElements();
      initializeBuffers(viewer);
    }
    if(viewer->isOpenGL_4_3())
    {
      QVector2D vp(viewer->width(), viewer->height());

      getEdgeContainer(0)->setViewport(vp);
      getEdgeContainer(0)->setWidth(GLfloat(d->line_Slider->value()));
    }
    getEdgeContainer(0)->setColor(this->color());
    getEdgeContainer(0)->draw(viewer, true);
  }
    if(d->draw_extremities)
    {
       Scene_group_item::drawEdges(viewer);
    }
}

void
Scene_polylines_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const {
  if(!visible())
    return;
  if(renderingMode() == Points)
  {
    if(!isInit(viewer))
      initGL(viewer);
    if ( getBuffersFilled() &&
         ! getBuffersInit(viewer))
    {
      initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }
    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->setGlPointSize(GLfloat(5));
    if(!getBuffersFilled())
    {
      computeElements();
      initializeBuffers(viewer);
    }

    getPointContainer(0)->setColor(this->color());
    getPointContainer(0)->draw(viewer, true);
    viewer->setGlPointSize(point_size);
  }
   if(d->draw_extremities)
   {
      Scene_group_item::drawPoints(viewer);
   }
}

QMenu* Scene_polylines_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_polylines_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // https://doc.qt.io/qt-5/qobject.html#property
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

        QMenu *container = new QMenu(tr("Line Width"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        connect(d->line_Slider, &QSlider::valueChanged, this, &Scene_polylines_item::itemChanged);

        sliderAction->setDefaultWidget(d->line_Slider);

        container->addAction(sliderAction);
        menu->addMenu(container);

        menu->setProperty(prop_name, true);
        container->menuAction()->setProperty("is_groupable", true);
        QAction* actionCreatePointSet =
            menu->addAction(tr("Create Point Set from Polyline"));
        actionCreatePointSet->setObjectName("actionCreatePointSet");
        connect(actionCreatePointSet, &QAction::triggered,
                this, &Scene_polylines_item::point_set_from_polyline);
    }
    return menu;
}

void Scene_polylines_item::invalidateOpenGLBuffers()
{
    setBuffersFilled(false);
    getEdgeContainer(0)->reset_vbos(ALL);
    getPointContainer(0)->reset_vbos(ALL);
    d->invalidate_stats();
    compute_bbox();


}

void Scene_polylines_item::change_corner_radii() {
    bool ok = true;
    double proposed_radius = std::sqrt(d->spheres_drawn_square_radius);
    if(proposed_radius == 0) {
        CGAL::Three::Scene_interface::Bbox b = bbox();
        proposed_radius = (std::max)(b.xmax() - b.xmin(),
                                     proposed_radius);
        proposed_radius = (std::max)(b.ymax() - b.ymin(),
                                     proposed_radius);
        proposed_radius = (std::max)(b.zmax() - b.zmin(),
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
        d->spheres_drawn_square_radius = r*r;
        d->draw_extremities = (r > 0);
        if(r>0 && !d->spheres)
        {
          d->spheres = new Scene_spheres_item(this, 0, false);
          d->spheres->setName("Corner spheres");
          d->spheres->setRenderingMode(Gouraud);
          connect(d->spheres, SIGNAL(destroyed()), this, SLOT(reset_spheres()));
          scene->addItem(d->spheres);
          scene->changeGroup(d->spheres, this);
          lockChild(d->spheres);
          Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
          {
            d->spheres->gl_initialization(qobject_cast<Viewer_interface*>(v));
          }
          d->spheres->invalidateOpenGLBuffers();
          d->computeSpheres();
          CGAL::Three::Three::mainWindow()->installEventFilter(d->spheres);
          connect(d->spheres, &Scene_spheres_item::destroyMe, this,
                  [this](){
            removeEventFilter(d->spheres);
            unlockChild(d->spheres);
            scene->erase(scene->item_id(d->spheres));
          });
        }
        else if(r>0 && d->spheres)
        {
          d->spheres->invalidateOpenGLBuffers();
          d->computeSpheres();
        }
        else if (r<=0 && d->spheres!=NULL)
        {
          unlockChild(d->spheres);
          scene->erase(scene->item_id(d->spheres));
        }
    Q_EMIT itemChanged();
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
  Q_EMIT itemChanged();
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
    invalidateOpenGLBuffers();
}

void Scene_polylines_item::reset_spheres()
{
  d->spheres = NULL;
}

void Scene_polylines_item::smooth(){
    for (Polylines_container::iterator pit=polylines.begin(),pit_end=polylines.end();pit!=pit_end;++pit)
        smooth(*pit);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

void Scene_polylines_item::smooth(std::vector<Point_3>& polyline){
    bool is_closed = polyline.front()==polyline.back();
    typedef K::Vector_3 Vector_3;

    std::size_t start = is_closed ? 0:1;
    std::size_t end   = polyline.size()-1;

    Vector_3 prev = (is_closed ? polyline[end-1] : polyline[0]) - CGAL::ORIGIN;

    for (std::size_t i=start; i!=end; ++i)
    {
        Vector_3 curr = polyline[i] - CGAL::ORIGIN;
        Vector_3 next = polyline[i+1] - CGAL::ORIGIN;

        polyline[i] = CGAL::ORIGIN+(prev+2*curr+next)/4;
        prev=curr;
    }

    if (is_closed) polyline[end]=polyline[0];
}

QString Scene_polylines_item::computeStats(int type)
{
  switch (type)
  {
  case NB_VERTICES:
    return QString::number(d->nb_vertices);
  case NB_EDGES:
    return QString::number(d->nb_edges);
  case MIN_LENGTH:
    return QString::number(d->min_length);
  case MAX_LENGTH:
    return QString::number(d->max_length);
  case MEAN_LENGTH:
    return QString::number(d->mean_length);
  default:
    return QString();
  }
}
CGAL::Three::Scene_item::Header_data Scene_polylines_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories
  data.categories.append(std::pair<QString,int>(QString("Properties"),5));


  //titles
  data.titles.append(QString("#Vertices"));
  data.titles.append(QString("#Segment Edges"));
  data.titles.append(QString("Shortest Segment Edge Length"));
  data.titles.append(QString("Longest Segment Edge Length"));
  data.titles.append(QString("Average Segment Edge Length"));
  return data;
}

void Scene_polylines_item::setWidth( int i)
{
  d->line_Slider->setValue(i);
  redraw();
}

void Scene_polylines_item::initializeBuffers(Viewer_interface *v) const
{
  d->initializeBuffers(v);
}

void Scene_polylines_item::computeElements() const
{
  d->computeElements();
}

void Scene_polylines_item::point_set_from_polyline()
{
  Scene_points_with_normal_item* new_ps_item
      = new Scene_points_with_normal_item();
  for(auto it = polylines.begin();
      it != polylines.end();
      ++it)
  {
    std::size_t length = it->size();
    for(std::size_t id = 0; id < length; ++id)
    {
      if(id != length - 1
         || it->at(id) != it->front())
        new_ps_item->point_set()->insert(it->at(id));
    }
  }
  new_ps_item->setName(tr("Points from %1").arg(name()));
  CGAL::Three::Three::scene()->addItem(new_ps_item);
  new_ps_item->invalidateOpenGLBuffers();
  new_ps_item->redraw();
}
