#include "Scene_polyhedron_shortest_path_item.h"
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Point_container.h>

#include "Scene_polylines_item.h"

#include <vector>
#include <Qt>
#include <QApplication>
#include <QKeySequence>
#include <fstream>

#include <CGAL/Surface_mesh_shortest_path/function_objects.h>
#include <CGAL/Three/Three.h>
#include <QString>
using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Point_container Pc;

Viewer_interface* (&getActiveViewer)() = Three::activeViewer;
typedef Scene_polyhedron_shortest_path_item It;
struct Scene_polyhedron_shortest_path_item_priv
{
  typedef CGAL::Three::Scene_interface::Bbox Bbox;

  typedef boost::property_map<Face_graph, CGAL::vertex_point_t>::type VertexPointMap;

  typedef boost::graph_traits<Face_graph> GraphTraits;
  typedef GraphTraits::face_descriptor face_descriptor;
  typedef GraphTraits::face_iterator face_iterator;

  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Face_graph>
  Surface_mesh_shortest_path_traits;
  typedef CGAL::Surface_mesh_shortest_path<Surface_mesh_shortest_path_traits>
  Surface_mesh_shortest_path;
  typedef Surface_mesh_shortest_path::Face_location Face_location;
  typedef CGAL::AABB_face_graph_triangle_primitive<Face_graph, VertexPointMap>
  AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Kernel, AABB_face_graph_primitive> AABB_face_graph_traits;
  typedef CGAL::AABB_tree<AABB_face_graph_traits> AABB_face_graph_tree;

  typedef Surface_mesh_shortest_path_traits::Barycentric_coordinates
  Barycentric_coordinates;
  typedef Surface_mesh_shortest_path_traits::Construct_barycentric_coordinates
  Construct_barycentric_coordinates;
  typedef Surface_mesh_shortest_path_traits::Ray_3 Ray_3;
  typedef Surface_mesh_shortest_path_traits::Point_3 Point_3;
  typedef Surface_mesh_shortest_path_traits::FT FT;
  Scene_polyhedron_shortest_path_item_priv(Scene_polyhedron_shortest_path_item *parent)
    : m_shortestPaths(NULL)
    , m_isTreeCached(false)
    , m_shiftHeld(false)
  {
    item = parent;
  }

  bool get_mouse_ray(QMouseEvent* mouseEvent, Kernel::Ray_3&);
  void recreate_shortest_path_object();
  void ensure_aabb_object();
  void ensure_shortest_paths_tree();

  bool run_point_select(const Kernel::Ray_3&);
  void remove_nearest_point(const Scene_polyhedron_shortest_path_item::Face_location& ray);
  void get_as_edge_point(Scene_polyhedron_shortest_path_item::Face_location& inOutLocation);
  void get_as_vertex_point(Scene_polyhedron_shortest_path_item::Face_location& inOutLocation);
  void compute_elements(void) const;
  void deinitialize()
  {
    if (m_shortestPaths)
    {
      delete m_shortestPaths;
      m_sceneInterface = NULL;
    }
  }

  Scene_polyhedron_shortest_path_item* item;
  Messages_interface* m_messages;
  QMainWindow* m_mainWindow;
  CGAL::Three::Scene_interface* m_sceneInterface;
  Scene_polyhedron_shortest_path_item::Surface_mesh_shortest_path* m_shortestPaths;
  Scene_polyhedron_shortest_path_item::AABB_face_graph_tree m_aabbTree;
  std::string m_deferredLoadFilename;
  Scene_polyhedron_shortest_path_item::Selection_mode m_selectionMode;
  Scene_polyhedron_shortest_path_item::Primitives_mode m_primitivesMode;
  bool m_isTreeCached;
  bool m_shiftHeld;

  mutable std::vector<float> vertices;
  mutable std::size_t nb_vertices;
};

void Scene_polyhedron_shortest_path_item::common_constructor()
{
  setPointContainer(0, new Pc(Vi::PROGRAM_NO_SELECTION, false));
}

Scene_polyhedron_shortest_path_item::Scene_polyhedron_shortest_path_item()
   :Scene_polyhedron_item_decorator(NULL, false)
{
  d = new Scene_polyhedron_shortest_path_item_priv(this);
  common_constructor();
}

Scene_polyhedron_shortest_path_item::Scene_polyhedron_shortest_path_item(
    Scene_face_graph_item* polyhedronItem,
    CGAL::Three::Scene_interface* sceneInterface,
    Messages_interface* messages,
    QMainWindow* mainWindow)
  :Scene_polyhedron_item_decorator(polyhedronItem, false)
{ d = new Scene_polyhedron_shortest_path_item_priv(this);
  common_constructor();
  initialize(polyhedronItem, sceneInterface, messages, mainWindow);
}

Scene_polyhedron_shortest_path_item::~Scene_polyhedron_shortest_path_item()
{
  deinitialize();
  delete d;
}

void Scene_polyhedron_shortest_path_item_priv::compute_elements() const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    vertices.resize(0);

     const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
     typedef Scene_polyhedron_shortest_path_item::Surface_mesh_shortest_path::
         Source_point_iterator Source_point_iterator;
    for( Source_point_iterator it = m_shortestPaths->source_points_begin();
         it != m_shortestPaths->source_points_end();
         ++it)
    {
      const Kernel::Point_3& p = m_shortestPaths->point(it->first, it->second);
      vertices.push_back(p.x() + offset.x);
      vertices.push_back(p.y() + offset.y);
      vertices.push_back(p.z() + offset.z);
    }
    QApplication::restoreOverrideCursor();
}

bool Scene_polyhedron_shortest_path_item::supportsRenderingMode(RenderingMode m) const
{
  switch (m)
  {
  case Points:
    return true;
  case PointsPlusNormals:
    return true;
  case Wireframe:
    return true;
  case Flat:
    return true;
  case FlatPlusEdges:
    return true;
  case Gouraud:
    return true;
  default:
    return true;
  }
}

void Scene_polyhedron_shortest_path_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
    if (supportsRenderingMode(renderingMode()))
    {
      drawPoints(viewer);
    }
}


void Scene_polyhedron_shortest_path_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
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

   viewer->setGlPointSize(4.0f);
   Pc* pc = getPointContainer(0);
   pc->setColor(QColor(Qt::green));
   pc->draw(viewer, true);
   viewer->setGlPointSize(1.0f);
}

Scene_polyhedron_shortest_path_item* Scene_polyhedron_shortest_path_item::clone() const
{
  return 0;
}

void Scene_polyhedron_shortest_path_item::set_selection_mode(Selection_mode mode)
{
  d->m_selectionMode = mode;
}

Scene_polyhedron_shortest_path_item::Selection_mode
Scene_polyhedron_shortest_path_item::get_selection_mode() const
{
  return d->m_selectionMode;
}

void Scene_polyhedron_shortest_path_item::set_primitives_mode(Primitives_mode mode)
{
  d->m_primitivesMode = mode;
}

Scene_polyhedron_shortest_path_item::Primitives_mode
Scene_polyhedron_shortest_path_item::get_primitives_mode() const
{
  return d->m_primitivesMode;
}

void Scene_polyhedron_shortest_path_item_priv::recreate_shortest_path_object()
{
  if (m_shortestPaths)
  {
    delete m_shortestPaths;
  }

  m_shortestPaths = new Scene_polyhedron_shortest_path_item::Surface_mesh_shortest_path(
        *(item->polyhedron()),
            CGAL::get(boost::vertex_index, *(item->polyhedron())),
            CGAL::get(CGAL::halfedge_index, *(item->polyhedron())),
            CGAL::get(CGAL::face_index, *(item->polyhedron())),
            CGAL::get(CGAL::vertex_point, *(item->polyhedron())));

  //m_shortestPaths->m_debugOutput = true;

  m_isTreeCached = false;
}

void Scene_polyhedron_shortest_path_item_priv::ensure_aabb_object()
{
  if (!m_isTreeCached)
  {
    m_shortestPaths->build_aabb_tree(m_aabbTree);
    m_isTreeCached = true;
  }
}

void Scene_polyhedron_shortest_path_item_priv::ensure_shortest_paths_tree()
{
  if (!m_shortestPaths->changed_since_last_build())
  {
    CGAL::Three::Three::information("Recomputing shortest paths tree...");
    m_shortestPaths->build_sequence_tree();
    CGAL::Three::Three::information("Done.");
  }
}

void Scene_polyhedron_shortest_path_item::poly_item_changed()
{
  d->recreate_shortest_path_object();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

void Scene_polyhedron_shortest_path_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  setBuffersFilled(false);
  getPointContainer(0)->reset_vbos(ALL);
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
  }
}

bool Scene_polyhedron_shortest_path_item_priv::get_mouse_ray(
    QMouseEvent* mouseEvent, Kernel::Ray_3& outRay)
{
  CGAL::QGLViewer* viewer = getActiveViewer();
  viewer->makeCurrent();
  const CGAL::qglviewer::Vec offset = viewer->offset();
  CGAL::qglviewer::Camera* camera = viewer->camera();
  bool found = false;
  CGAL::qglviewer::Vec point = camera->pointUnderPixel(mouseEvent->pos(), found) - offset;
  if(found)
  {
    const CGAL::qglviewer::Vec orig = camera->position() - offset;
    outRay = Ray_3(Point_3(orig.x, orig.y, orig.z), Point_3(point.x, point.y, point.z));
  }

  return found;
}

void Scene_polyhedron_shortest_path_item_priv::remove_nearest_point(
    const Scene_polyhedron_shortest_path_item::Face_location& faceLocation)
{
  Surface_mesh_shortest_path_traits::Compute_squared_distance_3 computeSquaredDistance3;

  const Point_3 pickLocation = m_shortestPaths->point(faceLocation.first,
                                                      faceLocation.second);

  Surface_mesh_shortest_path::Source_point_iterator found =
      m_shortestPaths->source_points_end();
  FT minDistance(0.0);

  for (Surface_mesh_shortest_path::Source_point_iterator it =
       m_shortestPaths->source_points_begin();
       it != m_shortestPaths->source_points_end(); ++it)
  {
    Point_3 sourceLocation = m_shortestPaths->point(it->first, it->second);
    FT distance = computeSquaredDistance3(sourceLocation, pickLocation);

    if (found == m_shortestPaths->source_points_end() || distance < minDistance)
    {
      found = it;
      minDistance = distance;
    }
  }

  if (found != m_shortestPaths->source_points_end())
  {
    m_shortestPaths->remove_source_point(found);
  }
}

void Scene_polyhedron_shortest_path_item_priv::get_as_edge_point(
    Scene_polyhedron_shortest_path_item::Face_location& inOutLocation)
{
  size_t minIndex = 0;
  FT minCoord(inOutLocation.second[0]);

  for (size_t i = 1; i < 3; ++i)
  {
    if (minCoord > inOutLocation.second[i])
    {
      minIndex = i;
      minCoord = inOutLocation.second[i];
    }
  }

  // The nearest edge is that of the two non-minimal barycentric coordinates
  size_t nearestEdge[2];
  size_t current = 0;

  for (size_t i = 0; i < 3; ++i)
  {
    if (i != minIndex)
    {
      nearestEdge[current] = i;
      ++current;
    }
  }

  Construct_barycentric_coordinates construct_barycentric_coordinates;

  Point_3 trianglePoints[3] = {
    m_shortestPaths->point(inOutLocation.first,
    construct_barycentric_coordinates(FT(1.0), FT(0.0), FT(0.0))),
    m_shortestPaths->point(inOutLocation.first,
    construct_barycentric_coordinates(FT(0.0), FT(1.0), FT(0.0))),
    m_shortestPaths->point(inOutLocation.first,
    construct_barycentric_coordinates(FT(0.0), FT(0.0), FT(1.0))),
  };

  CGAL::Surface_mesh_shortest_paths_3::Parametric_distance_along_segment_3<
      Surface_mesh_shortest_path_traits> parametricDistanceSegment3;

  Point_3 trianglePoint = m_shortestPaths->point(inOutLocation.first, inOutLocation.second);

  FT distanceAlongSegment = parametricDistanceSegment3(trianglePoints[nearestEdge[0]],
      trianglePoints[nearestEdge[1]], trianglePoint);

  FT coords[3] = { FT(0.0), FT(0.0), FT(0.0), };

  coords[nearestEdge[1]] = distanceAlongSegment;
  coords[nearestEdge[0]] = FT(1.0) - distanceAlongSegment;

  inOutLocation.second = construct_barycentric_coordinates(coords[0], coords[1], coords[2]);
}

void Scene_polyhedron_shortest_path_item_priv::get_as_vertex_point(
    Scene_polyhedron_shortest_path_item::Face_location& inOutLocation)
{
  size_t maxIndex = 0;
  FT maxCoord(inOutLocation.second[0]);

  for (size_t i = 1; i < 3; ++i)
  {
    if (inOutLocation.second[i] > maxCoord)
    {
      maxIndex = i;
      maxCoord = inOutLocation.second[i];
    }
  }

  FT coords[3] = { FT(0.0), FT(0.0), FT(0.0), };
  coords[maxIndex] = FT(1.0);

  Construct_barycentric_coordinates construct_barycentric_coordinates;
  inOutLocation.second = construct_barycentric_coordinates(coords[0], coords[1], coords[2]);
}

bool Scene_polyhedron_shortest_path_item_priv::run_point_select(const Ray_3& ray)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  ensure_aabb_object();

  Face_location faceLocation = m_shortestPaths->locate(ray, m_aabbTree);

  if (faceLocation.first == GraphTraits::null_face())
  {
    CGAL::Three::Three::information(QObject::tr("Shortest Paths: No face under cursor."));
    QApplication::restoreOverrideCursor();
    return false;
  }
  else
  {
    boost::property_map<Face_graph, CGAL::face_index_t>::type fimap
        = get(CGAL::face_index, *item->polyhedron());

    CGAL::Three::Three::information(QObject::tr("Shortest Paths: Selected Face: %1; Barycentric coordinates: %2 %3 %4")
                            .arg(get(fimap, faceLocation.first))
                            .arg(double(faceLocation.second[0]))
        .arg(double(faceLocation.second[1]))
        .arg(double(faceLocation.second[2])));
    switch (m_selectionMode)
    {
    case It::INSERT_POINTS_MODE:
      switch (m_primitivesMode)
      {
      case It::VERTEX_MODE:
        get_as_vertex_point(faceLocation);
        m_shortestPaths->add_source_point(faceLocation.first, faceLocation.second);
        break;
      case It::EDGE_MODE:
        get_as_edge_point(faceLocation);
        m_shortestPaths->add_source_point(faceLocation.first, faceLocation.second);
        break;
      case It::FACE_MODE:
        m_shortestPaths->add_source_point(faceLocation.first, faceLocation.second);
        break;
      }
      break;
    case It::REMOVE_POINTS_MODE:
      remove_nearest_point(faceLocation);
      break;
    case It::SHORTEST_PATH_MODE:
      switch (m_primitivesMode)
      {
      case It::VERTEX_MODE:
        get_as_vertex_point(faceLocation);
        break;
      case It::EDGE_MODE:
        get_as_edge_point(faceLocation);
        break;
      case It::FACE_MODE:
        break;
      }

      if (m_shortestPaths->number_of_source_points() > 0)
      {
        ensure_shortest_paths_tree();

        Scene_polylines_item* polylines = new Scene_polylines_item();

        polylines->polylines.push_back(Scene_polylines_item::Polyline());

        CGAL::Three::Three::information(QObject::tr("Computing shortest path polyline..."));

        QElapsedTimer time;
        time.start();
        //~ m_shortestPaths->m_debugOutput=true;
        m_shortestPaths->shortest_path_points_to_source_points(
              faceLocation.first,
              faceLocation.second,
              std::back_inserter(polylines->polylines.back()));
        std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
        if(!polylines->polylines.front().empty())
        {
          polylines->setName(
                QObject::tr("%1 (shortest path)").arg(item->polyhedron_item()->name()));
          polylines->setColor(Qt::red);
          this->m_sceneInterface->setSelectedItem(-1);
          this->m_sceneInterface->addItem(polylines);
          this->m_sceneInterface->changeGroup(polylines, item->parentGroup());
        }
        else
          delete polylines;
      }
      else
      {
        CGAL::Three::Three::warning(QObject::tr("No source points to compute shortest paths from."));
      }
      break;
    }
    item->invalidateOpenGLBuffers();
    item->redraw();
    QApplication::restoreOverrideCursor();
    return true;
  }
}



bool Scene_polyhedron_shortest_path_item::eventFilter(QObject*, QEvent* event)
{
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();
    d->m_shiftHeld = modifiers.testFlag(Qt::ShiftModifier);
  }

  if (event->type() == QEvent::MouseButtonPress && d->m_shiftHeld)
  {
    QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
    if(mouseEvent->button() == Qt::LeftButton)
    {
      Ray_3 mouseRay;

      if (d->get_mouse_ray(mouseEvent, mouseRay))
      {
        if (d->run_point_select(mouseRay))
        {
          return true;
        }
      }
    }
  }

  return false;
}

bool Scene_polyhedron_shortest_path_item::load(const std::string& file_name)
{
  d->m_deferredLoadFilename = file_name;
  return true;
}

bool Scene_polyhedron_shortest_path_item::deferred_load(
    Scene_face_graph_item* polyhedronItem,
    CGAL::Three::Scene_interface* sceneInterface,
    Messages_interface* messages,
    QMainWindow* mainWindow)
{
  initialize(polyhedronItem, sceneInterface, messages, mainWindow);

  std::ifstream inFile(d->m_deferredLoadFilename.c_str());

  if (!inFile)
  {
    return false;
  }

  d->m_shortestPaths->clear();

  std::vector<face_descriptor> listOfFaces;
  listOfFaces.reserve(CGAL::num_faces(*polyhedron()));
  face_iterator current, end;
  for (boost::tie(current, end) = CGAL::faces(*polyhedron()); current != end; ++current)
  {
    listOfFaces.push_back(*current);
  }

  std::string line;
  std::size_t faceId;
  Barycentric_coordinates location;
  Construct_barycentric_coordinates construct_barycentric_coordinates;

  while (std::getline(inFile, line))
  {
    std::istringstream lineStream(line);
    FT coords[3];
    lineStream >> faceId >> coords[0] >> coords[1] >> coords[2];

    location = construct_barycentric_coordinates(coords[0], coords[1], coords[2]);

    // std::cout << "Read in face: " << faceId << " , " << location << std::endl;

    d->m_shortestPaths->add_source_point(listOfFaces[faceId], location);
  }

  return true;
}

bool Scene_polyhedron_shortest_path_item::save(const std::string& file_name) const
{
  std::ofstream out(file_name.c_str());
  boost::property_map<Face_graph, CGAL::face_index_t>::type fimap
      = get(CGAL::face_index, *polyhedron());
  if (!out)
  {
    return false;
  }

  for(Surface_mesh_shortest_path::Source_point_iterator it =
      d->m_shortestPaths->source_points_begin();
      it != d->m_shortestPaths->source_points_end();
      ++it)
  {
    out << get(fimap, it->first) << " " << it->second[0] << " " << it->second[1]
        << " " << it->second[3] << std::endl;
  }

  return true;
}

void Scene_polyhedron_shortest_path_item::initialize(
    Scene_face_graph_item* polyhedronItem,
    CGAL::Three::Scene_interface* sceneInterface,
    Messages_interface* messages, QMainWindow* mainWindow)
{
  d->m_mainWindow = mainWindow;
  d->m_messages = messages;
  this->poly_item = polyhedronItem;
  d->m_sceneInterface = sceneInterface;
  connect(polyhedronItem, SIGNAL(item_is_about_to_be_changed()), this,
          SLOT(poly_item_changed()));
  Q_FOREACH(CGAL::QGLViewer* viewer, CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(this);
  connect(d->m_mainWindow, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
  d->m_mainWindow->installEventFilter(this);
  d->recreate_shortest_path_object();
}

void Scene_polyhedron_shortest_path_item::deinitialize()
{
  d->deinitialize();
  this->poly_item = NULL;
}

bool Scene_polyhedron_shortest_path_item::isFinite() const
{
  return true;
}

bool Scene_polyhedron_shortest_path_item::isEmpty() const
{
  return false;
}

void Scene_polyhedron_shortest_path_item::compute_bbox() const
{
  setBbox(polyhedron_item()->bbox());
}

QString Scene_polyhedron_shortest_path_item::toolTip() const
{
  return QString();
}

void Scene_polyhedron_shortest_path_item::computeElements() const
{
  d->compute_elements();
  getPointContainer(0)->allocate(
        Pc::Vertices,
        d->vertices.data(),
        static_cast<int>(d->vertices.size()*sizeof(float)));
  d->nb_vertices = d->vertices.size();
  setBuffersFilled(true);
}

void Scene_polyhedron_shortest_path_item::initializeBuffers(Viewer_interface *v) const
{
  getPointContainer(0)->initializeBuffers(v);
  getPointContainer(0)->setFlatDataSize(d->nb_vertices);
  d->vertices.clear();
  d->vertices.shrink_to_fit();
}

void Scene_polyhedron_shortest_path_item::connectNewViewer(QObject* o)
{
  o->installEventFilter(this);
}
