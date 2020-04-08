#define CGAL_data_type float
#define CGAL_GL_data_type GL_FLOAT
#include "Scene_points_with_normal_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/algorithm.h>

#include <QObject>
#include <QApplication>
#include <QMenu>
#include <QSlider>
#include <QWidgetAction>
#include <CGAL/Qt/manipulatedCameraFrame.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

#include <CGAL/boost/graph/properties_Surface_mesh.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#endif // CGAL_LINKED_WITH_TBB

const std::size_t limit_fast_drawing = 300000; //arbitraty large value

typedef CGAL::Three::Point_container Pc;
typedef CGAL::Three::Edge_container Ec;
typedef CGAL::Three::Viewer_interface VI;
typedef Scene_points_with_normal_item_priv Priv;
struct Scene_points_with_normal_item_priv
{
  enum Point_container_id{
    Points =0,
    Shaded_points,
    Selected_points,
    Selected_shaded_points
  };

  void init_values(Scene_points_with_normal_item* parent)
  {
    item = parent;
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
    is_point_slider_moving = false;
    normal_Slider = new QSlider(Qt::Horizontal);
    normal_Slider->setValue(CGAL::Three::Three::getDefaultNormalLength());
    point_Slider = new QSlider(Qt::Horizontal);
    point_Slider->setMinimum(1);
    point_Slider->setValue(CGAL::Three::Three::getDefaultPointSize());
    point_Slider->setMaximum(25);
    item->setPointContainer(Priv::Selected_shaded_points, new Pc(VI::PROGRAM_WITH_LIGHT,
                                                 false));
    item->setPointContainer(Priv::Selected_points, new Pc(VI::PROGRAM_NO_SELECTION,
                                                 false));
    item->setPointContainer(Priv::Shaded_points, new Pc(VI::PROGRAM_WITH_LIGHT,
                                                 false));
    item->setPointContainer(Priv::Points, new Pc(VI::PROGRAM_NO_SELECTION,
                                                 false));
    item->setEdgeContainer(0, new Ec(VI::PROGRAM_NO_SELECTION,
                                                 false));

  }
  Scene_points_with_normal_item_priv(Scene_points_with_normal_item* parent)
    :m_points(new Point_set)
  {
    init_values(parent);
  }
  Scene_points_with_normal_item_priv(const Scene_points_with_normal_item& toCopy, Scene_points_with_normal_item* parent)
    : m_points(new Point_set(*toCopy.d->m_points))
  {
    init_values(parent);
  }

  Scene_points_with_normal_item_priv(const SMesh& input_mesh, Scene_points_with_normal_item* parent)
    : m_points(new Point_set)
  {
   init_values(parent);
   boost::graph_traits<SMesh>::vertex_iterator v;
    m_points->add_normal_map();
    for (v = const_cast<SMesh&>(input_mesh).vertices_begin();
         v != const_cast<SMesh&>(input_mesh).vertices_end(); v++)
    {
      boost::graph_traits<SMesh>::vertex_descriptor vd(*v);
      const Kernel::Point_3& p = input_mesh.point(vd);
      Kernel::Vector_3 n =
        CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, input_mesh);
      m_points->insert(p,n);
    }
  }

  ~Scene_points_with_normal_item_priv()
  {
    if(m_points)
    {
      delete m_points;
      m_points = NULL;
    }
    delete normal_Slider;
    delete point_Slider;
  }
  bool isPointSliderMoving() { return is_point_slider_moving; }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices() const;

  Point_set* m_points;
  std::string m_comments;
  QAction* actionDeleteSelection;
  QAction* actionResetSelection;
  QAction* actionSelectDuplicatedPoints;
  QSlider* normal_Slider;
  QSlider* point_Slider;
  mutable bool is_point_slider_moving;
  mutable std::vector<CGAL_data_type> positions_lines;
  mutable std::vector<CGAL_data_type> normals;
  mutable std::vector<CGAL_data_type> positions_normals;
  mutable std::vector<CGAL_data_type> positions_selected_normals;
  mutable std::vector<CGAL_data_type> colors_points;
  mutable std::size_t nb_points;
  mutable std::size_t nb_selected_points;
  mutable std::size_t nb_lines;
  mutable QOpenGLShaderProgram *program;

  Scene_points_with_normal_item* item;
};

class Fill_buffers {

  Point_set* point_set;
  std::vector<Point_set::Index>& indices;
  std::vector<CGAL_data_type>& positions_lines;
  std::vector<CGAL_data_type>& positions_normals;
  bool has_normals;
  const CGAL::qglviewer::Vec offset;
  double length;
  std::size_t size_p;
  std::size_t offset_normal_indices;

public:
  Fill_buffers(Point_set* point_set,
               std::vector<Point_set::Index>& indices,
               std::vector<CGAL_data_type>& positions_lines,
               std::vector<CGAL_data_type>& positions_normals,
               bool has_normals,
               const CGAL::qglviewer::Vec offset,
               double length,
               std::size_t offset_normal_indices = 0)
    : point_set (point_set)
    , indices (indices)
    , positions_lines (positions_lines)
    , positions_normals (positions_normals)
    , has_normals (has_normals)
    , offset (offset)
    , length (length)
    , offset_normal_indices (offset_normal_indices)
  {
    if (has_normals)
      size_p = 6;
    else
      size_p = 3;
  }

#ifdef CGAL_LINKED_WITH_TBB
  void operator()(const tbb::blocked_range<std::size_t>& r) const
  {
    for( std::size_t i = r.begin(); i != r.end(); ++i)
      apply (i);
  }
#endif // CGAL_LINKED_WITH_TBB

  void apply (std::size_t i) const
  {
    const Point_set::Index& idx = indices[i];
    const Kernel::Point_3& p = point_set->point(idx);

    positions_lines[i * size_p    ] = p.x() + offset.x;
    positions_lines[i * size_p + 1] = p.y() + offset.y;
    positions_lines[i * size_p + 2] = p.z() + offset.z;

    if(has_normals)
    {
      const Kernel::Vector_3& n = point_set->normal(idx);
      Kernel::FT normalizer = 1.0/CGAL::sqrt(n.squared_length());
      Point_set_3<Kernel>::Point q = p + length * n * normalizer;
      positions_lines[i * size_p + 3] = q.x() + offset.x;
      positions_lines[i * size_p + 4] = q.y() + offset.y;
      positions_lines[i * size_p + 5] = q.z() + offset.z;

      positions_normals[(i - offset_normal_indices) * 3    ] = n.x();
      positions_normals[(i - offset_normal_indices) * 3 + 1] = n.y();
      positions_normals[(i - offset_normal_indices) * 3 + 2] = n.z();
    }
  }
};



Scene_points_with_normal_item::Scene_points_with_normal_item()
{
    setRenderingMode(Points);
    is_selected = true;
    d = new Scene_points_with_normal_item_priv(this);
}

// Copy constructor
Scene_points_with_normal_item::Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy)
{

  d = new Scene_points_with_normal_item_priv(toCopy, this);

  if (!has_normals())
  {
    setRenderingMode(Points);
    is_selected = true;
  }
  else{
    setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
    is_selected = true;
  }

  invalidateOpenGLBuffers();
}

// Converts polyhedron to point set

Scene_points_with_normal_item::Scene_points_with_normal_item(const SMesh& input_mesh)
{
  // Converts Polyhedron vertices to point set.
  // Computes vertices normal from connectivity.
  d = new Scene_points_with_normal_item_priv(input_mesh, this);
  setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
  is_selected = true;
  invalidateOpenGLBuffers();
}

Scene_points_with_normal_item::~Scene_points_with_normal_item()
{
  delete d;
}



void Scene_points_with_normal_item_priv::
initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  item->getEdgeContainer(0)->initializeBuffers(viewer);
  item->getPointContainer(Priv::Points)->initializeBuffers(viewer);
  item->getPointContainer(Priv::Shaded_points)->initializeBuffers(viewer);
  item->getPointContainer(Priv::Selected_points)->initializeBuffers(viewer);
  item->getPointContainer(Priv::Selected_shaded_points)->initializeBuffers(viewer);

  ////Clean-up
  item->getPointContainer(Priv::Points)->setFlatDataSize(nb_points - nb_selected_points);
  item->getPointContainer(Priv::Shaded_points)->setFlatDataSize(nb_points - nb_selected_points);
  item->getPointContainer(Priv::Selected_points)->setFlatDataSize(nb_selected_points);
  item->getPointContainer(Priv::Selected_shaded_points)->setFlatDataSize(nb_selected_points);
  item->getEdgeContainer(0)->setFlatDataSize(nb_lines);

  positions_lines             .resize(0);
  normals                     .resize(0);
  positions_normals           .resize(0);
  positions_selected_normals  .resize(0);
  colors_points               .resize(0);

  positions_lines             .shrink_to_fit();
  normals                     .shrink_to_fit();
  positions_normals           .shrink_to_fit();
  positions_selected_normals  .shrink_to_fit();
  colors_points               .shrink_to_fit();
}

void Scene_points_with_normal_item_priv::compute_normals_and_vertices() const
{
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first())->offset();
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_lines.resize(0);
    normals.resize(0);
    positions_normals.resize(0);
    positions_selected_normals.resize(0);
    normals.resize(0);
    colors_points.resize(0);

    //Shuffle container to allow quick display random points
    std::vector<Point_set::Index> indices;
    indices.reserve (m_points->size());
    std::copy (m_points->begin(), m_points->end(), std::back_inserter(indices));

    CGAL::cpp98::random_shuffle (indices.begin(), indices.end() - m_points->nb_selected_points());
    if (m_points->nb_selected_points() != 0)
      CGAL::cpp98::random_shuffle (indices.end() - m_points->nb_selected_points(), indices.end());

    //if item has normals, points will be one point out of two in the lines data.
    //else points will be lines and lines discarded.
    double average_spacing = 0;
    double normal_length =0;
    double length_factor =0;
    if (item->has_normals())
    {
      // Store normals
      Kernel::Sphere_3 region_of_interest = m_points->region_of_interest();
      positions_lines.resize(m_points->size() * 6);
      positions_normals.resize((m_points->size() - m_points->nb_selected_points()) * 3);
      positions_selected_normals.resize(m_points->nb_selected_points() * 3);

      // we can't afford computing real average spacing just for display, 0.5% of bbox will do
      average_spacing = 0.005 * item->diagonalBbox();
      normal_length = (std::min)(average_spacing, std::sqrt(
                                   region_of_interest.squared_radius() / 1000.));
      length_factor = 10.0/100*normal_Slider->value();
    }
    else
    {
      positions_lines.resize(m_points->size() * 3);
    }

    Fill_buffers fill_buffers (m_points, indices, positions_lines, positions_normals,
                               item->has_normals(), offset, normal_length * length_factor);
    Fill_buffers fill_buffers_2 (m_points, indices, positions_lines, positions_selected_normals,
                                 item->has_normals(), offset, normal_length * length_factor,
                                 m_points->first_selected() - m_points->begin());

#ifdef CGAL_LINKED_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0,
                                                 m_points->size() - m_points->nb_selected_points()),
                      fill_buffers);
    tbb::parallel_for(tbb::blocked_range<size_t>(m_points->size() - m_points->nb_selected_points(),
                                                 m_points->size()),
                      fill_buffers_2);
#else
    for (std::size_t i = 0; i < indices.size() - m_points->nb_selected_points(); ++ i)
      fill_buffers.apply (i);
    for (std::size_t i = indices.size() - m_points->nb_selected_points(); i < indices.size(); ++ i)
      fill_buffers_2.apply (i);
#endif

    //The colors
    if (m_points->has_colors())
    {
        colors_points.reserve((m_points->size() - m_points->nb_selected_points()) * 6);

        for (std::size_t i = 0; i < indices.size() - m_points->nb_selected_points(); ++ i)
        {
          colors_points.push_back (m_points->red(indices[i]));
          colors_points.push_back (m_points->green(indices[i]));
          colors_points.push_back (m_points->blue(indices[i]));
          colors_points.push_back (m_points->red(indices[i]));
          colors_points.push_back (m_points->green(indices[i]));
          colors_points.push_back (m_points->blue(indices[i]));
        }
    }

    nb_lines = positions_lines.size();
    if(item->has_normals())
      nb_points = positions_lines.size()/2;
    else
      nb_points = positions_lines.size();
    nb_selected_points = m_points->nb_selected_points() * 3;
    //edges
    if(item->has_normals())
    {
      item->getEdgeContainer(0)->allocate(Ec::Vertices, positions_lines.data(),
                                       static_cast<int>(nb_lines*sizeof(CGAL_data_type)));
      if (!(colors_points.empty()))
      {
        item->getEdgeContainer(0)->allocate(Ec::Colors, colors_points.data(),
                                            static_cast<int>(colors_points.size()
                                                             *sizeof(CGAL_data_type)));
      }
      //shaded points
      item->getPointContainer(Priv::Shaded_points)->setStride(Pc::Vertices,
                                           static_cast<int>(6*sizeof(CGAL_GL_data_type)));
      item->getPointContainer(Priv::Shaded_points)->allocate(Pc::Vertices,  positions_lines.data(),
                                       static_cast<int>(nb_lines*sizeof(CGAL_data_type)));
      item->getPointContainer(Priv::Shaded_points)->allocate(Pc::Normals, positions_normals.data(),
                                           static_cast<int>(positions_normals.size()
                                                            *sizeof(CGAL_data_type)));
    }
    //points
    if(!item->has_normals()) {
      item->getPointContainer(Priv::Points)->setStride(Pc::Vertices,
                                            0);
    }
    else{
      item->getPointContainer(Priv::Points)->setStride(Pc::Vertices,
                                            static_cast<int>(6*sizeof(CGAL_GL_data_type)));
    }
    item->getPointContainer(Priv::Points)->allocate(Pc::Vertices, positions_lines.data(),
                                         static_cast<int>(positions_lines.size()
                                                          *sizeof(CGAL_data_type)));
    if (!(colors_points.empty()))
    {
      item->getPointContainer(Priv::Points)->setStride(Pc::Colors,6*sizeof(CGAL_data_type));
      item->getPointContainer(Priv::Points)->allocate(Pc::Colors, colors_points.data(),
                                           static_cast<int>(colors_points.size()
                                                            *sizeof(CGAL_data_type)));
    }
    //selected points
    if(!item->has_normals()) {
      item->getPointContainer(Priv::Selected_points)->setStride(Pc::Vertices,
                                            0);
      item->getPointContainer(Priv::Selected_points)->setOffset(Pc::Vertices,
                                            static_cast<int>(
                                              3*
                                              (m_points->size()
                                               -m_points->nb_selected_points())
                                              *sizeof(CGAL_data_type) ));
    }
    else{
      item->getPointContainer(Priv::Selected_points)->setStride(Pc::Vertices,
                                            static_cast<int>(6*sizeof(CGAL_GL_data_type)));
      item->getPointContainer(Priv::Selected_points)->setOffset(Pc::Vertices,
                                            static_cast<int>(
                                              6*
                                              (m_points->size()
                                               -m_points->nb_selected_points())
                                              *sizeof(CGAL_data_type) ));
      item->getPointContainer(Priv::Selected_shaded_points)->setStride(Pc::Vertices,
                                            static_cast<int>(6*sizeof(CGAL_GL_data_type)));
      item->getPointContainer(Priv::Selected_shaded_points)->setOffset(Pc::Vertices,
                                            static_cast<int>(
                                              6*
                                              (m_points->size()
                                               -m_points->nb_selected_points())
                                              *sizeof(CGAL_data_type) ));
      item->getPointContainer(Priv::Selected_shaded_points)->allocate(Pc::Vertices,
                                                                      positions_lines.data(),
                                           static_cast<int>(positions_lines.size()
                                                            *sizeof(CGAL_data_type)));
      item->getPointContainer(Priv::Selected_shaded_points)->allocate(Pc::Normals,
                                                                      positions_selected_normals.data(),
                                           static_cast<int>(positions_selected_normals.size()
                                                            *sizeof(CGAL_data_type)));
    }
    item->getPointContainer(Priv::Selected_points)->allocate(Pc::Vertices, positions_lines.data(),
                                         static_cast<int>(positions_lines.size()
                                                          *sizeof(CGAL_data_type)));
    if (!(colors_points.empty()))
    {
      item->getPointContainer(Priv::Points)->setStride(Pc::Colors,6*sizeof(CGAL_data_type));
      item->getPointContainer(Priv::Points)->allocate(Pc::Colors, colors_points.data(),
                                           static_cast<int>(colors_points.size()
                                                            *sizeof(CGAL_data_type)));
    }
    QApplication::restoreOverrideCursor();
}

// Duplicates scene item
Scene_points_with_normal_item*
Scene_points_with_normal_item::clone() const
{
  return new Scene_points_with_normal_item(*this);
}

// Is selection empty?
bool Scene_points_with_normal_item::isSelectionEmpty() const
{
  return (d->m_points->nb_selected_points() == 0);
}

// Delete selection
void Scene_points_with_normal_item::deleteSelection()
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Delete " << d->m_points->nb_selected_points() << " points...";

  // Delete selected points
  d->m_points->delete_selection();

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Invert selection
void Scene_points_with_normal_item::invertSelection()
{
  d->m_points->invert_selection();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Select everything
void Scene_points_with_normal_item::selectAll()
{
  d->m_points->select_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
// Reset selection mark
void Scene_points_with_normal_item::resetSelection()
{
  // Un-select all points
  d->m_points->unselect_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
  //Select duplicated points
void Scene_points_with_normal_item::selectDuplicates()
{
  std::set<Kernel::Point_3> unique_points;
  std::vector<Point_set::Index> unselected, selected;
  for (Point_set::iterator ptit = d->m_points->begin(); ptit!= d->m_points->end(); ++ ptit)
    if ( !unique_points.insert(d->m_points->point(*ptit)).second)
      selected.push_back (*ptit);
    else
      unselected.push_back (*ptit);

  for (std::size_t i = 0; i < unselected.size(); ++ i)
    *(d->m_points->begin() + i) = unselected[i];
  for (std::size_t i = 0; i < selected.size(); ++ i)
    *(d->m_points->begin() + (unselected.size() + i)) = selected[i];

  if (selected.empty ())
  {
    d->m_points->unselect_all();
  }
  else
  {
    d->m_points->set_first_selected
      (d->m_points->begin() + unselected.size());
  }

  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

#ifdef CGAL_LINKED_WITH_LASLIB
// Loads point set from .LAS file
bool Scene_points_with_normal_item::read_las_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_las_point_set (stream, *(d->m_points)) &&
            !isEmpty();

  std::cerr << d->m_points->info();

  if (!d->m_points->has_normal_map())
  {
    setRenderingMode(Points);
  }
  else{
    setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
  }
  if (d->m_points->check_colors())
    std::cerr << "-> Point set has colors" << std::endl;

  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .LAS file
bool Scene_points_with_normal_item::write_las_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->reset_indices();

  return stream &&
    CGAL::write_las_point_set (stream, *(d->m_points));
}

#endif // LAS

// Loads point set from .PLY file
bool Scene_points_with_normal_item::read_ply_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_ply_point_set (stream, *(d->m_points), d->m_comments) &&
            !isEmpty();
    d->point_Slider->setValue(CGAL::Three::Three::getDefaultPointSize());
  std::cerr << d->m_points->info();

  if (!d->m_points->has_normal_map())
  {
    setRenderingMode(Points);
  }
  else{
    setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
  }
  if (d->m_points->check_colors())
    std::cerr << "-> Point set has colors" << std::endl;

  std::cerr << "[Comments from PLY input]" << std::endl << d->m_comments;

  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .PLY file
bool Scene_points_with_normal_item::write_ply_point_set(std::ostream& stream, bool binary) const
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->reset_indices();

  if (!stream)
    return false;

  if (binary)
    CGAL::set_binary_mode (stream);

  CGAL::write_ply_point_set (stream, *(d->m_points), d->m_comments);

  return true;
}

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
    CGAL::read_off_point_set(stream, *(d->m_points)) &&
            !isEmpty();
  d->point_Slider->setValue(CGAL::Three::Three::getDefaultPointSize());
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->reset_indices();

  return stream &&
    CGAL::write_off_point_set (stream, *(d->m_points));
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_xyz_point_set (stream, *(d->m_points)) &&
    !isEmpty();
  d->point_Slider->setValue(CGAL::Three::Three::getDefaultPointSize());
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .XYZ file
bool Scene_points_with_normal_item::write_xyz_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->reset_indices();

  return stream &&
    CGAL::write_xyz_point_set (stream, *(d->m_points));
}

QString
Scene_points_with_normal_item::toolTip() const
{
  Q_ASSERT(d->m_points != NULL);

  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point_set_3</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(d->m_points->size())
    .arg(color().name());
}

bool Scene_points_with_normal_item::supportsRenderingMode(RenderingMode m) const
{
  switch ( m )
  {
  case Points:
    return true;
  case ShadedPoints:
  case PointsPlusNormals:
    return has_normals();

  default:
    return false;
  }
}

void Scene_points_with_normal_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        (d->nb_lines/6 > limit_fast_drawing)) // arbitrary large value
      ratio_displayed = 6 * limit_fast_drawing / (double)(d->nb_lines);
    if(!isInit(viewer))
      initGL(viewer);
    if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
    {
      d->initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }
    if(!d->m_points->has_colors())
      getEdgeContainer(0)->setColor(color());
    std::size_t real_size =
        getEdgeContainer(0)->getFlatDataSize();
    getEdgeContainer(0)->setFlatDataSize(ratio_displayed * real_size);
    getEdgeContainer(0)->draw( viewer, !d->m_points->has_colors());
    getEdgeContainer(0)->setFlatDataSize(real_size);

}
void Scene_points_with_normal_item::
drawPoints(CGAL::Three::Viewer_interface* viewer) const
{

  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  viewer->setGlPointSize(GLfloat(d->point_Slider->value()));
  double ratio_displayed = 1.0;
  if ((viewer->inFastDrawing () || d->isPointSliderMoving())
      &&((d->nb_points )/3 > limit_fast_drawing)) // arbitrary large value
    ratio_displayed = 3 * limit_fast_drawing / (double)(d->nb_points);

  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
     ! getBuffersInit(viewer))
  {
    d->initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    d->initializeBuffers(viewer);
  }

  if(has_normals() && renderingMode() == ShadedPoints)
  {
    getPointContainer(Priv::Shaded_points)->setColor(color());
    std::size_t real_size =
        getPointContainer(Priv::Shaded_points)->getFlatDataSize();
    getPointContainer(Priv::Shaded_points)->setFlatDataSize(ratio_displayed * real_size);
    getPointContainer(Priv::Shaded_points)->draw( viewer, true);
    getPointContainer(Priv::Shaded_points)->setFlatDataSize(real_size);

    real_size =
        getPointContainer(Priv::Selected_shaded_points)->getFlatDataSize();
    getPointContainer(Priv::Selected_shaded_points)->setColor(QColor(Qt::red));
    getPointContainer(Priv::Selected_shaded_points)->setFlatDataSize(ratio_displayed * real_size);
    getPointContainer(Priv::Selected_shaded_points)->draw( viewer, true);
    getPointContainer(Priv::Selected_shaded_points)->setFlatDataSize(real_size);
  }
  else
  {
    if(!d->m_points->has_colors())
      getPointContainer(Priv::Points)->setColor(color());
    std::size_t real_size =
        getPointContainer(Priv::Points)->getFlatDataSize();
    getPointContainer(Priv::Points)->setFlatDataSize(ratio_displayed * real_size);
    getPointContainer(Priv::Points)->setFlatDataSize(real_size);
    getPointContainer(Priv::Points)->draw( viewer, !d->m_points->has_colors());

    real_size =
        getPointContainer(Priv::Selected_points)->getFlatDataSize();
    getPointContainer(Priv::Selected_points)->setColor(QColor(Qt::red));
    getPointContainer(Priv::Selected_points)->setFlatDataSize(ratio_displayed * real_size);
    getPointContainer(Priv::Selected_points)->setFlatDataSize(real_size);
    getPointContainer(Priv::Selected_points)->draw( viewer, true);
  }

  viewer->setGlPointSize(point_size);
}
// Gets wrapped point set
Point_set* Scene_points_with_normal_item::point_set()
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points;
}
const Point_set* Scene_points_with_normal_item::point_set() const
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points;
}

// Gets wrapped point set
std::string& Scene_points_with_normal_item::comments()
{
  return d->m_comments;
}
const std::string& Scene_points_with_normal_item::comments() const
{
  return d->m_comments;
}

bool
Scene_points_with_normal_item::isEmpty() const
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points->empty();
}

void
Scene_points_with_normal_item::compute_bbox()const
{
  Q_ASSERT(d->m_points != NULL);

  Kernel::Iso_cuboid_3 bbox = d->m_points->bounding_box();
  setBbox(Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax()));
}

void Scene_points_with_normal_item::computes_local_spacing(int k)
{
  typedef Kernel Geom_traits;

  typedef CGAL::Search_traits_3<Geom_traits> SearchTraits_3;
  typedef CGAL::Search_traits_adapter <Point_set::Index, Point_set::Point_map, SearchTraits_3> Search_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef Neighbor_search::Distance Distance;

  // build kdtree
  Tree tree(d->m_points->begin(),
            d->m_points->end(),
            Tree::Splitter(),
            Search_traits (d->m_points->point_map())
            );
  Distance tr_dist(d->m_points->point_map());

  if (!(d->m_points->has_property_map<double> ("radius")))
    d->m_points->add_radius();

  // Compute the radius of each point = (distance max to k nearest neighbors)/2.
  {
    int i=0;
    for (Point_set::iterator it=d->m_points->begin(); it!=d->m_points->end(); ++it, ++i)
    {
      Neighbor_search search(tree, d->m_points->point(*it), k+1, 0, true, tr_dist);
      double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
      d->m_points->radius(*it) = sqrt(maxdist2)/2.;
    }
  }

  d->m_points->set_radii_uptodate(true);
}

QMenu* Scene_points_with_normal_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_points_with_normal_item.";

    QMenu* menu = Scene_item::contextMenu();

    //add a slider to modify the normals length
    // Use dynamic properties:
    // https://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
      if(has_normals())
      {
        QMenu *container = new QMenu(tr("Normals Length"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        if((d->nb_points)/3 <= limit_fast_drawing)
        {
          connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
          connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::itemChanged);
        }
        else
        {
          connect(d->normal_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
          connect(d->normal_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::itemChanged);
        }
        sliderAction->setDefaultWidget(d->normal_Slider);
        container->menuAction()->setProperty("is_groupable", true);
        container->addAction(sliderAction);
        menu->addMenu(container);
      }
        QMenu *container = new QMenu(tr("Points Size"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        connect(d->point_Slider, &QSlider::sliderPressed, this, &Scene_points_with_normal_item::pointSliderPressed);
        connect(d->point_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::pointSliderReleased);
        connect(d->point_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::itemChanged);

        sliderAction->setDefaultWidget(d->point_Slider);
        container->menuAction()->setProperty("is_groupable", true);
        container->addAction(sliderAction);
        menu->addMenu(container);

        d->actionDeleteSelection = menu->addAction(tr("Delete Selection"));
        d->actionDeleteSelection->setObjectName("actionDeleteSelection");
        connect(d->actionDeleteSelection, SIGNAL(triggered()),this, SLOT(deleteSelection()));

        d->actionResetSelection = menu->addAction(tr("Reset Selection"));
        d->actionResetSelection->setObjectName("actionResetSelection");
        connect(d->actionResetSelection, SIGNAL(triggered()),this, SLOT(resetSelection()));

        d->actionSelectDuplicatedPoints = menu->addAction(tr("Select duplicated points"));
        d->actionSelectDuplicatedPoints->setObjectName("actionSelectDuplicatedPoints");
        connect(d->actionSelectDuplicatedPoints, SIGNAL(triggered()),this, SLOT(selectDuplicates()));
        QAction* resetColorsAction = menu->addAction(tr("Make Unicolor"));
        resetColorsAction->setObjectName("resetColorsAction");
        connect(resetColorsAction, &QAction::triggered, this, &Scene_points_with_normal_item::resetColors);
        menu->setProperty(prop_name, true);
    }
    QAction* actionColor = menu->findChild<QAction*>(tr("resetColorsAction"));
    actionColor->setVisible(d->m_points->has_colors());
    if (isSelectionEmpty())
    {
        d->actionDeleteSelection->setDisabled(true);
        d->actionResetSelection->setDisabled(true);
    }
    else
    {
        d->actionDeleteSelection->setDisabled(false);
        d->actionResetSelection->setDisabled(false);
    }

    return menu;
}

bool Scene_points_with_normal_item::has_normals() const { return d->m_points->has_normal_map(); }

void Scene_points_with_normal_item::invalidateOpenGLBuffers()
{
    setBuffersFilled(false);
    getPointContainer(Priv::Points)->reset_vbos(Scene_item_rendering_helper::ALL);
    getPointContainer(Priv::Selected_points)->reset_vbos(Scene_item_rendering_helper::ALL);
    getPointContainer(Priv::Selected_shaded_points)->reset_vbos(Scene_item_rendering_helper::ALL);
    getPointContainer(Priv::Shaded_points)->reset_vbos(Scene_item_rendering_helper::ALL);
    getEdgeContainer(0)->reset_vbos(Scene_item_rendering_helper::ALL);

    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      if(viewer == NULL)
        continue;
      setBuffersInit(viewer, false);
    }
}

void Scene_points_with_normal_item::pointSliderPressed()
{
  d->is_point_slider_moving = true;
}

void Scene_points_with_normal_item::pointSliderReleased()
{
  d->is_point_slider_moving = false;
}

void Scene_points_with_normal_item::copyProperties(Scene_item *item)
{
  Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(item);
  if(!point_item)
    return;
  int value = point_item->getPointSliderValue();
  setPointSize(value);
  if(has_normals())
    d->normal_Slider->setValue(point_item->getNormalSliderValue());
}

int Scene_points_with_normal_item::getNormalSliderValue()
{
  return d->normal_Slider->value();
}

int Scene_points_with_normal_item::getPointSliderValue()
{
  return d->point_Slider->value();
}

void Scene_points_with_normal_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);
  if(d && d->m_points && item == this)
  {
    delete d->m_points;
    d->m_points = NULL;
  }
}

void Scene_points_with_normal_item::
zoomToPosition(const QPoint &, CGAL::Three::Viewer_interface *viewer) const
{
  if (point_set()->nb_selected_points() == 0)
    return;
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  // fit plane to triangles
  Point_set points;
  Bbox selected_points_bbox;
  for(Point_set::const_iterator it = point_set()->first_selected();
      it != point_set()->end();
      ++it)
  {
    points.insert(point_set()->point(*it));
    selected_points_bbox += point_set()->point(*it).bbox();
  }
  Kernel::Plane_3 plane;
  Kernel::Point_3 center_of_mass;
  CGAL::linear_least_squares_fitting_3
      (points.points().begin(),points.points().end(),plane, center_of_mass,
       CGAL::Dimension_tag<0>());

  Kernel::Vector_3 plane_normal= plane.orthogonal_vector();
  plane_normal = plane_normal/(CGAL::sqrt(plane_normal.squared_length()));
  Kernel::Point_3 centroid(center_of_mass.x() + offset.x,
                           center_of_mass.y() + offset.y,
                           center_of_mass.z() + offset.z);

  CGAL::qglviewer::Quaternion new_orientation(CGAL::qglviewer::Vec(0,0,-1),
                                        CGAL::qglviewer::Vec(-plane_normal.x(), -plane_normal.y(), -plane_normal.z()));
  double max_side = (std::max)((std::max)(selected_points_bbox.xmax() - selected_points_bbox.xmin(),
                                          selected_points_bbox.ymax() - selected_points_bbox.ymin()),
                               selected_points_bbox.zmax() - selected_points_bbox.zmin());
  //put the camera in way we are sure the longest side is entirely visible on the screen
  //See openGL's frustum definition
  double factor = max_side/(tan(viewer->camera()->aspectRatio()/
                                  (viewer->camera()->fieldOfView()/2)));

  Kernel::Point_3 new_pos = centroid + factor*plane_normal ;
  viewer->camera()->setSceneCenter(CGAL::qglviewer::Vec(centroid.x(),
                                                  centroid.y(),
                                                  centroid.z()));
  viewer->moveCameraToCoordinates(QString("%1 %2 %3 %4 %5 %6 %7").arg(new_pos.x())
                                                                 .arg(new_pos.y())
                                                                 .arg(new_pos.z())
                                                                 .arg(new_orientation[0])
                                                                 .arg(new_orientation[1])
                                                                 .arg(new_orientation[2])
                                                                 .arg(new_orientation[3]));

}

void Scene_points_with_normal_item::setPointSize(int size)
{
  d->point_Slider->setValue(size);
}

void Scene_points_with_normal_item::setNormalSize(int size)
{
  d->normal_Slider->setValue(size);
}

void Scene_points_with_normal_item::resetColors()
{
  d->m_points->remove_colors();
  invalidateOpenGLBuffers();
  redraw();
}

void Scene_points_with_normal_item::computeElements()const
{
  d->compute_normals_and_vertices();
  setBuffersFilled(true);

}

void Scene_points_with_normal_item::initializeBuffers(Viewer_interface * v) const
{
  d->initializeBuffers(v);
}
