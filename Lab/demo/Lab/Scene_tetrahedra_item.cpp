#include <QLabel>
#include <QColor>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Three.h>
#include "Scene_tetrahedra_item.h"
#include "C3t3_type.h"
#include "CGAL_double_edit.h"

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tri;

struct tet_item_priv{
  mutable Scene_c3t3_item* c3t3_item;
  mutable std::vector<float> positions;
  mutable std::vector<float> filter_values[4];
  mutable std::size_t positions_size;
  mutable std::vector<float> normals;
  mutable std::vector<float> colors;
  mutable CGAL::qglviewer::Vec offset;
  float min_threshold;
  float max_threshold;
  float threshold_1;
  float threshold_2;
  int filter_index;
  double extrema[4][2];
  QLabel* minMinLabel;
  QLabel* minMaxLabel;
  QLabel* maxMinLabel;
  QLabel* maxMaxLabel;

  DoubleEdit* minEdit;
  DoubleEdit* maxEdit;

  void draw_triangle(const Tr::Bare_point& a,
                     const Tr::Bare_point& b,
                     const Tr::Bare_point& c,
                     const QColor& color)
  {
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    Geom_traits::Vector_3 n = cross_product(b - a, c - a);
    n = n / CGAL::sqrt(n*n);

    auto push_normal = [this](auto n) {
      normals.push_back(static_cast<float>(n.x()));
      normals.push_back(static_cast<float>(n.y()));
      normals.push_back(static_cast<float>(n.z()));
    };
    auto push_color = [this](auto c) {
      colors.push_back(static_cast<float>(c.redF()));
      colors.push_back(static_cast<float>(c.greenF()));
      colors.push_back(static_cast<float>(c.blueF()));
    };

    positions.push_back(static_cast<float>(a.x()+offset.x));
    positions.push_back(static_cast<float>(a.y()+offset.y));
    positions.push_back(static_cast<float>(a.z()+offset.z));

    positions.push_back(static_cast<float>(b.x()+offset.x));
    positions.push_back(static_cast<float>(b.y()+offset.y));
    positions.push_back(static_cast<float>(b.z()+offset.z));

    positions.push_back(static_cast<float>(c.x()+offset.x));
    positions.push_back(static_cast<float>(c.y()+offset.y));
    positions.push_back(static_cast<float>(c.z()+offset.z));

    for (int i = 0; i<3; i++)
    {
      push_normal(n);
      push_color(color);
    }
  }


};
Scene_tetrahedra_item::Scene_tetrahedra_item(Scene_c3t3_item* c3t3_item)
  : Scene_item_rendering_helper(),
    d(new tet_item_priv)
{
  d->c3t3_item = c3t3_item;
  d->positions_size = 0;
  d->min_threshold = 0.0f;
  d->max_threshold = 1.0f;
  d->threshold_1 = 0.0f;
  d->threshold_2 = 1.0f;
  d->filter_index = 0;
  d->extrema[0][0] = 360.0;
  d->extrema[0][1] = 360;
  d->extrema[1][0] = 0;
  d->extrema[1][1] = 0;
  d->extrema[2][0] = 1.0;
  d->extrema[2][1] = 0;
  d->extrema[3][0] = (std::numeric_limits<double>::max)();
  d->extrema[3][1] = 0.0;
  d->minMinLabel = nullptr;
  d->minMaxLabel = nullptr;
  d->maxMinLabel = nullptr;
  d->maxMaxLabel = nullptr;

  d->minEdit = nullptr;
  d->maxEdit = nullptr;
  setFlatMode();
  setTriangleContainer(0,
                       new Tri(Vi::PROGRAM_TETRA_FILTERING, false));

}

void Scene_tetrahedra_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer))
  {
    initGL(viewer);
  }
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
  getTriangleContainer(0)->setColor(this->color());
  getTriangleContainer(0)->getVao(viewer)->bind();
  getTriangleContainer(0)->getVao(viewer)->program->setUniformValue("min_threshold", d->min_threshold);
  getTriangleContainer(0)->getVao(viewer)->program->setUniformValue("max_threshold", d->max_threshold);
  getTriangleContainer(0)->getVao(viewer)->release();
  getTriangleContainer(0)->draw(viewer, true);

}


void Scene_tetrahedra_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  compute_bbox();

}

void Scene_tetrahedra_item::computeElements()const
{
  CGAL::Three::Three::CursorScopeGuard guard{QCursor(Qt::WaitCursor)};
  for(int i=0; i<4; ++i)
  {
    d->filter_values[i].clear();
    d->filter_values[i].shrink_to_fit();
  }
  d->offset = CGAL::Three::Three::mainViewer()->offset();
  const C3t3& c3t3 = d->c3t3_item->c3t3();
  Geom_traits::Construct_point_3 wp2p
      = c3t3.triangulation().geom_traits().construct_point_3_object();
  Geom_traits::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
      = c3t3.triangulation().geom_traits().compute_approximate_dihedral_angle_3_object();
  d->extrema[0][0] = 360.0;
  d->extrema[0][1] = 0;
  d->extrema[1][0] = 360;
  d->extrema[1][1] = 0;
  d->extrema[2][0] =  9999999999;
  d->extrema[2][1] =  0;
  d->extrema[3][0] = 1.0;
  d->extrema[3][1] = 0.0;
  std::vector<double> tmp_min, tmp_max, tmp_rrr, tmp_v;
  tmp_min.reserve(c3t3.number_of_cells_in_complex());
  tmp_max.reserve(c3t3.number_of_cells_in_complex());
  tmp_rrr.reserve(c3t3.number_of_cells_in_complex());
  tmp_v.reserve(c3t3.number_of_cells_in_complex());


  for(C3t3::Triangulation::Cell_iterator
      cit = c3t3.triangulation().finite_cells_begin(),
      end = c3t3.triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if(!c3t3.is_in_complex(cit))
    {
      continue;
    }
    const Tr::Bare_point& pa = wp2p(cit->vertex(0)->point());
    const Tr::Bare_point& pb = wp2p(cit->vertex(1)->point());
    const Tr::Bare_point& pc = wp2p(cit->vertex(2)->point());
    const Tr::Bare_point& pd = wp2p(cit->vertex(3)->point());
    const QColor color = d->c3t3_item->getSubdomainIndexColor(cit->subdomain_index());
    //0 - 1 : dihedral angle
    double min_dihedral_angle = 360.0;
    double max_dihedral_angle = -360.0;
    double a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pa, pb, pc, pd)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pa, pc, pb, pd)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pa, pd, pb, pc)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pb, pc, pa, pd)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pb, pd, pa, pc)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(pc, pd, pa, pb)));
    if(a < min_dihedral_angle) { min_dihedral_angle = static_cast<float>(a); }
    if(a > max_dihedral_angle) { max_dihedral_angle = static_cast<float>(a); }
    tmp_min.push_back(min_dihedral_angle);
    tmp_max.push_back(max_dihedral_angle);
    if(min_dihedral_angle < d->extrema[0][0]) { d->extrema[0][0]=min_dihedral_angle; }
    if(min_dihedral_angle > d->extrema[0][1]) { d->extrema[0][1]=min_dihedral_angle; }
    if(max_dihedral_angle < d->extrema[1][0]) { d->extrema[1][0]=max_dihedral_angle; }
    if(max_dihedral_angle > d->extrema[1][1]) { d->extrema[1][1]=max_dihedral_angle; }
    //3 : volume
    double v = std::abs(CGAL::volume(pa, pb, pc, pd));
    tmp_v.push_back(v);
    if(v < d->extrema[3][0]) { d->extrema[3][0] = static_cast<float>(v); }
    if(v > d->extrema[3][1]) { d->extrema[3][1] = static_cast<float>(v); }
    //2 : Radius-radius ratio
    double sumar = std::sqrt(CGAL::squared_area(pa,pb,pc))+std::sqrt(CGAL::squared_area(pb,pc,pd))+
        std::sqrt(CGAL::squared_area(pc,pd,pa)) + std::sqrt(CGAL::squared_area(pd,pb,pa));
    double inradius = 3*v/sumar;
    double circumradius = std::sqrt(CGAL::squared_radius(pa, pb, pc, pd));
    double rrr = inradius/circumradius*3; //*3 so that the perfect tet ratio is 1 instead of 1/3
    if(rrr < d->extrema[2][0]) { d->extrema[2][0] = static_cast<float>(rrr); }
    if(rrr > d->extrema[2][1]) { d->extrema[2][1] = static_cast<float>(rrr); }
    tmp_rrr.push_back(rrr);
    d->draw_triangle(pb, pa, pc, color);
    d->draw_triangle(pa, pb, pd, color);
    d->draw_triangle(pa, pd, pc, color);
    d->draw_triangle(pb, pc, pd, color);
  }
  d->positions_size = d->positions.size();
  for(std::size_t i = 0; i< tmp_min.size(); ++i)
  {
    float min = tmp_min[i]/(d->extrema[0][1]-d->extrema[0][0]) - d->extrema[0][0]/(d->extrema[0][1]-d->extrema[0][0]);
    float max = tmp_max[i]/(d->extrema[1][1]-d->extrema[1][0]) - d->extrema[1][0]/(d->extrema[1][1]-d->extrema[1][0]);
    float rrr = tmp_rrr[i]/(d->extrema[2][1]-d->extrema[2][0]) - d->extrema[2][0]/(d->extrema[2][1]-d->extrema[2][0]);
    float v = static_cast<float>(tmp_v[i]/(d->extrema[3][1]-d->extrema[3][0]) - d->extrema[3][0]/(d->extrema[3][1] - d->extrema[3][0]));
    for(int j = 0; j< 12; ++j)
    {
      d->filter_values[0].push_back(min);
      d->filter_values[1].push_back(max);
      d->filter_values[2].push_back(rrr);
      d->filter_values[3].push_back(v);
    }
  }

  getTriangleContainer(0)->allocate(
        Tri::Flat_vertices, d->positions.data(),
        static_cast<int>(d->positions.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::Flat_normals, d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));
  getTriangleContainer(0)->allocate(
        Tri::FColors, d->colors.data(),
        static_cast<int>(d->colors.size()*sizeof(float)));
  //Use the Distances vbo for the float filter values for convenience. They are probably not a distance, but the mechanism is already here so let's use it.
  setBuffersFilled(true);
  updateFilter();
}

Scene_item* Scene_tetrahedra_item::clone() const
{
  return new Scene_tetrahedra_item(d->c3t3_item);
}
void Scene_tetrahedra_item::compute_bbox() const
{
  if (isEmpty())
    setBbox(Bbox());
  else {
    bool bbox_init = false;
    CGAL::Bbox_3 result;
    const C3t3& c3t3 = d->c3t3_item->c3t3();
    for (Tr::Finite_vertices_iterator
         vit = c3t3.triangulation().finite_vertices_begin(),
         end = c3t3.triangulation().finite_vertices_end();
         vit != end; ++vit)
    {
      if(vit->in_dimension() == -1) continue;
      if (bbox_init)
        result = result + vit->point().bbox();
      else
      {
        result = vit->point().bbox();
        bbox_init = true;
      }
    }
    setBbox(Bbox(result.xmin(), result.ymin(), result.zmin(),
                 result.xmax(), result.ymax(), result.zmax()));
  }

}


void Scene_tetrahedra_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  getTriangleContainer(0)->initializeBuffers(viewer);
  getTriangleContainer(0)->setFlatDataSize(d->positions_size);

  d->positions.clear();
  d->positions.shrink_to_fit();
  d->normals.clear();
  d->normals.shrink_to_fit();
  d->colors.clear();
  d->colors.shrink_to_fit();
}

void Scene_tetrahedra_item::setMinThreshold(int i)
{
  d->threshold_1 = 0.001*i;
  updateThresholds();
}

void Scene_tetrahedra_item::setMaxThreshold(int i)
{
  d->threshold_2 = 0.001*i;
  updateThresholds();
}

void Scene_tetrahedra_item::setMinThreshold(void)
{
  double i = d->minEdit->text().toDouble();
  d->threshold_1 = 0.001*i;
  updateThresholds(false);
}

void Scene_tetrahedra_item::setMaxThreshold(void)
{
  double i = d->maxEdit->text().toDouble();
  d->threshold_2 = 0.001*i;
  updateThresholds(false);
}

void Scene_tetrahedra_item::setFilter(int i)
{
  d->filter_index = i;
  updateFilter();
  redraw();
}

void Scene_tetrahedra_item::updateFilter() const
{
  CGAL::Three::Vbo* vbo = getTriangleContainer(0)->getVbo(Tri::Distances);
  vbo->allocated = false;
  vbo->dataSize=0;
  getTriangleContainer(0)->allocate(Tri::Distances, d->filter_values[d->filter_index].data(),
      static_cast<int>(d->filter_values[d->filter_index].size()*sizeof(float)));
  d->minMinLabel->setText(QString("%1").arg(d->extrema[d->filter_index][0]));
  d->minMaxLabel->setText(QString("%1").arg(d->extrema[d->filter_index][1]));

  d->maxMinLabel->setText(QString("%1").arg(d->extrema[d->filter_index][0]));
  d->maxMaxLabel->setText(QString("%1").arg(d->extrema[d->filter_index][1]));

  double a = 1/(d->extrema[d->filter_index][1] - d->extrema[d->filter_index][0]);
  double b = - d->extrema[d->filter_index][0]/(d->extrema[d->filter_index][1] - d->extrema[d->filter_index][0]);

  d->minEdit->setValue((d->min_threshold - b)/a);
  d->maxEdit->setValue((d->max_threshold - b)/a);
  for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
  {
    setBuffersInit(static_cast<Vi*>(v), false);
  }
}

void Scene_tetrahedra_item::setMinMinLabelPointer(QLabel* ptr){ d->minMinLabel = ptr; }
void Scene_tetrahedra_item::setMinMaxLabelPointer(QLabel* ptr){ d->minMaxLabel = ptr; }
void Scene_tetrahedra_item::setMaxMinLabelPointer(QLabel* ptr){ d->maxMinLabel = ptr; }
void Scene_tetrahedra_item::setMaxMaxLabelPointer(QLabel* ptr){ d->maxMaxLabel = ptr; }
void Scene_tetrahedra_item::setMinEditPointer(DoubleEdit* ptr){ d->minEdit= ptr; }
void Scene_tetrahedra_item::setMaxEditPointer(DoubleEdit* ptr){ d->maxEdit = ptr; }

void Scene_tetrahedra_item::updateThresholds(bool update)
{
  d->max_threshold = (std::max)(d->threshold_1, d->threshold_2);
  d->min_threshold = (std::min)(d->threshold_1, d->threshold_2);
  if(update)
  {
    double a = 1/(d->extrema[d->filter_index][1] - d->extrema[d->filter_index][0]);
    double b = - d->extrema[d->filter_index][0]/(d->extrema[d->filter_index][1] - d->extrema[d->filter_index][0]);
    d->minEdit->setValue((d->min_threshold - b)/a);
    d->maxEdit->setValue((d->max_threshold - b)/a);
  }
  redraw();
}

Scene_c3t3_item* Scene_tetrahedra_item::c3t3_item()
{
  return d->c3t3_item;
}
