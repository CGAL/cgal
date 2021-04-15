#include <QApplication>
#include "Scene_facegraph_transform_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Edge_container.h>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Edge_container Ec;

struct Scene_facegraph_transform_item_priv
{
  Scene_facegraph_transform_item_priv(const CGAL::qglviewer::Vec& pos,FaceGraph* sm,
                                       const QString name, Scene_facegraph_transform_item *parent)
    : manipulable(false),
      frame(new CGAL::Three::Scene_item::ManipulatedFrame()),
      facegraph(sm),
      center_(pos),
      item_name(name)
  {
    item = parent;
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    frame->setPosition(pos+offset);
    nb_lines = 0;
  }
  ~Scene_facegraph_transform_item_priv()
{
  delete frame;
}
  void compute_elements() const;

  bool manipulable;
  CGAL::qglviewer::ManipulatedFrame* frame;
  FaceGraph* facegraph;
  CGAL::qglviewer::Vec center_;
  Scene_facegraph_transform_item *item;
  QMatrix4x4 f_matrix;
  const QString item_name;

  mutable std::vector<float> positions_lines;
  mutable std::size_t nb_lines;
};

Scene_facegraph_transform_item::Scene_facegraph_transform_item(const CGAL::qglviewer::Vec& pos, FaceGraph* sm,
                                                                 const QString name)
{
  d = new Scene_facegraph_transform_item_priv(pos,sm, name, this);
  setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  invalidateOpenGLBuffers();
}


void Scene_facegraph_transform_item_priv::compute_elements() const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_lines.resize(0);
    typedef Kernel::Point_3                        Point;
    typedef boost::graph_traits<FaceGraph>::edge_iterator        Edge_iterator;
    typedef boost::property_map<FaceGraph, CGAL::vertex_point_t>::type VPmap;
    VPmap vpmap = get(CGAL::vertex_point, *facegraph);
    Edge_iterator he;
    for(he = edges(*facegraph).begin();
        he != edges(*facegraph).end();
        he++)
    {
        const Point& a = get(vpmap, target(halfedge(*he, *facegraph), *facegraph));
        const Point& b = get(vpmap, target(opposite(halfedge(*he, *facegraph), *facegraph), *facegraph));
        positions_lines.push_back(a.x()-center_.x);
        positions_lines.push_back(a.y()-center_.y);
        positions_lines.push_back(a.z()-center_.z);

        positions_lines.push_back(b.x()-center_.x);
        positions_lines.push_back(b.y()-center_.y);
        positions_lines.push_back(b.z()-center_.z);

    }
    QApplication::restoreOverrideCursor();
}

void Scene_facegraph_transform_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
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
    Ec* ec = getEdgeContainer(0);
    ec->setColor(this->color());
    ec->setFrameMatrix(d->f_matrix);
    ec->draw(viewer, true);
}

QString Scene_facegraph_transform_item::toolTip() const {
    return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                       "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                       "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
            .arg(d->item_name);
}
bool Scene_facegraph_transform_item::keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_S){
    Q_EMIT stop();
        return true;
    }
    return false;
}

void
Scene_facegraph_transform_item::compute_bbox() const {
  typedef boost::property_map<FaceGraph,CGAL::vertex_point_t>::type VPmap;
  VPmap vpmap = get(CGAL::vertex_point,*d->facegraph);
    const Kernel::Point_3& p = get(vpmap, *vertices(*d->facegraph).begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(boost::graph_traits<FaceGraph>::vertex_iterator it = vertices(*d->facegraph).begin();
        it != vertices(*d->facegraph).end();
        ++it) {
      bbox = bbox + get(vpmap, *it).bbox();
    }
    CGAL::qglviewer::Vec vmin(bbox.xmin(),bbox.ymin(),bbox.zmin());
    CGAL::qglviewer::Vec vmax(bbox.xmax(),bbox.ymax(),bbox.zmax());
    setBbox(Bbox(vmin.x,vmin.y,vmin.z,
                 vmax.x,vmax.y,vmax.z));
}


void Scene_facegraph_transform_item::invalidateOpenGLBuffers()
{
  compute_bbox();
  setBuffersFilled(false);
  getEdgeContainer(0)->reset_vbos(ALL);
}

bool Scene_facegraph_transform_item::manipulatable() const { return d->manipulable; }
CGAL::Three::Scene_item::ManipulatedFrame* Scene_facegraph_transform_item::manipulatedFrame() { return d->frame; }
void Scene_facegraph_transform_item::setManipulatable(bool b = true) { d->manipulable = b;}
const CGAL::qglviewer::Vec& Scene_facegraph_transform_item::center() const { return d->center_; }
Scene_facegraph_transform_item::~Scene_facegraph_transform_item() { delete d; Q_EMIT killed(); }
void Scene_facegraph_transform_item::setFMatrix(double matrix[16])
{
  for (int i=0; i<16; ++i)
    d->f_matrix.data()[i] = (float)matrix[i];
}

FaceGraph *Scene_facegraph_transform_item::getFaceGraph()
{
 return d->facegraph;
}

void Scene_facegraph_transform_item::computeElements() const
{
  d->compute_elements();

  getEdgeContainer(0)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));
  d->nb_lines = d->positions_lines.size();
  setBuffersFilled(true);
}
void Scene_facegraph_transform_item::initializeBuffers(Viewer_interface *v) const
{
  getEdgeContainer(0)->initializeBuffers(v);
  getEdgeContainer(0)->setFlatDataSize(d->nb_lines);
  d->positions_lines.clear();
  d->positions_lines.shrink_to_fit();
}
