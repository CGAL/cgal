#include <QApplication>
#include "Scene_facegraph_transform_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Viewer_interface.h>

struct Scene_facegraph_transform_item_priv
{
  Scene_facegraph_transform_item_priv(const qglviewer::Vec& pos,FaceGraph* sm,
                                       const QString name, Scene_facegraph_transform_item *parent)
    : manipulable(false),
      frame(new CGAL::Three::Scene_item::ManipulatedFrame()),
      facegraph(sm),
      center_(pos),
      item_name(name)
  {
    item = parent;
    const qglviewer::Vec offset = static_cast<Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    frame->setPosition(pos+offset);
    nb_lines = 0;
  }
  ~Scene_facegraph_transform_item_priv()
{
  delete frame;
}
  void initialize_buffers(Viewer_interface *viewer) const;
  void compute_elements() const;
  enum VAOs {
      Edges=0,
      NbOfVaos
  };
  enum VBOs {
      Vertices = 0,
      NbOfVbos
  };

  bool manipulable;
  qglviewer::ManipulatedFrame* frame;
  FaceGraph* facegraph;
  qglviewer::Vec center_;
  Scene_facegraph_transform_item *item;
  QMatrix4x4 f_matrix;
  const QString item_name;

  mutable QOpenGLShaderProgram *program;
  mutable std::vector<float> positions_lines;
  mutable std::size_t nb_lines;
};

Scene_facegraph_transform_item::Scene_facegraph_transform_item(const qglviewer::Vec& pos, FaceGraph* sm,
                                                                 const QString name):
    Scene_item(Scene_facegraph_transform_item_priv::NbOfVbos,Scene_facegraph_transform_item_priv::NbOfVaos)
{
  d = new Scene_facegraph_transform_item_priv(pos,sm, name, this);
    invalidateOpenGLBuffers();
}


void Scene_facegraph_transform_item_priv::initialize_buffers(CGAL::Three::Viewer_interface *viewer =0) const
{
    //vao for the edges
    {
        program = item->getShaderProgram(Scene_facegraph_transform_item::PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        item->vaos[Edges]->bind();
        item->buffers[Vertices].bind();
        item->buffers[Vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        item->buffers[Vertices].release();
        item->vaos[Edges]->release();

        program->release();
    }
    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);

    item->are_buffers_filled = true;
}


void Scene_facegraph_transform_item_priv::compute_elements() const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_lines.resize(0);
    typedef Kernel::Point_3		        Point;
    typedef boost::graph_traits<FaceGraph>::edge_iterator	Edge_iterator;
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
    if(!are_buffers_filled)
        d->initialize_buffers(viewer);
    vaos[Scene_facegraph_transform_item_priv::Edges]->bind();
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer,PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    QColor color = this->color();
    d->program->setAttributeValue("colors",color);
    d->program->setUniformValue("f_matrix", d->f_matrix);
    d->program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_lines/3));
    vaos[Scene_facegraph_transform_item_priv::Edges]->release();
    d->program->release();

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
    qglviewer::Vec min(bbox.xmin(),bbox.ymin(),bbox.zmin());
    qglviewer::Vec max(bbox.xmax(),bbox.ymax(),bbox.zmax());
    _bbox = Bbox(min.x,min.y,min.z,
                 max.x,max.y,max.z);
}


void Scene_facegraph_transform_item::invalidateOpenGLBuffers()
{
    d->compute_elements();
    are_buffers_filled = false;
    compute_bbox();
}

bool Scene_facegraph_transform_item::manipulatable() const { return d->manipulable; }
CGAL::Three::Scene_item::ManipulatedFrame* Scene_facegraph_transform_item::manipulatedFrame() { return d->frame; }
void Scene_facegraph_transform_item::setManipulatable(bool b = true) { d->manipulable = b;}
const qglviewer::Vec& Scene_facegraph_transform_item::center() const { return d->center_; }
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
