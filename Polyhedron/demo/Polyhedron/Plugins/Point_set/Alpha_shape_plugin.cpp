#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QDockWidget>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include "Scene_points_with_normal_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene.h"

#include "ui_Alpha_shape_widget.h"


typedef double coord_type;


typedef Kernel K;

typedef K::Point_3      Point_3;
typedef K::Vector_3     Vector_3;
typedef K::Segment_3    Segment_3;
typedef K::Ray_3        Ray_3;
typedef K::Line_3       Line;
typedef K::Triangle_3   Triangle_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb1;
typedef CGAL::Alpha_shape_vertex_base_3<K, Vb1> Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>   Fb;

typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds> Triangulation_3;

typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

typedef Alpha_shape_3::Cell  Cell;
typedef Alpha_shape_3::Vertex Vertex;
typedef Alpha_shape_3::Edge Edge;
typedef Alpha_shape_3::Facet Facet;
typedef Alpha_shape_3::Cell_handle  Cell_handle;
typedef Alpha_shape_3::Vertex_handle Vertex_handle;

typedef Alpha_shape_3::Cell_circulator  Cell_circulator;

typedef Alpha_shape_3::Locate_type Locate_type;

typedef Alpha_shape_3::Cell_iterator  Cell_iterator;
typedef Alpha_shape_3::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_3::Edge_iterator  Edge_iterator;


typedef Alpha_shape_3::Coord_type Coord_type;
typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;


class Scene_alpha_shape_item : public CGAL::Three::Scene_item
{
  Q_OBJECT
public :
  Scene_alpha_shape_item(Scene_points_with_normal_item* , int alpha);
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE {
    return (m == Flat);
  }
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawPoints(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void invalidateOpenGLBuffers()Q_DECL_OVERRIDE;
  void computeElements() const;
  Scene_item* clone() const Q_DECL_OVERRIDE{return 0;}
  QString toolTip() const Q_DECL_OVERRIDE{return QString();}
  bool isEmpty() const Q_DECL_OVERRIDE{ return false;}
  bool isFinite() const Q_DECL_OVERRIDE{ return true;}
  Alpha_shape_3 alpha_shape;
  void createPolygonSoup(std::vector<Kernel::Point_3>& points,
                         std::vector<std::vector<std::size_t> >& polys) const;
  std::size_t getNbofAlphas()const { return alpha_shape.number_of_alphas(); }
public Q_SLOTS:
  void alpha_changed(int);
private:
  mutable std::vector<float> vertices;
  mutable std::vector<unsigned int> indices;
  mutable QOpenGLShaderProgram *program;
  mutable QOpenGLShaderProgram facet_program;
  using CGAL::Three::Scene_item::initializeBuffers;
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;
  Point_set point_set;
}; //end of class Scene_alpha_shape_item



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
class DockWidget :
    public QDockWidget,
    public Ui::AlphaShapesWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};

using namespace CGAL::Three;
class Q_DECL_EXPORT Alpha_shape_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)Q_DECL_OVERRIDE {
    this->scene = scene_interface;
    this->mw = mainWindow;
    QAction* actionAlphaShapes= new QAction("Alpha Shapes", mw);
    if(actionAlphaShapes) {
      connect(actionAlphaShapes, &QAction::triggered,
              this, &Alpha_shape_plugin::on_actionAlphaShapes_triggered);
      _actions << actionAlphaShapes;
    }

    dock_widget = new DockWidget("Alpha Shapes", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    dock_widget->as_itemPushButton->setEnabled(false);
    dock_widget->poly_itemPushButton->setEnabled(false);
    addDockWidget(dock_widget);
    as_item = NULL;

    connect(dock_widget->as_itemPushButton, &QPushButton::clicked,
            this, &Alpha_shape_plugin::on_as_itemPushButton_clicked);
    connect(dock_widget->poly_itemPushButton, &QPushButton::clicked,
            this, &Alpha_shape_plugin::on_poly_itemPushButton_clicked);

    connect(dock_widget->horizontalSlider, SIGNAL(valueChanged(int)),
            dock_widget->spinBox, SLOT(setValue(int)));

    connect(dock_widget->spinBox, SIGNAL(valueChanged(int)),
            dock_widget->horizontalSlider, SLOT(setValue(int)));

    connect(static_cast<Scene*>(scene), SIGNAL(itemIndexSelected(int)),
            this,SLOT(on_scene_selection_changed(int)));


  }
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return qobject_cast<Scene_points_with_normal_item*>( scene->item( scene->mainSelectionIndex() ) );
  }
  QList<QAction*> actions() const Q_DECL_OVERRIDE {
    return _actions;
  }
public Q_SLOTS:

  void on_itemDestroyed()
  {
    as_item = NULL;
  }
  void on_actionAlphaShapes_triggered()
  {
    if(dock_widget->isVisible()) { dock_widget->hide(); }
    else {
      dock_widget->show();
      Scene_points_with_normal_item* sel_item =
          qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
      if(sel_item)
        on_as_itemPushButton_clicked();
    }
  }

  void on_as_itemPushButton_clicked()
  {
    Scene_points_with_normal_item* ps_item =
        qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!ps_item)
      return;
    //Only one alpha_shape item at the same time
    if(as_item)
    {
      scene->erase(scene->item_id(as_item));
      as_item = NULL;
    }
    ps_item->setVisible(false);
    as_item = new Scene_alpha_shape_item(ps_item, dock_widget->spinBox->value());

    as_item->setName(QString("%1 (alpha shape)").arg(ps_item->name()));
    as_item->setFlatMode();

    connect(dock_widget->horizontalSlider, &QSlider::valueChanged,
            as_item, &Scene_alpha_shape_item::alpha_changed);
    connect(as_item, &Scene_alpha_shape_item::aboutToBeDestroyed,
            this, &Alpha_shape_plugin::on_itemDestroyed);
    scene->setSelectedItem(scene->addItem(as_item));
  }

  void on_poly_itemPushButton_clicked()
  {
    Scene_alpha_shape_item* as_item =
        qobject_cast<Scene_alpha_shape_item*>(scene->item(scene->mainSelectionIndex()));
    if(!as_item)
      return;

    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t> > polys;
    as_item->createPolygonSoup(points, polys);
    Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
    new_item->init_polygon_soup(points.size(), polys.size());
    BOOST_FOREACH(const Kernel::Point_3& p, points)
    {
      new_item->new_vertex(p.x(), p.y(), p.z());
    }
    BOOST_FOREACH(const std::vector<std::size_t>& poly, polys)
    {
      new_item->new_triangle(poly[0], poly[1], poly[2]);
    }
    QString name = as_item->name().left(as_item->name().size() - QString("(alpha shape)").size());
    name.append("(polygon soup)");
    new_item->setName(name);
    as_item->setVisible(false);
    scene->addItem(new_item);
  }

  void on_scene_selection_changed(int i)
  {
    dock_widget->as_itemPushButton->setEnabled(false);
    dock_widget->poly_itemPushButton->setEnabled(false);

    if( qobject_cast<Scene_points_with_normal_item*>( scene->item(i) ) &&
        !as_item)
      dock_widget->as_itemPushButton->setEnabled(true);
    else if( qobject_cast<Scene_alpha_shape_item*>( scene->item(i) ) )
    {
      dock_widget->horizontalSlider->setMaximum(static_cast<int>(qobject_cast<Scene_alpha_shape_item*>( scene->item(i) )->getNbofAlphas()));
      dock_widget->spinBox->setMaximum(static_cast<int>(qobject_cast<Scene_alpha_shape_item*>( scene->item(i) )->getNbofAlphas()));
      dock_widget->poly_itemPushButton->setEnabled(true);
    }
  }
  void closure() Q_DECL_OVERRIDE
  {
    dock_widget->hide();
  }
private:
  QList<QAction*> _actions;
  DockWidget* dock_widget;
  Scene_alpha_shape_item* as_item;
}; //end of class Alpha_shape_plugin




/***********************************
 ***********************************
 ** Item's functions declarations **
 ***********************************
 ***********************************/
#include <CGAL/Timer.h>
Scene_alpha_shape_item::Scene_alpha_shape_item(Scene_points_with_normal_item *point_set_item, int alpha)
  : Scene_item(1, 1), point_set(*point_set_item->point_set())
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  _bbox = point_set_item->bbox();
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
  vertices.reserve(point_set.size() * 3);
  CGAL::Timer timer;
  timer.start();
  alpha_shape.make_alpha_shape(point_set.points().begin(), point_set.points().end());
  double res = timer.time();
  timer.stop();
  qDebug()<<"Alpha shape done : "<<res<<" s.";
  //set ids
  unsigned int i=0;
  for(Alpha_shape_3::Vertex_iterator it = alpha_shape.all_vertices_begin();
      it != alpha_shape.all_vertices_end();
      ++it)
  {
    it->info() = i++;
    vertices.push_back(it->point().x()+offset.x);
    vertices.push_back(it->point().y()+offset.y);
    vertices.push_back(it->point().z()+offset.z);
  }

  const char vertex_source[] =
  {
    "#version 120 \n"
    "attribute highp vec4 vertex;\n"
    "attribute highp vec3 colors;\n"
    "uniform highp mat4 mvp_matrix;\n"
    "uniform highp mat4 mv_matrix; \n"
    "varying highp vec4 fP; \n"
    "varying highp vec4 color;; \n"
    "void main(void)\n"
    "{\n"
    "   color = vec4(colors, 1.0); \n"
    "   fP = mv_matrix * vertex; \n"
    "   gl_Position = mvp_matrix * vertex;\n"
    "}"
  };

  facet_program.addShaderFromSourceCode(QOpenGLShader::Vertex, vertex_source);
  facet_program.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/cgal/Polyhedron_3/resources/shader_flat.f");
  facet_program.link();
  facet_program.bind();
  invalidateOpenGLBuffers();
  alpha_changed(alpha);
  QApplication::restoreOverrideCursor();
}
void Scene_alpha_shape_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
  {
    computeElements();
    initializeBuffers(viewer);
  }
  QMatrix4x4 mvp_mat;
  QMatrix4x4 mv_mat;
  GLdouble d_mat[16];
  viewer->camera()->getModelViewMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mv_mat.data()[i] = GLfloat(d_mat[i]);
  viewer->camera()->getModelViewProjectionMatrix(d_mat);
  for (int i=0; i<16; ++i)
    mvp_mat.data()[i] = GLfloat(d_mat[i]);
  QVector4D position(0.0f,0.0f,1.0f, 1.0f );
  QVector4D ambient(0.4f, 0.4f, 0.4f, 0.4f);
  // Diffuse
  QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
  // Specular
  QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);
  int is_both_sides;
  viewer->glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE,
                        &is_both_sides);
  program = &facet_program;
  program->bind();
  program->setUniformValue("mvp_matrix", mvp_mat);
  program->setUniformValue("mv_matrix", mv_mat);
  program->setUniformValue("light_pos", position);
  program->setUniformValue("light_diff",diffuse);
  program->setUniformValue("light_spec", specular);
  program->setUniformValue("light_amb", ambient);
  program->setUniformValue("spec_power", 51.8f);
  program->setUniformValue("is_two_side", is_both_sides);
  program->setUniformValue("is_selected", false);

  vaos[0]->bind();

  program->setAttributeValue("colors", this->color());
  viewer->glDrawElements(GL_TRIANGLES, static_cast<GLuint>(indices.size()),
                         GL_UNSIGNED_INT, indices.data());
  program->release();
  vaos[0]->release();

  drawPoints(viewer);
}

void Scene_alpha_shape_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
  {
    computeElements();
    initializeBuffers(viewer);
  }

  vaos[0]->bind();
  attribBuffers(viewer, PROGRAM_NO_SELECTION);
  program = getShaderProgram(PROGRAM_NO_SELECTION);
  program->bind();
  program->setAttributeValue("colors", QColor(255,0,0));
  viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(point_set.size()));
  program->release();
  vaos[0]->release();
}

void Scene_alpha_shape_item::invalidateOpenGLBuffers()
{
  are_buffers_filled = false;
}

void Scene_alpha_shape_item::computeElements() const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  indices.resize(0);

  std::list<Facet> facets;
  alpha_shape.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);

  for(std::list<Facet>::iterator fit = facets.begin();
      fit != facets.end();
      ++fit) {
    const Cell_handle& ch = fit->first;
    const int index = fit->second;

    const unsigned int& a = ch->vertex((index+1)&3)->info();
    const unsigned int& b = ch->vertex((index+2)&3)->info();
    const unsigned int& c = ch->vertex((index+3)&3)->info();

    indices.push_back(a); indices.push_back(b); indices.push_back(c);
  }
  QApplication::restoreOverrideCursor();
}

void Scene_alpha_shape_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
{
  program = &facet_program;
  program->bind();
  vaos[0]->bind();
  buffers[0].bind();
  buffers[0].allocate(vertices.data(),
                      static_cast<GLsizei>(vertices.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  buffers[0].release();
  vaos[0]->release();
  program->release();

  program = viewer->getShaderProgram(PROGRAM_NO_SELECTION);
  program->bind();
  vaos[0]->bind();
  buffers[0].bind();
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  buffers[0].release();
  vaos[0]->release();
  program->release();
  are_buffers_filled = true;
}

void Scene_alpha_shape_item::alpha_changed(int i)
{
  std::size_t n = i;
  if (alpha_shape.number_of_alphas() > 0
      && n > 0){
    if(n < alpha_shape.number_of_alphas()){
      alpha_shape.set_alpha(alpha_shape.get_nth_alpha(static_cast<int>(n)));
    } else {
      Alpha_iterator alpha_end_it = alpha_shape.alpha_end();
      alpha_shape.set_alpha((*(--alpha_end_it))+1);
    }
  } else {
    alpha_shape.set_alpha(0);
  }
  invalidateOpenGLBuffers();
  itemChanged();
}
void Scene_alpha_shape_item::createPolygonSoup(std::vector<Kernel::Point_3>&points,
                                               std::vector<std::vector<std::size_t> > &polys)const
{
  //fill points
  bool is_first = true;
  for(Alpha_shape_3::Vertex_iterator it = alpha_shape.all_vertices_begin();
      it != alpha_shape.all_vertices_end();
      ++it)
  {
    if(is_first)
      is_first = false;
    else
      points.push_back(it->point());
  }
  //fill polygons

  for(std::size_t i=0; i<indices.size(); i+=3)
  {
    std::vector<std::size_t> poly;
    poly.resize(3);
    for(int j=0; j<3; ++j)
    {
      poly[j] = indices[i+j] - 1;
    }
    polys.push_back(poly);
  }
}
#include "Alpha_shape_plugin.moc"
