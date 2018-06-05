#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_spheres_item.h"

#include <CGAL/compute_average_spacing.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QSlider>
#include <QHBoxLayout>
#include <QMenu>
#include <QWidgetAction>
#include <QLabel>
#include <QPushButton>
#include <QColorDialog>

#include "ui_Color_spheres_widget.h"

using namespace CGAL::Three;
class Colored_spheres_item
    :public Scene_spheres_item
{
  Q_OBJECT
public :
  Colored_spheres_item(Scene_group_item *parent, double radius);
  ~Colored_spheres_item()
  {
    delete radius_slider;
  }
  // Indicates if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE {
    return (m == Wireframe || m == Gouraud);
  }
  Scene_item* clone() const Q_DECL_OVERRIDE {return 0;}
  QString toolTip() const Q_DECL_OVERRIDE {return QString();}
  QMenu* contextMenu() Q_DECL_OVERRIDE;   
  
public Q_SLOTS:
  void resize_spheres();
  
private:
  QSlider* radius_slider;
  double radius;
};


class DockWidget :
    public QDockWidget,
    public Ui::ColorSpheresWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};


class Colored_spheres_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionBubbles;
  }
  
  bool applicable(QAction*) const {
    Scene_surface_mesh_item* smitem = qobject_cast<Scene_surface_mesh_item*>(
          scene->item(scene->mainSelectionIndex()));
    if(smitem)
      return true;
    return false;
  }
  void closure()
  {
    //  dock_widget->hide();
  }
public Q_SLOTS:
  void pop();
private:
  DockWidget* dock_widget;
  QAction* actionBubbles;
  mutable std::vector<SMesh::Vertex_index> vertices;
  mutable std::set<QColor> colors;
  //ColoredSpheresWidget* dock_widget;
}; // end Colored_spheres_plugin


void Colored_spheres_plugin::init(QMainWindow *mainWindow, 
                                  Scene_interface *scene_interface, Messages_interface *)
{
  scene = scene_interface;
  mw = mainWindow;
  actionBubbles = new QAction(tr("Put Some Colors Inside Your Eyes"), mainWindow);
  connect(actionBubbles, SIGNAL(triggered()),
          this, SLOT(pop()));
  dock_widget = new ColoredSpheresWidget("Bubbles Party", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  addDockWidget(dock_widget);
}

void Colored_spheres_plugin::pop()
{
  Scene_surface_mesh_item* smitem = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->mainSelectionIndex()));
  if(!smitem)
    return;
  SMesh& mesh = *smitem->face_graph();
  double radius = 0.15 * CGAL::compute_average_spacing<CGAL::Sequential_tag>(mesh.points(), 3);
  Colored_spheres_item* spheres = new Colored_spheres_item(NULL, radius);
  QColor c;
  
  BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
  {
    std::size_t deg = smitem->face_graph()->degree(vi);
    if(deg != 4)
    {
      c = spheres->color();
      c = QColor::fromHsv((c.hue()+75*deg)%360, c.saturation(),c.lightness(), c.alpha());
      colors.insert(c);
      const CGAL::qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
      const EPICK::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);
      EPICK::Point_3 center(mesh.point(vi) + offset);
      
      typedef unsigned char UC;
      spheres->add_sphere(EPICK::Sphere_3(center, radius),
                          CGAL::Color(UC(c.red()), UC(c.green()), UC(c.blue())));
      vertices.push_back(vi);
    }
  }
  spheres->invalidateOpenGLBuffers();
  spheres->setName("Bubbles");
  spheres->setRenderingMode(Gouraud);
  //add colors to dock widget
  BOOST_FOREACH(QColor col, colors)
  {
    QLabel * text = new QLabel(QString("degree: "));
    QPushButton* color_button = new QPushButton();
    
  }
  dock_widget->show();
  scene->addItem(spheres);
}

Colored_spheres_item::Colored_spheres_item(Scene_group_item *parent, double radius)
  :Scene_spheres_item(parent),
    radius(radius)
{
  radius_slider = new QSlider(Qt::Horizontal);
  radius_slider->setMinimum(1);
  radius_slider->setValue(100);
  radius_slider->setMaximum(200);
}

QMenu* Colored_spheres_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_points_with_normal_item.";
  
  QMenu* menu = Scene_item::contextMenu();
  
  //add a slider to modify the normals length
  // Use dynamic properties:
  // http://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();
  
  if(!menuChanged) {
    QMenu *container = new QMenu(tr("Bubbles Size"));
    QWidgetAction *sliderAction = new QWidgetAction(0);
    connect(radius_slider, &QSlider::sliderReleased, this, &Colored_spheres_item::resize_spheres);
    
    sliderAction->setDefaultWidget(radius_slider);
    
    container->addAction(sliderAction);
    menu->addMenu(container);
    
    menu->setProperty(prop_name, true);
  }
  
  return menu;
}

void Colored_spheres_item::resize_spheres()
{
  for(std::size_t i = 0; i< numberOfSpheres(); ++i)
  {
    setRadius(radius * (radius_slider->value()/100.0),i);
    invalidateOpenGLBuffers();
    redraw();
  }
}

#include "Colored_spheres_plugin.moc"
