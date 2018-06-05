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
#include <QInputDialog>


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
    dock_widget->hide();
  }
public Q_SLOTS:
  void pop();
  void addColor();
  void invalidateDisplay();
private:
  DockWidget* dock_widget;
  QAction* actionBubbles;
  mutable std::vector<SMesh::Vertex_index> vertices;
  std::vector<QColor> colors;
  Colored_spheres_item* spheres;
}; // end Colored_spheres_plugin


void Colored_spheres_plugin::init(QMainWindow *mainWindow, 
                                  Scene_interface *scene_interface, Messages_interface *)
{
  scene = scene_interface;
  mw = mainWindow;
  actionBubbles = new QAction(tr("Put Some Colors Inside Your Eyes"), mainWindow);
  connect(actionBubbles, SIGNAL(triggered()),
          this, SLOT(pop()));
  dock_widget = new DockWidget("Bubbles Party", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  addDockWidget(dock_widget);
  spheres = NULL;
}

void Colored_spheres_plugin::invalidateDisplay()
{
  QGridLayout* layout = dock_widget->gridLayout;
  //clear existing layout
  if ( layout != NULL )
  {
    QLayoutItem* item;
    while ( ( item = layout->takeAt( 0 ) ) != NULL )
    {
      delete item->widget();
      delete item;
    }
  }
  //fill it up again
  for(std::size_t i = 0; i < colors.size(); ++i)
  {
    if(colors[i].isValid())
    {
      QLabel* text = new QLabel(dock_widget);
      text->setText(QString("degree: %1").arg(i));
      QPushButton* color_button = new QPushButton(dock_widget);
      QPalette palette;
      palette.setColor(QPalette::Button,colors[i]);
      color_button->setPalette(palette);
      layout->addWidget(color_button, i, 0);
      layout->addWidget(text, i, 1);
    }
  }
  layout->setColumnStretch(0,0);
  layout->setColumnStretch(1,1);
  
  dock_widget->update();
}

void Colored_spheres_plugin::addColor()
{
  int new_deg = QInputDialog::getInt(mw, "Add a new degree", "Degree:", 3);
  if(new_deg < static_cast<int>(colors.size()))
  {
    if(colors[new_deg].isValid())
    {
      QMessageBox::warning(mw, "Warning", "This degree is already there.");
      return;
    }
  }
  else
  {
    colors.resize(new_deg+1, QColor());
  }
  QColor color = spheres->color();
  color = QColor::fromHsv((color.hue()+55*new_deg)%360, (new_deg*25)%155+100,color.lightness(), color.alpha());
  colors[new_deg] = color;
  invalidateDisplay();
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
  std::size_t max_deg = 0;
  BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
  {
    std::size_t deg = smitem->face_graph()->degree(vi);
    if(deg > max_deg)
      max_deg = deg;
  }
  colors.resize(max_deg + 1);
  BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
  {
    std::size_t deg = smitem->face_graph()->degree(vi);
    if(deg != 4)
    {
      c = spheres->color();
      c = QColor::fromHsv((c.hue()+ 55*deg)%360, (deg*25)%155+100,c.lightness(), c.alpha());
      colors[deg] = c;
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
  invalidateDisplay();
  connect(dock_widget->pushButton, &QPushButton::clicked,
          this, &Colored_spheres_plugin::addColor);
  dock_widget->show();
  scene->addItem(spheres);
  this->spheres = spheres;
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
