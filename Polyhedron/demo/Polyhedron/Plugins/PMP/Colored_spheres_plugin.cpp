#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "SMesh_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_spheres_item.h"
#include "Scene_facegraph_item_k_ring_selection.h"

#include <CGAL/compute_average_spacing.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QMessageBox>
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
std::size_t INVALID_ID = static_cast<std::size_t>(-1);

//   I t e m - C l a s s
class Colored_spheres_item
    :public Scene_spheres_item
    
{
  Q_OBJECT
public :
  Colored_spheres_item(Scene_surface_mesh_item *parent, double radius);
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
  void select(double orig_x, double orig_y, double orig_z, 
              double dir_x, double dir_y, double dir_z) Q_DECL_OVERRIDE;
  
  
public Q_SLOTS:
  void resize_spheres();
  
private:
  QSlider* radius_slider;
  Scene_surface_mesh_item* parent;
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

//   P l u g i n - C l a s s
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
  void invalidateDisplay(bool update_spheres = false);
  void selected(const std::set<fg_vertex_descriptor> &m);
  void changeColor();
private:
  QColor changeDegree(size_t vid, bool &ok);
  DockWidget* dock_widget;
  QAction* actionBubbles;
  std::vector<std::size_t> vertex_sphere_map;
  std::vector<bool> vertex_has_sphere_map;
  fg_vertex_descriptor sel_v;
  fg_vertex_descriptor last_v;
  std::vector<QColor> colors;
  QColor selection_color;
  Colored_spheres_item* spheres;
  Scene_facegraph_item_k_ring_selection k_ring_selector;
  std::vector<std::size_t> vertex_deg_map;
}; // end Colored_spheres_plugin

//   P l u g i n - D e f i n i t i o n s
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
  selection_color = QColor(255,255,255);
}

void Colored_spheres_plugin::invalidateDisplay(bool update_spheres)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
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
  QLabel* text = new QLabel(dock_widget);
  text->setText(QString("selection"));
  QPushButton* color_button = new QPushButton(dock_widget);
  connect(color_button, &QPushButton::clicked,
          this, &Colored_spheres_plugin::changeColor);
  QPalette palette;
  palette.setColor(QPalette::Button,selection_color);
  color_button->setPalette(palette);
  layout->addWidget(color_button, 0, 0);
  layout->addWidget(text, 0, 1);
  for(std::size_t i = 0; i < colors.size(); ++i)
  {
    if(colors[i].isValid())
    {
      QLabel* text = new QLabel(dock_widget);
      text->setText(QString("degree: %1").arg(i));
      color_button = new QPushButton(dock_widget);
      connect(color_button, &QPushButton::clicked,
              this, &Colored_spheres_plugin::changeColor);
      QPalette palette;
      palette.setColor(QPalette::Button,colors[i]);
      color_button->setPalette(palette);
      layout->addWidget(color_button, i+1, 0);
      layout->addWidget(text, i+1, 1);
    }
  }
  layout->setColumnStretch(0,0);
  layout->setColumnStretch(1,1);
  if(update_spheres)
  {
    SMesh& mesh = *k_ring_selector.poly_item->face_graph();
    BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
    {
      if(vertex_has_sphere_map[std::size_t(vi)])
      {
        std::size_t deg = vertex_deg_map[std::size_t(vi)];
        QColor c = colors[deg];
        typedef unsigned char UC;
        spheres->setColor(CGAL::Color(UC(c.red()), 
                                      UC(c.green()), 
                                      UC(c.blue())), 
                          vertex_sphere_map[std::size_t(vi)]);
      }
    }
    if(std::size_t(sel_v) < vertex_has_sphere_map.size()
       && vertex_has_sphere_map[std::size_t(sel_v)])
    {
      typedef unsigned char UC;
      spheres->setColor(CGAL::Color(UC(selection_color.red()), 
                                    UC(selection_color.green()), 
                                    UC(selection_color.blue())), 
                        vertex_sphere_map[std::size_t(sel_v)]);
    }   
    spheres->invalidateOpenGLBuffers();
    spheres->redraw();
  }
  
  dock_widget->update();
  QApplication::restoreOverrideCursor();
}

QColor Colored_spheres_plugin::changeDegree(std::size_t vid, bool& ok)
{
  int new_deg = QInputDialog::getInt(mw, "Add a new degree", "Degree:", 3, 1, 100, 1, &ok);
  if(!ok)
    return QColor();
  if(new_deg >= static_cast<int>(colors.size()))
  {
    colors.resize(new_deg+1, QColor());
  }
  QColor color = spheres->color();
  if(colors[new_deg].isValid())
  {
    color = colors[new_deg];
  }
  else
  {
    color = QColor::fromHsv((color.hue()+55*new_deg)%360, (new_deg*25)%155+100,color.lightness(), color.alpha());
    colors[new_deg] = color;
    color = colors[new_deg];
  }
  vertex_deg_map[vid]=new_deg;
  invalidateDisplay();
  ok = true;
  return color;
}

void Colored_spheres_plugin::pop()
{
  Scene_surface_mesh_item* smitem = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->mainSelectionIndex()));
  if(!smitem)
    return;
  SMesh& mesh = *smitem->face_graph();
  vertex_sphere_map.resize(vertices(mesh).size());
  vertex_deg_map.resize(vertex_sphere_map.size());
  vertex_has_sphere_map.resize(vertex_sphere_map.size(), false);
  double radius = 0.3 * CGAL::compute_average_spacing<CGAL::Sequential_tag>(mesh.points(), 3);
  Colored_spheres_item* spheres = new Colored_spheres_item(smitem, radius);
  QColor c;
  std::size_t max_deg = 0;
  BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
  {
    std::size_t deg = smitem->face_graph()->degree(vi);
    vertex_deg_map[std::size_t(vi)] = deg;
    if(deg > max_deg)
      max_deg = deg;
  }
  colors.resize(max_deg + 1);
  BOOST_FOREACH(SMesh::Vertex_index vi, mesh.vertices())
  {
    std::size_t deg = vertex_deg_map[std::size_t(vi)];
    if(deg != 4)
    {
      c = spheres->color();
      c = QColor::fromHsv((c.hue()+ 55*deg)%360, (deg*25)%155+100,c.lightness(), c.alpha());
      colors[deg] = c;
      vertex_sphere_map[std::size_t(vi)] = spheres->numberOfSpheres();
      vertex_has_sphere_map[std::size_t(vi)] = true;
      const CGAL::qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
      const EPICK::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);
      EPICK::Point_3 center(mesh.point(vi) + offset);
      
      typedef unsigned char UC;
      spheres->add_sphere(EPICK::Sphere_3(center, radius),
                          CGAL::Color(UC(c.red()), UC(c.green()), UC(c.blue())));
    }
  }
  spheres->invalidateOpenGLBuffers();
  spheres->setName("Bubbles");
  spheres->setRenderingMode(Gouraud);
  //add colors to dock widget
  invalidateDisplay();
  dock_widget->show();
  scene->addItem(spheres);
  this->spheres = spheres;
  
  connect(&k_ring_selector, SIGNAL(selected(const std::set<fg_vertex_descriptor>&)), this,
          SLOT(selected(const std::set<fg_vertex_descriptor>&)));
  k_ring_selector.init(smitem, mw, Scene_facegraph_item_k_ring_selection::Active_handle::VERTEX, -1);
  k_ring_selector.k_ring = 0;
}

void Colored_spheres_plugin::selected(const std::set<fg_vertex_descriptor> &m)
{
  sel_v = *m.begin();
  static bool is_one_currently_selected = false;
  static std::size_t id = INVALID_ID;
  if(!is_one_currently_selected)
  {
    if(!vertex_has_sphere_map[std::size_t(sel_v)])
      return;
    else
      id = vertex_sphere_map[std::size_t(sel_v)];
    typedef unsigned char UC;
    spheres->setColor(CGAL::Color(UC(selection_color.red()), 
                                  UC(selection_color.green()), 
                                  UC(selection_color.blue())), 
                      id);
    last_v = sel_v;
    is_one_currently_selected = true;
  }
  else
  {
    if(sel_v != last_v 
       && vertex_has_sphere_map[std::size_t(sel_v)])
    {
      QMessageBox::warning(mw, "Warning", "Cannot overwrite an existing sphere.");
      return;
    }
    bool ok;
    QColor c = changeDegree(std::size_t(sel_v), ok);
    if(!ok)
      return;
    spheres->setCenter(k_ring_selector.poly_item->polyhedron()->point(SMesh::Vertex_index(sel_v)), id);
    vertex_sphere_map[std::size_t(last_v)] = INVALID_ID;
    vertex_has_sphere_map[std::size_t(last_v)] = false;
    vertex_sphere_map[std::size_t(sel_v)] = id;
    vertex_has_sphere_map[std::size_t(sel_v)] = true;
    typedef unsigned char UC;
    spheres->setColor(CGAL::Color(UC(c.red()),UC(c.green()),UC(c.blue())), id);
    is_one_currently_selected = false;
    sel_v = SMesh::Vertex_index(INVALID_ID);
  }
  spheres->invalidateOpenGLBuffers();
  spheres->redraw();
}
void Colored_spheres_plugin::changeColor()
{
  QPushButton* sender_button = qobject_cast<QPushButton*>(sender());
  QGridLayout* layout = dock_widget->gridLayout;
  QWidget* button;
  bool found(false);
  int r;
  for(r = 0; r < layout->rowCount(); ++r)
  {
    QLayoutItem* litem = layout->itemAtPosition(r,0);
    if(litem){
      button = litem->widget();
      if( button == sender_button)
      {
        found = true;
        break;
      }
    }
  }
  if(!found)
    return;
  //case selection color
  QColor new_color = QColorDialog::getColor();
  if(!new_color.isValid())
    return;
  if(r == 0)
  {
    selection_color = new_color;
  }
  else
  {
    QString text = qobject_cast<QLabel*>(layout->itemAtPosition(r,1)->widget())->text();
    QRegExp rx("(\\d+)");
    QString res;
    int start = text.indexOf(rx);
    for(; start<text.length(); ++start)
    {
      res.append(text.at(start));
    }
    int deg = res.toInt();
    colors[deg] = new_color;
  }
  invalidateDisplay(true);
}

//   I t e m - D e f i n i t i o n s
Colored_spheres_item::Colored_spheres_item(Scene_surface_mesh_item *parent, double radius)
  :Scene_spheres_item(NULL),
    parent(parent),
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

void Colored_spheres_item::select(double orig_x, double orig_y, double orig_z, double dir_x, double dir_y, double dir_z)
{
  parent->select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z); 
}
#include "Colored_spheres_plugin.moc"
