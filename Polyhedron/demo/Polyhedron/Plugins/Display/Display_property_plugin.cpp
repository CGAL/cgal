
#include "ui_Display_property.h"

#include "Color_map.h"
#include "Color_ramp.h"
#include "id_printing.h"
#include "Messages_interface.h"

#include "Scene.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "triangulate_primitive.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Buffer_for_vao.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <QAbstractItemView>
#include <QAction>
#include <QApplication>
#include <QColor>
#include <QColorDialog>
#include <QInputDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <QSlider>
#include <QObject>
#include <QPalette>
#include <QStyleFactory>

#include <boost/algorithm/clamp.hpp>
#include <boost/range/value_type.hpp>

#include <type_traits>
#include <unordered_map>
#include <vector>

#define ARBITRARY_DBL_MIN 1.0E-17
#define ARBITRARY_DBL_MAX 1.0E+17

namespace PMP = CGAL::Polygon_mesh_processing;

using namespace CGAL::Three;

Viewer_interface* (&getActiveViewer)() = Three::activeViewer;

class DockWidget
  : public QDockWidget,
    public Ui::DisplayPropertyWidget
{
public:
  DockWidget(const QString& name, QWidget *parent)
    : QDockWidget(name, parent)
  {
    setupUi(this);
  }
};

class Display_property_plugin
  : public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionDisplayProperties;

  DockWidget* dock_widget;

  // coloring choice and legend
  double rm = 1.;
  double gm = 0.;
  double bm = 0.;
  double rM = 0.;
  double gM = 1.;
  double bM = 0.;
  double rI = 0.;
  double gI = 1.;
  double bI = 0.;

  double expand_radius = 0.;
  double maxEdgeLength = -1.;

  Color_ramp color_ramp;
  std::vector<QColor> color_map;
  QPixmap legend;

  // tracks whether property 'i' (aka i-th in propertyBox) applies to vertices or faces
  enum class Property_simplex_type { VERTEX, FACE };
  std::vector<Property_simplex_type> property_simplex_types;

  enum Extremum
  {
    MIN_VALUE,
    MAX_VALUE
  };
  enum CurvatureType
  {
    MEAN_CURVATURE,
    GAUSSIAN_CURVATURE,
  };

public:
  bool applicable(QAction*) const override
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return false;

    return qobject_cast<Scene_surface_mesh_item*>(item) ||
           qobject_cast<Scene_points_with_normal_item*>(item);
  }

  QList<QAction*> actions() const override
  {
    return QList<QAction*>() << actionDisplayProperties;
  }

  void init(QMainWindow* mw,
            Scene_interface* sc,
            Messages_interface*) override
  {
    this->scene = sc;
    this->mw = mw;

    // Main action
    actionDisplayProperties = new QAction(QString("Display Properties"), mw);
    actionDisplayProperties->setProperty("submenuName", "Color");

    connect(actionDisplayProperties, SIGNAL(triggered()),
            this, SLOT(openDialog()));

    Scene* scene_item = static_cast<Scene*>(scene);
    connect(scene_item, SIGNAL(itemIndexSelected(int)),
            this, SLOT(onItemIndexSelected(int)));

    // Dock Widget
    dock_widget = new DockWidget("Property Display", mw);
    addDockWidget(dock_widget);

    dock_widget->setVisible(false);
    dock_widget->setEnabled(false);

    connect(dock_widget->propertyBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(onNewPropertySelected(int)));

    QPalette palette(Qt::red);
    dock_widget->minColorButton->setPalette(palette);
    dock_widget->minColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->minColorButton->update();

    palette = QPalette(Qt::green);
    dock_widget->maxColorButton->setPalette(palette);
    dock_widget->maxColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->maxColorButton->update();

    palette = QPalette(Qt::green);
    dock_widget->initColorButton->setPalette(palette);
    dock_widget->initColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->initColorButton->update();

    displayRampLegend();

    // lambda to generate the three connect(), for each color button (min, max, map)
    auto connect_color_buttons = [this](QPushButton* colorButton,
                                        double& r, double& g, double& b)
                                 {
                                    connect(colorButton, &QPushButton::pressed,
                                            this, [this, colorButton, &r, &g, &b]()
                                                  {
                                                    QColor color = QColorDialog::getColor();
                                                    if(!color.isValid())
                                                      return;

                                                    r = color.redF();
                                                    g = color.greenF();
                                                    b = color.blueF();

                                                    QPalette palette(color);
                                                    colorButton->setPalette(palette);
                                                    colorButton->update();

                                                    displayLegend();
                                                  });
                                 };

    connect_color_buttons(dock_widget->minColorButton, rm, gm, bm);
    connect_color_buttons(dock_widget->maxColorButton, rM, gM, bM);
    connect_color_buttons(dock_widget->initColorButton, rI, gI, bI);

    displayLegend();

    connect(dock_widget->colorizeButton, SIGNAL(clicked(bool)),
            this, SLOT(colorize()));

    connect(dock_widget->zoomToMinButton, &QPushButton::pressed,
            this, &Display_property_plugin::on_zoomToMinButton_pressed);
    connect(dock_widget->zoomToMaxButton, &QPushButton::pressed,
            this, &Display_property_plugin::on_zoomToMaxButton_pressed);

    connect(dock_widget->expandingRadiusSlider, &QSlider::valueChanged,
            this, &Display_property_plugin::setExpandingRadius);
  }

private Q_SLOTS:
  void openDialog()
  {
    if(!dock_widget->isVisible())
      dock_widget->show();
    dock_widget->raise();
  }

  void closure() override
  {
    dock_widget->hide();
  }

private:
    void disableExtremeValues()
  {
    dock_widget->extremeValuesGroup->setEnabled(false);
    dock_widget->zoomToMinButton->setEnabled(false);
    dock_widget->zoomToMaxButton->setEnabled(false);
  }

  void enableExtremeValues()
  {
    dock_widget->extremeValuesGroup->setEnabled(true);
    dock_widget->zoomToMinButton->setEnabled(true);
    dock_widget->zoomToMaxButton->setEnabled(true);
  }

  void resetExtremeValues()
  {
    // setup some dummy values such that the legend can be displayed
    const std::string& property_name = dock_widget->propertyBox->currentText().toStdString();

    if(property_name == "Smallest Angle Per Face" || property_name == "Largest Angle Per Face")
    {
      dock_widget->minBox->setRange(0, 360);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(0, 360);
      dock_widget->maxBox->setValue(0);
    }
    else if(property_name == "Scaled Jacobian")
    {
      dock_widget->minBox->setRange(-1000, 1000);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(-1000, 1000);
      dock_widget->maxBox->setValue(0);
    }
    else if(property_name == "Face Area")
    {
      dock_widget->minBox->setRange(-1000, 1000);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(-1000, 1000);
      dock_widget->maxBox->setValue(0);
    }
    else if (property_name == "Interpolated Corrected Mean Curvature")
    {
      dock_widget->minBox->setRange(-99999999, 99999999);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(-99999999, 99999999);
      dock_widget->maxBox->setValue(0);
    }
    else if (property_name == "Interpolated Corrected Gaussian Curvature")
    {
      dock_widget->minBox->setRange(-99999999, 99999999);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(-99999999, 99999999);
      dock_widget->maxBox->setValue(0);
    }
    else
    {
      dock_widget->minBox->setRange(-99999999, 99999999);
      dock_widget->minBox->setValue(0);
      dock_widget->maxBox->setRange(-99999999, 99999999);
      dock_widget->maxBox->setValue(0);
    }
  }

  void displayRampLegend()
  {
    color_ramp = Color_ramp(rm, rM, gm, gM, bm, bM);

    const int height = 256;
    const int width = 140;
    const int cell_width = width / 3;
    const int top_margin = 5;
    const int left_margin = 5;
    const int drawing_height = height - 2*top_margin;
    const int text_height = 20;

    legend = QPixmap(width, height + text_height);
    legend.fill(QColor(200, 200, 200));

    QPainter painter(&legend);
    painter.setPen(Qt::black);
    painter.setBrush(QColor(200, 200, 200));

    const double min_value = dock_widget->minBox->value();
    const double max_value = dock_widget->maxBox->value();

    // Build legend data
    std::vector<double> graduations(100);
    for(int i=0; i<100; ++i)
      graduations[i] = i / 100.0;

    int i = 0;
    for(std::vector<double>::iterator it = graduations.begin(), end = graduations.end(); it != end; ++it, i+=2)
    {
      QColor color(255 * color_ramp.r(*it),
                   255 * color_ramp.g(*it),
                   255 * color_ramp.b(*it));
      painter.fillRect(left_margin, drawing_height - top_margin - i, cell_width, 2, color);
    }

    // draw right vertical line
    painter.setPen(Qt::blue);
    painter.drawLine(QPoint(left_margin + cell_width + 10,
                            drawing_height - top_margin + 2),
                     QPoint(left_margin + cell_width + 10,
                            drawing_height - top_margin - static_cast<int>(graduations.size())*2 + 2));

    // draw min value and max value
    painter.setPen(Qt::blue);
    QRect min_text_rect(left_margin + cell_width + 10,
                        drawing_height - top_margin, 100, text_height);
    painter.drawText(min_text_rect, Qt::AlignCenter, QObject::tr("%1").arg(min_value, 0, 'f', 3));

    QRect max_text_rect(left_margin + cell_width + 10,
                        drawing_height - top_margin - 200, 100, text_height);
    painter.drawText(max_text_rect, Qt::AlignCenter, QObject::tr("%1").arg(max_value, 0, 'f', 3));

    dock_widget->legendLabel->setPixmap(legend);
  }

  template<typename ValueType>
  void displayMapLegend(const std::vector<ValueType>& values)
  {
    const std::size_t size = (std::min)(color_map.size(), std::size_t(1024));

    const int text_height = 20;
    const int height = text_height * static_cast<int>(size) + text_height;
    const int width = 140;
    const int cell_width = width / 3;
    const int top_margin = 15;
    const int left_margin = 5;
    const int drawing_height = height - text_height + top_margin;

    legend = QPixmap(width, height);
    legend.fill(QColor(200, 200, 200));

    QPainter painter(&legend);
    painter.setPen(Qt::black);
    painter.setBrush(QColor(200, 200, 200));

    int j = 0;
    int tick_height = text_height;
    for(std::size_t i=0; i<size; ++i, j+=tick_height)
    {
      QColor color(color_map[i].red(),
                   color_map[i].green(),
                   color_map[i].blue());

      painter.fillRect(left_margin,
                       drawing_height - top_margin - j,
                       cell_width,
                       tick_height,
                       color);

      QRect text_rect(left_margin + cell_width + 10, drawing_height - top_margin - j, 50, text_height);
      painter.drawText(text_rect, Qt::AlignCenter, QObject::tr("%1").arg(values[i], 0, 'f', 3, QLatin1Char(' ')));
    }

    if(color_map.size() > size)
    {
      QRect text_rect(left_margin + cell_width + 10, 0, 50, text_height);
      painter.drawText(text_rect, Qt::AlignCenter, QObject::tr("[...]"));
    }

    // draw right vertical line
    painter.setPen(Qt::blue);
    painter.drawLine(QPoint(left_margin + cell_width + 10,
                            drawing_height - top_margin + tick_height),
                     QPoint(left_margin + cell_width + 10,
                            drawing_height - top_margin  - static_cast<int>(size)*tick_height + tick_height));

    dock_widget->legendLabel->setPixmap(legend);
  }

  void displayLegend()
  {
    if(dock_widget->colorRampRadioButton->isChecked())
      displayRampLegend();
    else
      color_map.clear();
  }

private:
  template<typename Simplex>
  bool isSMPropertyScalar(const std::string& name,
                          const SMesh& mesh) const;

  void detectSMScalarProperties(SMesh& mesh)
  {
    std::vector<std::string> vprop = mesh.properties<vertex_descriptor>();
    for(const std::string& s : vprop)
    {
      if(isSMPropertyScalar<vertex_descriptor>(s, mesh))
      {
        dock_widget->propertyBox->addItem(s.c_str());
        property_simplex_types.push_back(Property_simplex_type::VERTEX);
      }
    }

    std::vector<std::string> fprop = mesh.properties<face_descriptor>();
    for(const std::string& s : fprop)
    {
      if(isSMPropertyScalar<face_descriptor>(s, mesh))
      {
        dock_widget->propertyBox->addItem(s.c_str());
        property_simplex_types.push_back(Property_simplex_type::FACE);
      }
    }
  }

  bool isPSPropertyScalar(const std::string& name,
                          const Point_set& ps) const;

  void detectPSScalarProperties(const Point_set& ps)
  {
    for(const auto& s : ps.properties())
      if(isPSPropertyScalar(s, ps))
        dock_widget->propertyBox->addItem(s.c_str());
  }

private:
  // This fills "dock_widget->propertyBox" with the properties that can be displayed for the item at position 'item_index'
  void detectScalarProperties(int item_index)
  {
    dock_widget->propertyBox->clear(); // calls onNewPropertySelected(-1)
    property_simplex_types.clear();

    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(item_index));
    Scene_points_with_normal_item* ps_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(item_index));

    if(sm_item)
    {
      dock_widget->propertyBox->addItems({"Smallest Angle Per Face",
                                          "Largest Angle Per Face",
                                          "Scaled Jacobian",
                                          "Face Area",
                                          "Interpolated Corrected Mean Curvature",
                                          "Interpolated Corrected Gaussian Curvature"});
      property_simplex_types = { Property_simplex_type::FACE,
                                 Property_simplex_type::FACE,
                                 Property_simplex_type::FACE,
                                 Property_simplex_type::FACE,
                                 Property_simplex_type::VERTEX,
                                 Property_simplex_type::VERTEX };
      detectSMScalarProperties(*(sm_item->face_graph()));
    }
    else if(ps_item)
    {
      detectPSScalarProperties(*(ps_item->point_set()));
    }

    int width = dock_widget->propertyBox->minimumSizeHint().width();
    dock_widget->propertyBox->view()->setMinimumWidth(width);
  }

private Q_SLOTS:
  // Called when a new geometric object is selected in the scene
  void onItemIndexSelected(int item_index)
  {
    // try to keep the same selected property, if possible
    const QString selected_property = dock_widget->propertyBox->currentText();

    detectScalarProperties(item_index);

    if(dock_widget->propertyBox->count() == 0)
      return;

    const int property_index = dock_widget->propertyBox->findText(selected_property);
    if(property_index == -1)
      dock_widget->propertyBox->setCurrentIndex(0);
    else
      dock_widget->propertyBox->setCurrentIndex(property_index);
  }

  // Called when a new property is selected in the combo box.
  // This function is only called if the index actually changed (doesn't trigger
  // if you click again the item)
  void onNewPropertySelected(int property_index)
  {
    resetExtremeValues(); // reset extreme value before the legend to get the proper values
    displayLegend();

    if(property_index >= 0 && property_index < dock_widget->propertyBox->count()) // valid property
    {
      dock_widget->setEnabled(true);
      disableExtremeValues(); // only available after coloring

      // Curvature property-specific slider
      const std::string& property_name = dock_widget->propertyBox->currentText().toStdString();
      const bool is_curvature_property = (property_name == "Interpolated Corrected Mean Curvature" ||
                                          property_name == "Interpolated Corrected Gaussian Curvature");
      dock_widget->expandingRadiusLabel->setVisible(is_curvature_property);
      dock_widget->expandingRadiusSlider->setVisible(is_curvature_property);
      dock_widget->expandingRadiusLabel->setEnabled(is_curvature_property);
      dock_widget->expandingRadiusSlider->setEnabled(is_curvature_property);
    }
    else // no or broken property
    {
      dock_widget->setEnabled(false);
      dock_widget->propertyBox->setEnabled(true);
    }
  }

private:
  void colorizePS(Scene_points_with_normal_item* ps_item)
  {
    ps_item->point_set()->add_colors();
    if(!displayPSProperty(dock_widget->propertyBox->currentText().toStdString(),
                          *(ps_item->point_set())))
    {
      return;
    }

    ps_item->invalidateOpenGLBuffers();
    ps_item->setRenderingMode(Points);
    ps_item->redraw();
  }

  void colorizeSM(Scene_surface_mesh_item* sm_item)
  {
    CGAL_assertion(static_cast<std::size_t>(dock_widget->propertyBox->count()) == property_simplex_types.size());

    // leave it flat if it was, otherwise set to flat+edges
    if(sm_item->renderingMode() != Flat && sm_item->renderingMode() != FlatPlusEdges)
      sm_item->setRenderingMode(FlatPlusEdges);

    const std::string& property_name = dock_widget->propertyBox->currentText().toStdString();
    if(property_name == "Smallest Angle Per Face")
    {
      displayExtremumAnglePerFace(sm_item, MIN_VALUE);
    }
    else if(property_name == "Largest Angle Per Face")
    {
      displayExtremumAnglePerFace(sm_item, MAX_VALUE);
    }
    else if(property_name == "Scaled Jacobian")
    {
      displayScaledJacobian(sm_item);
    }
    else if(property_name == "Face Area")
    {
      displayArea(sm_item);
    }
    else if(property_name == "Interpolated Corrected Mean Curvature")
    {
      displayInterpolatedCurvatureMeasure(sm_item, MEAN_CURVATURE);
      sm_item->setRenderingMode(Gouraud);
    }
    else if(property_name == "Interpolated Corrected Gaussian Curvature")
    {
      displayInterpolatedCurvatureMeasure(sm_item, GAUSSIAN_CURVATURE);
      sm_item->setRenderingMode(Gouraud);
    }
    else
    {
      const int property_index = dock_widget->propertyBox->currentIndex();
      if(property_simplex_types.at(property_index) == Property_simplex_type::VERTEX)
      {
        if(!displaySMProperty<vertex_descriptor>(dock_widget->propertyBox->currentText().toStdString(),
                                                 *(sm_item->face_graph())))
        {
          return;
        }
        sm_item->setRenderingMode(GouraudPlusEdges);
      }
      else if(property_simplex_types.at(property_index) == Property_simplex_type::FACE)
      {
        if(!displaySMProperty<face_descriptor>(dock_widget->propertyBox->currentText().toStdString(),
                                               *(sm_item->face_graph())))
        {
          return;
        }
      }
    }

    sm_item->invalidateOpenGLBuffers();
    sm_item->redraw();
  }

private Q_SLOTS:
  void colorize()
  {
    const int property_index = dock_widget->propertyBox->currentIndex();
    if(property_index < 0 || property_index >= dock_widget->propertyBox->count())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    enableExtremeValues();

    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_points_with_normal_item* ps_item = qobject_cast<Scene_points_with_normal_item*>(item);
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);

    if(sm_item)
      colorizeSM(sm_item);
    else if(ps_item)
      colorizePS(ps_item);

    // @todo emit a new SIGNAL on successful coloring, something like "colorChanged()"
    // itemChanged is too strong and would conflict with the connection below

    QApplication::restoreOverrideCursor();

    // below is a hackish way to call the connection only once.
    // This is required because the itemChanged signal is currently emitted when the colors are reset...
    //
    // @todo do not emit itemChanged when the colors are reset
    // @todo with qt6, single connection can be performed with `static_cast<Qt::ConnectionType>(Qt::SingleShotConnection)`
    // see https://www.kdab.com/single-shot-connections/
    auto connection = std::make_shared<QMetaObject::Connection>();
    *connection = connect(item, &Scene_surface_mesh_item::itemChanged,
                          this, [this, item, ps_item, sm_item, connection]()
                                {
                                  QObject::disconnect(*connection);

                                  this->removeDisplayPluginProperties(item); // meaningful only for sm_item

                                  // @todo Scene_item doesn't have resetColors()...
                                  if(ps_item)
                                    ps_item->resetColors();
                                  else if(sm_item)
                                    sm_item->resetColors();

                                  if(item == scene->item(scene->mainSelectionIndex()))
                                    onItemIndexSelected(scene->item_id(item));
                                });
  }

private:
  void removeDisplayPluginProperty(Scene_item* item,
                                   const std::string& property_name)
  {
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(!sm_item)
      return;

    SMesh* sm = sm_item->face_graph();
    if(sm == nullptr)
      return;

    // Here we only target the property maps added by this plugin, so 'double' is fine
    SMesh::Property_map<face_descriptor, double> property;
    bool found;
    std::tie(property, found) = sm->property_map<face_descriptor, double>(property_name);
    if(found)
      sm->remove_property_map(property);
  }

  void removeDisplayPluginProperties(Scene_item* item)
  {
    removeDisplayPluginProperty(item, "f:display_plugin_smallest_angle");
    removeDisplayPluginProperty(item, "f:display_plugin_largest_angle");
    removeDisplayPluginProperty(item, "f:display_plugin_scaled_jacobian");
    removeDisplayPluginProperty(item, "f:display_plugin_area");
    removeDisplayPluginProperty(item, "v:display_plugin_interpolated_corrected_mean_curvature");
    removeDisplayPluginProperty(item, "v:display_plugin_interpolated_corrected_Gaussian_curvature");
  }

  void displayExtremumAnglePerFace(Scene_surface_mesh_item* sm_item,
                                   const Extremum extremum)
  {
    SMesh* sm = sm_item->face_graph();
    if(sm == nullptr)
      return;

    // if the face is already a triangle, do not extract it from the mesh
    auto triangular_face_sector_angle = [](halfedge_descriptor h,
                                           const SMesh& mesh) -> double
                                        {
                                          auto vpm = get(boost::vertex_point, mesh);
                                          return CGAL::approximate_angle(get(vpm, source(h, mesh)),
                                                                        get(vpm, target(h, mesh)),
                                                                        get(vpm, target(next(h, mesh), mesh)));
                                        };

    // for non-triangular faces, extract a one-face mesh and triangulate it
    auto single_face_sector_angle = [](halfedge_descriptor h,
                                       const SMesh& mesh) -> double
                                    {
                                      CGAL_precondition(!is_border(h, mesh) && is_border(opposite(h, mesh), mesh));

                                      auto vpm = get(boost::vertex_point, mesh);

                                      double sector_angle = 0;
                                      do
                                      {
                                        sector_angle += CGAL::approximate_angle(get(vpm, source(h, mesh)),
                                                                                get(vpm, target(h, mesh)),
                                                                                get(vpm, target(next(h, mesh), mesh)));
                                        h = opposite(next(h, mesh), mesh);
                                      }
                                      while(!is_border(h, mesh));

                                      return sector_angle;
                                    };

    bool not_initialized;
    SMesh::Property_map<face_descriptor, double> fangle;

    if(extremum == MIN_VALUE)
      std::tie(fangle, not_initialized) = sm->add_property_map<face_descriptor, double>("f:display_plugin_smallest_angle", 0);
    else
      std::tie(fangle, not_initialized) = sm->add_property_map<face_descriptor, double>("f:display_plugin_largest_angle", 0);

    SMesh& mesh = *sm;
    auto vpm = get(boost::vertex_point, mesh);

    if(not_initialized)
    {
      for(face_descriptor f : faces(mesh))
      {
        if(CGAL::is_triangle(halfedge(f, mesh), mesh))
        {
          if(extremum == MIN_VALUE)
          {
            fangle[f] = (std::min)({triangular_face_sector_angle(halfedge(f, mesh), mesh),
                                    triangular_face_sector_angle(next(halfedge(f, mesh), mesh), mesh),
                                    triangular_face_sector_angle(prev(halfedge(f, mesh), mesh), mesh)});
          }
          else
          {
            fangle[f] = (std::max)({triangular_face_sector_angle(halfedge(f, mesh), mesh),
                                    triangular_face_sector_angle(next(halfedge(f, mesh), mesh), mesh),
                                    triangular_face_sector_angle(prev(halfedge(f, mesh), mesh), mesh)});
          }
        }
        else
        {
          SMesh local_smesh;
          auto local_vpm = get(boost::vertex_point, local_smesh);
          std::vector<vertex_descriptor> local_vertices;

          for(halfedge_descriptor h : halfedges_around_face(halfedge(f, mesh), mesh))
          {
            local_vertices.push_back(CGAL::add_vertex(local_smesh));
            put(local_vpm, local_vertices.back(), get(vpm, target(h, mesh)));
          }

          face_descriptor local_f = CGAL::Euler::add_face(local_vertices, local_smesh);

          // walk the border of the face to walk the halfedge of the pre-triangulation face
          halfedge_descriptor local_border_h = opposite(halfedge(local_f, local_smesh), local_smesh);
          CGAL_assertion(is_border(local_border_h, local_smesh));

          PMP::triangulate_faces(local_smesh);

          double extremum_angle_in_face = ARBITRARY_DBL_MAX;
          halfedge_descriptor local_border_end_h = local_border_h;
          do
          {
            double angle = single_face_sector_angle(opposite(local_border_h, local_smesh), local_smesh);
            if(extremum == MIN_VALUE)
              extremum_angle_in_face = (std::min)(extremum_angle_in_face, angle);
            else
              extremum_angle_in_face = (std::max)(extremum_angle_in_face, angle);

            local_border_h = next(local_border_h, local_smesh);
          }
          while(local_border_h != local_border_end_h);

          fangle[f] = extremum_angle_in_face;
        }
      }
    }

    if(extremum == MIN_VALUE)
      displaySMProperty<face_descriptor>("f:display_plugin_smallest_angle", mesh);
    else
      displaySMProperty<face_descriptor>("f:display_plugin_largest_angle", mesh);
  }

  double scaled_jacobian(const face_descriptor f,
                         const SMesh& mesh) const;

  void displayScaledJacobian(Scene_surface_mesh_item* sm_item)
  {
    SMesh* sm = sm_item->face_graph();
    if(sm == nullptr)
      return;

    bool not_initialized;
    SMesh::Property_map<face_descriptor, double> fjacobian;
    std::tie(fjacobian, not_initialized) = sm->add_property_map<face_descriptor, double>("f:display_plugin_scaled_jacobian", 0);

    if(not_initialized)
    {
      for(face_descriptor f : faces(*sm))
        fjacobian[f] = scaled_jacobian(f, *sm);
    }

    displaySMProperty<face_descriptor>("f:display_plugin_scaled_jacobian", *sm);
  }

  double area(const face_descriptor f,
              const SMesh& mesh) const
  {
    if(CGAL::is_triangle(halfedge(f, mesh), mesh))
      return PMP::face_area(f, mesh);

    auto vpm = get(boost::vertex_point, mesh);

    // create a local version of the mesh, triangulate it, sum the triangle areas
    SMesh local_smesh;
    auto local_vpm = get(boost::vertex_point, local_smesh);
    std::vector<vertex_descriptor> local_vertices;

    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, mesh), mesh))
    {
      local_vertices.push_back(CGAL::add_vertex(local_smesh));
      put(local_vpm, local_vertices.back(), get(vpm, target(h, mesh)));
    }

    CGAL::Euler::add_face(local_vertices, local_smesh);
    PMP::triangulate_faces(local_smesh);
    return PMP::area(local_smesh);
  }

  void displayArea(Scene_surface_mesh_item* sm_item)
  {
    SMesh* sm = sm_item->face_graph();
    if(sm == nullptr)
      return;

    bool not_initialized;
    SMesh::Property_map<face_descriptor, double> farea;
    std::tie(farea, not_initialized) = sm->add_property_map<face_descriptor, double>("f:display_plugin_area", 0);

    if(not_initialized)
    {
      for(face_descriptor f : faces(*sm))
        farea[f] = area(f, *sm);
    }

    displaySMProperty<face_descriptor>("f:display_plugin_area", *sm);
  }

private Q_SLOTS:
  void setExpandingRadius()
  {
    double sliderMin = dock_widget->expandingRadiusSlider->minimum();
    double sliderMax = dock_widget->expandingRadiusSlider->maximum() - sliderMin;
    double val = dock_widget->expandingRadiusSlider->value() - sliderMin;
    sliderMin = 0;

    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(sm_item == nullptr)
      return;

    SMesh& smesh = *(sm_item->face_graph());

    auto vpm = get(CGAL::vertex_point, smesh);

    // @todo use the upcoming PMP::longest_edge
    if(maxEdgeLength < 0)
    {
      auto edge_range = CGAL::edges(smesh);

      if(num_edges(smesh) == 0)
      {
        expand_radius = 0;
        dock_widget->expandingRadiusLabel->setText(tr("Expanding Radius: %1").arg(expand_radius));
        return;
      }

      auto eit = std::max_element(edge_range.begin(), edge_range.end(),
                                  [&, vpm, smesh](auto l, auto r)
                                  {
                                    auto res = EPICK().compare_squared_distance_3_object()(
                                        get(vpm, source((l), smesh)),
                                        get(vpm, target((l), smesh)),
                                        get(vpm, source((r), smesh)),
                                        get(vpm, target((r), smesh)));
                                    return (res == CGAL::SMALLER);
                                });

      CGAL_assertion(eit != edge_range.end());

      maxEdgeLength = PMP::edge_length(*eit, smesh);
    }

    double outMax = 5 * maxEdgeLength, base = 1.2;

    expand_radius = (pow(base, val) - 1) * outMax / (pow(base, sliderMax) - 1);
    dock_widget->expandingRadiusLabel->setText(tr("Expanding Radius: %1").arg(expand_radius));
  }

private:
  void displayInterpolatedCurvatureMeasure(Scene_surface_mesh_item* item,
                                           CurvatureType mu_index)
  {
    if(mu_index != MEAN_CURVATURE && mu_index != GAUSSIAN_CURVATURE)
      return;

    std::string tied_string = (mu_index == MEAN_CURVATURE) ? "v:display_plugin_interpolated_corrected_mean_curvature"
                                                           : "v:display_plugin_interpolated_corrected_Gaussian_curvature";

    SMesh& smesh = *item->face_graph();

    const auto vnm = smesh.property_map<vertex_descriptor, EPICK::Vector_3>("v:normal_before_perturbation").first;
    const bool vnm_exists = smesh.property_map<vertex_descriptor, EPICK::Vector_3>("v:normal_before_perturbation").second;

    // compute once and store the value per vertex
    bool non_init;
    SMesh::Property_map<vertex_descriptor, double> mu_i_map;
    std::tie(mu_i_map, non_init) = smesh.add_property_map<vertex_descriptor, double>(tied_string, 0);
    if(non_init)
    {
      if(vnm_exists)
      {
        if(mu_index == MEAN_CURVATURE)
        {
          PMP::interpolated_corrected_curvatures(smesh,
                                                 CGAL::parameters::vertex_mean_curvature_map(mu_i_map)
                                                                  .ball_radius(expand_radius)
                                                                  .vertex_normal_map(vnm));
        }
        else
        {
          PMP::interpolated_corrected_curvatures(smesh,
                                                 CGAL::parameters::vertex_Gaussian_curvature_map(mu_i_map)
                                                                  .ball_radius(expand_radius)
                                                                  .vertex_normal_map(vnm));
        }
      }
      else
      {
        if(mu_index == MEAN_CURVATURE)
        {
          PMP::interpolated_corrected_curvatures(smesh,
                                                  CGAL::parameters::vertex_mean_curvature_map(mu_i_map)
                                                                   .ball_radius(expand_radius));
        }
        else
        {
          PMP::interpolated_corrected_curvatures(smesh,
                                                 CGAL::parameters::vertex_Gaussian_curvature_map(mu_i_map)
                                                                  .ball_radius(expand_radius));
        }
      }
    }

    displaySMProperty<vertex_descriptor>(tied_string, smesh);
  }

private:
  template<typename Functor>
  bool call_on_PS_property(const std::string& name,
                           const Point_set& ps,
                           const Functor& functor) const;

  template<typename Simplex, typename Functor>
  bool call_on_SM_property(const std::string& name,
                           const SMesh& sm,
                           const Functor& functor) const;

private:
  template<typename PM>
  bool displayPSProperty(Point_set& ps,
                         PM pm)
  {
    PSDisplayer<PM> display_property(ps, pm, this);
    return display_property();
  }

  bool displayPSProperty(const std::string& name,
                         Point_set& ps);

  // -
  template<typename PM>
  bool displaySMProperty(SMesh& mesh,
                         PM pm,
                         vertex_descriptor)
  {
    SMVertexDisplayer<PM> display_property(mesh, pm, this);
    return display_property();
  }

  template<typename PM>
  bool displaySMProperty(SMesh& mesh,
                         PM pm,
                         face_descriptor)
  {
    SMFaceDisplayer<PM> display_property(mesh, pm, this);
    return display_property();
  }

  template<typename Simplex>
  bool displaySMProperty(const std::string& name,
                         SMesh& mesh);

private:
  template <typename SimplexRange>
  auto SimplexWithPropertyExtremum(const SimplexRange& simplex_range,
                                   const SMesh& mesh,
                                   const std::string& property_name,
                                   const Extremum extremum) const
  {
    using Simplex = typename boost::range_value<SimplexRange>::type;

    // We don't know what's the type of the property map so we can't simply do
    //   mesh.property_map<Simplex, TYPE>(property_name),
    // we have to try all acceptable types.
    Simplex extremum_s;
    call_on_SM_property<Simplex>(property_name, mesh,
                                 [extremum, &extremum_s, &simplex_range](const auto pmap) -> bool
                                 {
                                   double extremum_value = (extremum == MIN_VALUE) ? ARBITRARY_DBL_MAX : - ARBITRARY_DBL_MAX;

                                   for(Simplex s : simplex_range)
                                   {
                                     if((extremum == MIN_VALUE && get(pmap, s) < extremum_value) ||
                                       (extremum == MAX_VALUE && get(pmap, s) > extremum_value))
                                     {
                                       extremum_value = get(pmap, s);
                                       extremum_s = s;
                                     }
                                   }

                                    return true;
                                  });

    CGAL_assertion(extremum_s != Simplex());
    return extremum_s;
  }

  template <typename SimplexRange>
  void zoomToSimplexWithPropertyExtremum(const SimplexRange& simplex_range,
                                         const SMesh& mesh,
                                         const std::string& property_name,
                                         const Extremum extremum) const
  {
    using Simplex = typename boost::range_value<SimplexRange>::type;

    const Simplex extremum_s = SimplexWithPropertyExtremum(simplex_range, mesh, property_name, extremum);

    QString sid;
    if(std::is_same<Simplex, vertex_descriptor>::value)
      sid = QString("v%1").arg(extremum_s);
    else
      sid = QString("f%1").arg(extremum_s);

    face_descriptor unused_fd;
    Point_3 unused_p;
    ::zoomToId(mesh, sid,
               getActiveViewer(),
               unused_fd, unused_p);
  };

  void zoomToPointWithPropertyExtremum(const Point_set& ps,
                                       const std::string& property_name,
                                       const Extremum extremum) const
  {
    Point_set::Index extremum_i = -1;
    call_on_PS_property(property_name, ps,
                        [extremum, &extremum_i, &ps](const auto pmap) -> bool
                        {
                          double extremum_value = (extremum == MIN_VALUE) ? ARBITRARY_DBL_MAX : - ARBITRARY_DBL_MAX;

                          for(Point_set::Index i : ps)
                          {
                            if((extremum == MIN_VALUE && get(pmap, i) < extremum_value) ||
                              (extremum == MAX_VALUE && get(pmap, i) > extremum_value))
                            {
                              extremum_value = get(pmap, i);
                              extremum_i = i;
                            }
                          }

                          return true;
                        });

    CGAL_assertion(extremum_i != Point_set::Index(-1));
    Point_3 unused_p;
    ::zoomToPoint(ps, extremum_i,
                  getActiveViewer(),
                  unused_p);
  }

  void on_zoomToButton_pressed(const Extremum extremum)
  {
    const int property_index = dock_widget->propertyBox->currentIndex();
    if(property_index < 0 || property_index >= dock_widget->propertyBox->count())
      return;

    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(sm_item)
    {
      const SMesh& mesh = *sm_item->face_graph();

      const std::string& property_name = dock_widget->propertyBox->currentText().toStdString();
      if(property_name == "Smallest Angle Per Face")
        zoomToSimplexWithPropertyExtremum(faces(mesh), mesh, "f:display_plugin_smallest_angle", extremum);
      else if(property_name == "Largest Angle Per Face")
        zoomToSimplexWithPropertyExtremum(faces(mesh), mesh, "f:display_plugin_largest_angle", extremum);
      else if(property_name == "Scaled Jacobian")
        zoomToSimplexWithPropertyExtremum(faces(mesh), mesh, "f:display_plugin_scaled_jacobian", extremum);
      else if(property_name == "Face Area")
        zoomToSimplexWithPropertyExtremum(faces(mesh), mesh, "f:display_plugin_area", extremum);
      else if(property_name == "Interpolated Corrected Mean Curvature")
        zoomToSimplexWithPropertyExtremum(vertices(mesh), mesh, "v:display_plugin_interpolated_corrected_mean_curvature", extremum);
      else if(property_name == "Interpolated Corrected Gaussian Curvature")
        zoomToSimplexWithPropertyExtremum(vertices(mesh), mesh, "v:display_plugin_interpolated_corrected_Gaussian_curvature", extremum);
      else if(property_simplex_types.at(property_index) == Property_simplex_type::VERTEX)
        zoomToSimplexWithPropertyExtremum(vertices(mesh), mesh, property_name, extremum);
      else if(property_simplex_types.at(property_index) == Property_simplex_type::FACE)
        zoomToSimplexWithPropertyExtremum(faces(mesh), mesh, property_name, extremum);
    }

    Scene_points_with_normal_item* ps_item = qobject_cast<Scene_points_with_normal_item*>(item);
    if(ps_item)
    {
      const Point_set& ps = *ps_item->point_set();

      const std::string& property_name = dock_widget->propertyBox->currentText().toStdString();
      zoomToPointWithPropertyExtremum(ps, property_name, extremum);
    }
  }

private Q_SLOTS:
  void on_zoomToMinButton_pressed()
  {
    on_zoomToButton_pressed(MIN_VALUE);
  }

  void on_zoomToMaxButton_pressed()
  {
    on_zoomToButton_pressed(MAX_VALUE);
  }

private:
  template<typename, typename, typename>
  friend class PropertyDisplayer;

  // CRTP used to display properties of surface meshes (vertex and face) and point sets.
  template <class DataSet, class PM, class CRTP_derived_class>
  struct PropertyDisplayer
  {
    using value_type = typename PM::value_type;

    DataSet& dataset;
    PM property_map;
    std::vector<value_type> values;
    Display_property_plugin* parent;

    PropertyDisplayer(DataSet& ds,
                      PM pm,
                      Display_property_plugin* parent)
      : dataset(ds), property_map(pm), parent(parent)
    { }

    bool operator()()
    {
      static_cast<CRTP_derived_class*>(this)->fill_values();
      std::sort(values.begin(), values.end());
      auto end = std::unique(values.begin(), values.end());

      parent->dock_widget->minBox->setValue(*values.begin());
      parent->dock_widget->maxBox->setValue(*(std::prev(end)));

      // fill color pmap
      if(parent->dock_widget->colorRampRadioButton->isChecked())
      {
        // scale a color ramp between min and max
        parent->displayRampLegend();
        static_cast<CRTP_derived_class*>(this)->color_with_ramp();
      }
      else
      {
        CGAL_assertion(parent->dock_widget->randomColorsRadioButton->isChecked());

        // generate color map
        parent->color_map.clear();
        compute_color_map(QColor(255 * parent->rI, 255 * parent->gI, 255 * parent->bI),
                          std::distance(values.begin(), end),
                          std::back_inserter(parent->color_map));

        // fill map
        std::unordered_map<value_type, std::size_t> value_index_map;
        std::size_t counter = 0;
        for(auto it=values.begin(); it!=end; ++it)
          value_index_map[*it] = counter++;

        static_cast<CRTP_derived_class*>(this)->color_with_map(value_index_map);
        parent->displayMapLegend(values);
      }

      return true;
    }
  }; // struct PropertyDisplayer

  template <class PM>
  struct PSDisplayer
    : public PropertyDisplayer<Point_set, PM, PSDisplayer<PM> >
  {
    using Base = PropertyDisplayer<Point_set, PM, PSDisplayer<PM> >;
    using value_type = typename PM::value_type;

    PSDisplayer(Point_set& ps,
                PM pm,
                Display_property_plugin* parent)
      : Base(ps, pm, parent)
    {}

    void fill_values()
    {
      for(const auto& p : this->dataset)
        this->values.push_back(this->property_map[p]);
    }

    void color_with_map(std::unordered_map<value_type, std::size_t>& value_index_map)
    {
      for(const auto& p : this->dataset)
      {
        CGAL::IO::Color color(this->parent->color_map[value_index_map[this->property_map[p]]].red(),
                              this->parent->color_map[value_index_map[this->property_map[p]]].green(),
                              this->parent->color_map[value_index_map[this->property_map[p]]].blue());
        this->dataset.set_color(p, color.red(), color.green(), color.blue());
      }
    }

    void color_with_ramp() const
    {
      double min = this->parent->dock_widget->minBox->value();
      double max = this->parent->dock_widget->maxBox->value();
      for(const auto& p : this->dataset)
      {
        if(min == max)
          min -= 1.;

        double val = (this->property_map[p] - min) / (max - min);
        val = boost::algorithm::clamp<double>(val, 0., 1.);

        CGAL::IO::Color color(255 * this->parent->color_ramp.r(val),
                              255 * this->parent->color_ramp.g(val),
                              255 * this->parent->color_ramp.b(val));
        this->dataset.set_color(p, color.red(), color.green(), color.blue());
      }
    }
  }; // PSDisplayer

  template <class PM>
  struct SMVertexDisplayer
    : public PropertyDisplayer<SMesh, PM, SMVertexDisplayer<PM> >
  {
    using Base = PropertyDisplayer<SMesh, PM, SMVertexDisplayer<PM> > ;
    using value_type = typename PM::value_type;

    SMVertexDisplayer(SMesh& mesh,
                      PM pm,
                      Display_property_plugin* parent)
      : Base(mesh, pm, parent)
    {}

    void fill_values()
    {
      for(vertex_descriptor v : vertices(this->dataset))
        this->values.push_back(this->property_map[v]);
    }

    void color_with_map(std::unordered_map<value_type, std::size_t>& value_index_map) const
    {
      SMesh::Property_map<vertex_descriptor, CGAL::IO::Color> vcolors =
          this->dataset.template add_property_map<vertex_descriptor, CGAL::IO::Color >("v:color", CGAL::IO::Color()).first;
      for(vertex_descriptor v : vertices(this->dataset))
      {
        CGAL::IO::Color color(this->parent->color_map[value_index_map[this->property_map[v]]].red(),
                              this->parent->color_map[value_index_map[this->property_map[v]]].green(),
                              this->parent->color_map[value_index_map[this->property_map[v]]].blue());
        vcolors[v] = color;
      }
    }

    void color_with_ramp() const
    {
      SMesh::Property_map<vertex_descriptor, CGAL::IO::Color> vcolors =
          this->dataset.template add_property_map<vertex_descriptor, CGAL::IO::Color >("v:color", CGAL::IO::Color()).first;

      double min = this->parent->dock_widget->minBox->value();
      double max = this->parent->dock_widget->maxBox->value();

      if(min == max)
        min -= 1.;

      for(vertex_descriptor v : vertices(this->dataset))
      {
        double val = (this->property_map[v] - min) / (max - min);
        val = boost::algorithm::clamp<double>(val, 0., 1.);

        CGAL::IO::Color color(255 * this->parent->color_ramp.r(val),
                              255 * this->parent->color_ramp.g(val),
                              255 * this->parent->color_ramp.b(val));
        vcolors[v] = color;
      }
    }
  }; // struct SMVertexDisplayer

  template <class PM>
  struct SMFaceDisplayer
    : public PropertyDisplayer<SMesh, PM, SMFaceDisplayer<PM> >
  {
    using Base = PropertyDisplayer<SMesh, PM, SMFaceDisplayer<PM> >;
    using value_type = typename PM::value_type;

    SMFaceDisplayer(SMesh& mesh,
                    PM pm,
                    Display_property_plugin* parent)
      : Base(mesh, pm, parent)
    { }

    void fill_values()
    {
      for(face_descriptor f : faces(this->dataset))
        this->values.push_back(this->property_map[f]);
    }

    void color_with_map(std::unordered_map<value_type, std::size_t>& value_index_map) const
    {
      SMesh::Property_map<face_descriptor, CGAL::IO::Color> fcolors =
          this->dataset.template add_property_map<face_descriptor, CGAL::IO::Color >("f:color", CGAL::IO::Color()).first;
      for(face_descriptor f : faces(this->dataset))
      {
        CGAL::IO::Color color(this->parent->color_map[value_index_map[this->property_map[f]]].red(),
                              this->parent->color_map[value_index_map[this->property_map[f]]].green(),
                              this->parent->color_map[value_index_map[this->property_map[f]]].blue());
        fcolors[f] = color;
      }
    }

    void color_with_ramp() const
    {
      SMesh::Property_map<face_descriptor, CGAL::IO::Color> fcolors =
          this->dataset.template add_property_map<face_descriptor, CGAL::IO::Color >("f:color", CGAL::IO::Color()).first;

      double min = this->parent->dock_widget->minBox->value();
      double max = this->parent->dock_widget->maxBox->value();

      if(min == max)
        min -= 1.;

      for(face_descriptor f : faces(this->dataset))
      {
        double val = (this->property_map[f] - min) / (max - min);
        val = boost::algorithm::clamp<double>(val, 0., 1.);

        CGAL::IO::Color color(255 * this->parent->color_ramp.r(val),
                              255 * this->parent->color_ramp.g(val),
                              255 * this->parent->color_ramp.b(val));
        fcolors[f] = color;
      }
    }
  }; // struct SMFaceDisplayer
};

  /// Code based on the verdict module of vtk

/*=========================================================================
  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or https://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

double
Display_property_plugin::
scaled_jacobian(const face_descriptor f,
                const SMesh& mesh) const
{
  boost::property_map<SMesh, boost::vertex_point_t>::type vpm = get(boost::vertex_point, mesh);
  std::vector<double> corner_areas(degree(f, mesh));

  std::vector<EPICK::Vector_3> edges;
  for(halfedge_descriptor hd : CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
  {
    edges.emplace_back(get(vpm, source(hd, mesh)),
                       get(vpm, target(hd, mesh)));
  }

  std::vector<EPICK::Vector_3> corner_normals;
  for(std::size_t i=0; i<edges.size(); ++i)
    corner_normals.push_back(CGAL::cross_product(edges[i], edges[(i+1)%(edges.size())]));

  EPICK::Vector_3 unit_center_normal = PMP::compute_face_normal(f, mesh);

  for(std::size_t i=0; i<corner_areas.size(); ++i)
    corner_areas[i] =  unit_center_normal*corner_normals[i];

  std::vector<double> length;
  for(std::size_t i=0; i<edges.size(); ++i)
  {
    length.push_back(CGAL::approximate_sqrt(edges[i].squared_length()));
    if(length[i] < ARBITRARY_DBL_MIN)
      return 0.;
  }

  double min_scaled_jac = ARBITRARY_DBL_MAX;
  for(std::size_t i=0; i<edges.size(); ++i)
  {
    double scaled_jac = corner_areas[i] / (length[i] * length[(i+edges.size()-1)%(edges.size())]);
    min_scaled_jac = (std::min)(scaled_jac, min_scaled_jac);
  }

  if(min_scaled_jac > 0)
    return (std::min)(min_scaled_jac, ARBITRARY_DBL_MAX);

  return (std::max)(min_scaled_jac, -ARBITRARY_DBL_MAX);
}

bool
Display_property_plugin::
isPSPropertyScalar(const std::string& name,
                   const Point_set& ps) const
{
  if(name == "red" || name == "green" || name == "blue")
    return false;

  // the dispatch function does the filtering we want: if it founds a property
  // with which it can call the functor, then it already has a property we want
  return call_on_PS_property(name, ps, [](auto) -> bool { return true; });
}

// Shenanigans to deal with the fact that the value type is not known
//
// find the property map with the correct value, and call the functor
template<typename Simplex>
bool
Display_property_plugin::
isSMPropertyScalar(const std::string& name,
                   const SMesh& mesh) const
{
  // do not detect this plugin's properties
  if(name == "f:display_plugin_smallest_angle" ||
     name == "f:display_plugin_largest_angle" ||
     name == "f:display_plugin_scaled_jacobian" ||
     name == "f:display_plugin_area" ||
     name == "v:display_plugin_interpolated_corrected_mean_curvature" ||
     name == "v:display_plugin_interpolated_corrected_Gaussian_curvature")
    return false;

  // the dispatch function does the filtering we want: if it founds a property
  // with which it can call the functor, then it already has a property we want
  return call_on_SM_property<Simplex>(name, mesh, [](auto) -> bool { return true; });
}

bool
Display_property_plugin::
displayPSProperty(const std::string& name,
                  Point_set& ps)
{
  return call_on_PS_property(name, ps,
                             [this, &ps](auto pmap) -> bool
                             {
                                return this->displayPSProperty(ps, pmap);
                             });
}

template<typename Simplex>
bool
Display_property_plugin::
displaySMProperty(const std::string& name,
                  SMesh& mesh)
{
  return call_on_SM_property<Simplex>(name, mesh,
                                      [this, &mesh](auto pmap) -> bool
                                      {
                                        return this->displaySMProperty(mesh, pmap, Simplex());
                                      });
}

template<typename Functor>
bool
Display_property_plugin::
call_on_PS_property(const std::string& name,
                    const Point_set& ps,
                    const Functor& functor) const
{
  if(ps.template property_map<std::int8_t>(name).second)
    return functor(ps.template property_map<std::int8_t>(name).first);
  else if(ps.template property_map<std::uint8_t>(name).second)
    return functor(ps.template property_map<std::uint8_t>(name).first);
  else if(ps.template property_map<std::int16_t>(name).second)
    return functor(ps.template property_map<std::int16_t>(name).first);
  else if(ps.template property_map<std::uint16_t>(name).second)
    return functor(ps.template property_map<std::uint16_t>(name).first);
  else if(ps.template property_map<std::int32_t>(name).second)
    return functor(ps.template property_map<std::int32_t>(name).first);
  else if(ps.template property_map<std::uint32_t>(name).second)
    return functor(ps.template property_map<std::uint32_t>(name).first);
  else if(ps.template property_map<std::int64_t>(name).second)
    return functor(ps.template property_map<std::int64_t>(name).first);
  else if(ps.template property_map<std::uint64_t>(name).second)
    return functor(ps.template property_map<std::uint64_t>(name).first);
  else if(ps.template property_map<float>(name).second)
    return functor(ps.template property_map<float>(name).first);
  else if(ps.template property_map<double>(name).second)
    return functor(ps.template property_map<double>(name).first);

  return false;
}

template <typename Simplex, typename Functor>
bool
Display_property_plugin::
call_on_SM_property(const std::string& name,
                    const SMesh& mesh,
                    const Functor& functor) const
{
  if(mesh.template property_map<Simplex, std::int8_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::int8_t>(name).first);
  else if(mesh.template property_map<Simplex, std::uint8_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::uint8_t>(name).first);
  else if(mesh.template property_map<Simplex, std::int16_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::int16_t>(name).first);
  else if(mesh.template property_map<Simplex, std::uint16_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::uint16_t>(name).first);
  else if(mesh.template property_map<Simplex, std::int32_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::int32_t>(name).first);
  else if(mesh.template property_map<Simplex, std::uint32_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::uint32_t>(name).first);
  else if(mesh.template property_map<Simplex, std::int64_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::int64_t>(name).first);
  else if(mesh.template property_map<Simplex, std::uint64_t>(name).second)
    return functor(mesh.template property_map<Simplex, std::uint64_t>(name).first);
  else if(mesh.template property_map<Simplex, float>(name).second)
    return functor(mesh.template property_map<Simplex, float>(name).first);
  else if(mesh.template property_map<Simplex, double>(name).second)
    return functor(mesh.template property_map<Simplex, double>(name).first);

  return false;
}

#include "Display_property_plugin.moc"
