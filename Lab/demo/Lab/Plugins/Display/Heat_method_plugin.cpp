#include "ui_Heat_method.h"

#include "Color_ramp.h"
#include "id_printing.h"
#include "Messages_interface.h"
#include "triangulate_primitive.h"

#include "Scene.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_surface_mesh_item.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Buffer_for_vao.h>

#include <CGAL/boost/bimap.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <QAbstractItemView>
#include <QAction>
#include <QApplication>
#include <QColor>
#include <QColorDialog>
#include <QInputDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QPalette>
#include <QStyleFactory>

// @fixme multiple selection items are broken, so don't create a selection item if there is already one

using namespace CGAL::Three;

Viewer_interface* (&getActiveViewer)() = Three::activeViewer;

class Scene_heat_item
  : public Scene_item_rendering_helper
{
  Q_OBJECT

private:
  Scene_surface_mesh_item* parent;
  SMesh* sm;

  mutable std::vector<float> verts;
  mutable std::vector<float> normals;
  mutable std::vector<unsigned int> idx;
  mutable std::vector<float> colors;
  mutable std::vector<float> heat_values;
  mutable std::size_t nb_idx;

public:
  Scene_heat_item(Scene_surface_mesh_item* item)
    : parent(item), sm(item->face_graph())
  {
    CGAL_precondition(is_triangle_mesh(*sm));

    setTriangleContainer(0, new Triangle_container(Viewer_interface::PROGRAM_HEAT_INTENSITY, true));
    setRenderingMode(Gouraud);
  }

  ~Scene_heat_item(){}

  Scene_item* clone() const override { return nullptr; }
  QString toolTip() const override{ return QString(); }
  Scene_surface_mesh_item* getParent() { return parent; }
  bool isEmpty() const override { return false; }

  SMesh* face_graph() { return sm; }

  void initializeBuffers(Viewer_interface *viewer) const override
  {
    getTriangleContainer(0)->initializeBuffers(viewer);
    getTriangleContainer(0)->setIdxSize(nb_idx);

    verts.resize(0);
    verts.shrink_to_fit();
    normals.resize(0);
    normals.shrink_to_fit();
    colors.resize(0);
    colors.shrink_to_fit();
    idx.clear();
    idx.shrink_to_fit();
  }

  void draw(Viewer_interface *viewer) const override
  {
    if(!visible())
      return;

    if(!isInit(viewer))
      initGL(viewer);

    if(getBuffersFilled() && !getBuffersInit(viewer))
    {
      initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }

    if(!getBuffersFilled())
    {
      computeElements();
      initializeBuffers(viewer);
    }

    getTriangleContainer(0)->setAlpha(1.0f);
    getTriangleContainer(0)->draw(viewer, false);
  }

  void compute_bbox() const override
  {
    setBbox(parent->bbox());
  }

  virtual bool supportsRenderingMode(RenderingMode m) const override
  {
    return (m == Gouraud);
  }

  virtual void invalidateOpenGLBuffers() override
  {
    setBuffersFilled(false);
    compute_bbox();
    getTriangleContainer(0)->reset_vbos(NOT_INSTANCED);
  }

  void computeElements() const override
  {
    typedef CGAL::Buffer_for_vao CPF;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    auto vpm = CGAL::get(CGAL::vertex_point, *sm);
    auto vnormals = sm->property_map<vertex_descriptor, EPICK::Vector_3>("v:normal").value();

    auto vcolors = sm->property_map<vertex_descriptor, CGAL::IO::Color>("v:color").value();
    auto vdist = sm->property_map<vertex_descriptor, double>("v:HM_Plugin_heat_intensity").value();

    verts.clear();
    normals.clear();
    idx.clear();
    colors.clear();

    idx.reserve(3 * num_faces(*sm));
    for(face_descriptor fd : faces(*sm))
    {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *sm),*sm))
        idx.push_back(source(hd, *sm));
    }

    const CGAL::qglviewer::Vec& o = static_cast<Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    EPICK::Vector_3 offset(o.x, o.y, o.z);

    for(vertex_descriptor vd : vertices(*sm))
    {
      const CGAL::IO::Color& c = vcolors[vd];
      colors.push_back(float(c.red()) / 255);
      colors.push_back(float(c.green()) / 255);
      colors.push_back(float(c.blue()) / 255);

      const EPICK::Point_3& p = get(vpm, vd) + offset;
      CPF::add_point_in_buffer(p, verts);

      const EPICK::Vector_3& n = vnormals[vd];
      CPF::add_normal_in_buffer(n, normals);

      heat_values.push_back(vdist[vd]);
    }

    nb_idx = idx.size();
    getTriangleContainer(0)->allocate(Triangle_container::Vertex_indices, idx.data(),
                                      static_cast<int>(idx.size()*sizeof(unsigned int)));
    getTriangleContainer(0)->allocate(Triangle_container::Smooth_vertices, verts.data(),
                                      static_cast<int>(num_vertices(*sm)*3*sizeof(float)));
    getTriangleContainer(0)->allocate(Triangle_container::Smooth_normals, normals.data(),
                                      static_cast<int>(num_vertices(*sm)*3*sizeof(float)));
    getTriangleContainer(0)->allocate(Triangle_container::VColors, colors.data(),
                                      static_cast<int>(colors.size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Triangle_container::Distances, heat_values.data(),
                                      static_cast<int>(heat_values.size()*sizeof(float)));

    compute_bbox();
    setBuffersFilled(true);

    QApplication::restoreOverrideCursor();
  }
}; // class Scene_heat_item

class DockWidget
  : public QDockWidget,
    public Ui::HeatMethodWidget
{
public:
  DockWidget(const QString& name, QWidget *parent)
    : QDockWidget(name, parent)
  {
    setupUi(this);
  }
};

class Heat_method_plugin
  : public QObject,
    public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

  using Vertex_distance_map = SMesh::Property_map<vertex_descriptor, double>;
  using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<SMesh>;
  using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<SMesh, CGAL::Heat_method_3::Intrinsic_Delaunay>;

private:
  QAction* actionHeatMethod;

  DockWidget* dock_widget;

  // coloring choice and legend
  double rm = 1.;
  double gm = 0.;
  double bm = 0.;
  double rM = 0.;
  double gM = 1.;
  double bM = 0.;

  Color_ramp color_ramp;
  QPixmap legend;

  // tracking which scene items have which sources, and which heat method builders
  boost::bimap<boost::bimaps::set_of<Scene_surface_mesh_item*>,
               boost::bimaps::set_of<Scene_polyhedron_selection_item*> > item_source_vertices;

  // the point of storing this is that in the Heat Method(s), a number of computations are only
  // dependent on the mesh, and not the sources, and such can be performed just once.
  std::unordered_map<Scene_surface_mesh_item*, Heat_method*> heat_methods;
  std::unordered_map<Scene_surface_mesh_item*, Heat_method_idt*> idt_heat_methods;

public:
  bool applicable(QAction*) const override
  {
    // Single item => it must be a mesh and the selection item will be created through the plugin's button
    if(scene->selectionIndices().size() == 1)
    {
      Scene_item* item = scene->item(scene->mainSelectionIndex());
      return qobject_cast<Scene_surface_mesh_item*>(item);
    }
    // Two items => it must be a surface mesh and a selection item (in any order)
    else if(scene->selectionIndices().size() == 2)
    {
      Scene_item* item1 = scene->item(scene->selectionIndices().front());
      Scene_item* item2 = scene->item(scene->selectionIndices().back());
      return ((qobject_cast<Scene_surface_mesh_item*>(item1) &&
               qobject_cast<Scene_polyhedron_selection_item*>(item2)) ||
              (qobject_cast<Scene_polyhedron_selection_item*>(item1) &&
               qobject_cast<Scene_surface_mesh_item*>(item2)));
    }

    return false;
  }

  QList<QAction*> actions() const override
  {
    return QList<QAction*>() << actionHeatMethod;
  }

  void init(QMainWindow* mw,
            Scene_interface* sc,
            Messages_interface*) override
  {
    this->scene = sc;
    this->mw = mw;

    actionHeatMethod = new QAction(QString("Heat Method"), mw);
    actionHeatMethod->setProperty("submenuName", "Color");

    connect(actionHeatMethod, SIGNAL(triggered()),
            this, SLOT(openDialog()));

    Scene* scene_item = static_cast<Scene*>(scene);
    connect(scene_item, SIGNAL(itemIndicesSelected(QList<int>)),
            this, SLOT(onItemIndicesSelected(QList<int>)));

    // Dock Widget
    dock_widget = new DockWidget("Heat Method", mw);
    addDockWidget(dock_widget);
    dock_widget->setVisible(false);

    connect(dock_widget->methodBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(onNewMethodSelected(int)));

    dock_widget->methodBox->addItems({"Heat Method",
                                      "Heat Method (Intrinsic Delaunay)"});

    connect(dock_widget->createSourceVerticesButton, SIGNAL(clicked()),
            this, SLOT(createSourceVerticesSelectionItem()));

    QPalette palette(Qt::red);
    dock_widget->minColorButton->setPalette(palette);
    dock_widget->minColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->minColorButton->update();

    palette = QPalette(Qt::green);
    dock_widget->maxColorButton->setPalette(palette);
    dock_widget->maxColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->maxColorButton->update();

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

                                                    displayRampLegend();
                                                  });
                                 };

    connect_color_buttons(dock_widget->minColorButton, rm, gm, bm);
    connect_color_buttons(dock_widget->maxColorButton, rM, gM, bM);

    // Main action connection
    connect(dock_widget->estimateDistancesButton, SIGNAL(clicked(bool)),
            this, SLOT(estimateDistances()));

    // Post coloring connection
    connect(dock_widget->zoomToMaxButton, &QPushButton::pressed,
            this, &Heat_method_plugin::on_zoomToMaxButton_pressed);
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
    void disableExtremeValue()
  {
    dock_widget->extremeValueGroup->setEnabled(false);
    dock_widget->zoomToMaxButton->setEnabled(false);
  }

  void enableExtremeValue()
  {
    dock_widget->extremeValueGroup->setEnabled(true);
    dock_widget->zoomToMaxButton->setEnabled(true);
  }

  void resetExtremeValue()
  {
    dock_widget->maxBox->setRange(0, 99999999);
    dock_widget->maxBox->setValue(0);
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

    const double min_value = 0;
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

private Q_SLOTS:
  // Called when new geometric objects are selected in the scene
  void onItemIndicesSelected(QList<int> selected_items)
  {
    resetExtremeValue();
    dock_widget->setEnabled(false);

    Scene_surface_mesh_item* sm_item = nullptr;
    Scene_polyhedron_selection_item* source_vertices = nullptr;

    if(selected_items.size() == 1)
    {
      Scene_item* item = scene->item(scene->mainSelectionIndex());
      source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item);
      if(source_vertices)
      {
        // While selecting a selection item, enable coloring if the selection item is linked to a sm_item
        if(item_source_vertices.right.count(source_vertices) == 0)
          return;

        dock_widget->setEnabled(true);
        dock_widget->createSourceVerticesButton->setEnabled(false);
        dock_widget->estimateDistancesButton->setEnabled(true);
        disableExtremeValue();
        return;
      }
      else if(qobject_cast<Scene_heat_item*>(item))
      {
        dock_widget->setEnabled(true);
        dock_widget->createSourceVerticesButton->setEnabled(false);
        dock_widget->estimateDistancesButton->setEnabled(false);
        disableExtremeValue();
        return;
      }

      sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    }
    else if(selected_items.size() == 2)
    {
      Scene_item* item1 = scene->item(selected_items.front());
      Scene_item* item2 = scene->item(selected_items.back());
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item1);
      source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item2);
      if(!sm_item)
      {
        sm_item = qobject_cast<Scene_surface_mesh_item*>(item2);
        source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item1);
      }
    }

    if(!sm_item)
      return;

    dock_widget->setEnabled(true);
    const bool has_sources = (source_vertices || item_source_vertices.left.count(sm_item) != 0);
    dock_widget->estimateDistancesButton->setEnabled(has_sources);
    dock_widget->createSourceVerticesButton->setEnabled(!has_sources);
    disableExtremeValue();
  }

  // This function is only called if the index actually changed (doesn't trigger
  // if you click again the item)
  void onNewMethodSelected(int method_index)
  {
    resetExtremeValue(); // reset extreme value before the legend to get the proper values
    displayRampLegend();

    if(method_index >= 0 && method_index < dock_widget->methodBox->count()) // valid method
    {
      dock_widget->setEnabled(true);
      disableExtremeValue(); // only available after displaying geodesic distances
    }
    else // no or broken method?
    {
      dock_widget->setEnabled(false);
      dock_widget->methodBox->setEnabled(true);
    }
  }

private:
  bool displayHeatIntensity(Scene_surface_mesh_item* sm_item,
                            Scene_polyhedron_selection_item* source_vertices,
                            const bool use_iDT = false)
  {
    SMesh& mesh = *sm_item->face_graph();

    SMesh::Property_map<vertex_descriptor, double> heat_intensity =
      mesh.add_property_map<vertex_descriptor, double>("v:HM_Plugin_heat_intensity", 0).first;

    auto initialize_hm_map = [this, sm_item, &mesh] (auto*& hm_ptr, auto& hm_map) -> void
                             {
                               using Method = std::decay_t<decltype(*hm_ptr)>;

                               auto it = hm_map.find(sm_item);
                               if(it != hm_map.end()) // method already exists
                               {
                                 hm_ptr = it->second;

                                 for(vertex_descriptor v : vertices(mesh))
                                    hm_ptr->remove_source(v);
                               }
                               else
                               {
                                 hm_ptr = new Method(mesh);
                                 hm_map[sm_item] = hm_ptr;
                               }

                               connect(sm_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
                                       this, [this, sm_item, &hm_map]()
                                             {
                                               item_source_vertices.left.erase(sm_item);

                                               auto it = hm_map.find(sm_item);
                                               if(it == hm_map.end())
                                                 return;
                                               delete it->second;
                                               hm_map.erase(it);
                                             });
                             };

    Heat_method* hm = nullptr;
    Heat_method_idt* hm_idt = nullptr;

    if(use_iDT)
      initialize_hm_map(hm_idt, idt_heat_methods);
    else
      initialize_hm_map(hm, heat_methods);

    for(auto v : source_vertices->selected_vertices)
    {
      if(use_iDT)
        hm_idt->add_source(v);
      else
        hm->add_source(v);
    }

    if(use_iDT)
      hm_idt->estimate_geodesic_distances(heat_intensity);
    else
      hm->estimate_geodesic_distances(heat_intensity);

    // Post treatment
    double max = 0;
    for(vertex_descriptor v : vertices(mesh))
    {
      double hi = heat_intensity[v];
      if(hi > max)
        max = hi;
    }

    displayRampLegend();

    auto [vcolors, vcolors_added] = mesh.add_property_map<vertex_descriptor, CGAL::IO::Color >("v:color", CGAL::IO::Color());
    for(vertex_descriptor v : vertices(mesh))
    {
      double h = heat_intensity[v] / max;
      CGAL::IO::Color color(255 * color_ramp.r(h),
                            255 * color_ramp.g(h),
                            255 * color_ramp.b(h));
      vcolors[v] = color;
    }

    // Create the colored item
    Scene_heat_item* heat_item = new Scene_heat_item(sm_item);
    heat_item->setName(tr("%1 (distance isolines)").arg(sm_item->name()));
    heat_item->setVisible(false);

    sm_item->invalidateOpenGLBuffers();
    sm_item->setRenderingMode(GouraudPlusEdges);
    sm_item->redraw();

    auto heat_item_id = scene->addItem(heat_item);
    scene->setSelectedItem(heat_item_id);

    // any change of sm_item destroys everything

    // @todo do not emit itemChanged when the colors are reset
    // @todo with qt6, single connection can be performed with `static_cast<Qt::ConnectionType>(Qt::SingleShotConnection)`
    // see https://www.kdab.com/single-shot-connections/
    auto connection = std::make_shared<QMetaObject::Connection>();
    *connection = connect(sm_item, &Scene_surface_mesh_item::itemChanged,
                          this, [this, sm_item, heat_item, connection]()
                                {
                                  QObject::disconnect(*connection);

                                  sm_item->resetColors();
                                  removePluginProperties(sm_item);
                                  scene->erase(scene->item_id(heat_item));
                                  onItemIndicesSelected(scene->selectionIndices());
                                });

    connect(sm_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
            this, [this, heat_item]()
                  {
                    scene->erase(scene->item_id(heat_item));
                    onItemIndicesSelected(scene->selectionIndices());
                  });

    // here because if it's put above, the setSelectedItem() might reset the value
    dock_widget->maxBox->setValue(max);

    return true;
  }

private Q_SLOTS:
  void estimateDistances()
  {
    // Get the mesh and source vertices items
    Scene_surface_mesh_item* sm_item = nullptr;
    Scene_polyhedron_selection_item* source_vertices = nullptr;

    if(scene->selectionIndices().size() == 1)
    {
      Scene_item* item = scene->item(scene->mainSelectionIndex());
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
      if(sm_item)
      {
        // a surface mesh item is selected, an existing associated selection item must exist
        source_vertices = item_source_vertices.left.at(sm_item);
      }
      else
      {
        // a selection item is selected, an existing associated mesh item must exist
        source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item);
        if(source_vertices)
          sm_item = item_source_vertices.right.at(source_vertices);
      }
    }
    else if(scene->selectionIndices().size() == 2)
    {
      // two items, for (possibly unlinked) sm_item and its associated selection
      Scene_item* item1 = scene->item(scene->selectionIndices().front());
      Scene_item* item2 = scene->item(scene->selectionIndices().back());
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item1);
      source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item2);
      if(!sm_item)
      {
        sm_item = qobject_cast<Scene_surface_mesh_item*>(item2);
        source_vertices = qobject_cast<Scene_polyhedron_selection_item*>(item1);
      }

      link_mesh_and_selection(sm_item, source_vertices);
    }
    else
    {
      QMessageBox::critical(mw, "Error","Unsupported selection of items.");
      return;
    }

    CGAL_assertion(sm_item && source_vertices);

    if(!is_triangle_mesh(*sm_item->face_graph()))
    {
      QApplication::restoreOverrideCursor();
      QMessageBox::critical(mw, "Error","The mesh must be triangulated.");
      return;
    }

    if(source_vertices->selected_vertices.empty())
    {
      QApplication::restoreOverrideCursor();
      QMessageBox::critical(mw, "Error","At least one source vertex is required.");
      return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);

    enableExtremeValue();

    const std::string& method_name = dock_widget->methodBox->currentText().toStdString();
    if(method_name == "Heat Method")
      displayHeatIntensity(sm_item, source_vertices);
    else if(method_name == "Heat Method (Intrinsic Delaunay)")
      displayHeatIntensity(sm_item, source_vertices, true /*use IDT*/);

    // @todo emit a new SIGNAL on successful coloring, something like "colorChanged()"
    // itemChanged is too strong and would conflict with the connection below
    sm_item->invalidateOpenGLBuffers();
    sm_item->redraw();

    QApplication::restoreOverrideCursor();
  }

private:
  void removePluginProperty(Scene_item* item,
                            const std::string& property_name)
  {
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(!sm_item)
      return;

    SMesh* sm = sm_item->face_graph();
    if(sm == nullptr)
      return;

    // Here we only target the property maps added by this plugin, so 'double' is fine
    std::optional<SMesh::Property_map<face_descriptor, double>> property = sm->property_map<face_descriptor, double>(property_name);
    if(property.has_value())
      sm->remove_property_map(property.value());
  }

  void removePluginProperties(Scene_item* item)
  {
    removePluginProperty(item, "v:HM_Plugin_heat_intensity");
  }

private Q_SLOTS:
  // deletion of the selection item removes the pair from the map
  void link_mesh_and_selection(Scene_surface_mesh_item* sm_item,
                               Scene_polyhedron_selection_item* source_vertices)
  {
    item_source_vertices.left.insert(std::make_pair(sm_item, source_vertices));

    connect(source_vertices, &Scene_polyhedron_selection_item::aboutToBeDestroyed,
            this, [this, sm_item]()
                  {
                    item_source_vertices.left.erase(sm_item);
                    onItemIndicesSelected(scene->selectionIndices());
                  });
  }

  void createSourceVerticesSelectionItem()
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(!sm_item)
    {
      QMessageBox::warning(mw, "Warning", "Select a surface mesh to add source vertices");
      dock_widget->createSourceVerticesButton->setChecked(false);
      return;
    }

    CGAL_assertion(item_source_vertices.left.count(sm_item) == 0);

    Scene_polyhedron_selection_item* source_vertices = new Scene_polyhedron_selection_item(sm_item, mw);
    source_vertices->setName(tr("%1 (source vertices)").arg(sm_item->name()));
    auto source_vertices_id = scene->addItem(source_vertices);
    scene->setSelectedItem(source_vertices_id);

    link_mesh_and_selection(sm_item, source_vertices);

    dock_widget->createSourceVerticesButton->setEnabled(false);
    dock_widget->estimateDistancesButton->setEnabled(true);
  }

private Q_SLOTS:
  void on_zoomToMaxButton_pressed()
  {
    Scene_surface_mesh_item* sm_item = nullptr;
    if(scene->selectionIndices().size() == 1)
    {
      Scene_item* item = scene->item(scene->mainSelectionIndex());
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    }
    else if(scene->selectionIndices().size() == 2)
    {
      Scene_item* item1 = scene->item(scene->selectionIndices().front());
      Scene_item* item2 = scene->item(scene->selectionIndices().back());
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item1);
      if(!sm_item)
        sm_item = qobject_cast<Scene_surface_mesh_item*>(item2);
    }

    const SMesh& mesh = *(sm_item->face_graph());

    auto heat_intensity = mesh.property_map<vertex_descriptor, double>("v:HM_Plugin_heat_intensity").value();

    double max = 0;
    vertex_descriptor max_v = boost::graph_traits<SMesh>::null_vertex();
    for(vertex_descriptor v : vertices(mesh))
    {
      if(heat_intensity[v] > max)
      {
        max = heat_intensity[v];
        max_v = v;
      }
    }

    CGAL_assertion(max_v != boost::graph_traits<SMesh>::null_vertex());

    face_descriptor unused_fd;
    Point_3 unused_p;
    ::zoomToId(mesh,
               QString("v%1").arg(static_cast<int>(max_v)),
               getActiveViewer(),
               unused_fd, unused_p);
  }
};

#include "Heat_method_plugin.moc"
