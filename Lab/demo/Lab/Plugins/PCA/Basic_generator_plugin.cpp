#include "ui_Basic_generator_widget.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/subdivision_method_3.h>

#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QMessageBox>
#include <QBitmap>
#include <QTabBar>

#include <array>
#include <iterator>
#include <list>
#include <vector>

using namespace CGAL::Three;

using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;

struct Face
  : public std::array<int,3>
{
  Face(int i, int j, int k)
  {
    (*this)[0] = i;
    (*this)[1] = j;
    (*this)[2] = k;
  }
};

class GeneratorWidget
  : public QDockWidget,
    public Ui::BasicGenerator
{
  Q_OBJECT

  public:
  GeneratorWidget(const QString& name, QWidget *parent)
    : QDockWidget(name, parent)
  {
    setupUi(this);
  }
};

class Q_DECL_EXPORT Basic_generator_plugin
  : public QObject,
    public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "basic_generator_plugin.json")

  enum Basic_objects
  {
    POINT_SET = 0,
    POLYLINE,
    GRID,
    SPHERE,
    TETRAHEDRON,
    HEXAHEDRON,
    PRISM,
    PYRAMID
  };

  static constexpr int object_types_n = 8;
  int instances[object_types_n];

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    for (int i=0; i<object_types_n; ++i)
      instances[i] = 0;

    QMenu* menuFile = mw->findChild<QMenu*>("menuFile");

    QMenu* menu = menuFile->findChild<QMenu*>("menuGenerateObject");
    if (!menu) {
      QAction* actionLoad = mw->findChild<QAction*>("actionLoadPlugin");
      menu = new QMenu(tr("Generate &Object"), menuFile);
      menu->setObjectName("menuGenerateObject");
      menuFile->insertMenu(actionLoad, menu);
    }

    QAction* actionPointSet    = new QAction("Po&int Set", mw);
    QAction* actionPolyline    = new QAction("Po&lyline", mw);
    QAction* actionGrid        = new QAction("&Grid", mw);
    QAction* actionSphere      = new QAction("&Sphere", mw);
    QAction* actionTetrahedron = new QAction("&Tetrahedron", mw);
    QAction* actionHexahedron  = new QAction("&Hexahedron", mw);
    QAction* actionPrism       = new QAction("P&rism", mw);
    QAction* actionPyramid     = new QAction("Py&ramid", mw);

    connect(actionPointSet, SIGNAL(triggered()),
            this, SLOT(on_actionPointSet_triggered()));
    _actions << actionPointSet;

    connect(actionPolyline, SIGNAL(triggered()),
            this, SLOT(on_actionPolyline_triggered()));
    _actions << actionPolyline;

    connect(actionGrid, SIGNAL(triggered()),
            this, SLOT(on_actionGrid_triggered()));
    _actions << actionGrid;

    connect(actionSphere, SIGNAL(triggered()),
            this, SLOT(on_actionSphere_triggered()));
    _actions << actionSphere;

    connect(actionTetrahedron, SIGNAL(triggered()),
            this, SLOT(on_actionTetrahedron_triggered()));
    _actions << actionTetrahedron;

    connect(actionHexahedron, SIGNAL(triggered()),
            this, SLOT(on_actionHexahedron_triggered()));
    _actions << actionHexahedron;

    connect(actionPrism, SIGNAL(triggered()),
            this, SLOT(on_actionPrism_triggered()));
    _actions << actionPrism;

    connect(actionPyramid, SIGNAL(triggered()),
            this, SLOT(on_actionPyramid_triggered()));
    _actions << actionPyramid;

    for (QAction* action : _actions) {
      menu->addAction(action);
    }

    dock_widget = new GeneratorWidget("Basic Objects", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);

    connect(dock_widget->generateButton, &QAbstractButton::clicked,
            this, &Basic_generator_plugin::on_generate_clicked);
    connect(dock_widget->selector_tabWidget, &QTabWidget::currentChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget, &GeneratorWidget::visibilityChanged,
            this, &Basic_generator_plugin::on_tab_changed);
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
    connect(dock_widget->prismCheckBox, &QCheckBox::checkStateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget->pyramidCheckBox, &QCheckBox::checkStateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
#else
    connect(dock_widget->prismCheckBox, &QCheckBox::stateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
    connect(dock_widget->pyramidCheckBox, &QCheckBox::stateChanged,
            this, &Basic_generator_plugin::on_tab_changed);
#endif
    connect(dock_widget->polygon_checkBox, SIGNAL(toggled(bool)),
            dock_widget->fill_checkBox, SLOT(setEnabled(bool)));
  }

  bool applicable(QAction*) const
  {
    // not in the Operations menu
    return false;
  }

  QList<QAction*> actions() const {
    return QList<QAction*>();
  }

public Q_SLOTS:
  void on_actionPointSet_triggered();
  void on_actionPolyline_triggered();
  void on_actionGrid_triggered();
  void on_actionSphere_triggered();
  void on_actionTetrahedron_triggered();
  void on_actionHexahedron_triggered();
  void on_actionPrism_triggered();
  void on_actionPyramid_triggered();

  void on_generate_clicked();
  void on_tab_changed();

  void closure() { dock_widget->hide(); }

private:
  QList<QAction*> _actions;
  GeneratorWidget* dock_widget;

  void generatePoints();
  void generateLines();
  void generateGrid();
  void generateSphere();
  void generateTetrahedron();
  void generateHexahedron();
  void generatePrism();
  void generatePyramid();
}; // class Basic_generator_plugin

// show the widget
void Basic_generator_plugin::on_actionPointSet_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(POINT_SET);
}
void Basic_generator_plugin::on_actionPolyline_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(POLYLINE);
}
void Basic_generator_plugin::on_actionGrid_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(GRID);
}
void Basic_generator_plugin::on_actionSphere_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(SPHERE);
}
void Basic_generator_plugin::on_actionTetrahedron_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(TETRAHEDRON);
}
void Basic_generator_plugin::on_actionHexahedron_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(HEXAHEDRON);
}
void Basic_generator_plugin::on_actionPrism_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(PRISM);
}
void Basic_generator_plugin::on_actionPyramid_triggered()
{
  dock_widget->show();
  dock_widget->raise();
  dock_widget->selector_tabWidget->tabBar()->setCurrentIndex(PYRAMID);
}

void Basic_generator_plugin::on_tab_changed()
{
  QString name;
  int nb = 0;
  switch(dock_widget->selector_tabWidget->currentIndex())
  {
    case POINT_SET:
    {
      name = QString("Point_set");
      nb = instances[POINT_SET];
      break;
    }
    case POLYLINE:
    {
      name = QString("Polyline");
      nb = instances[POLYLINE];
      break;
    }
    case GRID:
    {
      name = QString("Grid");
      nb = instances[GRID];
      break;
    }
    case SPHERE:
    {
      name = QString("Sphere");
      nb = instances[SPHERE];
      break;
    }
    case TETRAHEDRON:
    {
      name = QString("Tetrahedron");
      nb = instances[TETRAHEDRON];
      break;
    }
    case HEXAHEDRON:
    {
      name = QString("Hexahedron");
      nb = instances[HEXAHEDRON];
      break;
    }
    case PRISM:
    {
      name = QString("Prism");
      nb = instances[PRISM];
      QPixmap pic;
      if (dock_widget->prismCheckBox->isChecked())
        pic = QPixmap(":/cgal/Lab/resources/prism.png");
      else
        pic = QPixmap(":/cgal/Lab/resources/prism-open.png");
      dock_widget->prism_picLabel->setPixmap(pic);
      dock_widget->prism_picLabel->show();
      break;
    }
    case PYRAMID:
    {
      name = QString("Pyramid");
      nb = instances[PYRAMID];
      QPixmap pic;
      if (dock_widget->pyramidCheckBox->isChecked())
        pic = QPixmap(":/cgal/Lab/resources/pyramid.png");
      else
        pic = QPixmap(":/cgal/Lab/resources/pyramid-open.png");
      dock_widget->pyramid_picLabel->setPixmap(pic);
      dock_widget->pyramid_picLabel->show();
      break;
    }
    default:
      break;
  };
  dock_widget->name_lineEdit->setText(name + "_" + QString::number(nb));
}
//generate
void Basic_generator_plugin::on_generate_clicked()
{
  switch(dock_widget->selector_tabWidget->currentIndex())
  {
  case POINT_SET:
    generatePoints();
    break;
  case POLYLINE:
    generateLines();
    break;
  case GRID:
    generateGrid();
    break;
  case SPHERE:
    generateSphere();
    break;
  case TETRAHEDRON:
    generateTetrahedron();
    break;
  case HEXAHEDRON:
    generateHexahedron();
    break;
  case PRISM:
    generatePrism();
    break;
  case PYRAMID:
    generatePyramid();
    break;
  default:
    break;
  };
  on_tab_changed();
}

void Basic_generator_plugin::generatePoints()
{
  QString text = dock_widget->point_textEdit->toPlainText();
  text.replace(',',' ');

  QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty())
    return;

  const bool is_2d = dock_widget->dim2_points_checkBox->isChecked();

  if (is_2d) {
    if (list.size() % 2 != 0) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of pairs.");
      msgBox->exec();
      return;
    }
  } else {
    if (list.size() % 3 != 0) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of triplets.");
      msgBox->exec();
      return;
    }
  }

  Scene_points_with_normal_item* item = new Scene_points_with_normal_item();

  double coord[3];
  int counter = 0;

  for (const QString& s : list) {
    if (!s.isEmpty()) {
      bool ok;
      double res = s.toDouble(&ok);
      if (!ok) {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR: Coordinates are invalid: " + s);
        msgBox->exec();
        return;
      } else {
        coord[counter] = res;
        ++counter;
      }
    }

    if (is_2d) {
      if (counter == 2) {
        const Point p(coord[0], coord[1], 0.0);
        item->point_set()->insert(p);
        counter = 0;
      }
    } else {
      if (counter == 3) {
        const Point p(coord[0], coord[1], coord[2]);
        item->point_set()->insert(p);
        counter = 0;
      }
    }
  }

  dock_widget->point_textEdit->clear();
  item->point_set()->unselect_all();
  item->setName(dock_widget->name_lineEdit->text());
  item->setColor(Qt::black);
  item->invalidateOpenGLBuffers();
  Scene_interface::Item_id id = scene->addItem(item);
  scene->setSelectedItem(id);

  ++instances[POINT_SET];
}

void Basic_generator_plugin::generateLines()
{
  QString text = dock_widget->line_textEdit->toPlainText();
  text.replace(',',' ');

  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;

  auto read_polyline = [&polylines](const QStringList& list,
                                    const bool is_2d,
                                    const bool is_closed) -> bool
  {
    int counter = -1;
    double coord[3];
    bool ok = true;

    if (!is_2d && (list.size() % 3 != 0)) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of triplets.");
      msgBox->exec();
      return false;
    } else if (is_2d && (list.size() % 2 != 0)) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of pairs.");
      msgBox->exec();
      return false;
    }

    polylines.back().reserve(list.size() + (is_closed ? 1 : 0));

    for (const QString& s : list) {
      if (!s.isEmpty()) {
        double res = s.toDouble(&ok);
        if (!ok) {
          QMessageBox *msgBox = new QMessageBox;
          msgBox->setWindowTitle("Error");
          msgBox->setText("ERROR: Coordinates are invalid: " + s);
          msgBox->exec();
          return false;
        } else {
          coord[++counter] = res;
        }
      }

      if (!is_2d && counter == 2) {
        Scene_polylines_item::Point_3 p(coord[0], coord[1], coord[2]);
        polylines.back().push_back(p);
        counter = -1;
      } else if (is_2d && counter == 1) {
        Scene_polylines_item::Point_3 p(coord[0], coord[1], 0);
        polylines.back().push_back(p);
        counter=-1;
      }
    }

    // close if not already closed
    if (is_closed)
    {
      if (!polylines.back().empty() && polylines.back().back()!=polylines.back().front())
        polylines.back().push_back(polylines.back().front());
    }

    return true;
  };

  const bool is_2d = dock_widget->dim2_polylines_checkBox->isChecked();
  const bool is_closed = dock_widget->polygon_checkBox->isChecked();
  const bool shall_fill = is_closed && dock_widget->fill_checkBox->isChecked();
  const bool oneperline = dock_widget->oneperline_checkBox->isChecked();

  bool ok = true;
  if (oneperline) {
    QStringList poly_list = text.split(QRegularExpression("\\n"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (poly_list.empty())
      return;

    for (const QString& qs : poly_list) {
      QStringList list = qs.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
      if (list.isEmpty())
        continue;
      polylines.emplace_back();
      ok = read_polyline(list, is_2d, is_closed);
      if (!ok)
        return;
    }
  } else {
    QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (list.isEmpty())
      return;
    polylines.emplace_back();
    ok = read_polyline(list, is_2d, is_closed);
    if (!ok)
      return;
  }

  if (ok) {
    dock_widget->line_textEdit->clear();
    if (shall_fill) {
      CGAL::Three::Three::CursorScopeGuard guard(Qt::WaitCursor);
      QApplication::processEvents();

      SMesh* poly = new SMesh;
      for (const auto& polyline : polylines) {
        if (polyline.size() < 4) { // no triangle, skip it (needs at least 3 + 1 repeat)
          QMessageBox::warning(mw, "Warning", "Needs at least 3 points to triangulate. Aborting.");
          delete poly;
          return;
        }

        std::vector<Face> patch;
        CGAL::Polygon_mesh_processing::triangulate_hole_polyline(polyline,
                                                                 std::back_inserter(patch),
                                                                 CGAL::parameters::use_delaunay_triangulation(true));

        if (patch.empty()) {
          QMessageBox::warning(mw, "Warning", "Triangulation failed.");
          delete poly;
          return;
        }

        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(polyline, patch, *poly);
      }

      Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item(poly);
      poly_item->setName(dock_widget->name_lineEdit->text());
      poly_item->setRenderingMode(FlatPlusEdges);
      scene->setSelectedItem(scene->addItem(poly_item));

      ++instances[POLYLINE];
    }
    else
    {
      Scene_polylines_item* item = new Scene_polylines_item();
      item->polylines = polylines;
      item->invalidateOpenGLBuffers();
      item->setName(dock_widget->name_lineEdit->text());
      item->setColor(Qt::black);
      QStringList polylines_metadata;
      item->setProperty("polylines metadata", polylines_metadata);
      Scene_interface::Item_id id = scene->addItem(item);
      scene->setSelectedItem(id);

      ++instances[POLYLINE];
    }
  }
}

void Basic_generator_plugin::generateGrid()
{
  struct Point_generator
  {
    std::size_t w, h;
    Point origin;
    Vector u, v;

    Point_generator(const std::size_t& w,
                    const std::size_t& h,
                    const Point& origin,
                    const Vector& u,
                    const Vector& v)
      : w(w), h(h), origin(origin), u(u), v(v)
    {}

    Point operator()(std::size_t i, std::size_t j) const
    {
      double alpha = (w == 1) ? 0.0 : double(i) / double(w - 1);
      double beta  = (h == 1) ? 0.0 : double(j) / double(h - 1);
      return origin + alpha * u + beta * v;
    }
  };

  typedef Scene_surface_mesh_item::Face_graph Face_graph;

  QString points_text;
  points_text = dock_widget->grid_lineEdit->text();
  points_text.replace(',',' ');

  QStringList list = points_text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty())
    return;

  if (list.size() != 9) {
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR: Input should consist of 9 values (origin and two points for axes).");
    msgBox->exec();
    return;
  }

  double coords[9];
  for (int j=0; j<9; ++j) {
    bool ok;
    coords[j] = list.at(j).toDouble(&ok);
    if (!ok) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Coordinates are invalid: " + list.at(j));
      msgBox->exec();
      return;
    }
  }

  Point pts[3];
  pts[0] = Point(coords[0], coords[1], coords[2]);
  pts[1] = Point(coords[3], coords[4], coords[5]);
  pts[2] = Point(coords[6], coords[7], coords[8]);

  const Vector u = pts[1] - pts[0];
  const Vector v = pts[2] - pts[0];

  using size_type = typename boost::graph_traits<Face_graph>::vertices_size_type;

  size_type nb_cells[2];
  nb_cells[0] = static_cast<size_type>(dock_widget->gridX_spinBox->value());
  nb_cells[1] = static_cast<size_type>(dock_widget->gridY_spinBox->value());

  if (nb_cells[0] < 1 || nb_cells[1] < 1) {
    QMessageBox* msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR: Number of grid cells must be at least 1 in each direction.");
    msgBox->exec();
    return;
  }

  // nb_points = nb_cells + 1
  Point_generator point_gen(nb_cells[0]+1, nb_cells[1]+1, pts[0], u, v);
  const bool triangulated = dock_widget->grid_checkBox->isChecked();
  Face_graph grid;
  CGAL::make_grid(nb_cells[0], nb_cells[1], grid, point_gen, triangulated);

  Scene_surface_mesh_item* grid_item = new Scene_surface_mesh_item(grid);
  grid_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(grid_item);
  scene->setSelectedItem(id);

  ++instances[GRID];
}

void Basic_generator_plugin::generateSphere()
{
  QString text = dock_widget->center_radius_lineEdit->text();
  text.replace(',',' ');
  QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty())
    return;

  if (list.size() != 4) {
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR: Input should consist of 4 values.");
    msgBox->exec();
    return;
  }

  bool ok = true;
  double coords[4];
  for (int i=0; i<4; ++i) {
    coords[i] = list.at(i).toDouble(&ok);
    if (!ok) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Coordinates are invalid: " + list.at(i));
      msgBox->exec();
      return;
    }
  }

  const double radius = coords[3];
  const Point center(coords[0], coords[1], coords[2]);
  Scene_surface_mesh_item::Face_graph sphere;
  make_icosahedron(sphere, center, radius);

  int precision = dock_widget->SphereSpinBox->value();
  if (precision != 0) {
    CGAL::Subdivision_method_3::Sqrt3_subdivision(sphere, CGAL::parameters::number_of_iterations(precision));
  }

  // project the points onto the sphere
  auto vpm = get(CGAL::vertex_point, sphere);
  for (auto vd : vertices(sphere)) {
    Vector vec(center, get(vpm, vd));
    vec = radius * vec / CGAL::sqrt(vec.squared_length());
    put(vpm, vd, Point(center.x() + vec.x(),
                       center.y() + vec.y(),
                       center.z() + vec.z()));
  }

  Scene_surface_mesh_item* sphere_item = new Scene_surface_mesh_item(sphere);
  sphere_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(sphere_item);
  scene->setSelectedItem(id);

  ++instances[SPHERE];
}

void Basic_generator_plugin::generateTetrahedron()
{
  Scene_surface_mesh_item::Face_graph tetrahedron;

  if (dock_widget->tabWidget->currentIndex() == 0) {
    QString point_texts[4];
    Point points[4];
    point_texts[0] = dock_widget->tetP0->text();
    point_texts[1] = dock_widget->tetP1->text();
    point_texts[2] = dock_widget->tetP2->text();
    point_texts[3] = dock_widget->tetP3->text();

    for (int i = 0; i < 4; ++i) {
      point_texts[i].replace(',',' ');
      QStringList list = point_texts[i].split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
      if (list.isEmpty())
        return;

      if (list.size() != 3) {
        QMessageBox* msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR: Input should consist of 3 values.");
        msgBox->exec();
        return;
      }

      double coords[3];
      for (int j=0; j<3; ++j) {
        bool ok;
        coords[j] = list.at(j).toDouble(&ok);
        if (!ok) {
          QMessageBox* msgBox = new QMessageBox;
          msgBox->setWindowTitle("Error");
          msgBox->setText("ERROR: Coordinates are invalid: " + list.at(j));
          msgBox->exec();
          return;
        }
      }
      points[i] = Point(coords[0], coords[1], coords[2]);
    }

    CGAL::make_tetrahedron(points[0], points[1], points[2], points[3], tetrahedron);
  }
  else
  {
    QString text = dock_widget->point_textEdit_2->toPlainText();
    text.replace(',',' ');

    QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (list.isEmpty())
      return;

    if (list.size() != 12) {
      QMessageBox* msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of 12 values.");
      msgBox->exec();
      return;
    }

    double coords[12];
    for (int i=0; i<12; ++i) {
      bool ok;
      coords[i] = list.at(i).toDouble(&ok);
      if (!ok) {
        QMessageBox* msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR: Coordinates are invalid: " + list.at(i));
        msgBox->exec();
        return;
      }
    }

    CGAL::make_tetrahedron(Point(coords[0], coords[1], coords[2]),
                           Point(coords[3], coords[4], coords[5]),
                           Point(coords[6], coords[7], coords[8]),
                           Point(coords[9], coords[10], coords[11]),
                           tetrahedron);
  }

  Scene_surface_mesh_item* tet_item = new Scene_surface_mesh_item(tetrahedron);
  tet_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(tet_item);
  scene->setSelectedItem(id);

  ++instances[TETRAHEDRON];
}

void Basic_generator_plugin::generateHexahedron()
{
  Scene_surface_mesh_item::Face_graph cube;

  if (dock_widget->tabWidget_2->currentIndex() == 0) {
    QString point_texts[8];
    point_texts[0] = dock_widget->cubeP0->text();
    point_texts[1] = dock_widget->cubeP1->text();
    point_texts[2] = dock_widget->cubeP2->text();
    point_texts[3] = dock_widget->cubeP3->text();
    point_texts[4] = dock_widget->cubeP4->text();
    point_texts[5] = dock_widget->cubeP5->text();
    point_texts[6] = dock_widget->cubeP6->text();
    point_texts[7] = dock_widget->cubeP7->text();

    Point points[8];
    for (int i=0; i<8; ++i) {
      point_texts[i].replace(',',' ');

      QStringList list = point_texts[i].split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
      if (list.isEmpty())
        return;

      if (list.size() != 3) {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR: Input should consist of 3 values.");
        msgBox->exec();
        return;
      }

      double coords[3];
      for (int j=0; j<3; ++j) {
        bool ok;
        coords[j] = list.at(j).toDouble(&ok);
        if (!ok) {
          QMessageBox *msgBox = new QMessageBox;
          msgBox->setWindowTitle("Error");
          msgBox->setText("ERROR: Coordinates are invalid: " + list.at(j));
          msgBox->exec();
          return;
        }
      }
      points[i] = Point(coords[0], coords[1], coords[2]);
    }

    CGAL::make_hexahedron(points[0], points[1], points[2], points[3],
                          points[4], points[5], points[6], points[7], cube);
  }
  else
  {
    QString text = dock_widget->extremaEdit->text();
    text.replace(',',' ');

    QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
    if (list.isEmpty())
      return;

    if (list.size() != 6) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Input should consist of 6 values.");
      msgBox->exec();
      return;
    }

    double coords[6];
    for (int i=0; i<6; ++i) {
      bool ok;
      coords[i] = list.at(i).toDouble(&ok);
      if (!ok) {
        QMessageBox *msgBox = new QMessageBox;
        msgBox->setWindowTitle("Error");
        msgBox->setText("ERROR: Coordinates are invalid: " + list.at(i));
        msgBox->exec();
        return;
      }
    }

    CGAL::make_hexahedron(Point(coords[0], coords[1], coords[5]),
                          Point(coords[3], coords[1], coords[5]),
                          Point(coords[3], coords[1], coords[2]),
                          Point(coords[0], coords[1], coords[2]),
                          Point(coords[0], coords[4], coords[2]),
                          Point(coords[0], coords[4], coords[5]),
                          Point(coords[3], coords[4], coords[5]),
                          Point(coords[3], coords[4], coords[2]),
                          cube);
  }

  Scene_surface_mesh_item* hex_item = new Scene_surface_mesh_item(cube);
  hex_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(hex_item);
  scene->setSelectedItem(id);

  ++instances[HEXAHEDRON];
}

void Basic_generator_plugin::generatePrism()
{
  QString text = dock_widget->prism_lineEdit->text();
  text.replace(',',' ');

  QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty())
    return;

  if (list.size() != 3) {
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR: Input should consist of 3 values.");
    msgBox->exec();
    return;
  }

  double coords[3];
  for (int i=0; i<3; ++i) {
    bool ok;
    coords[i] = list.at(i).toDouble(&ok);
    if (!ok) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Coordinates are invalid: " + list.at(i));
      msgBox->exec();
      return;
    }
  }

  const int nb_vertices = dock_widget->prismSpinBox->value();
  const double height = dock_widget->prismHeightSpinBox->value();
  const double radius = dock_widget->prismBaseSpinBox->value();
  const bool is_closed = dock_widget->prismCheckBox->isChecked();

  Scene_surface_mesh_item::Face_graph prism;
  make_regular_prism(nb_vertices,
                     prism,
                     Point(coords[0], coords[1], coords[2]),
                     height,
                     radius,
                     is_closed);

  Scene_surface_mesh_item* prism_item = new Scene_surface_mesh_item(prism);
  prism_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(prism_item);
  scene->setSelectedItem(id);

  ++instances[PRISM];
}

void Basic_generator_plugin::generatePyramid()
{
  QString text = dock_widget->pyramid_lineEdit->text();
  text.replace(',',' ');

  QStringList list = text.split(QRegularExpression("\\s+"), CGAL_QT_SKIP_EMPTY_PARTS);
  if (list.isEmpty())
    return;

  if (list.size() != 3) {
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR: Input should consist of 3 values.");
    msgBox->exec();
    return;
  }

  double coords[3];
  for (int i=0; i<3; ++i) {
    bool ok;
    coords[i] = list.at(i).toDouble(&ok);
    if (!ok) {
      QMessageBox *msgBox = new QMessageBox;
      msgBox->setWindowTitle("Error");
      msgBox->setText("ERROR: Coordinates are invalid: " + list.at(i));
      msgBox->exec();
      return;
    }
  }

  const int nb_vertices = dock_widget->pyramidSpinBox->value();
  const double height = dock_widget->pyramidHeightSpinBox->value();
  const double radius = dock_widget->pyramidBaseSpinBox->value();
  const bool is_closed = dock_widget->pyramidCheckBox->isChecked();

  Scene_surface_mesh_item::Face_graph pyramid;
  make_pyramid(nb_vertices,
               pyramid,
               Point(coords[0], coords[1], coords[2]),
               height,
               radius,
               is_closed);

  Scene_surface_mesh_item* pyramid_item = new Scene_surface_mesh_item(pyramid);
  pyramid_item->setName(dock_widget->name_lineEdit->text());
  Scene_interface::Item_id id = scene->addItem(pyramid_item);
  scene->setSelectedItem(id);

  ++instances[PYRAMID];
}

#include "Basic_generator_plugin.moc"
