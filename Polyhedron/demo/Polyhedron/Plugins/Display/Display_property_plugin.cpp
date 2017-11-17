#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QColorDialog>
#include <QPalette>
#include <QColor>
#include "Messages_interface.h"
#include "Scene_surface_mesh_item.h"
#include "Color_ramp.h"
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "ui_Display_property.h"
#include "id_printing.h"
#include "Scene.h"

#define ARBITRARY_DBL_MIN 1.0E-30
#define ARBITRARY_DBL_MAX 1.0E+30

class DockWidget :
    public QDockWidget,
    public Ui::DisplayPropertyWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};

typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

class DisplayPropertyPlugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:

  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    CGAL::Three::Scene_item* item = scene->item(scene->mainSelectionIndex());
    return qobject_cast<Scene_surface_mesh_item*>(item);
  }

  QList<QAction*> actions() const Q_DECL_OVERRIDE
  {
    return _actions;
  }

  QColor textColor(const QColor& color)
  {
    QColor text_color (255, 255, 255);
    if (color.red() * 0.299 + color.green() * 0.587 + color.blue() * 0.114 > 128)
      text_color = QColor (0, 0, 0);
    return text_color;
  }

  void init(QMainWindow* mw, CGAL::Three::Scene_interface* sc, Messages_interface*) Q_DECL_OVERRIDE
  {
    this->scene = sc;
    this->mw = mw;

    QAction *actionDisplayAngles= new QAction(QString("Display Properties"), mw);

    rm = 1.0;
    rM = 0.0;
    gm = 0.0;
    gM = 1.0;
    bm = 0.0;
    bM = 0.0;
    actionDisplayAngles->setProperty("submenuName", "Color");

    if(actionDisplayAngles) {
      connect(actionDisplayAngles, SIGNAL(triggered()),
              this, SLOT(openDialog()));
      _actions << actionDisplayAngles;

    }
    dock_widget = new DockWidget("Property Displaying", mw);
    dock_widget->setVisible(false);
    addDockWidget(dock_widget);
    QPalette palette(Qt::red);
    dock_widget->minColorButton->setPalette(palette);
    dock_widget->minColorButton->update();

    palette = QPalette(Qt::green);
    dock_widget->maxColorButton->setPalette(palette);
    dock_widget->maxColorButton->update();
    connect(dock_widget->colorizeButton, SIGNAL(clicked(bool)),
            this, SLOT(colorize()));

    connect(dock_widget->rampButton, SIGNAL(clicked(bool)),
            this, SLOT(replaceRamp()));

    connect(dock_widget->propertyBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(on_propertyBox_currentIndexChanged(int)));
    connect(dock_widget->zoomToMinButton, &QPushButton::pressed,
            this, &DisplayPropertyPlugin::on_zoomToMinButton_pressed);
    connect(dock_widget->zoomToMaxButton, &QPushButton::pressed,
            this, &DisplayPropertyPlugin::on_zoomToMaxButton_pressed);
    connect(dock_widget->minColorButton, &QPushButton::pressed,
            this, [this]()
    {
      QColor minColor = QColorDialog::getColor();
      rm = minColor.redF();
      gm = minColor.greenF();
      bm = minColor.blueF();
      QPalette palette(minColor);
      dock_widget->minColorButton->setAutoFillBackground(true);
      dock_widget->minColorButton->setPalette(palette);
      dock_widget->minColorButton->update();
    });
    connect(dock_widget->maxColorButton, &QPushButton::pressed,
            this, [this]()
    {
      QColor maxColor = QColorDialog::getColor();
      QPalette palette(maxColor);
      rM = maxColor.redF();
      gM = maxColor.greenF();
      bM = maxColor.blueF();

      dock_widget->maxColorButton->setAutoFillBackground(true);
      dock_widget->maxColorButton->setPalette(palette);
      dock_widget->maxColorButton->update();
    });
    dock_widget->zoomToMaxButton->setEnabled(false);
    dock_widget->zoomToMinButton->setEnabled(false);
    Scene* scene_obj =static_cast<Scene*>(scene);
    connect(scene_obj, &Scene::itemIndexSelected,
            this, &DisplayPropertyPlugin::enableButtons);
    on_propertyBox_currentIndexChanged(0);

  }
private Q_SLOTS:
  void openDialog()
  {
    if(dock_widget->isVisible()) { dock_widget->hide(); }
    else                         { replaceRamp(); dock_widget->show(); }
  }

  void colorize()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    item->face_graph()->collect_garbage();
    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
      displayAngles(item);
      break;
    default:
      displayScaledJacobian(item);
      break;
    }
    QApplication::restoreOverrideCursor();
    item->invalidateOpenGLBuffers();
    item->redraw();
    dock_widget->zoomToMinButton->setEnabled(true);
    dock_widget->zoomToMaxButton->setEnabled(true);
  }

  void enableButtons()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(! item )
    {
      dock_widget->zoomToMinButton->setEnabled(false);
      dock_widget->zoomToMaxButton->setEnabled(false);
    }

    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
      dock_widget->zoomToMinButton->setEnabled(angles_max.count(item)>0 );
      dock_widget->zoomToMaxButton->setEnabled(angles_max.count(item)>0 );
      break;
    case 1:
      dock_widget->zoomToMinButton->setEnabled(jacobian_max.count(item)>0);
      dock_widget->zoomToMaxButton->setEnabled(jacobian_max.count(item)>0);
      break;
    default:
      break;
    }
  }

  void resetProperty()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(sender());
    if(!item)
      return;
    SMesh& smesh = *item->face_graph();
    SMesh::Property_map<face_descriptor, double> jacobians;
    bool found;
    boost::tie(jacobians, found) = smesh.property_map<face_descriptor,double>("f:jacobian");
    if(found)
    {
      smesh.remove_property_map(jacobians);
    }
    SMesh::Property_map<face_descriptor, double> angles;
    boost::tie(angles, found) = smesh.property_map<face_descriptor,double>("f:angle");
    if(found)
    {
      smesh.remove_property_map(angles);
    }
  }

  void displayScaledJacobian(Scene_surface_mesh_item* item)
  {

    SMesh& smesh = *item->face_graph();
    //compute and store the jacobian per face
    bool non_init;
    SMesh::Property_map<face_descriptor, double> fjacobian;
    boost::tie(fjacobian, non_init) = smesh.add_property_map<face_descriptor, double>("f:jacobian", 0);
    if(non_init)
    {
      double res_min = ARBITRARY_DBL_MAX,
          res_max = -ARBITRARY_DBL_MAX;
      SMesh::Face_index min_index, max_index;
      for(boost::graph_traits<SMesh>::face_iterator fit = faces(smesh).begin();
          fit != faces(smesh).end();
          ++fit)
      {
        fjacobian[*fit] = scaled_jacobian(*fit, smesh);
        if(fjacobian[*fit] > res_max)
        {
          res_max = fjacobian[*fit];
          max_index = *fit;
        }
        if(fjacobian[*fit] < res_min)
        {
          res_min = fjacobian[*fit];
          min_index = *fit;
        }
      }
      jacobian_min.erase(item);
      jacobian_min.insert(std::make_pair(item, std::make_pair(res_min, min_index)));
      jacobian_max.erase(item);
      jacobian_max.insert(std::make_pair(item, std::make_pair(res_max, max_index)));
      connect(item, &Scene_surface_mesh_item::itemChanged,
              this, &DisplayPropertyPlugin::resetProperty);
    }
    //scale a color ramp between min and max
    double max = maxBox;
    double min = minBox;
    //fill f:color pmap
    SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
        smesh.add_property_map<face_descriptor, CGAL::Color >("f:color", CGAL::Color()).first;
    for(boost::graph_traits<SMesh>::face_iterator fit = faces(smesh).begin();
        fit != faces(smesh).end();
        ++fit)
    {
      if(min == max)
        --min;
      double f = (fjacobian[*fit]-min)/(max-min);
      if(f<min)
        f = min;
      if(f>max)
        f = max;
      CGAL::Color color(
            255*color_ramp.r(f),
            255*color_ramp.g(f),
            255*color_ramp.b(f));
      fcolors[*fit] = color;
    }
    dock_widget->minBox->setValue(jacobian_min[item].first-0.01);
    dock_widget->maxBox->setValue(jacobian_max[item].first);
  }

  void displayAngles(Scene_surface_mesh_item* item)
  {
    SMesh& smesh = *item->face_graph();
    typedef boost::property_map<SMesh, boost::vertex_point_t>::type PMap;
    PMap pmap = get(boost::vertex_point, smesh);
    //compute and store smallest angle per face
    bool non_init;
    SMesh::Property_map<face_descriptor, double> fangle;
    boost::tie(fangle, non_init) = smesh.add_property_map<face_descriptor, double>("f:angle", 0);
    if(non_init)
    {
      double res_min = ARBITRARY_DBL_MAX,
          res_max = -ARBITRARY_DBL_MAX;
      SMesh::Face_index index_min, index_max;
      for(boost::graph_traits<SMesh>::face_iterator fit = faces(smesh).begin();
          fit != faces(smesh).end();
          ++fit)
      {
        bool is_face_triangle = is_triangle(halfedge(*fit, smesh), smesh);
        bool normal_is_ok = true;
        EPICK::Vector_3 normal(0,0,0);

        EPICK::Orientation orientation = CGAL::POSITIVE;
        if(!is_face_triangle)
        {
          face_descriptor f = *fit;
          CGAL::Halfedge_around_face_circulator<SMesh>
              he(halfedge(f, smesh), smesh),
              he_end(he);
          do{
            normal_is_ok = true;

            //Initializes the facet orientation

            EPICK::Point_3 S,T;
            T = get(pmap, source(*he, smesh));
            S = get(pmap, target(*he, smesh));
            EPICK::Vector_3 V1((T-S).x(), (T-S).y(), (T-S).z());
            S = get(pmap,source(next(*he,smesh), smesh));
            T = get(pmap, target(next(*he,smesh), smesh));
            EPICK::Vector_3 V2((T-S).x(), (T-S).y(), (T-S).z());

            if(normal == EPICK::Vector_3(0,0,0))
              normal_is_ok = false;
            {
              normal = CGAL::cross_product(V1, V2);
            }
            if(normal_is_ok)
            {
              orientation = EPICK::Orientation_3()(V1, V2, normal);
              if( orientation == CGAL::COPLANAR )
                normal_is_ok = false;
            }
          }while( ++he != he_end && !normal_is_ok);
        }

        std::vector<float> local_angles;
        local_angles.reserve(degree(*fit, smesh));
        BOOST_FOREACH(halfedge_descriptor hd,
                      halfedges_around_face(halfedge(*fit, smesh),smesh))
        {
          halfedge_descriptor hdn = next(hd, smesh);
          EPICK::Vector_3 v1(get(pmap, source(hd, smesh)), get(pmap, target(hd, smesh))),
              v2(get(pmap, target(hdn, smesh)), get(pmap, source(hdn, smesh)));
          float norm1(CGAL::approximate_sqrt(v1.squared_length())), norm2(CGAL::approximate_sqrt(v2.squared_length()));
          float dot_prod = v1*v2;
          float angle = std::acos(dot_prod/(norm1*norm2));
          if(is_face_triangle || !normal_is_ok)
            local_angles.push_back(angle * 180/CGAL_PI);
          else
          {
            bool is_convex = true;
            EPICK::Orientation res = EPICK::Orientation_3()(v1, v2, normal) ;
            if(res!= orientation && res != CGAL::ZERO)
              is_convex = false;
            local_angles.push_back(is_convex ? angle * 180/CGAL_PI : 360 - angle * 180/CGAL_PI );
          }
        }
        std::sort(local_angles.begin(), local_angles.end());
        fangle[*fit]=local_angles.front();

        if(fangle[*fit] > res_max)
        {
          res_max = fangle[*fit];
          index_max = *fit;
        }
        if(fangle[*fit] < res_min)
        {
          res_min = fangle[*fit];
          index_min = *fit;
        }
      }
      angles_min.erase(item);
      angles_min.insert(std::make_pair(item, std::make_pair(res_min, index_min)));
      angles_max.erase(item);
      angles_max.insert(std::make_pair(item, std::make_pair(res_max, index_max)));

      connect(item, &Scene_surface_mesh_item::itemChanged,
              this, &DisplayPropertyPlugin::resetProperty);
    }
    //scale a color ramp between min and max

    float max = maxBox;
    float min = minBox;

    //fill f:color pmap
    SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
        smesh.add_property_map<face_descriptor, CGAL::Color >("f:color", CGAL::Color()).first;
    for(boost::graph_traits<SMesh>::face_iterator fit = faces(smesh).begin();
        fit != faces(smesh).end();
        ++fit)
    {
      if(min == max)
        --min;
      float f = (fangle[*fit]-min)/(max-min);
      if(f<0)
        f = 0;
      if(f>1)
        f = 1;
      CGAL::Color color(
            255*color_ramp.r(f),
            255*color_ramp.g(f),
            255*color_ramp.b(f));
      fcolors[*fit] = color;
    }
    dock_widget->minBox->setValue(angles_min[item].first);
    dock_widget->maxBox->setValue(angles_max[item].first);
  }

  void replaceRamp()
  {
    color_ramp = Color_ramp(rm, rM, gm, gM, bm, bM);
    displayLegend();
    minBox = dock_widget->minBox->value();
    maxBox = dock_widget->maxBox->value();
  }

  void on_propertyBox_currentIndexChanged(int)
  {
    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
    {
      dock_widget->minBox->setMinimum(0);
      dock_widget->minBox->setMaximum(360);
      dock_widget->minBox->setValue(0);

      dock_widget->maxBox->setMinimum(0);
      dock_widget->maxBox->setMaximum(360);
      Scene_surface_mesh_item* item =
          qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
      if(! item )
        dock_widget->maxBox->setValue(180);
      else if(is_triangle_mesh(*item->face_graph()))
        dock_widget->maxBox->setValue(60);
      else if(is_quad_mesh(*item->face_graph()))
        dock_widget->maxBox->setValue(90);
      break;
    }
    default:
      dock_widget->minBox->setMinimum(-1000);
      dock_widget->minBox->setMaximum(1000);
      dock_widget->minBox->setValue(0);

      dock_widget->maxBox->setMinimum(-1000);
      dock_widget->maxBox->setMaximum(1000);
      dock_widget->maxBox->setValue(2);
      break;
    }
    replaceRamp();
    enableButtons();
  }

  void closure()Q_DECL_OVERRIDE
  {
    dock_widget->hide();
  }

  void on_zoomToMinButton_pressed()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    face_descriptor dummy_fd;
    Point_3 dummy_p;
    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
    {
      ::zoomToId(*item->face_graph(),
                     QString("f%1").arg(angles_min[item].second),
                     qobject_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first()),
                     dummy_fd,
                     dummy_p);
    }
      break;
    case 1:
    {
      ::zoomToId(*item->face_graph(),
                     QString("f%1").arg(jacobian_min[item].second),
                     qobject_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first()),
                     dummy_fd,
                     dummy_p);
    }
      break;
    default:
      break;
    }
  }

  void on_zoomToMaxButton_pressed()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    face_descriptor dummy_fd;
    Point_3 dummy_p;
    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
    {
      ::zoomToId(*item->face_graph(),
                 QString("f%1").arg(angles_max[item].second),
                 qobject_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first()),
                 dummy_fd,
                 dummy_p);
    }
      break;
    case 1:
    {
      ::zoomToId(*item->face_graph(),
                 QString("f%1").arg(jacobian_max[item].second),
                 qobject_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first()),
                 dummy_fd,
                 dummy_p);
    }
      break;
    default:
      break;
    }
  }

private:
  void displayLegend()
  {
    // Create an legend_ and display it
    const int height = 256;
    const int width = 90;
    const int cell_width = width/3;
    const int top_margin = 5;
    const int left_margin = 5;
    const int drawing_height = height - top_margin * 2;
    const int text_height = 20;

    legend_ = QPixmap(width, height + text_height);
    legend_.fill(QColor(200, 200, 200));

    QPainter painter(&legend_);
    painter.setPen(Qt::black);
    painter.setBrush(QColor(200, 200, 200));

    // Build legend_ data
    double min_value(dock_widget->minBox->value()),
        max_value(dock_widget->maxBox->value());
    std::vector<double> graduations(100);
    for(int i=0; i<100; ++i)
      graduations[i] = i/100.0;

    // draw
    int i=0;
    for (std::vector<double>::iterator it = graduations.begin(), end = graduations.end();
         it != end; ++it, i+=2)
    {
      QColor color(255*color_ramp.r(*it),
                   255*color_ramp.g(*it),
                   255*color_ramp.b(*it));
      painter.fillRect(left_margin,
                       drawing_height - top_margin - i,
                       cell_width,
                       2,
                       color);
    }

    // draw right vertical line
    painter.setPen(Qt::blue);

    painter.drawLine(QPoint(left_margin + cell_width+10, drawing_height - top_margin),
                     QPoint(left_margin + cell_width+10,
                            drawing_height - top_margin - static_cast<int>(graduations.size())*2));


    // draw min value and max value
    painter.setPen(Qt::blue);
    QRect min_text_rect(left_margin + cell_width+10,drawing_height - top_margin,
                        50, text_height);
    painter.drawText(min_text_rect, Qt::AlignCenter, tr("%1").arg(min_value, 0, 'f', 1));

    QRect max_text_rect(left_margin + cell_width+10, drawing_height - top_margin - 200,
                        50, text_height);
    painter.drawText(max_text_rect, Qt::AlignCenter, tr("%1").arg(max_value, 0, 'f', 1));

    dock_widget->legendLabel->setPixmap(legend_);
  }
  double scaled_jacobian(const face_descriptor& f , const SMesh &mesh);
  QList<QAction*> _actions;
  Color_ramp color_ramp;
  DockWidget* dock_widget;
  double rm;
  double rM;
  double gm;
  double gM;
  double bm;
  double bM;
  std::map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > jacobian_min;
  std::map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > jacobian_max;

  std::map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > angles_min;
  std::map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > angles_max;

  double minBox;
  double maxBox;
  QPixmap legend_;
};

  /// Code based on the verdict module of vtk

  /*=========================================================================
  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

  double DisplayPropertyPlugin::scaled_jacobian( const face_descriptor& f , const SMesh& mesh)
  {
    boost::property_map<SMesh, boost::vertex_point_t>::type
        pmap = get(boost::vertex_point, mesh);
    std::vector<double> corner_areas(degree(f, mesh));
    std::vector<EPICK::Vector_3> edges;
    BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
    {
      edges.push_back(EPICK::Vector_3(get(pmap, source(hd, mesh)), get(pmap, target(hd, mesh))));
    }
    std::vector<EPICK::Vector_3> corner_normals;
    for(std::size_t i = 0; i < edges.size(); ++i)
    {
      corner_normals.push_back(CGAL::cross_product(edges[i], edges[(i+1)%(edges.size())]));
    }


    EPICK::Vector_3 unit_center_normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, mesh);
    unit_center_normal *= 1.0/CGAL::approximate_sqrt(unit_center_normal.squared_length());

    for(std::size_t i = 0; i < corner_areas.size(); ++i)
    {
      corner_areas[i] =  unit_center_normal*corner_normals[i];
    }
    std::vector<double> length;
    for(std::size_t i=0; i<edges.size(); ++i)
    {
      length.push_back(CGAL::approximate_sqrt(edges[i].squared_length()));
      if( length[i] < ARBITRARY_DBL_MIN)
        return 0.0;
    }
    double min_scaled_jac = ARBITRARY_DBL_MAX;
    for(std::size_t i=0; i<edges.size(); ++i)
    {
      double scaled_jac = corner_areas[i] / (length[i] * length[(i+edges.size()-1)%(edges.size())]);
      min_scaled_jac = (std::min)( scaled_jac, min_scaled_jac );
    }

    if( min_scaled_jac > 0 )
      return (double) (std::min)( min_scaled_jac, ARBITRARY_DBL_MAX );
    return (double) (std::max)( min_scaled_jac, -ARBITRARY_DBL_MAX );

  }


#include "Display_property_plugin.moc"
