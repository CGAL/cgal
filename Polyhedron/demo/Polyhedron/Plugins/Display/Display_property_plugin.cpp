#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QColorDialog>
#include <QPalette>
#include <QColor>
#include <QStyleFactory>
#include <QMessageBox>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include "Scene_points_with_normal_item.h"

#include "Messages_interface.h"
#include "Scene_surface_mesh_item.h"
#include "Color_ramp.h"
#include <boost/unordered_map.hpp>
#include "ui_Display_property.h"
#include "id_printing.h"
#include "Scene.h"
#include "triangulate_primitive.h"
#include <CGAL/Buffer_for_vao.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Dynamic_property_map.h>

#define ARBITRARY_DBL_MIN 1.0E-30
#define ARBITRARY_DBL_MAX 1.0E+30


//Item for heat values
typedef CGAL::Three::Triangle_container Tri;
typedef CGAL::Three::Viewer_interface VI;

class Scene_heat_item
    : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
  
public: 
  Scene_heat_item(Scene_surface_mesh_item* item)
    :sm(item->face_graph()), parent(item)
  {
    setTriangleContainer(0, new Triangle_container(VI::PROGRAM_HEAT_INTENSITY,
                                                   true));
    setRenderingMode(Gouraud);
  }
  Scene_item* clone() const{return nullptr;}
  QString toolTip() const{return QString(); }\
  void select(double orig_x,
             double orig_y,
             double orig_z,
             double dir_x,
             double dir_y,
             double dir_z)
  {
    parent->select( orig_x, orig_y, orig_z, 
                    dir_x, dir_y, dir_z);
  }
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const
  {
    getTriangleContainer(0)->initializeBuffers(viewer); 
    getTriangleContainer(0)->setIdxSize(idx.size());
    verts.resize(0);
    normals .resize(0);
    colors.resize(0);
    idx.clear();
    idx.shrink_to_fit();
    colors.shrink_to_fit();
    verts.shrink_to_fit();
    normals.shrink_to_fit();
    
    are_buffers_filled = true;
  }
  
  void draw(CGAL::Three::Viewer_interface *viewer) const
  {
    
    if(!isInit() && viewer->context()->isValid())
      initGL();
    if(!are_buffers_filled)
      initialize_buffers(viewer);
    getTriangleContainer(0)->setAlpha(1.0f);
    getTriangleContainer(0)->draw( viewer, false);
  }
  void compute_bbox() const
  {
    SMesh::Property_map<vertex_descriptor, Point_3> pprop = sm->points();
    CGAL::Bbox_3 bbox ;
    
    BOOST_FOREACH(vertex_descriptor vd,vertices(*sm))
    {
      bbox = bbox + pprop[vd].bbox();
    }
    _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                 bbox.xmax(),bbox.ymax(),bbox.zmax());
    is_bbox_computed = true;
  }
  Scene_item::Bbox bbox() const {
    if(!is_bbox_computed)
      compute_bbox();
    is_bbox_computed = true;
    return _bbox;
  }

  ~Scene_heat_item(){}
  virtual bool supportsRenderingMode(RenderingMode m) const { return m==Gouraud; }
  virtual void invalidateOpenGLBuffers()
  {
    computeElements();
    is_bbox_computed = false;
    are_buffers_filled = false;
  }
  void triangulate_convex_facet(face_descriptor fd,
                                boost::property_map< SMesh, boost::vertex_index_t >::type *im) const
  {
    const CGAL::qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    EPICK::Vector_3 offset = EPICK::Vector_3(v_offset.x, v_offset.y, v_offset.z);
    
    EPICK::Point_3 p0,p1,p2;
    SMesh::Halfedge_around_face_circulator he(halfedge(fd, *sm), *sm);
    SMesh::Halfedge_around_face_circulator he_end = he;
    
    while(next(*he, *sm) != prev(*he_end, *sm))
    {
      ++he;
      vertex_descriptor v0(target(*he_end, *sm)),
          v1(target(*he, *sm)),
          v2(target(next(*he, *sm), *sm));
      p0 = sm->point(v0) + offset;
      p1 = sm->point(v1) + offset;
      p2 = sm->point(v2) + offset;
      idx.push_back((*im)[v0]);
      idx.push_back((*im)[v1]);
      idx.push_back((*im)[v2]);
    }
  }
  void triangulate_facet(face_descriptor fd,
                         SMesh::Property_map<face_descriptor, EPICK::Vector_3> *fnormals,
                         boost::property_map< SMesh, boost::vertex_index_t >::type *im) const
  {
    //Computes the normal of the facet
    EPICK::Vector_3 normal = get(*fnormals, fd);
    
    //check if normal contains NaN values
    if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
    {
      qDebug()<<"Warning : normal is not valid. Facet not displayed";
      return;
    }
    
    typedef FacetTriangulator<SMesh, EPICK, boost::graph_traits<SMesh>::vertex_descriptor> FT;
    const CGAL::qglviewer::Vec off = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    EPICK::Vector_3 offset(off.x,off.y,off.z);
    FT triangulation(fd,normal,sm, offset);
    //iterates on the internal faces
    for(FT::CDT::Finite_faces_iterator
        ffit = triangulation.cdt->finite_faces_begin(),
        end = triangulation.cdt->finite_faces_end();
        ffit != end; ++ffit)
    {
      if(ffit->info().is_external)
        continue;
      //add the vertices to the positions
      //adds the vertices, normals and colors to the appropriate vectors
      //adds the indices to the appropriate vector
      idx.push_back((*im)[triangulation.v2v[ffit->vertex(0)]]);
      idx.push_back((*im)[triangulation.v2v[ffit->vertex(1)]]);
      idx.push_back((*im)[triangulation.v2v[ffit->vertex(2)]]);
    }
  }
  
  void computeElements() const
  {
    typedef EPICK::Point_3 Point;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const CGAL::qglviewer::Vec o = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    EPICK::Vector_3 offset(o.x, o.y, o.z);
    SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
        sm->points();
    SMesh::Property_map<vertex_descriptor, EPICK::Vector_3 > vnormals =
        sm->property_map<vertex_descriptor, EPICK::Vector_3 >("v:normal").first;
    SMesh::Property_map<face_descriptor, EPICK::Vector_3 > fnormals =
        sm->property_map<face_descriptor, EPICK::Vector_3 >("f:normal").first;
    typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
    typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
    SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
        sm->property_map<vertex_descriptor, CGAL::Color >("v:color").first;
    SMesh::Property_map<vertex_descriptor, float> vdist=
        sm->property_map<vertex_descriptor, float >("v:dist").first;    
    typedef CGAL::Buffer_for_vao<float, unsigned int> CPF;
    verts.clear();
    normals.clear();
    idx.clear();
    colors.clear();
    boost::property_map< SMesh, boost::vertex_index_t >::type
        im = get(boost::vertex_index, *sm);
    
    idx.reserve(num_faces(*sm) * 3);
    BOOST_FOREACH(face_descriptor fd, faces(*sm))
    {
      if(is_triangle(halfedge(fd,*sm),*sm))
      {
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *sm),*sm))
        {
          idx.push_back(source(hd, *sm));
        }
      }
      else
      {
        std::vector<Point> facet_points;
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, *sm),*sm))
        {
          facet_points.push_back(positions[target(hd, *sm)]);
        }
        bool is_convex = CPF::is_facet_convex(facet_points, fnormals[fd]);
        
        if(is_convex && is_quad(halfedge(fd,*sm),*sm) )
        {
          halfedge_descriptor hd = halfedge(fd,*sm);
          //1st half
          idx.push_back(source(hd, *sm));
          idx.push_back(source(next(hd, *sm), *sm));
          idx.push_back(source(next(next(hd, *sm), *sm), *sm));
          
          //2nd half
          idx.push_back(source(hd, *sm));
          idx.push_back(source(next(next(hd, *sm), *sm), *sm));
          idx.push_back(source(prev(hd, *sm), *sm));
        }    
        else if(is_convex)
        {
          triangulate_convex_facet(fd, &im);
        }
        else
        {
          triangulate_facet(fd, &fnormals, &im);
        }
      }
    }
    BOOST_FOREACH(vertex_descriptor vd, vertices(*sm))
    {
      CGAL::Color c = vcolors[vd];
      colors.push_back((float)c.red()/255);
      colors.push_back((float)c.green()/255);
      colors.push_back((float)c.blue()/255);
      
      
      Point p = positions[vd] + offset;
      CPF::add_point_in_buffer(p, verts);
      EPICK::Vector_3 n = vnormals[vd];
      CPF::add_normal_in_buffer(n, normals);
      heat_values.push_back(vdist[vd]);
    }
    
    getTriangleContainer(0)->allocate(Tri::Vertex_indices, idx.data(),
                                      static_cast<int>(idx.size()*sizeof(unsigned int)));
    getTriangleContainer(0)->allocate(Tri::Smooth_vertices, verts.data(),
                                      static_cast<int>(num_vertices(*sm)*3*sizeof(float)));
    
    getTriangleContainer(0)->allocate(Tri::Smooth_normals, normals.data(),
                                      static_cast<int>(num_vertices(*sm)*3*sizeof(float)));
    getTriangleContainer(0)->allocate(Tri::VColors, colors.data(),
                                      static_cast<int>(colors.size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Tri::Distances, heat_values.data(),
                                      static_cast<int>(heat_values.size()*sizeof(float)));
    compute_bbox();
     QApplication::restoreOverrideCursor();
  }
  
  bool isEmpty() const {return false;}
  SMesh *face_graph() { return sm;}
  Scene_surface_mesh_item* getParent() { return parent; }

private:
  SMesh* sm;
  Scene_surface_mesh_item* parent;
  mutable std::vector<float> normals;
  mutable std::vector<unsigned int> idx;
  mutable std::vector<float> verts;
  mutable std::vector<float> colors;
  mutable std::vector<float> heat_values;
}; // end class Scene_heat_item

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
  typedef SMesh::Property_map<boost::graph_traits<SMesh>::vertex_descriptor, double> Vertex_distance_map;
  typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<SMesh> Heat_method;
  typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<SMesh, CGAL::Heat_method_3::Intrinsic_Delaunay> Heat_method_idt;
  typedef CGAL::dynamic_vertex_property_t<bool>                        Vertex_source_tag;
  typedef boost::property_map<SMesh, Vertex_source_tag>::type Vertex_source_map;
  
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
    this->current_item = NULL;

    QAction *actionDisplayAngles= new QAction(QString("Display Properties"), mw);
    QAction *actionHeatMethod= new QAction(QString("Heat Method"), mw);
    actionHeatMethod->setProperty("submenuName", "Color");

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
      if(actionHeatMethod)
      {
        connect(actionHeatMethod, &QAction::triggered,
                this, [this](){
          this->dock_widget->propertyBox->setCurrentIndex(2);
          this->dock_widget->show();
        });
      }
      _actions << actionDisplayAngles;
      _actions << actionHeatMethod;

    }
    dock_widget = new DockWidget("Property Displaying", mw);
    dock_widget->setVisible(false);
    addDockWidget(dock_widget);
    QPalette palette(Qt::red);
    dock_widget->minColorButton->setPalette(palette);
    dock_widget->minColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->minColorButton->update();

    palette = QPalette(Qt::green);
    dock_widget->maxColorButton->setPalette(palette);
    dock_widget->maxColorButton->setStyle(QStyleFactory::create("Fusion"));
    dock_widget->maxColorButton->update();
    connect(dock_widget->colorizeButton, SIGNAL(clicked(bool)),
            this, SLOT(colorize()));

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
      if (!minColor.isValid())
      {
        return;
      }
      
      rm = minColor.redF();
      gm = minColor.greenF();
      bm = minColor.blueF();
      QPalette palette(minColor);
      dock_widget->minColorButton->setPalette(palette);
      dock_widget->minColorButton->update();
      replaceRamp();
    });
    connect(dock_widget->maxColorButton, &QPushButton::pressed,
            this, [this]()
    {
      QColor maxColor = QColorDialog::getColor();
      if(!maxColor.isValid())
        return;
      QPalette palette(maxColor);
      rM = maxColor.redF();
      gM = maxColor.greenF();
      bM = maxColor.blueF();

      dock_widget->maxColorButton->setPalette(palette);
      dock_widget->maxColorButton->update();
      replaceRamp();
    });

    connect(dock_widget->sourcePointsButton, SIGNAL(toggled(bool)),
            this, SLOT(on_sourcePointsButton_toggled(bool)));

    connect(dock_widget->resetButton, &QPushButton::pressed,
            this, &DisplayPropertyPlugin::resetRampExtremas);

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
    else{
      replaceRamp(); 
      dock_widget->show();
      dock_widget->raise(); }
  }

  void resetRampExtremas()
  {
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    item->face_graph()->collect_garbage();
    bool ok;
    switch(dock_widget->propertyBox->currentIndex())
    {
    case 0:
      ok = resetAngles(item);
      break;
    default:
      ok = resetScaledJacobian(item);
      break;
    }
    QApplication::restoreOverrideCursor();
    if(!ok)
      QMessageBox::warning(mw, "Error", "You must first run colorize once to initialize the values.");
  }
  
  void colorize()
  {
    Scene_heat_item* h_item = nullptr;
    Scene_surface_mesh_item* item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
    {
      h_item = qobject_cast<Scene_heat_item*>(scene->item(scene->mainSelectionIndex()));
      if(!h_item)
        return;
      item = h_item->getParent();
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);

    replaceRamp();
    item->face_graph()->collect_garbage();

    switch(dock_widget->propertyBox->currentIndex()){
    case 0:
      displayAngles(item);
      break;
      case 1:
        displayScaledJacobian(item);
        break;
    case 2:
      if(!displayHeatIntensity(item))
        return;
      item->setRenderingMode(Gouraud);
      break;
    default:  // Heat Method (Intrinsic Delaunay)
      if(!displayHeatIntensity(item, true))
        return;
      item->setRenderingMode(Gouraud);
      break;
    }

    connect(item, &Scene_surface_mesh_item::itemChanged,
            this, [item](){
      bool does_exist;
      SMesh::Property_map<face_descriptor, double> pmap;
      boost::tie(pmap, does_exist) = 
          item->face_graph()->property_map<face_descriptor,double>("f:jacobian");
      if(does_exist)
        item->face_graph()->remove_property_map(pmap);
      boost::tie(pmap, does_exist) = 
          item->face_graph()->property_map<face_descriptor,double>("f:angle");
      if(does_exist)
        item->face_graph()->remove_property_map(pmap);
    });
    QApplication::restoreOverrideCursor();
    item->invalidateOpenGLBuffers();
    item->redraw();
    if(dock_widget->propertyBox->currentIndex() != 2){
      dock_widget->zoomToMinButton->setEnabled(true);
      dock_widget->zoomToMaxButton->setEnabled(true);}
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
  }

  bool resetScaledJacobian(Scene_surface_mesh_item* item)
  {
    SMesh& smesh = *item->face_graph();
    if(!smesh.property_map<face_descriptor, double>("f:jacobian").second)
    {
      return false;
    }
    dock_widget->minBox->setValue(jacobian_min[item].first-0.01);
    dock_widget->maxBox->setValue(jacobian_max[item].first);
    return true;
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
  }

  bool resetAngles(Scene_surface_mesh_item* item)
  {
    SMesh& smesh = *item->face_graph();
    if(!smesh.property_map<face_descriptor, double>("f:angle").second)
    {
      return false;
    }
    dock_widget->minBox->setValue(angles_min[item].first);
    dock_widget->maxBox->setValue(angles_max[item].first);
    return true;
  }

  // AF: This function gets called when we click on the button "Colorize"
  bool displayHeatIntensity(Scene_surface_mesh_item* item, bool iDT = false)
  {
    SMesh& mesh = *item->face_graph();
    bool found = is_source.find(item) != is_source.end();
    if(!found
       || ! source_points
       || source_points->point_set()->is_empty())
    {
      QApplication::restoreOverrideCursor();
      QMessageBox::warning(mw, "Warning","Source vertices are needed for this property.");
      return false;
    }
    if(!is_triangle_mesh(mesh))
    {
      QApplication::restoreOverrideCursor();
      QMessageBox::warning(mw,"Error","The mesh must be triangulated.");
      return false;
    }
    Heat_method * hm = NULL;
    Heat_method_idt * hm_idt = NULL;
    SMesh::Property_map<vertex_descriptor, double> heat_intensity =
      mesh.add_property_map<vertex_descriptor, double>("v:heat_intensity", 0).first;
    if(! iDT){
      if(mesh_heat_method_map.find(item) != mesh_heat_method_map.end()){
        hm = mesh_heat_method_map[item];
      }else {
        hm = new Heat_method(mesh);
        mesh_heat_method_map[item] = hm;
      }
      connect(item, &Scene_surface_mesh_item::aboutToBeDestroyed,
              [this,item](){
                auto it =  mesh_heat_method_map.find(item);
                delete it->second;
                mesh_heat_method_map.erase(it);
              }
              );
    } else {
      if(mesh_heat_method_idt_map.find(item) != mesh_heat_method_idt_map.end()){
        hm_idt = mesh_heat_method_idt_map[item];
      }else {
        hm_idt = new Heat_method_idt(mesh);
        mesh_heat_method_idt_map[item] = hm_idt;
      }
      connect(item, &Scene_surface_mesh_item::aboutToBeDestroyed,
              [this,item](){
                auto it = mesh_heat_method_idt_map.find(item);
                if(it == mesh_heat_method_idt_map.end())
                  return;
                Heat_method_idt *hm_idt = it->second;
                delete hm_idt;
                mesh_heat_method_idt_map.erase(it);
              }
              );
    }

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      if(get(is_source[item], vd)){
        if(iDT){
          hm_idt->add_source(vd);
        } else
          hm->add_source(vd);
      }
      else
      {
        if(iDT){
          hm_idt->remove_source(vd);
        } else
          hm->remove_source(vd);
      }
    }

    if(iDT){
      hm_idt->estimate_geodesic_distances(heat_intensity);
    }else{
      hm->estimate_geodesic_distances(heat_intensity);
    }

    double max = 0;
    double min = (std::numeric_limits<double>::max)();

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      double hi = heat_intensity[vd];
      if(hi < min)
        min = hi;
      if(hi > max)
        max = hi;
    }
    color_ramp = Color_ramp(rm, rM, gm, gM, bm, bM);
    dock_widget->minBox->setValue(min);
    dock_widget->maxBox->setValue(max);

    //}
    SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
        mesh.add_property_map<vertex_descriptor, CGAL::Color >("v:color", CGAL::Color()).first;
    SMesh::Property_map<vertex_descriptor, float> vdist=
        mesh.add_property_map<vertex_descriptor, float >("v:dist", 0.0).first;
    for(boost::graph_traits<SMesh>::vertex_iterator vit = vertices(mesh).begin();
        vit != vertices(mesh).end();
        ++vit)
    {
      double h =(heat_intensity[*vit]-min)/(max-min);
      CGAL::Color color(
            255*color_ramp.r(h),
            255*color_ramp.g(h),
            255*color_ramp.b(h));
      vcolors[*vit] = color;
      vdist[*vit]=h;
    }
    Scene_group_item* group;
    if(mesh_heat_item_map.find(item) != mesh_heat_item_map.end())
    {
      group = mesh_heat_item_map[item]->parentGroup();
      group->unlockChild(mesh_heat_item_map[item]);
      scene->erase(scene->item_id(mesh_heat_item_map[item]));
    }
    else
    {
      group = new Scene_group_item("Heat Visualization");
      scene->addItem(group);
      scene->changeGroup(item, group);
      scene->changeGroup(source_points, group);
      group->lockChild(item);
      group->lockChild(source_points);
    }
    mesh_heat_item_map[item] = new Scene_heat_item(item);
    mesh_heat_item_map[item]->setName(tr("%1 heat").arg(item->name()));
    scene->addItem(mesh_heat_item_map[item]);
    scene->changeGroup(mesh_heat_item_map[item], group);
    group->lockChild(mesh_heat_item_map[item]);
    item->setVisible(false);
    displayLegend();
    if(dock_widget->sourcePointsButton->isChecked())
      dock_widget->sourcePointsButton->toggle();
    return true;
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
      dock_widget->groupBox->  setEnabled(true);
      dock_widget->groupBox_3->setEnabled(true);

      dock_widget->sourcePointsButton->setEnabled(false);

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
    case 1:
      dock_widget->groupBox->  setEnabled(true);
      dock_widget->groupBox_3->setEnabled(true);
      dock_widget->sourcePointsButton->setEnabled(false);

      dock_widget->minBox->setMinimum(-1000);
      dock_widget->minBox->setMaximum(1000);
      dock_widget->minBox->setValue(0);

      dock_widget->maxBox->setMinimum(-1000);
      dock_widget->maxBox->setMaximum(1000);
      dock_widget->maxBox->setValue(2);
      break;
    default:
      dock_widget->maxBox->setMinimum(0);
      dock_widget->maxBox->setMaximum(99999999);
      dock_widget->groupBox->  setEnabled(false);
      dock_widget->groupBox_3->setEnabled(false);
      dock_widget->sourcePointsButton->setEnabled(true);

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
                     qobject_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first()),
                     dummy_fd,
                     dummy_p);
    }
      break;
    case 1:
    {
      ::zoomToId(*item->face_graph(),
                     QString("f%1").arg(jacobian_min[item].second),
                     qobject_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first()),
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
                 qobject_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first()),
                 dummy_fd,
                 dummy_p);
    }
      break;
    case 1:
    {
      ::zoomToId(*item->face_graph(),
                 QString("f%1").arg(jacobian_max[item].second),
                 qobject_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first()),
                 dummy_fd,
                 dummy_p);
    }
      break;
    default:
      break;
    }
  }

  void on_sourcePointsButton_toggled(bool b)
  {
    if(b)
    {
      Scene_heat_item* h_item = nullptr;
      Scene_surface_mesh_item* item =
          qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
      if(!item)
      {
        h_item = qobject_cast<Scene_heat_item*>(scene->item(scene->mainSelectionIndex()));
        if(h_item)
          item = h_item->getParent();
      }
      if(!item)
      {
        QMessageBox::warning(mw, "Warning", "You must select a Surface_mesh_item to make this work. Aborting.");
        dock_widget->sourcePointsButton->setChecked(false);
        return;
      }
      current_item = item;
      connect(current_item, &Scene_surface_mesh_item::aboutToBeDestroyed,
              this, [this]()
      {
        dock_widget->sourcePointsButton->setChecked(false);
      });
      if(mesh_sources_map.find(item) == mesh_sources_map.end())
      {
        source_points = new Scene_points_with_normal_item();
        source_points->setName(QString("Source vertices for %1").arg(current_item->name()));
        source_points->setColor(QColor(Qt::red));
        source_points->setPointSize(5);
        scene->addItem(source_points);
        connect(source_points, &Scene_points_with_normal_item::aboutToBeDestroyed,
                [this](){
          boost::unordered_map<Scene_surface_mesh_item*, Scene_points_with_normal_item*>::iterator it;
          for(it = mesh_sources_map.begin();
              it != mesh_sources_map.end();
              ++it)
          {
            if(it->second == source_points)
            {
              mesh_sources_map.erase(it);
              break;
            }
          }
        });
      mesh_sources_map[current_item] = source_points;
      }
      else
      {
        source_points=mesh_sources_map[current_item];
      }
      connect(item, SIGNAL(selected_vertex(void*)), this, SLOT(on_vertex_selected(void*)));
      bool non_init = is_source.find(item) == is_source.end();
      if(non_init)
      {
        Vertex_source_map map = get(Vertex_source_tag(), *item->face_graph());
        is_source.insert(std::make_pair(item, map));
        connect(item, &Scene_surface_mesh_item::itemChanged,
                this, &DisplayPropertyPlugin::resetProperty);
        connect(item, &Scene_surface_mesh_item::aboutToBeDestroyed,
                [this, item](){
          if(is_source.find(item) != is_source.end())
          {
            is_source.erase(item);
          }
        });
      }
    }
    else
    {
      if(!current_item)
        return;
      disconnect(current_item, SIGNAL(selected_vertex(void*)), this, SLOT(on_vertex_selected(void*)));
      current_item = NULL;
    }
  }

  void on_vertex_selected(void* void_ptr)
  {
    typedef boost::graph_traits<SMesh>::vertices_size_type size_type;
    size_type h = static_cast<size_type>(reinterpret_cast<std::size_t>(void_ptr));
    vertex_descriptor vd = static_cast<vertex_descriptor>(h) ;
    bool found = is_source.find(current_item) != is_source.end();
    if(found)
    {
      if(!get(is_source[current_item], vd))
      {
        put(is_source[current_item], vd, true);
        source_points->point_set()->insert(current_item->face_graph()->point(vd));
      }
      else
      {
        put(is_source[current_item], vd, false);
        Point_set::iterator it;
        for(it = source_points->point_set()->begin(); it != source_points->point_set()->end(); ++it)
          if(source_points->point_set()->point(*it) == current_item->face_graph()->point(vd))
          {
            source_points->point_set()->remove(it);
            source_points->point_set()->collect_garbage();
            break;
          }
      }
    }

   source_points->invalidateOpenGLBuffers();
   source_points->itemChanged();
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
  boost::unordered_map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > jacobian_min;
  boost::unordered_map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > jacobian_max;

  boost::unordered_map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > angles_min;
  boost::unordered_map<Scene_surface_mesh_item*, std::pair<double, SMesh::Face_index> > angles_max;
  boost::unordered_map<Scene_surface_mesh_item*, Vertex_source_map> is_source;


  double minBox;
  double maxBox;
  QPixmap legend_;

  Scene_surface_mesh_item* current_item;
  Scene_points_with_normal_item* source_points;
  boost::unordered_map<Scene_surface_mesh_item*, Scene_points_with_normal_item*> mesh_sources_map;
  boost::unordered_map<Scene_surface_mesh_item*, Scene_heat_item*> mesh_heat_item_map;

  boost::unordered_map<Scene_surface_mesh_item*, Heat_method*> mesh_heat_method_map;
  boost::unordered_map<Scene_surface_mesh_item*, Heat_method_idt*> mesh_heat_method_idt_map;
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
