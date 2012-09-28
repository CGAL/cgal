#include "Scene_combinatorial_map_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_interface.h"

#include <QObject>
#include <QMenu>
#include <QAction>
#include <QtDebug>
 #include <QKeyEvent>
#include <CGAL/corefinement_operations.h>

Scene_combinatorial_map_item::Scene_combinatorial_map_item(Scene_interface* scene,void* address):last_known_scene(scene),volume_to_display(0),exportSelectedVolume(NULL),address_of_A(address){m_combinatorial_map=NULL;}
Scene_combinatorial_map_item::~Scene_combinatorial_map_item(){if (m_combinatorial_map!=NULL) delete m_combinatorial_map;}

Scene_combinatorial_map_item* Scene_combinatorial_map_item::clone() const{return NULL;}

Kernel::Vector_3 Scene_combinatorial_map_item::compute_face_normal(Combinatorial_map_3::Dart_const_handle adart) const
{
  typedef Combinatorial_map_3::Dart_of_orbit_const_range<1> Dart_in_facet_range;
  typedef Kernel::Vector_3 Vector_3;
  Vector_3 normal = CGAL::NULL_VECTOR;
  
  Dart_in_facet_range vertices=combinatorial_map().darts_of_orbit<1>(adart);
  Kernel::Point_3 points[3];
  int index=0;
  Dart_in_facet_range::const_iterator pit=vertices.begin();
  for (;pit!=vertices.end() && index!=3;++pit,++index ){
    points[index]=pit->attribute<0>()->point();
  }
  
  if (index!=3) return normal;
  
  do{
    Vector_3 n = CGAL::cross_product(points[2]-points[1],points[0]-points[1]);
    if (n != Vector_3(0,0,0) )
      normal = normal + (n / std::sqrt(n*n));
    points[0]=points[1];
    points[1]=points[2];
    if ( pit==vertices.end() ) break;
    
    points[2]=pit->attribute<0>()->point();
    ++pit;
  }while(true);
  
  return normal == Vector_3(0,0,0)? normal : normal / std::sqrt(normal * normal);
}  

void Scene_combinatorial_map_item::set_next_volume(){
  ++volume_to_display;
  volume_to_display=volume_to_display%(combinatorial_map().attributes<3>().size()+1);
  emit itemChanged();

  if (exportSelectedVolume!=NULL && ( volume_to_display==1 || volume_to_display==0 ) )
    exportSelectedVolume->setEnabled(!exportSelectedVolume->isEnabled());
}


template <class Predicate> 
void Scene_combinatorial_map_item::export_as_polyhedron(Predicate pred,const QString& name) const {
  typedef Combinatorial_map_3::Dart_const_handle Dart_handle;
  typedef Combinatorial_map_3::One_dart_per_cell_const_range<3> One_dart_per_vol_range;
  typedef CGAL::internal::Import_volume_as_polyhedron<Polyhedron::HalfedgeDS> Volume_import_modifier;
  
  std::vector<Dart_handle> darts;
  One_dart_per_vol_range cell_range=combinatorial_map().template one_dart_per_cell<3>();

  
  for (One_dart_per_vol_range::const_iterator it = cell_range.begin();it!= cell_range.end() ; ++it )
    if ( pred(it) ){
      darts.push_back(it);
      if (Predicate::only_one_run) break;
    }

  if (!darts.empty())
  {
    Volume_import_modifier modifier=Predicate::swap_orientation?
      Volume_import_modifier(combinatorial_map(),darts.begin(),darts.end(),Predicate::swap_orientation):
      Volume_import_modifier(combinatorial_map(),darts.begin(),darts.end());

    Polyhedron* new_poly=new Polyhedron();
    new_poly->delegate(modifier);
    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(new_poly);
    new_item->setName(name);
    last_known_scene->addItem(new_item);  
  }  
}

struct Select_volume{
  static const bool only_one_run=true;
  static const bool swap_orientation=false;
  Select_volume(std::size_t i):volume_to_select(i),index(0){}
  template <class Dart_handle>
  bool operator() (Dart_handle){
    return ++index==volume_to_select;
  }
private:
  std::size_t volume_to_select;
  std::size_t index;
};

void Scene_combinatorial_map_item::export_current_volume_as_polyhedron() const {
  if (volume_to_display==0) return; //no volume selected
  
  Select_volume predicate(volume_to_display);
  export_as_polyhedron(predicate,QString("%1_%2").arg(this->name()).arg(volume_to_display-1));
}

struct Select_union{
  static const bool only_one_run=false;
  static const bool swap_orientation=true;
  template <class Dart_handle>
  bool operator() (Dart_handle d){ return d->template attribute<3>()->info().outside.size()==2; }
};

struct Select_inter{
  static const bool only_one_run=false;
  static const bool swap_orientation=false;
  template <class Dart_handle>
  bool operator() (Dart_handle d){ return d->template attribute<3>()->info().inside.size()==2; }
};

struct Select_A_minus_B{
  static const bool only_one_run=false;
  static const bool swap_orientation=false;
  Select_A_minus_B(void* address):address_of_A(address){}
  template <class Dart_handle>
  bool operator() (Dart_handle d){ 
    return d->template attribute<3>()->info().inside.size()==1 &&
           static_cast<void*>(*d->template attribute<3>()->info().inside.begin())==address_of_A;  
  }
private:
  void* address_of_A;
};

struct Select_B_minus_A{
  static const bool only_one_run=false;
  static const bool swap_orientation=false;
  Select_B_minus_A(void* address):address_of_A(address){}
  template <class Dart_handle>
  bool operator() (Dart_handle d){ 
    return d->template attribute<3>()->info().inside.size()==1 &&
           static_cast<void*>(*d->template attribute<3>()->info().inside.begin())!=address_of_A;  
  }
private:
  void* address_of_A;
};

void Scene_combinatorial_map_item::export_union_as_polyhedron() const {
  export_as_polyhedron(Select_union(),QString("%1_union_%2").arg("A").arg("B"));
}
void Scene_combinatorial_map_item::export_intersection_as_polyhedron() const{
  export_as_polyhedron(Select_inter(),QString("%1_inter_%2").arg("A").arg("B"));  
}
void Scene_combinatorial_map_item::export_A_minus_B_as_polyhedron() const{
  Select_A_minus_B predicate(address_of_A);
  export_as_polyhedron(predicate,QString("%1_minus_%2").arg("A").arg("B"));
}
void Scene_combinatorial_map_item::export_B_minus_A_as_polyhedron() const{
  Select_B_minus_A predicate(address_of_A);
  export_as_polyhedron(predicate,QString("%1_minus_%2").arg("B").arg("A"));
}

QMenu* Scene_combinatorial_map_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_combinatorial_map_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // http://doc.trolltech.com/lastest/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {
    QAction* actionSelectNextVolume = 
      menu->addAction(tr("Iterate over volumes"));
    actionSelectNextVolume->setObjectName("actionSelectNextVolume");
    connect(actionSelectNextVolume, SIGNAL(triggered()),this, SLOT(set_next_volume()));

    exportSelectedVolume = 
      menu->addAction(tr("Export current volume as polyhedron"));
    exportSelectedVolume->setObjectName("exportSelectedVolume");
    connect(exportSelectedVolume, SIGNAL(triggered()),this, SLOT(export_current_volume_as_polyhedron()));
    exportSelectedVolume->setEnabled(volume_to_display!=0);
    menu->setProperty(prop_name, true);
    
    if(is_from_corefinement()){
      //Export union as polyhedron
      QAction* exportUnion = 
        menu->addAction(tr("Export union as polyhedron"));
      exportUnion->setObjectName("exportUnion");
      connect(exportUnion, SIGNAL(triggered()),this, SLOT(export_union_as_polyhedron()));

      //Export intersection as polyhedron
      QAction* exportIntersection = 
        menu->addAction(tr("Export intersection as polyhedron"));
      exportIntersection->setObjectName("exportIntersection");
      connect(exportIntersection, SIGNAL(triggered()),this, SLOT(export_intersection_as_polyhedron()));

      //Export A minus B as polyhedron
      QAction* exportAMinusB = 
        menu->addAction(tr("Export A minus B as polyhedron"));
      exportAMinusB->setObjectName("exportAMinusB");
      connect(exportAMinusB, SIGNAL(triggered()),this, SLOT(export_A_minus_B_as_polyhedron()));

      //Export B minus A as polyhedron
      QAction* exportBMinusA = 
        menu->addAction(tr("Export B minus A as polyhedron"));
      exportBMinusA->setObjectName("exportBMinusA");
      connect(exportBMinusA, SIGNAL(triggered()),this, SLOT(export_B_minus_A_as_polyhedron()));
      
    }
  }
  return menu;
}

bool Scene_combinatorial_map_item::keyPressEvent(QKeyEvent* e){
  if (e->key()==Qt::Key_N){
    set_next_volume();
    return true;
  }
  return false;
}

void Scene_combinatorial_map_item::direct_draw() const {
  typedef Combinatorial_map_3::One_dart_per_cell_const_range<3> Volume_dart_range;
  typedef Combinatorial_map_3::One_dart_per_incident_cell_const_range<2,3> Facet_in_volume_drange;
  typedef Combinatorial_map_3::Dart_of_orbit_const_range<1> Dart_in_facet_range;
  Volume_dart_range dart_per_volume_range = combinatorial_map().one_dart_per_cell<3>();
  
  std::size_t index = 0;
  for (Volume_dart_range::const_iterator vit=dart_per_volume_range.begin();vit!=dart_per_volume_range.end();++vit)
  {
    if (++index!=volume_to_display && volume_to_display!=0) continue;
    Facet_in_volume_drange facet_range=combinatorial_map().one_dart_per_incident_cell<2,3>(vit);
    
    for(Facet_in_volume_drange::const_iterator fit=facet_range.begin();fit!=facet_range.end();++fit){
      Dart_in_facet_range vertices=combinatorial_map().darts_of_orbit<1>(fit);
      Kernel::Vector_3 normal = compute_face_normal(fit);
      
      ::glBegin(GL_POLYGON);  
      ::glNormal3d(normal.x(),normal.y(),normal.z());
    
      for (Dart_in_facet_range::const_iterator pit=vertices.begin();pit!=vertices.end();++pit ){
        const Kernel::Point_3& p= pit->attribute<0>()->point();
        ::glVertex3d(p.x(),p.y(),p.z());
      }      
      ::glEnd(); 
    }
  }
}

void Scene_combinatorial_map_item::direct_draw_edges() const {
  typedef Combinatorial_map_3::One_dart_per_cell_const_range<1> Edge_darts;
  Edge_darts darts=combinatorial_map().one_dart_per_cell<1>();
  ::glBegin(GL_LINES);
  for (Edge_darts::const_iterator dit=darts.begin();dit!=darts.end();++dit){
    CGAL_assertion(!dit->is_free(1));
    const Kernel::Point_3& a = dit->attribute<0>()->point();
    const Kernel::Point_3& b = dit->beta(1)->attribute<0>()->point();
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());  
  }
  ::glEnd();
}

void Scene_combinatorial_map_item::draw_points() const{
  typedef Combinatorial_map_3::Attribute_const_range<0>::type Point_range;
  const Point_range& points=combinatorial_map().attributes<0>();
  ::glBegin(GL_POINTS);
  for(Point_range::const_iterator pit=boost::next(points.begin());pit!=points.end();++pit){
    const Kernel::Point_3& p=pit->point();
    ::glVertex3d(p.x(),p.y(),p.z());
  }
  ::glEnd();
}

bool Scene_combinatorial_map_item::isEmpty() const {return combinatorial_map().number_of_darts()==0;}

Scene_combinatorial_map_item::Bbox 
Scene_combinatorial_map_item::bbox() const {
  typedef Combinatorial_map_3::Attribute_const_range<0>::type Point_range;
  const Point_range& points=combinatorial_map().attributes<0>();
  CGAL::Bbox_3 bbox=points.begin()->point().bbox();
  for(Point_range::const_iterator pit=boost::next(points.begin());pit!=points.end();++pit)
    bbox=bbox+pit->point().bbox();
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}


QString Scene_combinatorial_map_item::toolTip() const{ 
  if(!m_combinatorial_map)
    return QString();

  std::vector<unsigned int> cells(5);
  for (unsigned int i=0; i<=4; ++i)
    cells[i]=i;
  std::vector<unsigned int> res = combinatorial_map().count_cells(cells);
  if (volume_to_display==0)
    return QObject::tr("<p>Combinatorial_map_3 <b>%1</b> (mode: %8, color: %9)</p>"
                       "<p>Number of darts: %2<br />"
                       "Number of vertices: %3<br />"
                       "Number of edges: %4<br />"
                       "Number of facets: %5<br />"
                       "Number of volumes: %6<br />"
                       "Number of connected components: %7</p>")
      .arg(this->name())
      .arg(combinatorial_map().number_of_darts())
      .arg(res[0])
      .arg(res[1])
      .arg(res[2])
      .arg(res[3])
      .arg(res[4])
      .arg(this->renderingModeName())
      .arg(this->color().name());
  return QObject::tr("<p>Combinatorial_map_3 <b>%1</b> (mode: %8, color: %9)</p>"
                     "<p>Number of darts: %2<br />"
                     "Number of vertices: %3<br />"
                     "Number of edges: %4<br />"
                     "Number of facets: %5<br />"
                     "Number of volumes: %6<br />"
                     "Number of connected components: %7 <br />"
                     "Currently Displaying facets of volume: %10 </p>")
    .arg(this->name())
    .arg(combinatorial_map().number_of_darts())
    .arg(res[0])
    .arg(res[1])
    .arg(res[2])
    .arg(res[3])
    .arg(res[4])
    .arg(this->renderingModeName())
    .arg(this->color().name())
    .arg(volume_to_display-1);  
}


#include "Scene_combinatorial_map_item.moc"
