#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_transform_item.h"
#include "Polyhedron_type.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polylines_item.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

class Polyhedron_demo_transform_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:

  Polyhedron_demo_transform_polyhedron_plugin():started(false){}

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTransformPolyhedron;
  }

  bool applicable() const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_polyhedron_transform_item*>(scene->item(scene->mainSelectionIndex()));
  }
  
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionTransformPolyhedron = new QAction("Affine transformation of polyhedron", mw);
    if(actionTransformPolyhedron) {
      connect(actionTransformPolyhedron, SIGNAL(triggered()),this, SLOT(go()));
    }
  }

  void start(Scene_polyhedron_item*);
  void end();  
  
private:

  QAction*  actionTransformPolyhedron;
  Scene_polyhedron_transform_item* transform_item;
  Scene_interface::Item_id tr_item_index;
  bool started;

public slots:
  void go();
  void transformed_killed();
}; // end class Polyhedron_demo_transform_polyhedron_plugin

void Polyhedron_demo_transform_polyhedron_plugin::go(){
  if (!started){
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item) return;    
    
    started=true;
    actionTransformPolyhedron->setText("Apply affine transformation");
    start(poly_item);
  }
  else
    end();    
}

void Polyhedron_demo_transform_polyhedron_plugin::transformed_killed(){
    started=false;
    actionTransformPolyhedron->setText("Affine transformation of polyhedron");
}

void Polyhedron_demo_transform_polyhedron_plugin::start(Scene_polyhedron_item* poly_item){
  QApplication::setOverrideCursor(Qt::PointingHandCursor);
  
  Scene_polyhedron_item::Bbox bbox = poly_item->bbox();
  double x=(bbox.xmin+bbox.xmax)/2;
  double y=(bbox.ymin+bbox.ymax)/2;
  double z=(bbox.zmin+bbox.zmax)/2;
  
  transform_item = new Scene_polyhedron_transform_item(qglviewer::Vec(x,y,z),poly_item,scene);
  transform_item->setManipulatable(true);
  transform_item->setColor(Qt::green);
  transform_item->setName(tr("Affine transformation of polyhedron"));
  connect(transform_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  tr_item_index=scene->addItem(transform_item);
}


struct Modifier_transform_vertices : public CGAL::Modifier_base<Polyhedron::HalfedgeDS> {
  typedef Polyhedron::HalfedgeDS HDS;
  
  CGAL::Aff_transformation_3<Kernel> transform;
  Kernel::Vector_3 frame_center_translation;
  Modifier_transform_vertices(const GLdouble* m,const qglviewer::Vec& tr):
    transform(m[0],m[4], m[8],m[12],
              m[1],m[5], m[9],m[13],
              m[2],m[6],m[10],m[14],
              /*m[3],m[7],m[11],*/m[15]),
    frame_center_translation(-tr.x,-tr.y,-tr.z)
  {
    CGAL_assertion(m[3]==0);
    CGAL_assertion(m[7]==0);
    CGAL_assertion(m[11]==0);
  }
  
  void operator()(HDS& hds)
  {
    for (HDS::Vertex_iterator it=hds.vertices_begin(),
                                   endit=hds.vertices_end();endit!=it;++it)
    {
      it->point() = transform( it->point() + frame_center_translation );
    }
  }
};


void Polyhedron_demo_transform_polyhedron_plugin::end(){
  QApplication::restoreOverrideCursor();
  const GLdouble* matrix = transform_item->manipulatedFrame()->matrix();
  Modifier_transform_vertices modifier(matrix,transform_item->center());
  Polyhedron* new_poly=new Polyhedron(*transform_item->getBase()->polyhedron());
  new_poly->delegate(modifier);
  
  Scene_polyhedron_item* new_item=new Scene_polyhedron_item(new_poly);
  new_item->setName(tr("%1_transformed").arg(transform_item->getBase()->name()));
  
  scene->replaceItem(tr_item_index,new_item);
  delete transform_item;
}
  

Q_EXPORT_PLUGIN2(Polyhedron_demo_transform_polyhedron_plugin, Polyhedron_demo_transform_polyhedron_plugin)

#include "Polyhedron_demo_transform_polyhedron_plugin.moc"
