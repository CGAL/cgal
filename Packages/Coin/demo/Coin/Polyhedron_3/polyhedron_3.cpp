#include "cgal_types.h"
#include <stdlib.h> // exit()
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoCallback.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoPickStyle.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoTexture2.h>

#include <qapplication.h>
#include <qmainwindow.h>
#include <qfiledialog.h>
#include <qpopupmenu.h>
#include <qmenubar.h>

#include <CGAL/IO/So_node_polyhedron_3.h>

#include <list>
#include <map>
#include <stack>

SoQtExaminerViewer * viewer;
SoSeparator *root;
Polyhedron P, P2;
Point_3 pp1, pp2, pp3;
SbVec3f normal3f;
//SoPrimitiveVertex v1p, v2p, v3p;
bool should_pick = false; //true only once when the key was pressed
bool right_button_found_polygon = false;
Node_polyhedron_3<Polyhedron> *poly;
SoComplexity * complexity_node;
/*
void render_custom(void *, SoAction *){
  if(!should_pick)
    return;
  should_pick = false;
  float params[4];
  glPushMatrix();
  //glTranslatef(0, 0, 0.1);

  glGetFloatv(GL_CURRENT_COLOR, &params[0]);
  //printf("\n%f  %f  %f  %f\n", params[0], params[1], params[2], params[3]);
  glColor3f(0.0, 1.0, 0.0);
  float mat_diffuse[] = {0.0, 1.0, 0.0, 0.0};
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  //glColorMaterial(GL_FRONT, GL_DIFFUSE);
  glBegin(GL_TRIANGLES);    
    glNormal3f(normal3f[0], normal3f[1], normal3f[2]);
    glVertex3f(v1p.getPoint()[0], v1p.getPoint()[1], v1p.getPoint()[2]);
    glVertex3f(v2p.getPoint()[0], v2p.getPoint()[1], v2p.getPoint()[2]);
    glVertex3f(v3p.getPoint()[0], v3p.getPoint()[1], v3p.getPoint()[2]);
  glEnd();
  glColor4fv(&params[0]);
  glPopMatrix();
}
*/

void
mouse_moved(void * ud, SoEventCallback *n)
{  
  const SoEvent *mme = (SoEvent *)n->getEvent();
  if(!right_button_found_polygon)
    return;
  SoQtExaminerViewer * viewer = (SoQtExaminerViewer *)ud;
  SoRayPickAction rp(viewer->getViewportRegion());
  rp.setPoint(mme->getPosition());
  rp.apply(viewer->getSceneManager()->getSceneGraph());
  SoPickedPoint * point = rp.getPickedPoint();
  if (point == NULL) {    
    return;
  }
  const SoDetail* pickDetail = point->getDetail();
  if(pickDetail != NULL && pickDetail->getTypeId() == SoPolyhedronDetail<Polyhedron>::getClassTypeId()){
    SoPolyhedronDetail<Polyhedron> *poly_detail = (SoPolyhedronDetail<Polyhedron> *) pickDetail;
    Facet_handle fh = poly_detail->find_face();
    Halfedge_handle hh = (*fh).halfedge();
    //Halfedge_handle h(&(*hh));
    P.erase_facet(hh);
    viewer->render();
  } 

}

void
mouse_button_pressed(void * ud, SoEventCallback * n)
{
  const SoMouseButtonEvent * mbe = (SoMouseButtonEvent *)n->getEvent();

  if (mbe->getButton() == SoMouseButtonEvent::BUTTON1 &&
      mbe->getState() == SoButtonEvent::DOWN) {
    SoQtExaminerViewer * viewer = (SoQtExaminerViewer *)ud;

    SoRayPickAction rp(viewer->getViewportRegion());
    rp.setPoint(mbe->getPosition());
    //rp.setPickAll(true);
    rp.apply(viewer->getSceneManager()->getSceneGraph());        
    SoPickedPoint * point = rp.getPickedPoint();    
    if (point == NULL) {
      (void)fprintf(stderr, "\n** MISS LEFT BUTTON! **\n\n");
      (void)fprintf(stdout, "\n");
      return;
    }
    
    SbVec3f v = point->getPoint();
    SbVec3f nv = point->getNormal();

    normal3f = nv;

    const SoDetail* pickDetail = point->getDetail();
    if(pickDetail != NULL && pickDetail->getTypeId() == SoPolyhedronDetail<Polyhedron>::getClassTypeId()){
      SoPolyhedronDetail<Polyhedron> *poly_detail = (SoPolyhedronDetail<Polyhedron> *) pickDetail;
      //v1p.setPoint(poly_detail->get_vertex(0)->getPoint());
      //v1p.setNormal(poly_detail->get_vertex(0)->getNormal());
      //v2p.setPoint(poly_detail->get_vertex(1)->getPoint());
      //v2p.setNormal(poly_detail->get_vertex(1)->getNormal());
      //v3p.setPoint(poly_detail->get_vertex(2)->getPoint());
      //v3p.setNormal(poly_detail->get_vertex(2)->getNormal());

      Facet_handle fh = poly_detail->find_face();
      Halfedge_handle hh = (*fh).halfedge();
      P.erase_facet(hh);
      should_pick = true;      
    }    
    viewer->render();
  }
  if (mbe->getButton() == SoMouseButtonEvent::BUTTON2 &&
      mbe->getState() == SoButtonEvent::DOWN) {
    if(right_button_found_polygon)
        right_button_found_polygon = false;
      else
        right_button_found_polygon = true;
    viewer->render();
  }
}


SoSeparator* get_main_scene(){
  Node_polyhedron_3<Polyhedron>::initClass();


  //read the polyhedron
  const char* iname = "cin";
  std::istream*    p_in  = &std::cin;
  std::ifstream    in;
  in.open("data\\venus.off");
  p_in = &in;
  if ( !*p_in)
    std::cout << "error: cannot open file for reading." <<endl;
  //CGAL::set_ascii_mode(* p_in);
  (*p_in) >> P;
  in.close();

  SoSeparator * sep1 = new SoSeparator;  
  complexity_node = new SoComplexity;
  SoMaterial * material = new SoMaterial;

  poly = new Node_polyhedron_3<Polyhedron>(P);

  material->shininess.setValue(1.0f);
  material->ambientColor.setValue(0.0f, 0.0f, 0.0f);
  material->diffuseColor.setValue(0.8f, 0.2f, 0.0);
  complexity_node->value.setValue(0.8f);
  
  sep1->ref();
  sep1->addChild(complexity_node);
  sep1->addChild(material);
  sep1->addChild(poly);

  return sep1;
}

class MyWindow : public QMainWindow{
  Q_OBJECT
public:
  MyWindow(){
    SoQt::init(this);

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Polyhedron", this, SLOT(load_polyhedron()), CTRL+Key_L);
    file->insertItem("&Save Polyhedron", this, SLOT(save_polyhedron()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("&Print to ps", this, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // edit widget menu
    QPopupMenu * edit = new QPopupMenu( this );
    menuBar()->insertItem( "&Edit", edit );
    QPopupMenu * complexity = new QPopupMenu( edit );    
    edit->insertItem("&Change complexity", complexity);
    complexity->insertItem("Bounding box", this, SLOT(change_bbox()));
    complexity->insertItem("Faces", this, SLOT(change_face()));
    complexity->insertItem("Smooth", this, SLOT(change_smooth()));
    
    QPopupMenu * draw = new QPopupMenu(this);
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate Polyhedron", this, SLOT(generate_polyhedron()), CTRL+Key_G);

    SoEventCallback *myEventCB = new SoEventCallback;  
  
    root = get_main_scene();    
    root->addChild(myEventCB);
    // Set up the ExaminerViewer
    viewer = new SoQtExaminerViewer(this);
    viewer->setSceneGraph(root);    
    viewer->viewAll();  
    viewer->show();

    //viewer->setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_BBOX); 

    setCentralWidget(viewer->getBaseWidget());
    // Set up the event callback. We want to pass the root of the
    // entire scene graph (including the camera) as the userData,
    // so we get the scene manager's version of the scene graph
    // root.

    myEventCB->addEventCallback(
        SoMouseButtonEvent::getClassTypeId(),
        mouse_button_pressed,
        viewer);

    myEventCB->addEventCallback(
        SoEvent::getClassTypeId(),
        mouse_moved,
        viewer);

  }
public slots:
  void new_instance(){
  }
  void new_window(){
  }
  void load_polyhedron(){
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "GeomView files (*.off) (*.cgal)", this ) );
    if ( s.isEmpty() )
        return;
    //read the polyhedron
    P.clear();
    const char* iname = "cin";
    std::istream*    p_in  = &std::cin;
    std::ifstream    in;
    in.open(s);
    p_in = &in;
    if ( !*p_in)
      std::cout << "error: cannot open file for reading." <<endl;
    //CGAL::set_ascii_mode(* p_in);
    (*p_in) >> P;
    in.close();
    poly->compute_normals_for_faces();
    poly->compute_normals_for_vertices();
    viewer->viewAll();
  }
  void save_polyhedron(){
    QFileDialog qfd(this, "Save Polyhedron", true);
    qfd.setViewMode(QFileDialog::Detail);    
    qfd.addFilter("CGAL files (*.cgal)");
    qfd.addFilter("GeomView files (*.off)");
    qfd.addFilter("VRML files (*.wrl)");
    qfd.addFilter("Inventor files (*.iv)");
    qfd.addFilter("Geometry files (*.cgal *.off *.wrl *.iv)");
    qfd.setMode(QFileDialog::AnyFile);

    QString fileName;
    if ( qfd.exec() == QDialog::Accepted )
      fileName = qfd.selectedFile();

    if ( !fileName.isNull() ) {
      // got a file name
      if(fileName.endsWith(".cgal")){
        std::ofstream out(fileName);
        CGAL::set_ascii_mode(out);
        out << P;
      }
      else if(fileName.endsWith(".off")){
        std::ofstream out(fileName);
        CGAL::set_ascii_mode(out);
        out << P;
      } else if(fileName.endsWith(".iv")){
        CGAL::File_writer_inventor writter;
      } else if(fileName.endsWith(".wrl")){
        CGAL::File_writer_VRML_2 writter;
      }
      
    }

  }
  void print_to_ps(){
  }
  void generate_polyhedron(){
  }
  void change_bbox(){
    complexity_node->value.setValue(0);
    poly->touch();
  }

  void change_face(){
    complexity_node->value.setValue(0.7f);
    poly->touch();
  }

  void change_smooth(){
    complexity_node->value.setValue(1.0f);
    poly->touch();
  }

private:
};

#include "polyhedron_3.moc"

int
main (int argc, char ** argv)
{
  QApplication app(argc, argv);
  MyWindow *mainwin = new MyWindow();
  app.setMainWidget(mainwin);   
  mainwin->resize(400, 400);
  mainwin->setCaption("Polyhedron_3 Demo");
  mainwin->show();
  app.exec();
  return 0;
}
