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

#include <CGAL/IO/So_node_polyhedron_3.h>

#include <list>
#include <map>
#include <stack>

SoQtExaminerViewer * viewer;
SoSeparator *root;
Polyhedron P, P2;
Point_3 pp1, pp2, pp3;
SbVec3f normal3f;
SoPrimitiveVertex v1p, v2p, v3p;
bool should_pick = false; //true only once when the key was pressed
bool right_button_found_polygon = false;



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
      return;
    }
    

    (void)fprintf(stdout, "\n");

    SbVec3f v = point->getPoint();
    SbVec3f nv = point->getNormal();

    normal3f = nv;

    const SoDetail* pickDetail = point->getDetail();
    if(pickDetail != NULL && pickDetail->getTypeId() == SoPolyhedronDetail<Polyhedron>::getClassTypeId()){
      SoPolyhedronDetail<Polyhedron> *poly_detail = (SoPolyhedronDetail<Polyhedron> *) pickDetail;
      v1p.setPoint(poly_detail->get_vertex(0)->getPoint());
      v1p.setNormal(poly_detail->get_vertex(0)->getNormal());
      v2p.setPoint(poly_detail->get_vertex(1)->getPoint());
      v2p.setNormal(poly_detail->get_vertex(1)->getNormal());
      v3p.setPoint(poly_detail->get_vertex(2)->getPoint());
      v3p.setNormal(poly_detail->get_vertex(2)->getNormal());

      Facet_handle fh = poly_detail->find_face();      
      Halfedge_around_facet_circulator haf = (*fh).facet_begin ();
      do{
        std::cout << (*(*haf).vertex()).point().x() <<"\n"<< (*(*haf).vertex()).point().y() <<"\n"<< (*(*haf).vertex()).point().z() << "\n";
      }while(++haf!=(*fh).facet_begin());
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



int
main (int argc, char ** argv)
{
  // Initialize Coin, and return a main window to use
  // If unsuccessful, exit
  QWidget * window = SoQt::init(argv[0]); 
  if (window==NULL) exit(1);    
  Node_polyhedron_3<Polyhedron>::initClass();

  //read the polyhedron
  const char* iname = "cin";
  std::istream*    p_in  = &std::cin;
  std::ifstream    in;
  in.open("data\\venus.off");
  p_in = &in;  
  if ( !*p_in)
    std::cout << "error: cannot open file for reading." <<endl;  
  CGAL::set_ascii_mode(* p_in);      
  (*p_in) >> P;
  in.close();

  root = new SoSeparator;
  SoSeparator * sep1 = new SoSeparator;
  SoSeparator * sep2 = new SoSeparator;
  SoComplexity * complexity = new SoComplexity;
  SoMaterial * material = new SoMaterial;

  Node_polyhedron_3<Polyhedron> *poly = new Node_polyhedron_3<Polyhedron>(P);
  SoEventCallback *myEventCB = new SoEventCallback;  
  

  material->diffuseColor.setValue(0.8f, 0.2f, 0.0);
  complexity->value.setValue(0.8f);


  root->ref();                   //increments the reference counter
  root->addChild(myEventCB);
  
  sep1->ref();
  sep1->addChild(complexity);
  sep1->addChild(material);
  sep1->addChild(poly);

  root->addChild(sep1);

  // Set up the ExaminerViewer
  viewer = new SoQtExaminerViewer(window);
  viewer->setSceneGraph(root);
  viewer->setTitle("Draw Style");
  viewer->viewAll();  
  viewer->show();

  viewer->setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_BBOX);  

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

  SoQt::show(window); // display the main window
  SoQt::mainLoop();   // main Coin event loop
  delete viewer;      // free all the viewers resources
  root->unref();      // decrements the reference counter  
  return 0;
}
