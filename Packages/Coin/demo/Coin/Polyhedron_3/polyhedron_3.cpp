#include "cgal_types.h"
#include <stdlib.h> // exit()
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoCallback.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoPickStyle.h>
#include <Inventor/nodes/SoColorIndex.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/manips/SoTrackballManip.h>
#include <Inventor/manips/SoTransformBoxManip.h>
#include <Inventor/draggers/SoDragPointDragger.h>
#include <Inventor/draggers/SoTransformBoxDragger.h>

#include "So_node_polyhedron_3.h"
#include "So_node_triangulation_3.h"
#include "So_node_tetrahedron_3.h"
#include "So_node_terrain.h"

#include <list>
#include <map>
#include <stack>

Delaunay dt;
Triangulation	T;
SoCube *cube;
SoQtExaminerViewer * viewer;
SoSeparator *root;
Polyhedron P, P2;
Tetrahedron_3 Tetra;
std::list<Cell_handle>  list_of_cell_handles;
Plane_3 plane(Point_3(0, 0, 0), Point_3(0, 1, 0), Point_3(1, 0, 0));
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

void render_cells(void *, SoAction *){
  std::list<Cell_handle>::const_iterator it = list_of_cell_handles.begin();
  while(it!=list_of_cell_handles.end()){
    Point c0 = (**it).vertex(0)->point();
    Point c1 = (**it).vertex(1)->point();
    Point c2 = (**it).vertex(2)->point();
    Point c3 = (**it).vertex(3)->point();
    
    Point_3 p0(c0.x(), c0.y(), c0.z());
    Point_3 p1(c1.x(), c1.y(), c1.z());
    Point_3 p2(c2.x(), c2.y(), c2.z());
    Point_3 p3(c3.x(), c3.y(), c3.z());

    
    double sqnorm;
    Vector_3 v_n, normal;    

    glPushMatrix();      
      glColor3f(1.0, 0.0, 0.0);
      glBegin(GL_TRIANGLES);
        normal = CGAL::cross_product(p2 - p0, p1 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);
          glNormal3f(v_n.x(), v_n.y(), v_n.z());          
        }        
        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p2.x(),p2.y(),p2.z());

        normal = CGAL::cross_product(p3 - p0, p2 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }

        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p2.x(),p2.y(),p2.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

        normal = CGAL::cross_product(p1 - p0, p3 - p0);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }        
        glVertex3f(p0.x(),p0.y(),p0.z());
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

        normal = CGAL::cross_product(p2 - p1, p3 - p1);
        sqnorm = normal * normal;
        if(sqnorm != 0){
          v_n = normal / std::sqrt(sqnorm);          
          glNormal3f(v_n.x(), v_n.y(), v_n.z());
        }        
        glVertex3f(p1.x(),p1.y(),p1.z());
        glVertex3f(p2.x(),p2.y(),p2.z());
        glVertex3f(p3.x(),p3.y(),p3.z());

      glEnd();
    glPopMatrix();

    it++;
  }
}

SbBool
writePickedPath (SoNode *root, 
   const SbViewportRegion &viewport, 
   const SbVec2s &cursorPosition)
{
   SoRayPickAction myPickAction(viewport);

   // Set an 8-pixel wide region around the pixel
   myPickAction.setPoint(cursorPosition);
   myPickAction.setRadius(8.0);

   // Start a pick traversal
   myPickAction.apply(root);
   const SoPickedPoint *myPickedPoint = 
            myPickAction.getPickedPoint();
   if (myPickedPoint == NULL) return FALSE;

   // Write out the path to the picked object
   SoWriteAction myWriteAction;
   SoPath *path = myPickedPoint->getPath();
   path->write(&myWriteAction);

   return TRUE;
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


void headCB(void *data, SoDragger *dragger)
{
  SoTransformBoxDragger *drg = (SoTransformBoxDragger *)dragger;
  SbVec3f vec_rot;
  float angle;
  drg->rotation.getValue(vec_rot, angle);
  printf("Translation_X: %f\n",drg->translation.getValue()[0]);
  printf("Translation_Y: %f\n",drg->translation.getValue()[1]);
  printf("Translation_Z: %f\n",drg->translation.getValue()[2]);  
  printf("Rotation_angle: %f\n", angle);
  printf("Rotation_X axe: %f\n", vec_rot[0]);
  printf("Rotation_Y axe  %f\n", vec_rot[1]);
  printf("Rotation_Z axe  %f\n", vec_rot[2]);

}

void draggerCB(void *data, SoDragger *dragger)
{

  SoTransformBoxDragger *drg = (SoTransformBoxDragger *)dragger;
  SbVec3f vec_rot;
  float angle;
  drg->rotation.getValue(vec_rot, angle);
    
  Vector_3 v(vec_rot[0], vec_rot[1], vec_rot[2]);
  Vector_3 u = v/sqrt(v*v);
  double x = u.x();
  double y = u.y();
  double z = u.z();
  double cosinus = cos(angle);
  double sinus = sin(angle);

  double m11 = x * x + cosinus * (1 - x * x);
  double m12 = x * y - cosinus * x * y - sinus * z;
  double m13 = x * z - cosinus * x * z + sinus * y;
  double m21 = y * x - cosinus * y * x + sinus * z;
  double m22 = y * y + cosinus * (1 - y * y);
  double m23 = y * z - cosinus * y * z - sinus * x;
  double m31 = z * x - cosinus * z * x - sinus * y;
  double m32 = z * y - cosinus * z * y + sinus * x;
  double m33 = z * z + cosinus * (1 - z * z);
  
  
  Aff t(m11, m12, m13, drg->translation.getValue()[0],
        m21, m22, m23, drg->translation.getValue()[1],
        m31, m32, m33, drg->translation.getValue()[2]);
          
  Plane_3 transf_plane = t.transform(plane);
  pp1 = transf_plane.point();
  pp2 = pp1 + 5*transf_plane.base1();
  pp3 = pp1 + 5*transf_plane.base2();



  Vertex_handle vh = T.infinite_vertex();

  std::list<Cell_handle>  list_of_cells;
  typedef std::list<Cell_handle> LIST;
  std::back_insert_iterator<LIST> backiter(list_of_cells);
  T.incident_cells(vh, backiter);
  
  bool found = false;
  Cell_handle result_handle;
  std::list<Cell_handle>::const_iterator it = list_of_cells.begin();  
  while(!found && it!=list_of_cells.end()){
    int index = (**it).index(vh);  //get the index of the infinite_vertex
    Cell_handle start_cell = (**it).neighbor(index);
    bool intersection_found = false;
    bool changed_res;
    for(int i=0; i<4; i++){
      Point p1 = (*start_cell).vertex(i)->point();
      bool res = transf_plane.has_on_positive_side(Point_3(p1.x(), p1.y(), p1.z()));
      if(i == 0)
        changed_res = res;
      else
        if(changed_res != res)
          intersection_found = true;
    }
    if(intersection_found){
      result_handle = start_cell;
      found = true;
    }
    else
      it++;
  }
  
  list_of_cell_handles.clear();
  if(found){    
    //list_of_cell_handles.push_back(result_handle);

    std::stack<Cell_handle, std::list<Cell_handle> > cell_stack;            //used for DFS
    std::map<Cell_handle, Cell_handle>  cell_map; //used for DFS
    typedef std::pair<Cell_handle, Cell_handle> C_PAIR;

    cell_stack.push(result_handle);
    cell_map.insert(C_PAIR(result_handle, result_handle));

    Cell_handle ch;
    while(!cell_stack.empty()){
      ch = cell_stack.top();      
      cell_stack.pop(); //done with this cell
      for(int i=0; i<4; i++) //visit all the neighbors
      {
        if(cell_map.find((*ch).neighbor(i)) == cell_map.end())
          if(!T.is_infinite((*ch).neighbor(i))){
            //test intersection
            bool intersection_found = false;
            bool changed_res;
            for(int j=0; j<4; j++){
              Point p1 = (*ch).vertex(j)->point();
              bool res = transf_plane.has_on_positive_side(Point_3(p1.x(), p1.y(), p1.z()));
              if(j == 0)
                changed_res = res;
              else
                if(changed_res != res)
                  intersection_found = true;;
            }
            if(intersection_found){
              cell_stack.push((*ch).neighbor(i));              
            }
            cell_map.insert(C_PAIR((*ch).neighbor(i), (*ch).neighbor(i)));
          }
      }//endfor
      list_of_cell_handles.push_back(ch);
    }//endwhile  

  }//end if(found)


}


//SoSeparator* get_plane(SoSFVec3f * draggerfield, SoSFRotation * draggerfieldr, Triangulation &t){
SoSeparator* get_plane(Triangulation &t){
    
  SoSeparator * sub = new SoSeparator;

  sub->ref();
  SoMaterial * mat = new SoMaterial;
  sub->addChild(mat);
  mat->transparency = 0.5;
  mat->diffuseColor.setValue(0.0, 1.0, 0.0);

  SoPickStyle * pickstyle = new SoPickStyle;
  sub->addChild(pickstyle);


  SoTransformBoxManip * plane_manip = new SoTransformBoxManip;
  //plane_manip->translation.setValue(5, 0, 0);
  sub->addChild(plane_manip);

  plane_manip->getDragger()->enableValueChangedCallbacks(true);
  plane_manip->getDragger()->addValueChangedCallback(draggerCB, NULL);

    Finite_vertices_iterator vit;
    Traits::FT xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;    
    for (vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) {
      if((*vit).point().x() < xmin)
        xmin = (*vit).point().x();
      if((*vit).point().y() < ymin)
        ymin = (*vit).point().y();
      if((*vit).point().z() < zmin)
        zmin = (*vit).point().z();
      if((*vit).point().x() > xmax)
        xmax = (*vit).point().x();
      if((*vit).point().y() > ymax)
        ymax = (*vit).point().y();
      if((*vit).point().z() > zmax)
        zmax = (*vit).point().z();
      vit++;
    }

  cube = new SoCube;
  sub->addChild(cube);
  cube->height.setValue(max(ymax-ymin, xmax - xmin));
  cube->width.setValue(max(ymax-ymin, xmax - xmin));
  cube->depth.setValue(0.05f);

  return sub;


}

int
main (int argc, char ** argv)
{
  // Initialize Coin, and return a main window to use
  // If unsuccessful, exit
  QWidget * window = SoQt::init(argv[0]); 
  if (window==NULL) exit(1);
    
  Node_tetrahedron_3<Kernel>::initClass();
  Node_polyhedron_3<Polyhedron>::initClass();
  Node_triangulation_3<Triangulation>::initClass();
  Node_terrain<Delaunay>::initClass();

/*
  //read the polyhedron
  const char* iname = "cin";
  std::istream*    p_in  = &std::cin;
  std::ifstream    in;
  in.open("d:\\coinsource\\coin_vc6\\data\\venus.off");
  p_in = &in;  
  if ( !*p_in)
    std::cout << "error: cannot open file for reading." <<endl;  
  CGAL::set_ascii_mode(* p_in);      
  (*p_in) >> P;
  in.close();
*/
/*
  in.open("d:\\coinsource\\coin_vc6\\data\\terrain.cin");
  std::istream_iterator<TPoint_3> begin(in);
  std::istream_iterator<TPoint_3> end;

  
  dt.insert(begin, end);
	*/
  int size = 100;
  int divs = 12;
  int NumTerrainIndices  = divs * divs * 2 * 3;
  int num_divisions = (int)sqrt(NumTerrainIndices/6);
	float plane_size = (float)size;

	float delta = plane_size/(float)num_divisions;
	float tex_delta = 2.0f/(float)num_divisions;
	float x_dist   = 0.0f;
	float z_dist   = 0.0f;
  srand( (unsigned)time( NULL ) );

	for (int i=0;i<=num_divisions;i++)
	{
    for (int j=0;j<=num_divisions;j++)
		{
			x_dist = (-0.5f * plane_size) + ((float)j*delta);
			z_dist = (-0.5f * plane_size) + ((float)i*delta);
      int rand_n = rand()%2;
      if(rand_n)
			  dt.push_back(TPoint_3(x_dist, 0, z_dist));
      else
        dt.push_back(TPoint_3(x_dist, 20, z_dist));
    }
  }


  /*
  in.open("d:\\coinsource\\coin_vc6\\david.off");
  p_in = &in;  
  if ( !*p_in)
    std::cout << "error: cannot open file for reading." <<endl;  
  CGAL::set_ascii_mode(* p_in);      
  (*p_in) >> P2;
  in.close();
*/
  //generate a triangulation
  std::list<Point> L;	 
  CGAL::Random_points_in_sphere_3<Point> g(4);
  for(int count=0; count<100; count++) {
    L.push_front(*g++);
  }
  int n = T.insert(L.begin(), L.end());


  root = new SoSeparator;
  SoSeparator * sep1 = new SoSeparator;
  SoSeparator * sep2 = new SoSeparator;
  SoSeparator * sep3 = new SoSeparator;
  SoSeparator * sep4 = new SoSeparator;
  Node_polyhedron_3<Polyhedron> *poly = new Node_polyhedron_3<Polyhedron>(P);
  Node_polyhedron_3<Polyhedron> *poly2 = new Node_polyhedron_3<Polyhedron>(P2);
  Node_triangulation_3<Triangulation> *triangl = new Node_triangulation_3<Triangulation>(T);
  Node_terrain<Delaunay> *terrain = new Node_terrain<Delaunay>(dt);
  SoCube * cube = new SoCube;

  SoMaterial * material = new SoMaterial;
  SoMaterial * material2 = new SoMaterial;
  SoColorIndex * indexC = new SoColorIndex;
  SoDrawStyle * style = new SoDrawStyle;
  SoDrawStyle * triangulation_style = new SoDrawStyle;  
  SoLightModel * lightm = new SoLightModel;
  SoPointLight * lightp = new SoPointLight;
  SoPointLight * lightb = new SoPointLight;
  lightp->location.setValue(0, 0, 1);
  lightp->color.setValue(1, 1, 1);
  lightb->location.setValue(0, 0, -1);
  lightb->color.setValue(1, 1, 0);
  SoTrackballManip * polyhedron_trackball = new SoTrackballManip;
  SoTrackballManip * triangulation_trackball = new SoTrackballManip;    
  SoEventCallback *myEventCB = new SoEventCallback;  
  SoTranslation * transl = new SoTranslation;
  SoRotation * rot = new SoRotation;
  SoTransformBoxManip * plane_manip = new SoTransformBoxManip;
  SoComplexity * complexity = new SoComplexity;
  

  style->style.setValue(SoDrawStyle::LINES);
  style->lineWidth = 4;
  triangulation_style->pointSize = 5;
  triangulation_style->lineWidth = 2;
  material->diffuseColor.setValue(0.8f, 0.2f, 0.0);
  material2->diffuseColor.setValue(0.0f, 0.0f, 1.0f);
  transl->translation.setValue(0, 3.5, 0);
  rot->rotation.setValue(SbVec3f(-0.348137f, 0540522.0f, 0.265922f), 3.225028f);

  root->ref();                   //increments the reference counter

  root->addChild(myEventCB);
  
  sep3->ref();
  sep3->addChild(lightp);
  SoCallback *triangle_callback = new SoCallback;
  triangle_callback->setCallback(render_custom);
  sep3->addChild(triangle_callback);
  root->addChild(sep3);
  

  sep2->ref();
  //sep1->addChild(lightp);
  sep2->addChild(material2);
  sep2->addChild(terrain);
  root->addChild(sep2);

  sep1->ref();
  //SoPickStyle * pickstyle = new SoPickStyle;
  //sep1->addChild(pickstyle);
  //pickstyle->style = SoPickStyle::UNPICKABLE;
  //sep1->addChild(polyhedron_trackball);  
  sep1->addChild(material);
  //sep1->addChild(lightm);
  //SoPickStyle * pickstyle2 = new SoPickStyle;
  //sep1->addChild(pickstyle2);  
  //sep1->addChild(style);
  //complexity->type.setValue(SoComplexity::SCREEN_SPACE);

  complexity->value.setValue(0.8f);
  sep1->addChild(complexity);
  sep1->addChild(lightp);
  sep1->addChild(lightb);
  //sep1->addChild(texture());
  sep1->addChild(poly);
  sep1->addChild(transl);
  
  //SoTransformBoxManip * head_manip = new SoTransformBoxManip;
  //sep1->addChild(head_manip);

  //head_manip->getDragger()->enableValueChangedCallbacks(true);
  //head_manip->getDragger()->addValueChangedCallback(headCB, NULL);

  //sep1->addChild(rot);

  //sep1->addChild(poly2);
  //root->addChild(sep1);

/*
  sep2->ref();
  sep2->addChild(material2);  
  //sep2->addChild(triangulation_trackball);    
  //sep2->addChild(triangulation_style);  
  sep2->addChild(triangl);
  SoCallback *cell_callback = new SoCallback;
  cell_callback->setCallback(render_cells);
  sep2->addChild(cell_callback);  
  root->addChild(sep2);
*/


/*

  sep4->ref();  
  sep4->addChild(get_plane(T));  
  root->addChild(sep4);
*/


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
