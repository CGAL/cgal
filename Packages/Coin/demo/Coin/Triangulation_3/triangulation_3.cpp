#include "cgal_types.h"
#include <stdlib.h> // exit()
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/events/SoMouseButtonEvent.h>
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
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/manips/SoTrackballManip.h>
#include <Inventor/manips/SoTransformBoxManip.h>
#include <Inventor/draggers/SoDragPointDragger.h>
#include <Inventor/draggers/SoTransformBoxDragger.h>

#include <CGAL/IO/So_node_triangulation_3.h>

#include <list>
#include <map>
#include <stack>

Triangulation	T;
SoCube *cube;
SoQtExaminerViewer * viewer;
SoSeparator *root;
std::list<Cell_handle>  list_of_cell_handles;
Plane_3 plane(Point_3(0, 0, 0), Point_3(0, 1, 0), Point_3(1, 0, 0));
Point_3 pp1, pp2, pp3;
SbVec3f normal3f;





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
  cube->height.setValue(std::max(ymax-ymin, xmax - xmin));
  cube->width.setValue(std::max(ymax-ymin, xmax - xmin));
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
 
  Node_triangulation_3<Triangulation>::initClass();

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
  SoMaterial * material = new SoMaterial;
  SoCallback *cell_callback = new SoCallback;
  
  material->diffuseColor.setValue(0.0f, 0.0f, 1.0f);
  cell_callback->setCallback(render_cells);

  Node_triangulation_3<Triangulation> *triangl = new Node_triangulation_3<Triangulation>(T);  


  sep1->ref();
  sep1->addChild(material);  
  sep1->addChild(triangl);  
  sep1->addChild(cell_callback);  
  root->addChild(sep1);

  sep2->ref();  
  sep2->addChild(get_plane(T));  
  root->addChild(sep2);



  // Set up the ExaminerViewer
  viewer = new SoQtExaminerViewer(window);
  viewer->setSceneGraph(root);
  viewer->setTitle("Triangulation_3 demo");  
  viewer->viewAll();  
  viewer->show();

  viewer->setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_BBOX);  

  SoQt::show(window); // display the main window
  SoQt::mainLoop();   // main Coin event loop
  delete viewer;      // free all the viewers resources
  root->unref();      // decrements the reference counter  
  return 0;
}
