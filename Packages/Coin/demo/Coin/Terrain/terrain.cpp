#include "cgal_types.h"
#include "fBm.h"
#include <stdlib.h> // exit()
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoComplexity.h>

#include <CGAL/IO/So_node_terrain.h>


#include <list>
#include <map>
#include <stack>

Delaunay dt;
SoQtExaminerViewer * viewer;
SoSeparator *root;

int
main (int argc, char ** argv)
{
  // Initialize Coin, and return a main window to use
  // If unsuccessful, exit
  QWidget * window = SoQt::init(argv[0]); 
  if (window==NULL) exit(1);
    
  Node_terrain<Delaunay>::initClass();

  Vector p;
  double value;
  double roughness = 0.5;
  int frequency    = 70;
  int gridSize = 50;
  double landscape[100][100];
  bool initFractals = true;

    for ( int x = 0; x < gridSize; x++ )
      for ( int y = 0; y < gridSize; y++ )
        landscape[x][y] = 0;  

  p.x = p.y = p.z = 0;
  // Initialise fbm routine
  if ( initFractals ) {
    initFractals = false;
    value = fBm( p, roughness, 2.0, 8.0, 1 );
  }

  
    // Fractalize grid
    for ( int x = 0; x < gridSize; x++ ) {
      for ( int y = 0; y < gridSize; y++ ) {
        for(int count = 0; count < 6; count++){
          p.x = (double) x / (101 - frequency);
          p.y = (double) y / (101 - frequency);
	        p.z = (double) landscape[x][y] / (101 - frequency);
	        value = fBm(p, roughness, 2.0, 8.0, 0);
	        landscape[x][y] += value;
          if(landscape[x][y] < 0)
            landscape[x][y] = 0;
        }        
        dt.push_back(TPoint_3(x, y, landscape[x][y]));
	    }
    }
  


/*  
  int num_divisions = 50;
  float x_dist   = 0.0f;
	float y_dist   = 0.0f;
  srand( (unsigned)time( NULL ) );
  int delta = 100;
	for (int i=0;i<=num_divisions;i++)
	{
    x_dist = -num_divisions/2*delta + ((float)i*delta);
    for (int j=0;j<=num_divisions;j++)
		{			
			y_dist = -num_divisions/2*delta + ((float)j*delta);
      int rand_n = rand()%200;
      dt.push_back(TPoint_3(x_dist, y_dist, rand_n));
    }
  }

*/
  root = new SoSeparator;
  SoSeparator * sep1 = new SoSeparator;
  SoMaterial * material = new SoMaterial;
  SoRotation * rot = new SoRotation;
  SoComplexity * complexity = new SoComplexity;
  Node_terrain<Delaunay> *terrain = new Node_terrain<Delaunay>(dt);

  rot->rotation.setValue(SbVec3f(1.0f, 0.0f, 0.0f), 30.0f);
  material->diffuseColor.setValue(0.8f, 0.2f, 0.0);
  complexity->value.setValue(1.0f);


  root->ref();                   //increments the reference counter

  sep1->ref();
  sep1->addChild(rot);
  sep1->addChild(complexity);
  sep1->addChild(material);
  sep1->addChild(terrain);
  root->addChild(sep1);

  // Set up the ExaminerViewer
  viewer = new SoQtExaminerViewer(window);
  viewer->setSceneGraph(root);
  viewer->setTitle("Draw Style");
  viewer->viewAll();  
  viewer->show();

  //viewer->setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_AS_IS); 
  //viewer->setDrawStyle(SoQtViewer::STILL, SoQtViewer::VIEW_WIREFRAME_OVERLAY);

  SoQt::show(window); // display the main window
  SoQt::mainLoop();   // main Coin event loop
  delete viewer;      // free all the viewers resources
  root->unref();      // decrements the reference counter  
  return 0;


}