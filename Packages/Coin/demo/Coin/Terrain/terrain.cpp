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

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/So_node_terrain.h>


#include <qmainwindow.h>
#include <qapplication.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qfiledialog.h>
#include <qimage.h>
#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qplatinumstyle.h>
#include <qslider.h>
#include <qstatusbar.h>
#include <qtimer.h>
#include <qtoolbutton.h>
#include <qprogressdialog.h>

#include <list>
#include <map>
#include <stack>

Delaunay dt;
SoQtExaminerViewer * viewer;
SoSeparator *root;
Node_terrain<Delaunay> *terrain;

SoSeparator* get_main_scene(){

    dt.clear();
    Vector p;
    double value;
    double roughness = 0.5;
    int frequency    = 70;
    int gridSize = 50;
    double landscape[100][100];
    bool initFractals = true;

    QProgressDialog progress( "Generating terrain...", "Cancel generate", gridSize,
                            NULL, "progress", true );
    progress.setMinimumDuration(0);
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
      progress.setProgress( x );    
      for ( int y = 0; y < gridSize; y++ ) {
        if ( progress.wasCancelled() )
          exit(1);
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
    progress.setProgress( gridSize );


  SoSeparator * sep1 = new SoSeparator;
  SoMaterial * material = new SoMaterial;
  SoRotation * rot = new SoRotation;
  SoComplexity * complexity = new SoComplexity;
  terrain = new Node_terrain<Delaunay>(dt);

  rot->rotation.setValue(SbVec3f(1.0f, 0.0f, 0.0f), 30.0f);
  material->diffuseColor.setValue(0.8f, 0.2f, 0.0);
  complexity->value.setValue(0.7f);
  sep1->ref();
  sep1->addChild(rot);
  sep1->addChild(complexity);
  sep1->addChild(material);
  sep1->addChild(terrain);

  return sep1;

}

class Qt_layer_show_triangulation : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_triangulation(){};
private:
  void draw(){    
    *widget << CGAL::BLUE;
    Finite_edges_iterator it = dt.finite_edges_begin();
    while(it!=dt.finite_edges_end()){
      //*widget << Segment_2(Point_2(0, 0), Point_2(1, 1));
      TPoint_3 p1 =  (*(*((*it).first)).vertex(((*it).second + 1)%3)).point();
      TPoint_3 p2 =  (*(*((*it).first)).vertex(((*it).second + 2)%3)).point();
      Segment_2 s = Segment_2(Point_2(p1.x(), p1.y()), Point_2(p2.x(), p2.y()));
      *widget << s;
      it++;
    }
  }
};


class Layout_widget : public QWidget{
public:
  Layout_widget(QWidget *parent, const char *name=0):
      QWidget(parent, name) {
    SoQt::init(this);
    Node_terrain<Delaunay>::initClass();
    
    QSizePolicy p(QSizePolicy::Expanding, QSizePolicy::Expanding, true);
    widget = new CGAL::Qt_widget(this);
    widget->setMouseTracking(TRUE);
    widget->resize(300, 300);
    widget->setSizePolicy(p);
    QBoxLayout *topLayout1 = new QHBoxLayout( this);
    
    viewer = new SoQtExaminerViewer(this);
    SoSeparator *root = get_main_scene();
    viewer->setSceneGraph(root);    
    viewer->show();
    viewer->getBaseWidget()->resize(300, 300);

    topLayout1->addWidget(viewer->getBaseWidget());
    topLayout1->addWidget(widget);        
  }
  ~Layout_widget(){}
  CGAL::Qt_widget* get_qt_widget(){return widget;}
  SoQtExaminerViewer* get_viewer(){return viewer;}
private:
  CGAL::Qt_widget     *widget;  
  SoQtExaminerViewer  *viewer;
};

class MyWindow : public QMainWindow{
  Q_OBJECT
public:
  MyWindow(){
    Layout_widget *cwidget = new Layout_widget(this, "Main_layout");
    widget = cwidget->get_qt_widget();
    viewer = cwidget->get_viewer();
    setCentralWidget(cwidget);
    widget->attach(&ls);
    widget->set_window(0, 50, 0, 50);
    *widget << CGAL::BackgroundColor(CGAL::BLACK);


    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Terrain", this, SLOT(load_terrain()), CTRL+Key_L);
    file->insertSeparator();
    file->insertItem("&Print to ps", this, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );
    // draw widget menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate terrain", this,
		      SLOT(generate_terrain()), CTRL+Key_G );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);
    this->addToolBar(stoolbar->toolbar(), Top, FALSE);


  }
public slots:
  void generate_terrain(){
    dt.clear();
    Vector p;
    double value;
    double roughness = 0.5;
    int frequency    = 70;
    int gridSize = 50;
    double landscape[100][100];
    bool initFractals = true;

    QProgressDialog progress( "Generating terrain...", "Cancel generate", gridSize,
                            NULL, "progress", true );
    progress.setMinimumDuration(0);
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
      progress.setProgress( x );    
      for ( int y = 0; y < gridSize; y++ ) {
        if ( progress.wasCancelled() )
          exit(1);
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
    progress.setProgress( gridSize );
    terrain->lock();
    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
    terrain->unlock();
    terrain->touch();
    viewer->viewAll();
  }
  void load_terrain(){      
  }
  void print_to_ps(){
  }
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  
  void new_instance()
  {
  }

  void new_window(){
    MyWindow *ed = new MyWindow();
    ed->resize(600, 300);
    ed->setCaption("Layer");    
    ed->show();
    ed->set_window(0, 50, 0, 50);
  }

private:
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  Qt_layer_show_triangulation ls;
  SoQtExaminerViewer *viewer;
};

#include "terrain.moc"

int
main (int argc, char ** argv)
{
   QApplication app(argc, argv);   
   MyWindow *mainwin = new MyWindow();
   app.setMainWidget(mainwin);   
   mainwin->resize(600, 300);
   mainwin->setCaption("Terrain");
   mainwin->show();
   
   app.exec();
}