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
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/SoOffScreenRenderer.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/So_node_terrain_constrained.h>


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
SoComplexity * complexity_node;
SoPerspectiveCamera *cam1;


SoSeparator* get_main_scene(){

  dt.clear();
  SoSeparator * sep1 = new SoSeparator;
  SoMaterial * material = new SoMaterial;
  //SoRotation * rot = new SoRotation;
  complexity_node = new SoComplexity;
  terrain = new Node_terrain<Delaunay>(dt);

  //rot->rotation.setValue(SbVec3f(1.0f, 0.0f, 0.0f), 30.0f);
  material->diffuseColor.setValue(0.5f, 0.6f, 0.0f);
  complexity_node->value.setValue(0.7f);

  cam1 = new SoPerspectiveCamera;  

  sep1->ref();
  sep1->addChild(cam1);
  //sep1->addChild(rot);
  sep1->addChild(complexity_node);
  sep1->addChild(material);
  sep1->addChild(new SoDirectionalLight);
  sep1->addChild(terrain);
  

  return sep1;

}

class Qt_layer_show_triangulation : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_show_triangulation(){};
private:
  void draw(){    
    
    Finite_edges_iterator it = dt.finite_edges_begin();
    while(it!=dt.finite_edges_end()){      
      TPoint_3 p1 =  (*(*((*it).first)).vertex(((*it).second + 1)%3)).point();
      TPoint_3 p2 =  (*(*((*it).first)).vertex(((*it).second + 2)%3)).point();
      Segment_2 s = Segment_2(Point_2(p1.x(), p1.y()), Point_2(p2.x(), p2.y()));
      if(dt.is_constrained((*it)))
        *widget << CGAL::RED;
      else
        *widget << CGAL::BLUE;
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
    *widget << CGAL::BackgroundColor(CGAL::BLACK);
    widget->setMouseTracking(TRUE);
    widget->resize(300, 300);
    widget->setSizePolicy(p);
    QBoxLayout *topLayout1 = new QHBoxLayout( this);
    
    viewer = new SoQtExaminerViewer(this);
    SoSeparator *root = get_main_scene();
    viewer->setSceneGraph(root);
    viewer->setDrawStyle(SoQtViewer::STILL, SoQtViewer::VIEW_WIREFRAME_OVERLAY);    
    viewer->getBaseWidget()->resize(300, 300);

    topLayout1->addWidget(viewer->getBaseWidget());
    //topLayout1->addWidget(viewer->getRenderAreaWidget());
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
    cwidget = new Layout_widget(this, "Main_layout");
    widget = cwidget->get_qt_widget();
    viewer = cwidget->get_viewer();    
    setCentralWidget(cwidget);
    widget->attach(&ls);    
    widget->set_window(0, 50, 0, 50);

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Load &Terrain Constrained", this, 
      SLOT(load_terrain_constrained()), CTRL+Key_T);
    file->insertItem("&Load Terrain", this, SLOT(load_terrain()), CTRL+Key_L);
    file->insertItem("&Save Terrain", this, SLOT(save_terrain()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("&Load Points", this, SLOT(load_masspts()), CTRL+Key_I);
    file->insertItem("&Load Constraints", this, SLOT(load_breaklines()), CTRL+Key_B);
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
    // draw widget menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("&Generate terrain", this,
		      SLOT(generate_terrain()), CTRL+Key_G );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this);    
    this->addDockWindow(stoolbar->toolbar(), DockRight, FALSE);
    this->show();
  }
public slots:

void load_masspts()
{
  terrain->lock();
  dt.clear();
  int i, j;
  int first = 0;
  char c;

  QString s( QFileDialog::getOpenFileName( QString::null,
                       "Mass Points(*.dat)", this ) );
  if ( s.isEmpty() ){
    return;
  }
  Bbox b;


  std::ifstream ins(s);
  do {
    ins.get(c);
    if(c == '#') {
      ins.ignore(200,'\n');
    }
  }while(c == '#');
  ins.unget();

  bool done = false;
  while(! done) {
    TPoint_3 p;
    ins >> i >> j >> p;
    if(ins.fail()){
      done = true;
      continue;
    }
    assert(j == 0);

    dt.insert(p);
    if(first == 0){
      b = p.bbox();
      first++;
    } else {
      b = b + p.bbox();
    }
  }  

    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
    cam1->position.setValue(SbVec3f(b.xmin(), b.ymin(), 2*b.zmax()));    
    cam1->pointAt(SbVec3f(b.xmin() + (b.xmax()-b.xmin())/100, b.ymin() + (b.ymax()-b.ymin())/100, b.zmax()), SbVec3f(0, 0, 1));
    widget->clear_history();
    widget->set_window(b.xmin(), b.xmax(), b.ymin(), b.ymax());
    terrain->unlock();
    terrain->touch();
    //viewer->viewAll();
    viewer->render();
    widget->redraw();
}

void load_breaklines()
{
  QString s( QFileDialog::getOpenFileName( QString::null,
                       "Break Lines(*.dat)", this ) );
  if ( s.isEmpty() ){
    return;
  }
  bool first = true;
  std::ifstream ins(s);
  char c;
  do {
    ins.get(c);
    if(c == '#') {
      ins.ignore(200,'\n');
    }
  }while(c == '#');
  ins.unget();

  bool done = false;
  while(! done) {
    TPoint_3 p, q;
    Vertex_handle vp, vq;
    int iq, j;
    ins >> iq >> j >> q;
    if(ins.fail()){
      done = true;
      continue;
    }
    if(first) { // we start a new polyline
      p = q;
      vp = dt.insert(p);
      first = false;
    } else {
      vq = dt.insert(q);
      dt.insert_constraint(vp,vq);
      p = q;
      vp = vq;
    }
  }

    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
    widget->clear_history();    
    terrain->touch();
    //viewer->viewAll();
    viewer->render();
    widget->redraw();

}


  void generate_terrain(){
    dt.clear();
    terrain->lock();
    Vector p;
    double value;
    double roughness = 0.5;
    int frequency    = 70;
    int gridSize = 200;
    double landscape[300][300];
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
        for(int count = 0; count < 10; count++){
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
    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
        
    Finite_vertices_iterator it = dt.finite_vertices_begin();        
    double  xmin = (*it).point().x(), ymin = (*it).point().y(), zmin = (*it).point().z(),
            xmax = (*it).point().x(), ymax = (*it).point().y(), zmax = (*it).point().z();
    while(it != dt.finite_vertices_end()){
      if((*it).point().x() < xmin)
        xmin = (*it).point().x();
      if((*it).point().y() < ymin)
        ymin = (*it).point().y();
      if((*it).point().x() > xmax)
        xmax = (*it).point().x();
      if((*it).point().y() > ymax)
        ymax = (*it).point().y();
      if((*it).point().z() > zmax)
        zmax = (*it).point().z();
      it++;
    }    
    cam1->position.setValue(SbVec3f(xmin, ymin, 2*zmax));
    cam1->pointAt(SbVec3f(xmin + (xmax-xmin)/5, ymin + (ymax-ymin)/5, zmax), SbVec3f(0, 0, 1));
    widget->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
    terrain->unlock();
    terrain->touch();
    //viewer->viewAll();
    viewer->render();
    widget->redraw();
    
  }
  void change_bbox(){
    complexity_node->value.setValue(0);
    terrain->touch();
  }

  void change_face(){
    complexity_node->value.setValue(0.7f);
    terrain->touch();
  }

  void change_smooth(){
    complexity_node->value.setValue(1.0f);
    terrain->touch();
  }

  // read number of points then the points and insert them in the triangulation
  // read the number of polylines and for each polyline read number of points and
  // insert them as long as the segment formed by each 2 points as a constraint
  void load_terrain_constrained(){
    dt.clear();
    int i, j;
    int number_of_points;

    QString s( QFileDialog::getOpenFileName( QString::null,
                        "Points and Constraint File (*.pc);; All files (*.*)", this ) );
    if ( s.isEmpty() ){
      return;
    }
    Bbox b;

    std::ifstream ins(s);
    ins >> number_of_points;
    for(i = 0; i < number_of_points; i++) {
      TPoint_3 p;
      ins >> p;
      dt.insert(p);
      if(i == 0){
        b = p.bbox();
      } else {
        b = b + p.bbox();
      }
    }


    int  number_of_polylines;
    ins >> number_of_polylines;

    for(i = 0; i < number_of_polylines; i++){
      ins >> number_of_points;
      TPoint_3 p;
      ins >> p;
      b = b + p.bbox();
      Vertex_handle ph = dt.insert(p);
      for(j = 1; j < number_of_points; j++) {
        TPoint_3 q;
        ins >> q;
        b = b+ q.bbox();
        Vertex_handle qh = dt.insert(q);
        dt.insert_constraint(ph,qh);
        ph = qh;
        p = q;
      }
    }
    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
    cam1->position.setValue(SbVec3f(b.xmin() - (b.xmax()-b.xmin())/2, b.ymin() - (b.ymax()-b.ymin())/2, 4*b.zmax()));
    cam1->pointAt(SbVec3f(b.xmin() + (b.xmax()-b.xmin())/2, b.ymin() + (b.ymax()-b.ymin())/2, b.zmin()), SbVec3f(0, 0, 1));
    widget->clear_history();
    widget->set_window(b.xmin(), b.xmax(), b.ymin(), b.ymax());
    terrain->touch();
    //viewer->viewAll();
    viewer->render();
    widget->redraw();
    
  }

  //load terrain using:
  // - the input operator if .cgal extension
  // - input operator of points and inserting points with insert method if no extension
  void load_terrain(){
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal);;All files (*.*)", this ) );
    if ( s.isEmpty() )
        return;
      
    std::ifstream in(s);
    //CGAL::set_ascii_mode(out);
    dt.clear();
    
    if(s.right(5) == ".cgal")
    {
      in >> dt;
    } else {
      TPoint_3 p;
      while (in){
        in >> p;
        dt.insert(p);
      }
    }
    //in >> dt;
    terrain->compute_normals_for_faces();
    terrain->compute_normals_for_vertices();
    
    Finite_vertices_iterator it = dt.finite_vertices_begin();        
    double  xmin = (*it).point().x(), ymin = (*it).point().y(), zmin = (*it).point().z(),
            xmax = (*it).point().x(), ymax = (*it).point().y(), zmax = (*it).point().z();
    while(it != dt.finite_vertices_end()){
      if((*it).point().x() < xmin)
        xmin = (*it).point().x();
      if((*it).point().y() < ymin)
        ymin = (*it).point().y();
      if((*it).point().x() > xmax)
        xmax = (*it).point().x();
      if((*it).point().y() > ymax)
        ymax = (*it).point().y();
      if((*it).point().z() > zmax)
        zmax = (*it).point().z();
      it++;
    }    
    cam1->position.setValue(SbVec3f(xmin , ymin , 4*zmax));
    cam1->pointAt(SbVec3f(xmin + (xmax-xmin)/2, ymin + (ymax-ymin)/2, zmin), SbVec3f(0, 0, 1));

    widget->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
    terrain->touch();    
    //viewer->viewAll();
    viewer->render();
    widget->redraw();
  }

  //save terrain using the output operator of Delaunay Triangulation
  void save_terrain(){
    QString fileName = 
      QFileDialog::getSaveFileName( "triangulation.cgal", 
				  "Cgal files (*.cgal);;All files (*.*)", this ); 

    if ( !fileName.isNull() ) {
      // got a file name
      std::ofstream out(fileName);
      //CGAL::set_ascii_mode(out);      
      out << dt;
    }
  }
  void print_to_ps(){
    
    QFileDialog qfd(this, "Save Texture", true);
    qfd.setViewMode(QFileDialog::Detail);    
    qfd.addFilter("Postscript files(*.eps)");
    qfd.setMode(QFileDialog::AnyFile);
    qfd.setCaption("Save Texture");
    
    QString fileName;
    if ( qfd.exec() == QDialog::Accepted )
      fileName = qfd.selectedFile();

    if ( !fileName.isNull() ) {
      SbViewportRegion myViewport(800, 600);
      SoOffscreenRenderer *myRenderer = new SoOffscreenRenderer(myViewport);      
      myRenderer->setBackgroundColor(SbColor(1.0f, 1.0f, 1.0f));
      if(!myRenderer->render(root)){
        delete myRenderer;
        return;
      }
      //texture->image.setValue(SbVec2s(textureWidth, textureHeight), SoOffscreenRenderer::RGB, myRenderer->getBuffer());
      //myRenderer->writeToRGB(fileName);
      myRenderer->writeToPostScript(fileName);
      delete myRenderer;
      return;
    }
  }
  void set_window(double xmin, double xmax, double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }
  
  void new_instance()
  {
    dt.clear();
    terrain->touch();
    widget->redraw();
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
  Layout_widget *cwidget;
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
   
   app.exec();
}