#include "Mesh.h"
#include <qapplication.h>
#include <qpainter.h>
#include <qfiledialog.h>
#include <stdlib.h>
#include <fstream>
#include "viewer.h"


TrFrame::TrFrame() : 
  QMainWindow(0, "Triangulation Frame", WDestructiveClose) {
  resize(800, 500);
  int id;
  QPopupMenu *tr=new QPopupMenu(this);
  menuBar()->insertItem("&Triangulation", tr);
  id = tr->insertItem("&Clear", this, SLOT(clearTriangulation()), CTRL+Key_C);
  tr->setWhatsThis(id, "Clears the triangulation");
  id = tr->insertItem("&Open", this, SLOT(openTriangulation()), CTRL+Key_O);
  tr->setWhatsThis(id, "Opens a triangulation");
  id = tr->insertItem("&Save", this, SLOT(saveTriangulation()), CTRL+Key_S);
  tr->setWhatsThis(id, "Save a triangulation");
  id = tr->insertItem("&Mesh", this, SLOT(mesh()), CTRL+Key_M);
  tr->setWhatsThis(id, "Computes the mesh"); 
  tr->insertSeparator();
  tr->insertItem("&Quit", qApp, SLOT(closeAllWindows()), CTRL+Key_Q);
  
  scr = new QScrollView(this);
  setCentralWidget(scr);
  trv = new TrViewer(this, scr);
  trv->resize(1000, 800);
  scr->addChild(trv);

  QToolBar *toolBar = new QToolBar(this);
  QPushButton *pbMesh = new QPushButton("Mesh", toolBar);
  connect(pbMesh, SIGNAL(clicked()), trv, SLOT(onMesh()));
  //pbMesh->setFixedSize(100, 32);
//   QPushButton *pbAutoMesh = new QPushButton("Auto Mesh", toolBar);
//   connect(pbAutoMesh, SIGNAL(toggled(bool)), trv, SLOT(toggleAutoMesh(bool)));
//   pbAutoMesh->setToggleButton(true);
  addToolBar(toolBar, "Controls");
  
  QToolBar *tbTrProps = new QToolBar(this);
  new QLabel("Triangulation sizes:", tbTrProps);
  new QLabel("  Width =  ", tbTrProps);
  editW = new QLineEdit(tbTrProps);
  QString sW;
  sW.setNum(trv->width());
  editW->setText(sW);
  new QLabel("  Height = ", tbTrProps);
  editH = new QLineEdit(tbTrProps);
  sW.setNum(trv->height());
  editH->setText(sW);
  QPushButton *pbChange = new QPushButton("Change", tbTrProps);
  connect(pbChange, SIGNAL(clicked()), this, SLOT(onChangeSizes()));
  addToolBar(tbTrProps, "Toolbar Propertied", Top,  true);
  //tbTrProps->setHorizontalStretchable(true);
  
  lblStatus = new QLabel("Ready", statusBar());
  statusBar()->addWidget(lblStatus);
}

void TrFrame::clearTriangulation() {
  trv->points.clear();
  trv->lines.clear();
  trv->mesh->reset();
  trv->update();
}

void TrFrame::saveTriangulation() {
  // trv->points.clear();
  // trv->lines.clear();
  //trv->mesh->clear();
  //trv->update();
  QFileDialog fd(this, "Save", TRUE);
  fd.setCaption("Save");
  fd.setMode(QFileDialog::AnyFile);
  fd.exec();
  ofstream of(fd.selectedFile());
  //of<<trv->mesh;
  trv->mesh->write(of);
}

void TrFrame::openTriangulation() {
//   trv->points.clear();
//   trv->lines.clear();
//   trv->mesh->clear();
//   trv->update();
  QFileDialog fd(this, "Open", TRUE);
  fd.setCaption("Open");
  fd.setMode(QFileDialog::ExistingFile);
  fd.exec();
  ifstream f(fd.selectedFile());
  //f>>trv->mesh;
  trv->mesh->read(f);
  trv->update();
}

void TrFrame::onChangeSizes() {
  int W, H;
  W = atoi(editW->text());
  H = atoi(editH->text());
  trv->resize(W, H);
}

void TrFrame::mesh() {
  trv->onMesh();
  update();
}

void TrFrame::setStatus(int x, int y) {
  QString s;
  s+=QString("x = ");
  s+=QString::number(x);
  s+=QString("  y = ");
  s+=QString::number(y);
  lblStatus->setText(s);
}

TrViewer::TrViewer(TrFrame *f, QWidget *parent) : 
  QWidget(parent)//, auto_mesh(false)
{
  frame = f;
  mesh= new Msh();
  setCursor(crossCursor);
  setMouseTracking(TRUE);
  //setFocusPolicy(QWidget::StrongFocus);
  setFocus();
}

void TrViewer::paintEvent(QPaintEvent *pe) {
    QPainter p(this);
  p.fillRect(0, 0, width(), height(), QColor(0,0,0));
  // paint the triangulation
  
   p.setPen(QColor(0, 0, 255));
   p.setBrush(QBrush(NoBrush));
   Msh::Edge_iterator eit = mesh->edges_begin();
   while(eit != mesh->edges_end()) {
     Msh::Face_handle face = (*eit).first;
     int nedge = (*eit).second;
     double x1, x2, y1, y2;
     double D;
     x1 = face->vertex(mesh->cw(nedge))->point().x();
     y1 = face->vertex(mesh->cw(nedge))->point().y();
     x2 = face->vertex(mesh->ccw(nedge))->point().x();
     y2 = face->vertex(mesh->ccw(nedge))->point().y();
     D = ::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
     
//      p.setPen(QColor(85, 85, 85));
//      p.drawEllipse((int)rint((x1+x2)/2.0-D/2.0), (int)rint((y1+y2)/2.0-D/2.0),      (int)D, (int)D); 
     if(face->is_constrained(nedge)) {
       p.setPen(QColor(255, 0, 0));
     } else {
       p.setPen(QColor(0, 0, 255));
     }
     p.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
     eit++;
   }
   p.setPen(QColor(0, 0, 255));
//   p.setPen(QColor(0, 0, 255));
//   Msh::Face_iterator fit = mesh->faces_begin();
//   while(fit != mesh->faces_end()) {
//     double x1, x2, x3, y1, y2, y3;
//     x1 =  (*fit).vertex(0)->point().x()/*.to_double()*/;
//     x2 =  (*fit).vertex(1)->point().x()/*.to_double()*/;
//     x3 =  (*fit).vertex(2)->point().x()/*.to_double()*/;
//     y1 =  (*fit).vertex(0)->point().y()/*.to_double()*/;
//     y2 =  (*fit).vertex(1)->point().y()/*.to_double()*/;
//     y3 =  (*fit).vertex(2)->point().y()/*.to_double()*/;
//     p.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
//     p.drawLine((int)x1, (int)y1, (int)x3, (int)y3);
//     p.drawLine((int)x3, (int)y3, (int)x2, (int)y2); 
//     fit++;
//   }
   Msh::Vertex_iterator vit = mesh->vertices_begin();
   p.setPen(QColor(0, 255, 0));
   p.setBrush(QBrush(QColor(0,255,0)));
   while(vit != mesh->vertices_end()){
     p.drawEllipse((int)(*vit).point().x()/*.to_double()*/-1, 
		  (int)(*vit).point().y()/*.to_double()*/-1, 3, 3);
     vit++;
   } 
   
  // paint the constraints
   

}

void TrViewer::mouseMoveEvent(QMouseEvent *me) {
  QString s;
  int x, y;
  if(dragging) {
    QPainter p(this);
    p.setPen(QColor(255,255,255));
    p.setRasterOp(XorROP);
    p.drawLine(startx, starty, oldx, oldy);
  }
  x = me->x();
  y = me->y();
  frame->setStatus(x<0?0:(x<width()?x:width()-1), 
		   y<0?0:(y<height()?y:height()-1));
  if(dragging) {
    QPainter p(this);
    p.setPen(QColor(255,255,255));
    p.setRasterOp(XorROP);
    p.drawLine(startx, starty, x, y);
    oldx=x;
    oldy=y;
  }
}

void TrViewer::mousePressEvent(QMouseEvent *me) {
  int x, y;
  if(me->button() == LeftButton) {
    x = me->x();
    y = me->y();
    
    startx = x<0?0:(x<width()?x:width()-1);
    starty = y<0?0:(y<height()?y:height()-1);
    oldx=startx;
    oldy=starty;
    dragging = true;
  }
}

void TrViewer::mouseReleaseEvent(QMouseEvent *me) {
  int x, y;
  if(me->button() == LeftButton && dragging) {
    x = me->x();
    y = me->y();
    endx = x<0?0:(x<width()?x:width()-1);
    endy = y<0?0:(y<height()?y:height()-1);
    dragging = false;
    
    if(startx == endx && starty == endy) {
      points.push_back(Point(startx, starty));
      cout<<"point : "<<startx<<" "<<starty<<endl;
      Msh::Point p(endx, endy);
      cout<<"inserting point... "<<flush;
    //p.y() = itp->y;
      //Msh::Vertex v(p);
      //v.set_point(p);
      mesh->insert(p);
      cout<<"ok"<<endl;
    } else {
      points.push_back(Point(startx, starty));
      points.push_back(Point(endx, endy));
      lines.push_back(Line(Point(startx, starty), Point(endx, endy)));
      cout<<"line : "<<startx<<" "<<starty<<" : "<<endx<<" "<<endy<<endl;
      Msh::Point p1(startx, starty); //mesh->insert(p1);
      Msh::Point p2(endx, endy); //mesh->insert(p2);
      mesh->insert(p1, p2);
    }
    update();
    /*    if(auto_mesh)
	  mesh->refine_mesh();*/
    //    update();
  }
}

void TrViewer::keyPressEvent(QKeyEvent *ke) {
  if(ke->key()==Key_Escape) {
    dragging=false;
    cout<<"aborted"<<endl;
    update();
  }
}

void TrViewer::onMesh(){
  try {
    mesh->refine_mesh();
    update();
  } catch(...) {
    cerr<<"CGAL threw an exception"<<endl;
  }
}

// void TrViewer::toggleAutoMesh(bool b) {
//   auto_mesh=b;
// }

// moc_source_file: viewer.h
#include "viewer.moc"
