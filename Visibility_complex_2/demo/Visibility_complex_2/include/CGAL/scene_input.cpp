// argv[1] = File to store the disks. Defaults to stdout, in which case no 
// constraints are input
// argv[2] = File to store the constraints. If not given, no constraints
// are input.

// First input the disks. Then press enter. Then input the constraints.
// To begin a left_xx (respectively right_xx) constraint, click the source
// disk with the left (respectively right button). Then enter the target
// disk similarily.

// It is possible to pop the disks or constraints previously input
// by pressing backspace.

#include<CGAL/basic.h>

#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else

#include<qapplication.h>
#include<qmainwindow.h>
#include<iostream>
#include<fstream>

#include<algorithm>


#include<CGAL/Cartesian.h>
#include<CGAL/squared_distance_2.h>
#include<CGAL/Gmpz.h>
#include<CGAL/Point_2.h>
#include<CGAL/IO/Ostream_iterator.h>
#include<CGAL/IO/Qt_widget.h>


typedef CGAL::Gmpz N;

typedef CGAL::Cartesian<N> K;
typedef CGAL::Point_2<K> Point_2;

#ifdef VC_SCENE_INPUT_CIRCLE

#include<CGAL/Visibility_complex_2/Circle_traits.h>
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;


#include <CGAL/IO/Qt_widget_get_circle.h>
typedef CGAL::Circle_2<K> Cgal_disk;
typedef CGAL::Qt_widget_get_circle<K> Get_disk;

int find_disk(Point_2 p,std::vector<Gt::Disk>::iterator begin,
              std::vector<Gt::Disk>::iterator end) {
  for (std::vector<Gt::Disk>::iterator i=begin;i!=end;++i) {
    if (CGAL::squared_distance(p,i->center())<=i->squared_radius())
      return i-begin;
  }
  return -1;
}
#elif defined(VC_SCENE_INPUT_POLYGON)

#include<CGAL/Visibility_complex_2/Polygon_traits.h>
typedef CGAL::Visibility_complex_2_polygon_traits<K> Gt;


#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>
typedef Gt::Disk Cgal_disk;
typedef CGAL::Qt_widget_get_polygon<Gt::Disk> Get_disk;

int find_disk(Point_2 p,std::vector<Gt::Disk>::iterator begin,
              std::vector<Gt::Disk>::iterator end) {
  for (std::vector<Gt::Disk>::iterator i=begin;i!=end;++i) {
    if (i->has_on_bounded_side(p)) return i-begin;
  }
  return -1;
}
#elif defined(VC_SCENE_INPUT_SEGMENT)
#include<CGAL/Visibility_complex_2/Segment_traits.h>
typedef CGAL::Visibility_complex_2_segment_traits<K> Gt;


#include <CGAL/IO/Qt_widget_get_segment.h>
typedef Gt::Disk Cgal_disk;
typedef CGAL::Qt_widget_get_segment<K> Get_disk;

int find_disk(Point_2 p,std::vector<Gt::Disk>::iterator begin,
              std::vector<Gt::Disk>::iterator end) {
  int j=-1;
  N min_distance=999999999;
  for (std::vector<Gt::Disk>::iterator i=begin;i!=end;++i) {
    N d=CGAL::squared_distance(p,*i);
    if (d<min_distance) {
      min_distance=d;
      j=i-begin;
    }
  }
  return j;
}
#elif defined(VC_SCENE_INPUT_POINT)
#include<CGAL/Visibility_complex_2/Point_traits.h>
typedef CGAL::Visibility_complex_2_point_traits<K> Gt;


#include <CGAL/IO/Qt_widget_get_point.h>
typedef Gt::Disk Cgal_disk;
typedef CGAL::Qt_widget_get_point<K> Get_disk;

int find_disk(Point_2 p,std::vector<Gt::Disk>::iterator begin,
              std::vector<Gt::Disk>::iterator end) {
  int j=-1;
  N min_distance=999999999;
  for (std::vector<Gt::Disk>::iterator i=begin;i!=end;++i) {
    N d=CGAL::squared_distance(p,*i);
    if (d<min_distance) {
      min_distance=d;
      j=i-begin;
    }
  }
  return j;
}

#elif defined(VC_SCENE_INPUT_ELLIPSE)

#include <CGAL/IO/Qt_widget_get_conic.h>
#include<CGAL/Visibility_complex_2/Ellipse_traits.h>
typedef CGAL::Visibility_complex_2_ellipse_traits<K> Gt;

typedef CGAL::Conic_2<K> Cgal_disk;
typedef CGAL::Qt_widget_get_conic<K> Get_disk;

int find_disk(Point_2 p,std::vector<Gt::Disk>::iterator begin,
              std::vector<Gt::Disk>::iterator end) {
  for (std::vector<Gt::Disk>::iterator i=begin;i!=end;++i) {
    if (i->has_on_convex_side(p))
      return i-begin;
  }
  return -1;
}

#endif

typedef Gt::Segment_2 Segment_2;
typedef Gt::Bitangent_2 Bitangent_2;
typedef CGAL::Visibility_complex_2_details::Constraint_input Constraint_input;




class Scene_input : public QMainWindow {
  Q_OBJECT

  class Scene: public CGAL::Qt_widget_layer {
    Scene_input* super;
    CGAL::Color disk_color,constraint_color;
    int ignored;
  public:
    Scene(Scene_input* s): super(s),disk_color(CGAL::BLACK),
                     constraint_color(CGAL::GREEN), ignored(-1) {}
    void draw() {
      widget->lock();
      QColor c=widget->color();
      *widget<<disk_color;
      for (std::vector<Gt::Disk>::iterator i=super->disks.begin();
           i!=super->disks.end();++i)
        if (i-super->disks.begin()!=ignored) *widget<<*i;
      *widget<<constraint_color;
      std::copy(super->bitangents.begin(),super->bitangents.end(),
                CGAL::Ostream_iterator<Segment_2,CGAL::Qt_widget>
                (*widget));
      widget->setColor(c);
      widget->unlock();
    }
    void ignore(int i) {
      ignored=i;
    }
    
  };

  class Get_constraint : public CGAL::Qt_widget_layer {
    Scene_input* super;
    bool sourceL,targetL;
    int source,target;
  public:
    Get_constraint(Scene_input* s):super(s), source(-1),target(-1),
                             cursor_source(Qt::CrossCursor),
                             cursor_target(Qt::ArrowCursor),
                             source_color(CGAL::RED)
    {}
    bool abort() {
      super->scene.ignore(-1);
      if (source==-1) 
        return false; 
      else {
        source=-1; target=-1; return true;
      }
    }
  protected:
    void mousePressEvent(QMouseEvent *e) {
      bool * pos;
      int * obj;
      if (source==-1) {
        pos=&sourceL; obj=&source;
      } else {
        pos=&targetL; obj=&target;
      }
      if (e->button() == Qt::LeftButton) 
        *pos=true;
      else if (e->button() == Qt::RightButton)
        *pos=false;
      else return;
      N x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      Point_2 p(x,y);
      *obj=find_disk(p,super->disks.begin(),super->disks.end());
      if (target==-1) {
        if (source!=-1) {
          super->scene.ignore(source);
          widget->setCursor(cursor_target);
          widget->redraw();
        }
        return;
      }
      Bitangent_2::Type_util tu;
      Constraint_input c(tu(sourceL,targetL),source,target);
      source=-1;
      target=-1;
      widget->setCursor(cursor_source);
      super->scene.ignore(-1);
      widget->new_object(CGAL::make_object(c));
    }
    void activating()
    {
      oldcursor = widget->cursor();
      if (source==-1) 
        widget->setCursor(cursor_source);
      else if (target==-1)
        widget->setCursor(cursor_target);
    }
  
    void deactivating()
    {
      widget->setCursor(oldcursor);
    }
    void draw() {
      if (source!=-1) {
        QColor c=widget->color();
        *widget<<source_color;
        *widget<<super->disks.begin()[source];
        widget->setColor(c);
      }
    }
    QCursor cursor_source;
    QCursor cursor_target;
    QCursor oldcursor;
    CGAL::Color source_color;
  };



public:
  std::vector<Gt::Disk> disks;
  std::vector<Gt::Bitangent_2> bitangents;
  std::vector<Constraint_input> constraints;
  char * filed;
  char * filec;
  Scene_input(int x,int y,char * f1,char *f2) :
    scene(this), get_constraint(this) {
    stage=0;
    filed=f1;
    filec=f2;
    widget=new CGAL::Qt_widget(this,"atchoum");
    setCentralWidget(widget);
    resize(x,y);
    connect(widget,SIGNAL(s_keyPressEvent(QKeyEvent *)),
            this,SLOT(keyPressEvent(QKeyEvent *)));
    widget->set_window(0, x, 0, y);
    widget->attach(&scene);
    widget->attach(&get_disk);
    connect(widget,SIGNAL(new_cgal_object(CGAL::Object)),
            this,SLOT(get_new_object(CGAL::Object)));
  };
  CGAL::Qt_widget* widget;
  Get_disk get_disk;
  
  int stage;
  Scene scene;
  Get_constraint get_constraint;
  public slots:
  void get_new_object(CGAL::Object obj) {
    Cgal_disk C;
    if (CGAL::assign(C, obj)) {
      disks.push_back(C);
      widget->redraw();
      return;
    }
    Constraint_input c;
    if (CGAL::assign(c,obj)) {
      constraints.push_back(c);
      bitangents.push_back(Bitangent_2(c.type(),&(disks.begin()[c.source()]),
                                       &(disks.begin()[c.target()])));
      widget->redraw();
    }
  };
  void keyPressEvent(QKeyEvent *e) {
    switch (e->key()) {
    case Key_Return:      
      {
        if (stage==0&&filec) {
          stage++;
          widget->detach(&get_disk);
          widget->attach(&get_constraint);
          return;
        }
        std::ostream* o;
        if (filed)
          o=new std::ofstream(filed);
        else
          o=&std::cout;
        std::copy(disks.begin(),disks.end(),
                  std::ostream_iterator<Cgal_disk>(*o,"\n"));
        if (filed) delete o;
        o=0;
        if (stage>0) {
          o=new std::ofstream(filec);
          std::copy(constraints.begin(),constraints.end(),
                    std::ostream_iterator<Constraint_input>(*o,"\n"));
          delete o;
          o=0;
        }
        exit(0);
      }
      break;

    case Key_Backspace:
      if (stage==0) {
        if (!disks.empty()) disks.pop_back();        
      } else {
        if (!get_constraint.abort()&&!constraints.empty()) {
          constraints.pop_back();
          bitangents.pop_back();
        }
      }
      widget->redraw();
    }
  }
};





#include "scene_input.moc"

int main(int argc,char ** argv) {
  int ac=1;
  const char * ctitle="Input scene";
  char title[12];
  std::copy(ctitle,ctitle+12,title);
  char * av[1]={title};
  QApplication app(ac,av);
  Scene_input * pouloum = new Scene_input (800,600,(argc>1?argv[1]:0),
                               (argc>2?argv[2]:0));
  app.setMainWidget(pouloum);
  pouloum->show();
  app.exec();
}
#endif
