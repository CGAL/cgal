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
#include<CGAL/Point_2.h>
#include<CGAL/Cartesian.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <iostream>
#include <fstream>
#include <string>


#include<CGAL/IO/Qt_widget.h>
#include<CGAL/Visibility_complex_2.h>



#include<algorithm>
#include<cmath>
#include <CGAL/squared_distance_2.h>

typedef long N;
typedef CGAL::Cartesian<N> K;

#include <CGAL/Point_2.h>
typedef CGAL::Point_2<K> Point_2;

#ifdef VC_SCENE_DISPLAY_CIRCLE

#include<CGAL/Visibility_complex_2/Circle_traits.h>
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.center().x()-d.radius();
  xM=d.center().x()+d.radius();
  ym=d.center().y()-d.radius();
  yM=d.center().y()+d.radius();
}

Point_2 ref_point(Gt::Disk&d) {
  return d.center();
}


#elif defined(VC_SCENE_DISPLAY_POLYGON)

#include<CGAL/Visibility_complex_2/Polygon_traits.h>
typedef CGAL::Visibility_complex_2_polygon_traits<K> Gt;


#include <CGAL/IO/Qt_widget_Polygon_2.h>

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.left_vertex()->x();
  xM=d.right_vertex()->x();
  ym=d.bottom_vertex()->y();
  yM=d.top_vertex()->y();
}

Point_2 ref_point(Gt::Disk&d) {
  return *d.right_vertex();
}

#elif defined(VC_SCENE_DISPLAY_SEGMENT)
#include<CGAL/Visibility_complex_2/Segment_traits.h>
typedef CGAL::Visibility_complex_2_segment_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=std::min(d.source().x(),d.target().x());
  xM=std::max(d.source().x(),d.target().x());
  ym=std::min(d.source().y(),d.target().y());
  yM=std::max(d.source().y(),d.target().y());
}

Point_2 ref_point(Gt::Disk&d) {
  return CGAL::midpoint(d.source(),d.target());
}

#elif defined(VC_SCENE_DISPLAY_POINT)
#include<CGAL/Visibility_complex_2/Point_traits.h>
typedef CGAL::Visibility_complex_2_point_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.x();
  xM=d.x();
  ym=d.y();
  yM=d.y();
}

Point_2 ref_point(Gt::Disk&d) {
  return d;
}

#endif


typedef Gt::Bitangent_2 Bitangent_2;
typedef CGAL::Visibility_complex_2<Gt>::Constraint_input Constraint_input;


#include <CGAL/IO/Ostream_iterator.h>




typedef CGAL::Point_2<K> Point_2;
typedef Gt::Segment_2 Segment_2;

bool pausep;



int main(int argc,char ** argv) {
  std::istream * di=&std::cin;
  if (argc>1) di=new std::ifstream(argv[1]);
  std::istream_iterator<Gt::Disk> disk_it(*di),disk_end;
  std::istream_iterator<Constraint_input> 
    constraint_it,constraint_end;
  std::istream * ci=0;
  if (argc>2) {
    ci=new std::ifstream(argv[2]);
    constraint_it=
      std::istream_iterator<Constraint_input>(*ci);
  }
  std::vector<Gt::Disk> disks(disk_it,disk_end);
  std::vector<Bitangent_2> bitangents;
  for (;constraint_it!=constraint_end;++constraint_it) {
    Constraint_input c(*constraint_it);
    bitangents.push_back(Bitangent_2(c.type(),&(disks.begin()[c.source()]),
                         &(disks.begin()[c.target()])));
  }
  N xm=999999999,xM=-999999999,ym=999999999,yM=-999999999;
  for (std::vector<Gt::Disk>::iterator i=disks.begin();i!=disks.end();++i) {
    N xmin,xmax,ymin,ymax;
    bbox(*i,xmin,xmax,ymin,ymax);
    xm=std::min(xmin,xm);
    xM=std::max(xmax,xM);
    ym=std::min(ymin,ym);
    yM=std::max(ymax,yM);
  }
  double xmin=CGAL_NTS to_double(xm);
  double xmax=CGAL_NTS to_double(xM);
  double ymin=CGAL_NTS to_double(ym);
  double ymax=CGAL_NTS to_double(yM);
  double dx=xmax-xmin;
  double dy=ymax-ymin;
  int width=600;
  int height=(int)(width*dy/dx);

  int ac=1;
  char* av[1]={
   "Scene" 
  };
  QApplication app(ac,av);
  CGAL::Qt_widget* w;
  w = new CGAL::Qt_widget();
  app.setMainWidget( w );
  w->resize(width,height);
  w->set_window(xmin-0.1*dx,xmax+0.1*dx,ymin-0.1*dy,ymax+0.1*dy);
  w->show();
  w->lock();
  *w << CGAL::BLACK;
  for (std::vector<Gt::Disk>::iterator i=disks.begin();i!=disks.end();++i) {
    *w<<*i;
    std::ostringstream n;
    n<<(i-disks.begin());
    Point_2 p=ref_point(*i);
    w->get_painter().drawText(w->x_pixel(p.x()-10),
                              w->y_pixel(p.y()),
                              n.str());

  }
  *w << CGAL::GREEN;
  std::copy(bitangents.begin(),bitangents.end(),
            CGAL::Ostream_iterator<Segment_2,CGAL::Qt_widget>(*w));
  w->unlock();
  return app.exec();
}
