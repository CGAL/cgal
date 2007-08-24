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
#include <CGAL/Gmpz.h>


typedef CGAL::Gmpz N;
typedef CGAL::Cartesian<N> K;

#include <CGAL/Point_2.h>
typedef CGAL::Point_2<K> Point_2;

#ifdef VC_SCENE_COMPUTE_BIT_CIRCLE

#include<CGAL/Visibility_complex_2/Circle_traits.h>
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.center().x()-d.radius();
  xM=d.center().x()+d.radius();
  ym=d.center().y()-d.radius();
  yM=d.center().y()+d.radius();
}

#elif defined(VC_SCENE_COMPUTE_BIT_POLYGON)

#include<CGAL/Visibility_complex_2/Polygon_traits.h>
typedef CGAL::Visibility_complex_2_polygon_traits<K> Gt;


#include <CGAL/IO/Qt_widget_Polygon_2.h>

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.left_vertex()->x();
  xM=d.right_vertex()->x();
  ym=d.bottom_vertex()->y();
  yM=d.top_vertex()->y();
}


#elif defined(VC_SCENE_COMPUTE_BIT_SEGMENT)
#include<CGAL/Visibility_complex_2/Segment_traits.h>
typedef CGAL::Visibility_complex_2_segment_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=std::min(d.source().x(),d.target().x());
  xM=std::max(d.source().x(),d.target().x());
  ym=std::min(d.source().y(),d.target().y());
  yM=std::max(d.source().y(),d.target().y());
}


#elif defined(VC_SCENE_COMPUTE_BIT_POINT)
#include<CGAL/Visibility_complex_2/Point_traits.h>
typedef CGAL::Visibility_complex_2_point_traits<K> Gt;

void bbox(Gt::Disk&d,N&xm,N&xM,N&ym,N&yM) {
  xm=d.x();
  xM=d.x();
  ym=d.y();
  yM=d.y();
}
#endif


typedef Gt::Bitangent_2 Bitangent_2;
typedef CGAL::Compute_free_bitangents_2<Gt> VC;
typedef VC::Constraint_input Constraint_input;


#include <CGAL/IO/Ostream_iterator.h>




typedef CGAL::Point_2<K> Point_2;
typedef Gt::Segment_2 Segment_2;

bool pausep;



int main(int argc,char ** argv) {
  std::istream * di=&std::cin;
  if (argc>1) di=new std::ifstream(argv[1]);
  std::istream_iterator<Gt::Disk> disk_it(*di),disk_end;
  std::vector<Gt::Disk> disks(disk_it,disk_end);
  std::istream_iterator<Constraint_input> 
    constraint_it,constraint_end;
  std::istream * ci=0;
  if (argc>2) {
    ci=new std::ifstream(argv[2]);
    constraint_it=
      std::istream_iterator<Constraint_input>(*ci);
  }
  std::vector<Bitangent_2> vb;
  VC()(disks.begin(),disks.end(),
       constraint_it,constraint_end,std::back_inserter(vb));
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
   "Free bitangents" 
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
  std::copy(disks.begin(),disks.end(),
            CGAL::Ostream_iterator<Gt::Disk,CGAL::Qt_widget>(*w));
  *w << CGAL::RED;
  std::copy(vb.begin(),vb.end(),
            CGAL::Ostream_iterator<Segment_2,CGAL::Qt_widget>(*w));
  w->unlock();
  std::cout<<vb.size()<<"\n";
  return app.exec();
}
