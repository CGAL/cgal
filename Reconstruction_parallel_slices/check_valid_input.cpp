#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;

template<class Iterator>
std::ostream& write_polygon(Iterator begin,Iterator end,double elev,int cc,std::ostream& out){
  double coords[3];
  coords[cc]=elev;
  
  out << std::distance(begin,end)+1 << std::endl;
  for (Iterator it=begin;it!=end;++it){
    coords[(cc+1)%3]=it->x();
    coords[(cc+2)%3]=it->y();
    out << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
  }
  coords[(cc+1)%3]=begin->x();
  coords[(cc+2)%3]=begin->y();
  out << coords[0] << " " << coords[1] << " " << coords[2] << "\n\n";
  return out;
}

bool read_polygon(std::list<Point_2>& points,int cc,double& elev,std::istream& input){
  double nbpt;
  if (input >> nbpt && input){
    double coords[3];
    input >> coords[0] >> coords[1] >> coords[2];
    std::stringstream sstream;
    elev=coords[cc];
    
    Point_2 p(coords[(cc+1)%3],coords[(cc+2)%3]);
    points.push_back(p);
    Point_2 prev=p;
    for (int i=1;i<nbpt;++i){
      input >> coords[0] >> coords[1] >> coords[2];
      if (coords[cc]!=elev){
        std::cout << coords[cc] << " vs " << elev << std::endl;
        exit(EXIT_FAILURE);
      }
      p=Point_2(coords[(cc+1)%3],coords[(cc+2)%3]);
      if (p!=prev)
        points.push_back(p);
      prev=p;
    }
    
    if (points.front()!=points.back() ){
      std::cout << "Polygon is not closed" << std::endl;
    }
    
    while (points.front()==points.back() ){
      if (points.empty()) break;
      points.pop_back();
    }
    return true;
  }
  return false;
}

bool set_polygon(Polygon_2& poly,std::list<Point_2>& points,double elev,int cc,std::ostream& bad){
    if (points.size()<=2){
      std::cout << "Polygon made of one or two points found" << std::endl;
      return false;
    }
      
    poly=Polygon_2(points.begin(),points.end());
    
    if ( !poly.is_simple() ){
      std::cout << "Non-simple polygon found\n";
      write_polygon(poly.vertices_begin(),poly.vertices_end(),elev,cc,bad);
      return false;
    }
  return true;
}

#include <CGAL/intersections.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>

typedef std::vector<Segment_2>                               Segments;

template<class PolygonIterator>
struct My_Box : public CGAL::Box_intersection_d::Box_d< double, 2, CGAL::Box_intersection_d::ID_EXPLICIT>
{
  PolygonIterator it_poly;
  std::pair<const Point_2*,const Point_2*> segment;
  
  Segment_2 get_segment() const {
    return Segment_2(*segment.first,*segment.second);
  }
  
  typedef CGAL::Box_intersection_d::Box_d< double, 2, CGAL::Box_intersection_d::ID_EXPLICIT> Base;
  My_Box(const Point_2& p1,const Point_2& p2,PolygonIterator it)
    :Base(p1.bbox()+p2.bbox()),it_poly(it),segment(&p1,&p2)
  {}
};

template <class PolygonIterator>
struct Report_inters{
  std::set<Polygon_2*>& m_intersecting_polygons;
  typedef My_Box<PolygonIterator> Box;
  
  Report_inters(std::set<Polygon_2*>& intersecting_polygons):m_intersecting_polygons(intersecting_polygons){}
  
  void operator() ( const Box& a, const Box& b) {
    if (a.it_poly==b.it_poly) return;
    
    if( CGAL::do_intersect(a.get_segment(),b.get_segment()) ){
      m_intersecting_polygons.insert(&(*a.it_poly));
      m_intersecting_polygons.insert(&(*b.it_poly));
    }
  }
};

template <class PolygonIterator>
void check_intersection_in_slice(PolygonIterator begin,PolygonIterator end,double elev,int cc,std::ofstream& ok,std::ofstream& bad){
  typedef My_Box<PolygonIterator> Box;
  std::set<Polygon_2*> intersecting_polygons;
  Report_inters<PolygonIterator> report_inters(intersecting_polygons);
  
  std::vector<Box> boxes;
  for (PolygonIterator it_poly=begin;it_poly!=end;++it_poly){
    const Point_2* pt_ptr=&( *CGAL::cpp0x::prev(it_poly->vertices_end()) );
    for (Polygon_2::Vertex_iterator itv=it_poly->vertices_begin();itv!=it_poly->vertices_end();++itv){
      boxes.push_back( Box(*pt_ptr,*itv,it_poly) );
      pt_ptr=&(*itv);
    }
  }
  
  CGAL::box_self_intersection_d( boxes.begin(), boxes.end(), report_inters);
  
  
  if (intersecting_polygons.empty())
    for (std::list<Polygon_2>::iterator it_poly=begin;it_poly!=end;++it_poly)
      write_polygon(it_poly->vertices_begin(),it_poly->vertices_end(),elev,cc,ok);
  else{
    std::cout << "Polygons intersecting" << std::endl;
    for (std::list<Polygon_2>::iterator it_poly=begin;it_poly!=end;++it_poly){
      if ( intersecting_polygons.find(&(*it_poly))!=intersecting_polygons.end() )
        write_polygon(it_poly->vertices_begin(),it_poly->vertices_end(),elev,cc,bad);
      else
        write_polygon(it_poly->vertices_begin(),it_poly->vertices_end(),elev,cc,ok);
    }
  }
}

int main(int ,char** argv){
  std::ifstream input(argv[1]);
  int cc=atoi(argv[2]);
  
  std::ofstream ok("ok.cgal");
  std::ofstream bad("bad.cgal");
    
  bool read_ok=true;
  std::list<Point_2> points;
  double elev;
  std::list<Polygon_2> slice;
  slice.push_back(Polygon_2());
  do{
    points.clear();
    read_ok=read_polygon(points,cc,elev,input);
  }
  while (read_ok && !set_polygon(slice.back(),points,elev,cc,bad));
  
  do{
    Polygon_2 poly;
    double poly_elev;
    
    do{
      points.clear();
      read_ok=read_polygon(points,cc,poly_elev,input);
    }
    while (read_ok && !set_polygon(poly,points,elev,cc,bad));
    
    if (!read_ok || poly_elev!=elev){
      //check for self-intersection in the slice
      check_intersection_in_slice(slice.begin(),slice.end(),elev,cc,ok,bad);

      slice.clear();
      if (read_ok){
        slice.push_back(poly);
        elev=poly_elev;
      }
    }
    else
      slice.push_back(poly);
  }while(read_ok);
}
