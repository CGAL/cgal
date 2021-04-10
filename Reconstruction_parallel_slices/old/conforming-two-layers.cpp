#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/nearest_vertex.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Random.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <iostream>
#include <fstream>

//CGAL::Random random;

struct VertexInfo2
{
  int index;
  double z;
};

struct FaceInfo2
{
  bool in_domain;
  bool volume;
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K> Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>      TDS;

typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS> CDT;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

typedef K::Point_3 Point_3;


typedef CDT::Locate_type Locate_type;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef CDT::All_faces_iterator All_faces_iterator;
typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
typedef CDT::Finite_edges_iterator Finite_edges_iterator;
typedef CDT::Finite_faces_iterator Finite_faces_iterator;


typedef boost::tuple<Vertex_handle, Vertex_handle, Vertex_handle> Facet;

Facet
make_canonical(Facet facet)
{
  if(&*facet.get<0>()<&*facet.get<1>()) std::swap(facet.get<0>(), facet.get<1>());
  if(&*facet.get<1>()<&*facet.get<2>()) std::swap(facet.get<1>(), facet.get<2>());
  if(&*facet.get<0>()<&*facet.get<1>()) std::swap(facet.get<0>(), facet.get<1>());
  return facet;
}


// For triangle-vertex cells the type is TOP_TRIANGLE or BOTTOM_TRIANGLE
// and the cell is defined by f0 and vh 
// For edge-edge cells the cell is defined by Edge(f0,i0) and Edge(f1,i1)
struct CellInfo3 {
  enum Type { TOP_TRIANGLE, BOTTOM_TRIANGLE, EDGE_EDGE };
  Type type;
  
  int index;
  Face_handle f0, f1;
  int i0, i1;
  Vertex_handle vh;
  bool cc;

};

struct VertexInfo3
{
  Vertex_handle v;
  int index;
};



typedef CGAL::Triangulation_data_structure_3<CGAL::Triangulation_vertex_base_with_info_3<VertexInfo3,K>, CGAL::Triangulation_cell_base_with_info_3<CellInfo3, K> > Tds_3;

typedef Tds_3::Vertex_handle Vertex_handle_3;
typedef Tds_3::Cell_handle Cell_handle_3;
typedef Tds_3::Facet Facet_3;
typedef Tds_3::Cell_iterator Cell_iterator_3;
typedef Tds_3::Vertex_iterator Vertex_iterator_3;


Tds_3 tds_3;

typedef std::map<Vertex_handle, Vertex_handle_3> V2V3_map;



typedef boost::tuple<Vertex_handle_3, Vertex_handle_3, Vertex_handle_3> Vertex_3_tuple;

Vertex_3_tuple
make_canonical(Vertex_3_tuple facet)
{
  if(&*facet.get<0>()<&*facet.get<1>()) std::swap(facet.get<0>(), facet.get<1>());
  if(&*facet.get<1>()<&*facet.get<2>()) std::swap(facet.get<1>(), facet.get<2>());
  if(&*facet.get<0>()<&*facet.get<1>()) std::swap(facet.get<0>(), facet.get<1>());
  return facet;
}

Vertex_3_tuple
make_vertex_3_tuple(const Facet_3& f)
{
  return make_canonical(boost::make_tuple(f.first->vertex((f.second+1)%4),
                                          f.first->vertex((f.second+2)%4),
                                          f.first->vertex((f.second+3)%4)));
}
    

typedef std::map<Vertex_3_tuple, Facet_3> VertexTripleCellMap;

V2V3_map v2v3_map;


Point_2 project(const Point_3& p)
{
  return Point_2(p.x(), p.y());
}


void
read(std::istream& is, CDT& cdt)
{
  int n;
  Point_3 p3;
  Point_2 p2;
  Vertex_handle vp, vq;
  while(is >> n){
    is >> p3;
    p2 = project(p3);
    vp = cdt.insert(p2);
    vp->info().z = p3.z();
    for(int i=1; i < n; i++){
      is >> p3;
      p2 = project(p3);
      vq = cdt.insert(p2, vp->face());
      vq->info().z = p3.z();
      cdt.insert_constraint(vp, vq);
      vp = vq;
    }
  }
}

// Index the vertices for writing polyline files
void
index(const CDT& cdt)
{
  int i = 0;
  for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
    it->info().index = i++;
  }
}


// We currently only support non-nested polygons
void
mark_domains(const CDT& cdt)
{
  for(All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().in_domain = true;
  }
  std::list<Face_handle> faces;
  faces.push_back(cdt.infinite_vertex()->face());
  while(! faces.empty()){
    Face_handle fh = faces.front();
    faces.pop_front();
    if(fh->info().in_domain==false){
      continue;
    }
    fh->info().in_domain = false;
    for(int i=0; i < 3; i++){
      if((! fh->is_constrained(i)) && (fh->neighbor(i)->info().in_domain)){
        faces.push_front(fh->neighbor(i));
      }
    }
  }
}



void
vertices_to_add(const CDT& cdtA, const CDT& cdtB, std::list<Point_2>& points) 
{
  for(Finite_faces_iterator it = cdtA.finite_faces_begin(); it != cdtA.finite_faces_end(); ++it){
    if(! it->info().in_domain){ 
      Point_2 center = cdtA.circumcenter(it);
      Locate_type lt;
      int li;
      Face_handle loc = cdtB.locate(center, lt, li);
      if( (lt == CDT::FACE) && (loc->info().in_domain)){
        points.push_back(center);
      } 
    }
  }
}

void
add_vertices_inside(CDT& cdt, const std::list<Point_2>& points)
{
  double z = cdt.finite_vertices_begin()->info().z;
  for(std::list<Point_2>::const_iterator it = points.begin(); it != points.end(); ++it){
    Vertex_handle vh = cdt.insert(*it);
    vh->info().z = z;
  }
}

// After refining a CDT we have to set the z values for the new vertices 
void
update_z(CDT& cdt)
{
  double z = cdt.finite_vertices_begin()->info().z;
  for(Finite_vertices_iterator it = cdt.finite_vertices_begin();  it != cdt.finite_vertices_end(); ++it){
    it->info().z = z;
  }
}


Vertex_handle_3 v2v3(Vertex_handle v)
{
  Vertex_handle_3 v3;
  V2V3_map::iterator vvi = v2v3_map.find(v);
  if(vvi != v2v3_map.end()){
    v3 = vvi->second;
  } else {
    v3 = tds_3.create_vertex();
    v2v3_map.insert(std::make_pair(v, v3));
  }
  v3->info().v = v;
  v3->point() = Point_3(v->point().x(), v->point().y(), v->info().z);
  return v3;
}


void
create_cells(const CDT& cdtA, const CDT& cdtB, CellInfo3::Type type)
{
  for(Finite_faces_iterator it = cdtA.finite_faces_begin();it != cdtA.finite_faces_end(); ++it){
    Point_2 center = cdtA.circumcenter(it);
    Vertex_handle vh = nearest_vertex(cdtB, center, Face_handle());

    Vertex_handle_3 v0, v1, v2, v3;
      v0 = v2v3(it->vertex(0)); 
      v1 = v2v3(it->vertex(1));
      v2 = v2v3(it->vertex(2));
      v3 = v2v3(vh);   
    if(type == CellInfo3::TOP_TRIANGLE){
      std::swap(v1,v2);
    }
    assert(CGAL::orientation(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point()) == CGAL::POSITIVE);
    assert(CGAL::orientation(v0->point(), v1->point(), v2->point(), v3->point()) == CGAL::POSITIVE);
    Cell_handle_3 ch = tds_3.create_cell(v0, v1, v2, v3, Cell_handle_3(), Cell_handle_3(), Cell_handle_3(), Cell_handle_3());
    ch->info().type = type;
    ch->info().f0 = it;
    it->info().volume = true;

  }
}


bool dual(const CDT& cdt, CDT::Edge e, Segment_2& s, Ray_2& r)
{
  bool is_voronoi_segment = true;
  Face_handle fh, nh;
  int fi, ni;
  boost::tie(fh, fi) = e;
  nh = fh->neighbor(fi);
  ni = nh->index(fh);
  if( (! cdt.is_infinite(fh))
      && (! cdt.is_infinite(nh))){
    s = Segment_2(cdt.circumcenter(fh), cdt.circumcenter(nh));
  } else {
    is_voronoi_segment = false;
    if(! cdt.is_infinite(fh)){
      const Point_2& p = fh->vertex(CDT::cw(fi))->point();
      const Point_2& q = fh->vertex(CDT::ccw(fi))->point();
      Line_2 l = CGAL::bisector(p,q);
      r = Ray_2(cdt.circumcenter(fh), l);
    } else {
      const Point_2& p = nh->vertex(CDT::cw(ni))->point();
      const Point_2& q = nh->vertex(CDT::ccw(ni))->point();
      Line_2 l = CGAL::bisector(p,q);
      r = Ray_2(cdt.circumcenter(nh), l);
    }
  }
  return is_voronoi_segment;
}


// TODO:: Replace this quadratic algorithm
void
create_cells_22(const CDT& top, const CDT& bottom)
{
   for(Finite_edges_iterator tit = top.finite_edges_begin(); tit != top.finite_edges_end(); ++tit){
     Segment_2 ts;
     Ray_2 tr;
     bool top_voronoi_segment = dual(top, *tit, ts, tr);
     
     for(Finite_edges_iterator bit = bottom.finite_edges_begin(); bit != bottom.finite_edges_end(); ++bit){
       Segment_2 bs;
       Ray_2 br;
       bool bottom_voronoi_segment = dual(bottom, *bit, bs, br);
       
       bool intersect = false;
       if(top_voronoi_segment)
         if(bottom_voronoi_segment){
           // std::cerr << "tsbs" << std::endl;
           if((ts == bs) || (ts.opposite() == bs)){
             std::cerr << "ts == bs" << std::endl;
           }
           intersect = CGAL::do_intersect(ts, bs);
         }else{
           // std::cerr << "tsbr" << std::endl;
           intersect = CGAL::do_intersect(ts, br);
         }else 
         if(bottom_voronoi_segment){
           // std::cerr << "trbs" << std::endl;
           intersect = CGAL::do_intersect(tr, bs);
         }else{
           // std::cerr << "trbr" << std::endl;
           intersect = CGAL::do_intersect(tr, br);
         }
       if(intersect){
         //std::cout << "they intersect" << std::endl;
       
         Vertex_handle_3 v0, v1, v2, v3;
         v0 = v2v3(tit->first->vertex(CDT::cw(tit->second))); 
         v1 = v2v3(tit->first->vertex(CDT::ccw(tit->second))); 
         v2 = v2v3(bit->first->vertex(CDT::cw(bit->second))); 
         v3 = v2v3(bit->first->vertex(CDT::ccw(bit->second))); 
         if(CGAL::orientation(v0->point(), v1->point(), v2->point(), v3->point()) != CGAL::POSITIVE){        
           std::swap(v2, v3);
         }
         if(CGAL::orientation(v0->point(), v1->point(), v2->point(), v3->point()) == CGAL::POSITIVE){
           Cell_handle_3 ch = tds_3.create_cell(v0, v1, v2, v3, Cell_handle_3(), Cell_handle_3(), Cell_handle_3(), Cell_handle_3());
           ch->info().type = CellInfo3::EDGE_EDGE;
           ch->info().f0 = tit->first;
           ch->info().i0 = tit->second;
           ch->info().f1 = bit->first;
           ch->info().i1 = bit->second;
         }
       }
     }
   }
}


void
connect_cells(const CDT& top, const CDT& bottom)
{
  VertexTripleCellMap vtcm;
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ++it){
    // Cells don't have neighbors on top or bottom
    int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
    for(int i=0; i < n; i++){
      Vertex_3_tuple t = make_vertex_3_tuple(Facet_3(it,i));
      VertexTripleCellMap::iterator vtcmi = vtcm.find(t);
      if(vtcmi != vtcm.end()){
        Facet_3 nf = vtcmi->second;
        assert(nf.first->neighbor(nf.second) == Cell_handle_3());
        nf.first->set_neighbor(nf.second, it);
        it->set_neighbor(i, nf.first);
      } else {
        vtcm.insert(std::make_pair(t,Facet_3(it,i)));
      }
    }
  }
}

void
is_valid()
{
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ++it){
    int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
    for(int i=0; i < 4; i++){
      Cell_handle_3 nh = it->neighbor(i);
      if(nh != Cell_handle_3()){
        for(int j = 0; j < 4; j++){
          if(j != i){
            Cell_handle_3 ch = it->neighbor(j);
            if(ch != Cell_handle_3()){
              assert(ch != nh);
            }
          }
        }
        int ni = nh->index(it);
        assert(nh->neighbor(ni) == it);
      }
    }
  }
}

void
create_tetrahedrization(const CDT& top, const CDT& bottom)
{
  create_cells(top, bottom, CellInfo3::TOP_TRIANGLE);
  create_cells(bottom, top, CellInfo3::BOTTOM_TRIANGLE);

    create_cells_22(top, bottom);

  connect_cells(top, bottom);

  is_valid();
}


void remove_cell(Cell_handle_3 it)
{
  for(int i = 0; i<4; i++){
    Cell_handle_3 nh = it->neighbor(i);
    if(nh != Cell_handle_3()){
      int ni = nh->index(it);
      nh->set_neighbor(ni, Cell_handle_3());
    }
  }
  if(it->info().type != CellInfo3::EDGE_EDGE){
    it->info().f0->info().volume = false;
  }
  tds_3.delete_cell(it);
}


void remove_cells()
{
  // First pass: Remove cells with an edge or face outside the domain 
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ){
    // initialize the connected component field for the third pass
    it->info().cc = false;
    Cell_iterator_3 next = it;
    ++next;
    if(it->info().type == CellInfo3::EDGE_EDGE){
      if( ( (! it->info().f0->info().in_domain) && (! it->info().f0->neighbor(it->info().i0)->info().in_domain) ) || 
          ( (! it->info().f1->info().in_domain) && (! it->info().f1->neighbor(it->info().i1)->info().in_domain) ) ){
        remove_cell(it);
      }
    } else {
      if(! it->info().f0->info().in_domain){
        remove_cell(it);
      }
    }
    it = next;
  }

  // Second pass: Remove t12 cells, which are not 1-solid or not 2-solid, or both
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ){
    Cell_iterator_3 next = it;
    ++next;
    if(it->info().type == CellInfo3::EDGE_EDGE){
      if(!  ( (it->info().f0->info().volume || it->info().f0->neighbor(it->info().i0)->info().volume)  && 
              (it->info().f1->info().volume || it->info().f1->neighbor(it->info().i1)->info().volume) ) ) {
        remove_cell(it);
      }
    }
    it = next;
  }

  // Third pass: Remove connected components which do not have at least one t12 cell
  std::list<Cell_handle_3> cells;
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ++it ) {
    if(it->info().type != CellInfo3::EDGE_EDGE){
      if(! it->info().cc){
        // explore the component
        bool found_edge_edge_cell = false;
        std::list<Cell_handle_3> component;
        std::list<Cell_handle_3> queue;
        it->info().cc = true;
        component.push_back(it);
        queue.push_back(it);
        while(! queue.empty()){
          Cell_handle_3 ch = queue.front();
          queue.pop_front();
          for(int i=0; i < 4; i++){
            Cell_handle_3 nh = ch->neighbor(i);
            if(nh != Cell_handle_3()){
              if(! nh->info().cc){
                if(nh->info().type == CellInfo3::EDGE_EDGE){
                  found_edge_edge_cell = true;
                }
                nh->info().cc = true;
                component.push_back(nh);
                queue.push_back(nh);
              }
            }
          }
        }
        if(! found_edge_edge_cell){
          std::copy(component.begin(), component.end(), std::back_inserter(cells));
        }
      }
    }
  }
  for(std::list<Cell_handle_3>::iterator it = cells.begin(); it != cells.end(); ++it){
    remove_cell(*it);
  }
}


void write_cells(std::ostream& os)
{
  int fc = 0;
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ++it){
    int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
    for(int i=0; i < n; i++){
      if(it->neighbor(i) == Cell_handle_3()){
        fc++;
      }
    }
  }
  os << "OFF\n" << tds_3.number_of_vertices() << " " << fc << " 0\n";

  int vc = 0;
  for(Vertex_iterator_3 it = tds_3.vertices_begin(); it != tds_3.vertices_end(); ++it){
    it->info().index = vc++;
    os << it->point() << std::endl;
  }
  
  for(Cell_iterator_3 it = tds_3.cells_begin(); it != tds_3.cells_end(); ++it){    
    int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
    for(int i=0; i < n; i++){
      if(it->neighbor(i) == Cell_handle_3()){
        os << "3 " 
           << it->vertex((i+1)%4)->info().index <<  " "  
           << it->vertex((i+2)%4)->info().index <<  " "  
           << it->vertex((i+3)%4)->info().index <<  std::endl;
      }
    }
  }
}

void write_layer(std::ostream& os, const CDT& cdt)
{
  double z =  cdt.finite_vertices_begin()->info().z;
  int number_of_domain_faces = 0;
  for(All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    if(it->info().in_domain){
      ++number_of_domain_faces;
    }
  }
  os << "OFF\n" << cdt.number_of_vertices() << " " << number_of_domain_faces << " 0\n";
  for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
    os << it->point() << " " << z << std::endl;
  }
  for(All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    if(it->info().in_domain){
       os << "3 " << it->vertex(0)->info().index << " "  << it->vertex(1)->info().index << " "  << it->vertex(2)->info().index << std::endl;
    }
  }
}

void write_voronoi(std::ostream& os, const CDT& cdt)
{
  double z =  cdt.finite_vertices_begin()->info().z;
  for(Finite_edges_iterator it = cdt.finite_edges_begin(); it != cdt.finite_edges_end(); ++it){
    Face_handle fh, nh;
    int fi;
    boost::tie(fh,fi) = *it;
    nh = fh->neighbor(fi);
    
    if(fh->info().in_domain && nh->info().in_domain){
      os << "2 " << cdt.circumcenter(fh) << " " << z+1 <<  " "  << cdt.circumcenter(nh) << " " << z+1 << std::endl;;  
    }
  }  
}


int main()
{
  tds_3.set_dimension(3);
  CDT bottom, top;
  
  std::ifstream input_1("top.cgal");
  std::ifstream input_2("bottom.cgal");
  read(input_1, top);
  read(input_2, bottom);


  CGAL::make_conforming_Gabriel_2(top);
  CGAL::make_conforming_Gabriel_2(bottom);

  mark_domains(top);
  mark_domains(bottom);

  std::list<Point_2> add_to_top, add_to_bottom;
  vertices_to_add(top, bottom, add_to_bottom);
  vertices_to_add(bottom, top, add_to_top);

  add_vertices_inside(bottom, add_to_bottom);
  add_vertices_inside(top, add_to_top);

  CGAL::make_conforming_Gabriel_2(top);
  CGAL::make_conforming_Gabriel_2(bottom);

  mark_domains(top);
  mark_domains(bottom);

  update_z(top);
  update_z(bottom);

  create_tetrahedrization(top, bottom);

  mark_domains(top);
  mark_domains(bottom);

  remove_cells();
  mark_domains(top);
  mark_domains(bottom);
  
  index(top);
  index(bottom);
  
  std::ofstream output_1("data/bottom.off");
  std::ofstream output_2("data/top.off");
  write_layer(output_1,bottom);
  write_layer(output_2, top);

  output_1.close();
  output_2.close();
  output_1.open("bottomVoronoi.cgal");
  output_2.open("topVoronoi.cgal");
  
  write_voronoi(output_1,bottom);
  write_voronoi(output_2, top);

  std::ofstream output("graph.off");
  write_cells(output);

  std::cerr << "done" << std::endl;
  return 0;
}
