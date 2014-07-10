#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/nearest_vertex.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <iostream>
#include <fstream>

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

void
mark_domains(const CDT& cdt);

void
vertices_to_add(const CDT& cdtA, const CDT& cdtB, std::list<Point_2>& points);

void
add_vertices_inside(CDT& cdt, const std::list<Point_2>& points);

void
update_z(CDT& cdt);



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






struct Layer {
  CDT cdt;
  double z;
  Tds_3 tds;
  typedef std::map<Vertex_handle, Vertex_handle_3> V2V3_map;
  V2V3_map v2v3_map;

  Layer()
  {
    tds.set_dimension(3);
  }

  Vertex_handle_3 v2v3(Vertex_handle v)
  {
    Vertex_handle_3 v3;
    V2V3_map::iterator vvi = v2v3_map.find(v);
    if(vvi != v2v3_map.end()){
      v3 = vvi->second;
    } else {
      v3 = tds.create_vertex();
      v2v3_map.insert(std::make_pair(v, v3));
    }
    v3->info().v = v;
    v3->point() = Point_3(v->point().x(), v->point().y(), v->info().z);
    return v3;
  }
  
};



void
create_tetrahedrization(Layer& top, Layer& bottom);

void remove_cells(Layer& layer);

struct Body {

  std::list<Layer> layers;
  std::vector<std::list<Layer>::iterator> indexed_layers;
  std::map<double, std::list<Layer>::iterator> z_layers;

  Layer& layer(double d)
  {
    std::map<double, std::list<Layer>::iterator>::iterator it = z_layers.find(d);
    if( it == z_layers.end()){
      Layer l;
      l.z = d;
      std::list<Layer>::iterator lit = layers.insert(layers.end(),l);
      it = z_layers.insert(std::make_pair(d, lit)).first;
    }
    return *(it->second);
  }

  int index_layers()
  {
    indexed_layers.clear();
    indexed_layers.reserve(z_layers.size());
    for(std::map<double, std::list<Layer>::iterator>::iterator it = z_layers.begin();
        it != z_layers.end();
        ++it){
      indexed_layers.push_back(it->second);
    }
    return indexed_layers.size();
  }
  

  Layer& operator[](int i)
  {
    return *(indexed_layers[i]);
  }

  void
  make_conforming_Gabriel_2() 
  {
    for(std::size_t i =0; i < indexed_layers.size(); i++){
      CGAL::make_conforming_Gabriel_2(indexed_layers[i]->cdt);
      mark_domains(indexed_layers[i]->cdt);
    }
  }

  void add_vertices()
  {
    std::vector<std::list<Point_2> > points(indexed_layers.size());
    for(std::size_t i =0; i < indexed_layers.size(); i++){
      if(i > 0){
        vertices_to_add(indexed_layers[i]->cdt, indexed_layers[i-1]->cdt, points[i-1]);
      }
      if(i <  indexed_layers.size()-1){
        vertices_to_add(indexed_layers[i]->cdt, indexed_layers[i+1]->cdt, points[i+1]);
      }
    }
    for(std::size_t i =0; i < indexed_layers.size(); i++){
      add_vertices_inside(indexed_layers[i]->cdt, points[i]);
    }
  }

  void update_z()
  {
    for(std::size_t i=0; i < indexed_layers.size(); i++){
      ::update_z(indexed_layers[i]->cdt);
    }
  }

  void
  create_tetrahedrization()
  {
    for(std::size_t i=1; i < indexed_layers.size(); i++){
      ::create_tetrahedrization(*indexed_layers[i], *indexed_layers[i-1]);
      // We call remove imediately afterwards as neighboring layers would use 'volume'
      remove_cells(*indexed_layers[i-1]);
    }
  }
  /*
  // index the vertices of all 2D triangulations
  void index_vertices()
  {
    int n = 0;
    for(int i=0; i < indexed_layers.size(); i++){
      n = index(indexed_layers[i]->cdt, n);
    }
  }
  */
  int number_of_vertices()
  {
    int vc = 0;
    for(std::size_t i=0; i < indexed_layers.size(); i++){
      CDT& cdt = indexed_layers[i]->cdt;
      vc += cdt.number_of_vertices();
    }
    return vc;
  }

  int number_of_faces()
  {
    int fc = 0;
    for(std::size_t i=0; i < indexed_layers.size()-1; i++){
      Tds_3& tds = indexed_layers[i]->tds;
      for(Cell_iterator_3 it = tds.cells_begin(); it != tds.cells_end(); ++it){
        int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
        for(int i=0; i < n; i++){
          if(it->neighbor(i) == Cell_handle_3()){
            fc++;
          }
        }
      } 
    }
    return fc;
  }

  void write_surface(std::ostream& os)
  {
    os << "OFF\n" << number_of_vertices() << " " << number_of_faces() << " 0\n";

    for(std::size_t i=0; i < indexed_layers.size(); i++){
      CDT& cdt = indexed_layers[i]->cdt;
      int vc = 0;
      for(Finite_vertices_iterator it = cdt.vertices_begin(); it != cdt.vertices_end(); ++it){
        it->info().index = vc++;
        os << it->point() << std::endl;
      }
    }
  
    for(std::size_t i=0; i < indexed_layers.size()-1; i++){
      Tds_3 tds = indexed_layers[i]->tds;
      for(Cell_iterator_3 it = tds.cells_begin(); it != tds.cells_end(); ++it){    
        int n = (it->info().type == CellInfo3::EDGE_EDGE) ? 4 : 3;
        for(int i=0; i < n; i++){
          if(it->neighbor(i) == Cell_handle_3()){
            os << "3 " 
               << it->vertex((i+1)%4)->info().v->info().index <<  " "  
               << it->vertex((i+2)%4)->info().v->info().index <<  " "  
               << it->vertex((i+3)%4)->info().v->info().index <<  std::endl;
          }
        }
      }
    }
  }

};



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




Point_2 project(const Point_3& p)
{
  return Point_2(p.x(), p.y());
}


void
read(Body& body, std::istream& is)
{
  int n;
  Point_3 p3;
  Point_2 p2;
  Vertex_handle vp, vq;
  while(is >> n){
    is >> p3;
    Layer& layer = body.layer(p3.z());   
    p2 = project(p3);
    vp = layer.cdt.insert(p2);
    vp->info().z = p3.z();
    for(int i=1; i < n; i++){
      is >> p3;
      p2 = project(p3);
      vq = layer.cdt.insert(p2, vp->face());
      vq->info().z = p3.z();
      layer.cdt.insert_constraint(vp, vq);
      vp = vq;
    }
  }
}



  

// Index the vertices for writing polyline files
int
index(const CDT& cdt, int i=0 )
{
  for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
    it->info().index = i++;
  }
  return i;
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





void
create_cells(Layer& layerA, Layer& layerB, CellInfo3::Type type)
{
  CDT& cdtA = layerA.cdt;
  CDT& cdtB = layerB.cdt;
  Tds_3& tds = layerA.tds;
  for(Finite_faces_iterator it = cdtA.finite_faces_begin();it != cdtA.finite_faces_end(); ++it){
    Point_2 center = cdtA.circumcenter(it);
    Vertex_handle vh = nearest_vertex(cdtB, center, Face_handle());

    Vertex_handle_3 v0, v1, v2, v3;
      v0 = layerA.v2v3(it->vertex(0)); 
      v1 = layerA.v2v3(it->vertex(1));
      v2 = layerA.v2v3(it->vertex(2));
      v3 = layerA.v2v3(vh);   
    if(type == CellInfo3::TOP_TRIANGLE){
      std::swap(v1,v2);
    }

    assert(CGAL::orientation(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point()) == CGAL::POSITIVE);
    assert(CGAL::orientation(v0->point(), v1->point(), v2->point(), v3->point()) == CGAL::POSITIVE);
    Cell_handle_3 ch = tds.create_cell(v0, v1, v2, v3, Cell_handle_3(), Cell_handle_3(), Cell_handle_3(), Cell_handle_3());
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
create_cells_22(Layer& topLayer, Layer& bottomLayer)
{
  CDT& top = topLayer.cdt;
  CDT& bottom = bottomLayer.cdt;

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
         if(bottom_voronoi_segment)
           intersect = CGAL::do_intersect(ts, bs);
         else
           intersect = CGAL::do_intersect(ts, br);
       else 
         if(bottom_voronoi_segment)
           intersect = CGAL::do_intersect(tr, bs);
         else
           intersect = CGAL::do_intersect(tr, br);

       if(intersect){
         //std::cout << "they intersect" << std::endl;
       
         Vertex_handle_3 v0, v1, v2, v3;
         v0 = bottomLayer.v2v3(tit->first->vertex(CDT::cw(tit->second))); 
         v1 = bottomLayer.v2v3(tit->first->vertex(CDT::ccw(tit->second))); 
         v2 = bottomLayer.v2v3(bit->first->vertex(CDT::cw(bit->second))); 
         v3 = bottomLayer.v2v3(bit->first->vertex(CDT::ccw(bit->second))); 
         if(CGAL::orientation(v0->point(), v1->point(), v2->point(), v3->point()) != CGAL::POSITIVE){        
           std::swap(v2, v3);
         }
         Cell_handle_3 ch = bottomLayer.tds.create_cell(v0, v1, v2, v3, Cell_handle_3(), Cell_handle_3(), Cell_handle_3(), Cell_handle_3());
         ch->info().type = CellInfo3::EDGE_EDGE;
         ch->info().f0 = tit->first;
         ch->info().i0 = tit->second;
         ch->info().f1 = bit->first;
         ch->info().i1 = bit->second;
       }
     }
   }
}


void
connect_cells(Layer& topLayer, Layer& bottomLayer)
{
  std::cerr << "connect cells " << topLayer.cdt.vertices_begin()->info().z << "  " << bottomLayer.cdt.vertices_begin()->info().z << std::endl;
  CDT& top = topLayer.cdt;
  CDT& bottom = bottomLayer.cdt;
  VertexTripleCellMap vtcm;
  for(Cell_iterator_3 it = bottomLayer.tds.cells_begin(); it != bottomLayer.tds.cells_end(); ++it){
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
is_valid(const Layer& layer)
{
  for(Cell_iterator_3 it = layer.tds.cells_begin(); it != layer.tds.cells_end(); ++it){
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
is_valid()
{
  
}


void
create_tetrahedrization(Layer& top, Layer& bottom)
{
  create_cells(top, bottom, CellInfo3::TOP_TRIANGLE);
  create_cells(bottom, top, CellInfo3::BOTTOM_TRIANGLE);

  create_cells_22(top, bottom);

  connect_cells(top, bottom);

  is_valid();
}


void remove_cell(Layer& layer, Cell_handle_3 it)
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
  layer.tds.delete_cell(it);
}


void remove_cells(Layer& layer)
{
  // First pass: Remove cells with an edge or face outside the domain 
  for(Cell_iterator_3 it = layer.tds.cells_begin(); it != layer.tds.cells_end(); ){
    // initialize the connected component field for the third pass
    it->info().cc = false;
    Cell_iterator_3 next = it;
    ++next;
    if(it->info().type == CellInfo3::EDGE_EDGE){
      if( ( (! it->info().f0->info().in_domain) && (! it->info().f0->neighbor(it->info().i0)->info().in_domain) ) || 
          ( (! it->info().f1->info().in_domain) && (! it->info().f1->neighbor(it->info().i1)->info().in_domain) ) ){
        remove_cell(layer, it);
      }
    } else {
      if(! it->info().f0->info().in_domain){
        remove_cell(layer, it);
      }
    }
    it = next;
  }

  // Second pass: Remove t12 cells, which are not 1-solid or not 2-solid, or both
  for(Cell_iterator_3 it = layer.tds.cells_begin(); it != layer.tds.cells_end(); ){
    Cell_iterator_3 next = it;
    ++next;
    if(it->info().type == CellInfo3::EDGE_EDGE){
      if(!  ( (it->info().f0->info().volume || it->info().f0->neighbor(it->info().i0)->info().volume)  && 
              (it->info().f1->info().volume || it->info().f1->neighbor(it->info().i1)->info().volume) ) ) {
        remove_cell(layer, it);
      }
    }
    it = next;
  }

  // Third pass: Remove connected components which do not have at least one t12 cell
  std::list<Cell_handle_3> cells;
  for(Cell_iterator_3 it = layer.tds.cells_begin(); it != layer.tds.cells_end(); ++it ) {
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
    remove_cell(layer, *it);
  }
}

#if 0

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
#endif 

#if 0
void bottom_top()
{
  tds_3.set_dimension(3);
  CDT bottom, top;
  
  read(std::ifstream("top.cgal"), top);
  read(std::ifstream("bottom.cgal"), bottom);

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

  remove_cells();

  index(top);
  index(bottom);
  
  write_layer(std::ofstream("bottom.off"),bottom);
  write_layer(std::ofstream("top.off"), top);

  write_voronoi(std::ofstream("bottomVoronoi.cgal"),bottom);
  write_voronoi(std::ofstream("topVoronoi.cgal"), top);

  write_cells(std::ofstream("graph.off"));

}
#endif

int main()
{
  Body body;
  std::ifstream input("z-const.cgal");
  read(body,input);
  body.index_layers();
  body.make_conforming_Gabriel_2();
  body.add_vertices();
  body.make_conforming_Gabriel_2();
  body.update_z();
  body.create_tetrahedrization();
  std::ofstream output("body.off");
  body.write_surface(output);

  std::cerr << "done" << std::endl;
  return 0;
}
