#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/IO/WKT.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>



typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;


struct FaceInfo {

  FaceInfo()
  {}

  int nesting_level[2];

  bool in_domain(int i) const
  {
    return nesting_level[i] % 2 == 1;
  }

  template <typename Fct>
  bool in_domain(const Fct& fct) const
  {
    return fct(in_domain(0), in_domain(1));
  }

};

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2                                         Point_2;

typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo,K>   Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>             Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,Tds,Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
typedef CDTplus::Face_handle                              Face_handle;
typedef CDTplus::Constraint_id                            Constraint_id;
typedef CDTplus::Edge                                     Edge;
typedef CDTplus::Context                                  Context;




void
mark_domains(CDTplus& ct,
             Face_handle start,
             int index,
             std::list<Edge>& border,
             const std::set<Constraint_id>& cids,
             int aorb)
{
  if(start->info().nesting_level[aorb] != -1){
    return;
  }
  std::list<Face_handle> queue;
  queue.push_back(start);

  while(! queue.empty()){
    Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level[aorb] == -1){
      fh->info().nesting_level[aorb] = index;
      for(int i = 0; i < 3; i++){
        Edge e(fh,i);
        Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level[aorb] == -1){
          if(ct.is_constrained(e)){
            for(Context c : ct.contexts(e.first->vertex(ct.cw(e.second)),
                                        e.first->vertex(ct.ccw(e.second)))){
              if(cids.find(c.id()) != cids.end()){
                border.push_back(e);
              }
            }
          }else{
            queue.push_back(n);
          }
        }
      }
    }
  }
}




void
mark_domains(CDTplus& cdt, const std::set<Constraint_id>& cids, int aorb)
{
  for(Face_handle f : cdt.all_face_handles()){
    f->info().nesting_level[aorb] = -1;
  }

  std::list<Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border, cids, aorb);
  while(! border.empty()){
    Edge e = border.front();
    border.pop_front();
    Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level[aorb] == -1){
      mark_domains(cdt, n, e.first->info().nesting_level[aorb]+1, border, cids, aorb);
    }
  }
}

template <typename Fct>
void to_stl(const CDTplus& cdt, const Fct& fct)
{
  std::ofstream out("cdt.stl");
  out.precision(17);
  out << "solid outside\n";
  for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit){
    if(fit->info().in_domain(fct)){
      out << "facet normal 0 0 0\n"
                       << "outer loop\n";
        for(int i=0; i < 3; ++i){
          out << "vertex " << fit->vertex(i)->point() << " 0\n";
        }
        out << "endloop\n"
            << "endfacet\n";
    }
  }
  out << "endsolid" << std::endl;
  out.close();
}



int
main( )
{
  CDTplus cdt;

  Polygon_with_holes_2 pA, pB;

  std::set<Constraint_id> idA, idB;

  {
    std::istringstream is("POLYGON ((0 0,  20 0, 20 30, 0 30, 0 0), (1 1, 1 2, 2 2, 2 1, 1 1 ) )");
    CGAL::IO::read_polygon_WKT(is, pA);
  }


  {
    std::istringstream is("POLYGON ((10 1,  30 1, 30 2, 20 2, 20 4, 10 4))");
    CGAL::IO::read_polygon_WKT(is, pB);
  }

  Constraint_id cidA = cdt.insert_constraint(pA.outer_boundary().vertices_begin(), pA.outer_boundary().vertices_end(), true);
  idA.insert(cidA);
  for(auto const& hole : pA.holes()){
    cidA = cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);
    idA.insert(cidA);
  }


  Constraint_id cidB = cdt.insert_constraint(pB.outer_boundary().vertices_begin(), pB.outer_boundary().vertices_end(), true);
  idB.insert(cidB);
  for(auto const& hole : pB.holes()){
    cidB = cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);
    idB.insert(cidB);
  }


  mark_domains(cdt, idA, 0);
  mark_domains(cdt, idB, 1);

  to_stl(cdt, [](bool a, bool b){ return a || b;} );
  return 0;
}
