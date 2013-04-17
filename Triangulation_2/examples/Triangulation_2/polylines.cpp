#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>

#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel        K;
typedef CGAL::Polygon_2<K>                                       Polygon_2;
typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>             TDS;
typedef CGAL::Exact_intersections_tag                            Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag>  CDT;
typedef CGAL::Polyline_constrained_triangulation_2<CDT>          PCT;

typedef PCT::Point                                               Point;
typedef PCT::Constraint_id                                       Cid;
typedef PCT::Vertex_handle                                       Vertex_handle;

void 
print(const PCT& cdt, Cid cid)
{
  typedef PCT::Vertices_in_constraint                              Vertices_in_constraint;

  std::cout << "Polyline constraint:" << std::endl;
  for(Vertices_in_constraint it = cdt.vertices_in_constraint_begin(cid);
      it != cdt.vertices_in_constraint_end(cid);
      it++){
    Vertex_handle vh = *it;
    std::cout << vh->point() << std::endl;
  }
}


void 
contexts(const PCT& pct)
{
  PCT::Subconstraint_iterator beg, end;
  beg = pct.subconstraints_begin();
  end = pct.subconstraints_end();

  for(; beg!=end; ++beg){
    Vertex_handle vp = beg->first, vq = beg->second;

    if(pct.number_of_enclosing_constraints(vp, vq) == 2){
      PCT::Context_iterator cbeg = pct.contexts_begin(vp,vq);
      PCT::Context_iterator cend = pct.contexts_end(vp,vq);
      std::cout << "subconstraint " << vp->point() << " " << vq->point() 
                << " is on constraints starting at:\n";
      for(; cbeg !=  cend; ++cbeg){
        PCT::Context c = *cbeg;
        std::cout << (*(c.vertices_begin()))->point() << std::endl;
      }
    }
  }
}

int
main( )
{
  PCT pct;

  pct.insert_constraint(Point(0,0), Point(1,1));

  std::vector<Point> points;
  points.push_back(Point(1,1));
  points.push_back(Point(5,2));
  points.push_back(Point(6,0));
  points.push_back(Point(3,0));
  Cid id1 = pct.insert_constraint(points.begin(), points.end());

  print(pct, id1);

  Polygon_2 poly;
  poly.push_back(Point(2,3));
  poly.push_back(Point(4,0));
  poly.push_back(Point(5,0));
  poly.push_back(Point(6,2));

  Cid id2 = pct.insert_constraint(poly);

  print(pct, id1);
  print(pct, id2);

  contexts(pct);

  return 0;
}
