
#include <CGAL/circulator.h>

#include "SM_common.h"

// trivial circulator testing
struct test_circ_size : public Surface_fixture
{
  void operator()()
  {
    Sm::Vertex_around_target_circulator vvc(m.halfedge(v),m);
    Sm::Halfedge_around_target_circulator hvc(m.halfedge(v),m);
    Sm::Face_around_target_circulator fvc(m.halfedge(v),m);
    Sm::Vertex_around_face_circulator vaf(m.halfedge(f1),m);
    Sm::Halfedge_around_face_circulator hfc(m.halfedge(f1),m);
    
    assert(CGAL::circulator_size(vvc) == 4);
    assert(CGAL::circulator_size(hvc) == 4);
    assert(CGAL::circulator_size(fvc) == 4);
    assert(CGAL::circulator_size(vaf) == 3);
    assert(CGAL::circulator_size(hfc) == 3);
 }
};



struct test_circ_distance : public Surface_fixture
{
  void operator()()
  {
    // the distance against end
    Sm::Vertex_around_target_circulator vvc(m.halfedge(v),m), vvce(vvc); 
    assert(CGAL::circulator_distance(vvc, vvce) == 4);

    vvc =  vvce = Sm::Vertex_around_target_circulator(m.halfedge(u),m);
    assert(CGAL::circulator_distance(vvc, vvce) == 2);

    vvc = vvce = Sm::Vertex_around_target_circulator(m.halfedge(w),m);
    assert(CGAL::circulator_distance(vvc, vvce) == 3);

    vvc = vvce = Sm::Vertex_around_target_circulator(m.halfedge(x),m);
    assert(CGAL::circulator_distance(vvc, vvce) == 3);

    Sm::Halfedge_around_target_circulator hvc(m.halfedge(v),m), hvce(hvc); 
    assert(CGAL::circulator_distance(hvc, hvce) == 4);
  
    Sm::Face_around_target_circulator fvc(m.halfedge(v),m), fvce(fvc);
    assert(CGAL::circulator_distance(fvc, fvce) == 4);

    Sm::Vertex_around_face_circulator vfc(m.halfedge(f1),m), vfce(vfc);
    assert(CGAL::circulator_distance(vfc, vfce) == 3);

    Sm::Halfedge_around_face_circulator hfc(m.halfedge(f1),m), hfce(hfc); 
    assert(CGAL::circulator_distance(hfc, hfce) == 3);
  }
};

struct test_emptiness : public Surface_fixture
{
  void operator()()
  {
    // add an isolated vertex
    Sm::Vertex_index iv = m.add_vertex(Point_3(2,2,0));
    assert(m.is_isolated(iv));

    Sm::Vertex_around_target_range vr = m.vertices_around_target(m.halfedge(iv));
    assert(is_empty_range(boost::begin(vr), boost::end(vr)));

    Sm::Face_around_target_range fr = m.faces_around_target(m.halfedge(iv));
    assert(is_empty_range(boost::begin(fr), boost::end(fr)));

    Sm::Halfedge_around_target_range hr = m.halfedges_around_target(m.halfedge(iv));
    assert(is_empty_range(boost::begin(hr), boost::end(hr)));
    // not true for everything else
    m.remove_vertex(iv);
    assert(m.is_removed(iv));
    Sm::Vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = m.vertices(); vb != ve; ++vb) {
      Sm::Vertex_around_target_range vr = m.vertices_around_target(m.halfedge(*vb));
      assert(!is_empty_range(boost::begin(vr), boost::end(vr)));
    }
  }
};

int main()
{
  test_circ_size tcs;
  tcs();
  test_circ_distance tcd;
  tcd();
  test_emptiness te;
  te();
  std::cerr << "done"<< std::endl;
  return 0;
}
