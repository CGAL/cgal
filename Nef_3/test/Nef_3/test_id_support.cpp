#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/ID_support_handler.h>
#include <cassert>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::SNC_indexed_items Items;
typedef CGAL::SNC_structure<K, Items, bool> SNC_structure;
typedef CGAL::SNC_decorator<SNC_structure> SNC_decorator;
typedef CGAL::SM_decorator<SNC_structure::Sphere_map> SM_decorator;
typedef CGAL::ID_support_handler<Items, SNC_decorator> ID_support_handler;

typedef SNC_structure::Vertex_handle Vertex_handle;
typedef SNC_structure::SHalfedge_handle SHalfedge_handle;
typedef SNC_structure::SHalfloop_handle SHalfloop_handle;
typedef SNC_structure::Sphere_circle Sphere_circle;

void run_test_case(int fw_id1, int bw_id1, int fw_id2, int bw_id2, bool test_loop) {
    ID_support_handler handler;

    handler.initialize_hash(fw_id1);
    handler.initialize_hash(bw_id1);
    handler.initialize_hash(fw_id2);
    handler.initialize_hash(bw_id2);

    SNC_structure snc;
    Vertex_handle v = snc.new_vertex();
    SM_decorator SD(&*v);

    SHalfedge_handle se = SD.new_shalfedge_pair();
    se->circle() = Sphere_circle(1,0,0);
    se->twin()->circle() = se->circle().opposite();

    SHalfedge_handle se1 = SD.new_shalfedge_pair();
    se1->set_forward_index(fw_id1);
    se1->set_backward_index(bw_id1);
    se1->twin()->set_forward_index(bw_id1);
    se1->twin()->set_backward_index(fw_id1);
    se1->circle() = se->circle();
    se1->twin()->circle() = se->circle().opposite();

    if (test_loop) {

        SHalfloop_handle sl2 = SD.new_shalfloop_pair();
        sl2->circle() = se->circle();
        sl2->twin()->circle() = se->circle().opposite();
        sl2->set_index(fw_id2);
        sl2->twin()->set_index(bw_id2);

        handler.handle_support(se, se1, sl2);

        // fw_id1 is merged with fw_id2
        assert(handler.get_hash(fw_id1) == handler.get_hash(fw_id2));

        // bw_id1 is merged with bw_id2
        assert(handler.get_hash(bw_id1) == handler.get_hash(bw_id2));
    } else {

        SHalfedge_handle se2 = SD.new_shalfedge_pair();
        se2->set_forward_index(fw_id2);
        se2->set_backward_index(bw_id2);
        se2->twin()->set_forward_index(bw_id2);
        se2->twin()->set_backward_index(fw_id2);
        se2->circle() = se->circle();
        se2->twin()->circle() = se->circle().opposite();

        handler.handle_support(se, se1, se2);

        // fw_id1 is merged with fw_id2
        assert(handler.get_hash(fw_id1) == handler.get_hash(fw_id2));

        // bw_id1 is merged with bw_id2
        assert(handler.get_hash(bw_id1) == handler.get_hash(bw_id2));
    }
}

int main() {
    const int id1 = 10;
    const int id2 = 11;

    const int id3 = 20;
    const int id4 = 21;

    // edge-edge
    run_test_case(id1, id2, id3, id4, false);
    run_test_case(id3, id4, id1, id2, false);

    // edge-loop
    run_test_case(id1, id2, id3, id4, true);
    run_test_case(id3, id4, id1, id2, true);

    return 0;
}
