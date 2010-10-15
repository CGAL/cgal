#include <CGAL/basic.h>

#include <cassert>
#include <fstream>

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/Apollonius_graph_2/Delage_traits/Apollonius_graph_new_traits_2.h>

typedef CGAL::Sign Sign;

template <class GT, class GT_test>
class Check_traits : public GT
{
    GT_test test;
public:

    typedef typename GT::Site_2 Site_2;

    Check_traits (const GT &gt = GT(), const GT_test &gt_test = GT_test())
        : GT(gt), test (gt_test)
    {}

    template <class Predicate, class New_predicate>
    struct Checked_predicate : public Predicate
    {
        New_predicate newp;
        typedef typename Predicate::result_type result_type;

        Checked_predicate (const Predicate &p, const New_predicate &np)
            : Predicate (p), newp (np)
        {}

        template <class T1>
        result_type operator() (const T1 &t1) const
        {
            result_type r1 = Predicate::operator() (t1);
            result_type r2 = newp (t1);
            assert(r1 == r2);
            return r1;
        }
        template <class T1, class T2>
        result_type operator() (const T1 &t1, const T2 &t2) const
        {
            result_type r1 = Predicate::operator() (t1, t2);
            result_type r2 = newp (t1, t2);
            assert(r1 == r2);
            return r1;
        }
        template <class T1, class T2, class T3>
        result_type operator() (const T1 &t1, const T2 &t2, const T3 &t3) const
        {
            result_type r1 = Predicate::operator() (t1, t2, t3);
            result_type r2 = newp (t1, t2, t3);
            assert(r1 == r2);
            return r1;
        }
        template <class T1, class T2, class T3, class T4>
        result_type operator() (const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4) const
        {
            result_type r1 = Predicate::operator() (t1, t2, t3, t4);
            result_type r2 = newp (t1, t2, t3, t4);
            assert(r1 == r2);
            return r1;
        }
        template <class T1, class T2, class T3, class T4, class T5>
        result_type operator() (const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5) const
        {
            result_type r1 = Predicate::operator() (t1, t2, t3, t4, t5);
            result_type r2 = newp (t1, t2, t3, t4, t5);
            assert(r1 == r2);
            return r1;
        }
        template <class T1, class T2, class T3, class T4, class T5, class T6>
        result_type operator() (const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6) const
        {
            result_type r1 = Predicate::operator() (t1, t2, t3, t4, t5, t6);
            result_type r2 = newp (t1, t2, t3, t4, t5, t6);
            assert(r1 == r2);
            return r1;
        }
    };

    typedef Checked_predicate<typename GT::Vertex_conflict_2,
                              typename GT_test::Vertex_conflict_2>
                                  Vertex_conflict_2;

    typedef Checked_predicate<typename GT::Finite_edge_interior_conflict_2,
                              typename GT_test::Finite_edge_interior_conflict_2>
                                  Finite_edge_interior_conflict_2;

    typedef Checked_predicate<typename GT::Infinite_edge_interior_conflict_2,
                              typename GT_test::Infinite_edge_interior_conflict_2>
                                  Infinite_edge_interior_conflict_2;

    Vertex_conflict_2
    vertex_conflict_2_object() const
    {
        return Vertex_conflict_2 (
                GT::vertex_conflict_2_object(),
                test.vertex_conflict_2_object());
    }
    Finite_edge_interior_conflict_2
    finite_edge_interior_conflict_2_object() const
    {
        return Finite_edge_interior_conflict_2 (
                GT::finite_edge_interior_conflict_2_object(),
                test.finite_edge_interior_conflict_2_object());
    }
    Infinite_edge_interior_conflict_2
    infinite_edge_interior_conflict_2_object() const
    {
        return Infinite_edge_interior_conflict_2 (
                GT::infinite_edge_interior_conflict_2_object(),
                test.infinite_edge_interior_conflict_2_object());
    }
};

typedef CGAL::MP_Float NT;
typedef CGAL::Simple_cartesian<NT> K;
typedef CGAL::Apollonius_graph_traits_2<K> GT_old;
typedef CGAL::Apollonius_graph_new_traits_2<K> GT_new;
typedef Check_traits<GT_old,GT_new> GT;
typedef CGAL::Apollonius_graph_2<GT> AG;

typedef AG::Site_2 Site;
typedef AG::Finite_faces_iterator Finite_faces_iterator;
typedef AG::Finite_vertices_iterator Finite_vertices_iterator;
typedef AG::Vertex_handle Vertex_handle;

void test_orientation (const AG &ag,
        Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
{
    const Site &s0 = v0->site();
    const Site &s1 = v1->site();
    const Site &s2 = v2->site();

    GT gt;
    GT_new gt_new;

    GT::Orientation_2 orientation1 = gt.orientation_2_object();
    GT_new::Orientation_new_2 orientation2 = gt_new.orientation_new_2_object();

    for (Finite_vertices_iterator _v = ag.finite_vertices_begin();
            _v != ag.finite_vertices_end(); ++_v)
    {
        Vertex_handle v = _v;
        if (v == v0) continue;

        Sign o1 = orientation1 (s0, s1, s2, s0, v->site());
        Sign o2 = orientation2 (s0, s1, s2,     v->site().point());

        assert(o1 == o2);
    }
}

void test_file (const char *filename)
{
    std::ifstream is (filename);

    assert(is);

    std::cout << "File " << filename << ": construction... " << std::flush;
    AG ag;
    Site s; while (is >> s) ag.insert (s);

    std::cout << "validation... " << std::flush;
    assert(ag.is_valid());

    std::cout << "OK" << std::endl;

    std::cout << "Test orientation... " << std::flush;

    for (Finite_faces_iterator f = ag.finite_faces_begin();
            f != ag.finite_faces_end(); ++f)
    {
        test_orientation (ag, f->vertex(0), f->vertex(1), f->vertex(2));
        test_orientation (ag, f->vertex(1), f->vertex(2), f->vertex(0));
        test_orientation (ag, f->vertex(2), f->vertex(0), f->vertex(1));
    }

    std::cout << "OK" << std::endl;
}

int main (void)
{
    test_file ("data/traits.dat");
    test_file ("data/algo.dat");

    return 0;
}
