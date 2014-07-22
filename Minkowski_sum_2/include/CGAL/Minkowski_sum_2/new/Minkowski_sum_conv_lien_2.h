#ifndef CGAL_MINKOWSKI_SUM_REDUCED_CONV_H
#define CGAL_MINKOWSKI_SUM_REDUCED_CONV_H

#include <CGAL/Arrangement_with_history_2.h>
#include "aabb/AABB_Collision_detector.h"
#include "Arr_SegmentData_traits.h"

#include <queue>
#include <boost/unordered_set.hpp>

namespace CGAL {
namespace internal {

template <class Kernel, class Container>
class Minkowski_sum_by_convolution_lien_2 {

private:

    typedef CGAL::Polygon_2<Kernel, Container> Polygon_2;

    // Kernel types:
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Vector_2 Vector_2;
    typedef typename Kernel::Direction_2 Direction_2;

    // Traits-related types:
    typedef Arr_segment_traits_2<Kernel> Traits_2_A;
    typedef Arr_SegmentData_traits<Traits_2_A> Traits_2;

    typedef typename Traits_2_A::Segment_2 Base_Segment_2;
    typedef typename Traits_2::X_monotone_curve_2 Segment_2;
    typedef std::list<Segment_2> Segments_list;

    typedef Arr_default_dcel<Traits_2> Dcel;

    typedef CGAL::Arrangement_with_history_2<Traits_2, Dcel> Arrangement_history_2;
    typedef typename Arrangement_history_2::Halfedge Halfedge;
    typedef typename Arrangement_history_2::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement_history_2::Edge_iterator Edge_iterator;
    typedef typename Arrangement_history_2::Halfedge_handle Halfedge_handle;
    typedef typename Arrangement_history_2::Vertex_handle Vertex_handle;
    typedef typename Arrangement_history_2::Face_iterator Face_iterator;
    typedef typename Arrangement_history_2::Face_handle Face_handle;
    typedef typename Arrangement_history_2::Hole_iterator Hole_iterator;
    typedef typename Arrangement_history_2::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
    typedef typename Arrangement_history_2::Ccb_halfedge_circulator Ccb_halfedge_circulator;

    typedef typename Arrangement_history_2::Originating_curve_iterator Originating_curve_iterator;
    typedef std::pair<int, int> State;

    // Data members:
    typename Kernel::Construct_translated_point_2 f_add;
    typename Kernel::Construct_vector_2 f_vector;
    typename Kernel::Construct_direction_2 f_direction;
    typename Kernel::Orientation_2 f_orientation;
    typename Kernel::Compare_xy_2 f_compare_xy;
    typename Kernel::Counterclockwise_in_between_2 f_ccw_in_between;
    typename Kernel::Compute_x_2 f_compute_x;
    typename Kernel::Compute_y_2 f_compute_y;

    typename Traits_2::Compare_y_at_x_2 f_compare_y_at_x;
    typename Traits_2::Compare_x_2 f_compare_x;

    AABBCollisionDetector<Kernel, Container> *aabb_collision_detector;

public:

    Minkowski_sum_by_convolution_lien_2() {
        // Obtain kernel functors.
        Kernel ker;

        f_add = ker.construct_translated_point_2_object();
        f_vector = ker.construct_vector_2_object();
        f_direction = ker.construct_direction_2_object();
        f_orientation = ker.orientation_2_object();
        f_compare_xy = ker.compare_xy_2_object();
        f_ccw_in_between = ker.counterclockwise_in_between_2_object();
        f_compare_x = Traits_2().compare_x_2_object();
        f_compare_y_at_x = Traits_2().compare_y_at_x_2_object();
        f_compute_x = ker.compute_x_2_object();
        f_compute_y = ker.compute_y_2_object();
    }

    template <class OutputIterator>
    OutputIterator operator()(const Polygon_2 &pgn1, const Polygon_2 &pgn2,
                              Polygon_2 &sum_bound, OutputIterator sum_holes) {

        CGAL_precondition(pgn1.is_simple());
        CGAL_precondition(pgn2.is_simple());
        CGAL_precondition(pgn1.orientation() == CGAL::COUNTERCLOCKWISE);
        CGAL_precondition(pgn2.orientation() == CGAL::COUNTERCLOCKWISE);

        const Polygon_2 inversed_p1 = transform(Aff_transformation_2<Kernel>(SCALING, -1), pgn1);
        aabb_collision_detector = new AABBCollisionDetector<Kernel, Container>(pgn2, inversed_p1);

        // compute the reduced convolution
        Segments_list reduced_conv;
        build_reduced_convolution(pgn1, pgn2, reduced_conv);

        // split the segments at their intersection points
        Arrangement_history_2 arr;
        CGAL::insert(arr, reduced_conv.begin(), reduced_conv.end());

        // trace outer loop
        get_outer_loop(arr, sum_bound);

        // trace holes
        for (Face_iterator itr = arr.faces_begin(); itr != arr.faces_end(); ++itr) {
            handle_face(arr, itr, inversed_p1, pgn2, sum_holes);
        }

        delete aabb_collision_detector;

        return sum_holes;
    }

private:

    // Builds the reduced convolution using the fiber grid approach. For each
    // starting vertex, try to add two outgoing next states. If a visited
    // vertex is reached, then do not explore further. This is a BFS like
    // iteration beginning from each vertex in the first column of the fiber
    // grid.
    void build_reduced_convolution(const Polygon_2 &pgn1, const Polygon_2 &pgn2, Segments_list &reduced_conv) const {
        unsigned int n1 = pgn1.size();
        unsigned int n2 = pgn2.size();

        // Init the direcions of both polygons
        std::vector<Direction_2> p1_dirs = directions_of_polygon(pgn1);
        std::vector<Direction_2> p2_dirs = directions_of_polygon(pgn2);

        boost::unordered_set<State> visited_vertices;
        boost::unordered_map<std::pair<int, int>, Point_2> points_map;

        // Init the queue with vertices from the first column
        std::queue<State> state_queue;
        for (int i = n1-1; i >= 0; --i) {
            state_queue.push(State(i, 0));
        }

        while (state_queue.size() > 0) {
            State curr_state = state_queue.front();
            state_queue.pop();

            int i1 = curr_state.first;
            int i2 = curr_state.second;

            if (visited_vertices.count(curr_state) > 0) {
                continue;
            }
            visited_vertices.insert(curr_state);

            // add two outgoing edges:
            int next_i1 = (i1+1) % n1;
            int next_i2 = (i2+1) % n2;
            int prev_i1 = (n1+i1-1) % n1;
            int prev_i2 = (n2+i2-1) % n2;

            // add geometric entites of the transition from state (i,j) to (i+1,j) and (i,j+1), if they are in the reduced convolution.

            bool is_end_coincide;
            bool is_start_coincide;

            // Add an edge from pgn2
            if (p1_dirs[prev_i1] == p2_dirs[i2] || f_ccw_in_between(p2_dirs[i2], p1_dirs[prev_i1], p1_dirs[i1])) {
                state_queue.push(State(i1, next_i2));

                if (check_convex(pgn1[prev_i1], pgn1[i1], pgn1[next_i1])) {
                    Point_2 start_point = get_point(i1, i2, points_map, pgn1, pgn2);
                    Point_2 end_point = get_point(i1, next_i2, points_map, pgn1, pgn2);

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point, end_point), Segment_Data_Label(state(i1, i2), state(i1, next_i2), cres, 1));

                    reduced_conv.push_back(conv_seg);
                }
            }

            // Add an edge from pgn1
            if (p2_dirs[i2] == p1_dirs[i1] || f_ccw_in_between(p1_dirs[i1], p2_dirs[prev_i2], p2_dirs[i2])) {
                state_queue.push(State(next_i1, i2));

                if (check_convex(pgn2[prev_i2], pgn2[i2], pgn2[next_i2])) {
                    Point_2 start_point = get_point(i1, i2, points_map, pgn1, pgn2);
                    Point_2 end_point = get_point(next_i1, i2, points_map, pgn1, pgn2);

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point, end_point), Segment_Data_Label(state(i1, i2), state(next_i1, i2), cres, 0));
                    reduced_conv.push_back(conv_seg);
                }
            }
        }
    }

    std::vector<Direction_2> directions_of_polygon(const Polygon_2 &pgn1) const {
        std::vector<Direction_2> directions;
        unsigned int n = pgn1.size();

        for (int i = 0; i < n-1; ++i) {
            directions.push_back(f_direction(f_vector(pgn1[i], pgn1[i+1])));
        }
        directions.push_back(f_direction(f_vector(pgn1[n-1], pgn1[0])));

        return directions;
    }

    bool check_convex(const Point_2 &prev, const Point_2 &curr, const Point_2 &next) const {
        return f_orientation(prev, curr, next) == LEFT_TURN;
    }

    // Gets point corresponding to a state (i,j) if exists, creates this point if asked for first time.
    Point_2 get_point(int i1, int i2, boost::unordered_map<std::pair<int, int>, Point_2> &points_map, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Point_2 result;

        if (points_map.count(State(i1, i2)) == 0) {
            result = f_add(pgn1[i1], Vector_2(Point_2(ORIGIN), pgn2[i2]));
            points_map[State(i1, i2)] = result;
        } else {
            result = points_map[State(i1, i2)];
        }

        return result;
    }

    // Put the outside loop of the arrangement in 'out_bound'
    void get_outer_loop(Arrangement_history_2 &arr, Polygon_2 &out_bound) {
        Ccb_halfedge_circulator circ_start = *(arr.unbounded_face()->holes_begin());
        Ccb_halfedge_circulator circ = circ_start;

        do {
            out_bound.push_back(circ->source()->point());
        } while (--circ != circ_start);
    }

    // Check whether the face is on the M-sum's border. Add it to 'holes' if it is.
    template <class OutputIterator>
    void handle_face(Arrangement_history_2 &arr, Face_handle itr, const Polygon_2 &reverse_pgn1, const Polygon_2 &pgn2, OutputIterator holes) {

        // If the face contains holes, it can't be on the Minkowski sum's border
        if (itr->holes_begin() != itr->holes_end()) {
            return;
        }

        Ccb_halfedge_circulator start = itr->outer_ccb();
        Ccb_halfedge_circulator circ = start;

        // The face needs to be orientable
        do {
            if (!check_originating_edge_has_same_direction(arr, circ)) {
                return;
            }
        } while (++circ != start);

        // Check whether this is a false hole (TODO: explain)
        if (checkCollisionDetection(arr, start, reverse_pgn1, pgn2)) {
            return;
        }

        // mark as hole
        circ = start;
        Polygon_2 pgn_hole;

        do {
            pgn_hole.push_back(circ->source()->point());
        } while (--circ != start);

        *holes = pgn_hole;
        ++holes;
    }

    // Check whether the originating edge(s) had the same direction as the current half edge
    bool check_originating_edge_has_same_direction(Arrangement_history_2 &arr, Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            if (segment_itr->label()._orientation == (CGAL::Comparison_result)he->direction()) {
                return false;
            }
        }

        return true;
    }

    /*
    This version assumes poly1 is reflected through origin.
    */
    bool checkCollisionDetection(Arrangement_history_2 &arr, Halfedge_handle &handle, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Point_2 mid_point = handle->source()->point();
        Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, mid_point)), pgn1);
        aabb_collision_detector->setTranslationPoint(mid_point);
        return aabb_collision_detector->checkCollision(t_pgn1, pgn2);
    }
};

} // namespace internal
} // namespace CGAL

#endif
