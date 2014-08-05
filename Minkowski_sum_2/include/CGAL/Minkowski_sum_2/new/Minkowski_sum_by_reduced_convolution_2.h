#ifndef CGAL_MINKOWSKI_SUM_BY_REDUCED_CONVOLUTION_2_H
#define CGAL_MINKOWSKI_SUM_BY_REDUCED_CONVOLUTION_2_H

#include <CGAL/basic.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Timer.h> // TODO: remove when optimization is done

#include <CGAL/Minkowski_sum_2/new/aabb/AABB_collision_detector_2.h>

#include <iostream> // TODO: remove when optimization is done
#include <queue>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {

template <class Kernel_, class Container_>
class Minkowski_sum_by_reduced_convolution_2 {

private:

    typedef Kernel_ Kernel;
    typedef Container_ Container;

    // Basic types:
    typedef CGAL::Polygon_2<Kernel, Container> Polygon_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Vector_2 Vector_2;
    typedef typename Kernel::Direction_2 Direction_2;
    typedef typename Kernel::Triangle_2 Triangle_2;
    typedef typename Kernel::FT FT;

    // Segment-related types:
    typedef Arr_segment_traits_2<Kernel> Traits_2;
    typedef typename Traits_2::X_monotone_curve_2 Segment_2;
    typedef std::list<Segment_2> Segment_list;
    typedef Arr_default_dcel<Traits_2> Dcel;
    typedef std::pair<int, int> State;

    // Arrangement-related types:
    typedef Arrangement_with_history_2<Traits_2, Dcel> Arrangement_history_2;
    typedef typename Arrangement_history_2::Halfedge_handle Halfedge_handle;
    typedef typename Arrangement_history_2::Face_iterator Face_iterator;
    typedef typename Arrangement_history_2::Face_handle Face_handle;
    typedef typename Arrangement_history_2::Ccb_halfedge_circulator Ccb_halfedge_circulator;
    typedef typename Arrangement_history_2::Originating_curve_iterator Originating_curve_iterator;

    // Function object types:
    typename Kernel::Construct_translated_point_2 f_add;
    typename Kernel::Construct_vector_2 f_vector;
    typename Kernel::Construct_direction_2 f_direction;
    typename Kernel::Orientation_2 f_orientation;
    typename Kernel::Compare_xy_2 f_compare_xy;
    typename Kernel::Counterclockwise_in_between_2 f_ccw_in_between;

public:

    Minkowski_sum_by_reduced_convolution_2() {
        // Obtain kernel functors
        Kernel ker;
        f_add = ker.construct_translated_point_2_object();
        f_vector = ker.construct_vector_2_object();
        f_direction = ker.construct_direction_2_object();
        f_orientation = ker.orientation_2_object();
        f_compare_xy = ker.compare_xy_2_object();
        f_ccw_in_between = ker.counterclockwise_in_between_2_object();
    }

    template <class OutputIterator>
    void operator()(const Polygon_2 &pgn1, const Polygon_2 &pgn2,
                    Polygon_2 &outer_boundary, OutputIterator holes) const {

            Timer timer; // TODO: remove when optimization is done
            timer.start();

        CGAL_precondition(pgn1.is_simple());
        CGAL_precondition(pgn2.is_simple());
        CGAL_precondition(pgn1.orientation() == COUNTERCLOCKWISE);
        CGAL_precondition(pgn2.orientation() == COUNTERCLOCKWISE);

            timer.stop();
            std::cout << timer.time() << " s: Preconditions" << std::endl;
            timer.reset();
            timer.start();

        // Initialize collision detector. It operates on pgn2 and on the inversed pgn1:
        const Polygon_2 inversed_pgn1 = transform(Aff_transformation_2<Kernel>(SCALING, -1), pgn1);
        AABB_collision_detector_2<Kernel, Container> collision_detector(pgn2, inversed_pgn1);

            timer.stop();
            std::cout << timer.time() << " s: AABB init" << std::endl;
            timer.reset();
            timer.start();

        // Compute the reduced convolution
        Segment_list reduced_convolution;
        build_reduced_convolution(pgn1, pgn2, reduced_convolution);

            timer.stop();
            std::cout << timer.time() << " s: Convolution" << std::endl;
            timer.reset();
            timer.start();

        // Insert the segments into an arrangement
        Arrangement_history_2 arr;
        insert(arr, reduced_convolution.begin(), reduced_convolution.end());

            timer.stop();
            std::cout << timer.time() << " s: Arrangement" << std::endl;
            timer.reset();
            timer.start();

        // Trace the outer loop and put it in 'outer_boundary'
        get_outer_loop(arr, outer_boundary);

            timer.stop();
            std::cout << timer.time() << " s: Outer Loop" << std::endl;
            timer.reset();
            timer.start();

        // Check for each face whether it is a hole in the M-sum. If it is, add it to 'holes'.
        for (Face_iterator face = arr.faces_begin(); face != arr.faces_end(); ++face) {
            handle_face(arr, face, holes, collision_detector);
        }

            timer.stop();
            std::cout << timer.time() << " s: Holes" << std::endl;
    }

private:

    // Builds the reduced convolution using a fiber grid approach. For each
    // starting vertex, try to add two outgoing next states. If a visited
    // vertex is reached, then do not explore further. This is a BFS-like
    // iteration beginning from each vertex in the first column of the fiber
    // grid.
    void build_reduced_convolution(const Polygon_2 &pgn1, const Polygon_2 &pgn2, Segment_list &reduced_convolution) const {
        unsigned int n1 = pgn1.size();
        unsigned int n2 = pgn2.size();

        // Init the direcions of both polygons
        std::vector<Direction_2> p1_dirs = directions_of_polygon(pgn1);
        std::vector<Direction_2> p2_dirs = directions_of_polygon(pgn2);

        // Contains states that were already visited
        boost::unordered_set<State> visited_states;

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

            // If this state was already visited, skip it
            if (visited_states.count(curr_state) > 0) {
                continue;
            }
            visited_states.insert(curr_state);

            int next_i1 = (i1+1) % n1;
            int next_i2 = (i2+1) % n2;
            int prev_i1 = (n1+i1-1) % n1;
            int prev_i2 = (n2+i2-1) % n2;

            // Try two transitions: From (i,j) to (i+1,j) and to (i,j+1). Add
            // the respective segments, if they are in the reduced convolution.
            for(int step_in_pgn1 = 0; step_in_pgn1 <= 1; step_in_pgn1++) {
                int new_i1, new_i2;
                if (step_in_pgn1) {
                    new_i1 = next_i1;
                    new_i2 = i2;
                } else {
                    new_i1 = i1;
                    new_i2 = next_i2;
                }

                // If the segment's direction lies counterclockwise in between
                // the other polygon's vertex' ingoing and outgoing directions,
                // the segment belongs to the full convolution.
                bool belongs_to_convolution;
                if (step_in_pgn1) {
                    belongs_to_convolution = f_ccw_in_between(p1_dirs[i1], p2_dirs[prev_i2], p2_dirs[i2]) ||
                                             p1_dirs[i1] == p2_dirs[i2];
                } else {
                    belongs_to_convolution = f_ccw_in_between(p2_dirs[i2], p1_dirs[prev_i1], p1_dirs[i1]) ||
                                             p2_dirs[i2] == p1_dirs[prev_i1];
                }

                if (belongs_to_convolution) {
                    state_queue.push(State(new_i1, new_i2));

                    // Only edges added to convex vertices can be on the M-sum's boundary.
                    // This filter only leaves the *reduced* convolution.
                    bool convex;
                    if (step_in_pgn1) {
                        convex = is_convex(pgn2[prev_i2], pgn2[i2], pgn2[next_i2]);
                    } else {
                        convex = is_convex(pgn1[prev_i1], pgn1[i1], pgn1[next_i1]);
                    }

                    if (convex) {
                        Point_2 start_point = get_point(i1, i2, pgn1, pgn2);
                        Point_2 end_point = get_point(new_i1, new_i2, pgn1, pgn2);
                        reduced_convolution.push_back(Segment_2(start_point, end_point));
                    }
                }
            }
        }
    }

    // Returns a sorted list of the polygon's edges
    std::vector<Direction_2> directions_of_polygon(const Polygon_2 &p) const {
        std::vector<Direction_2> directions;
        unsigned int n = p.size();

        for (int i = 0; i < n-1; ++i) {
            directions.push_back(f_direction(f_vector(p[i], p[i+1])));
        }
        directions.push_back(f_direction(f_vector(p[n-1], p[0])));

        return directions;
    }

    bool is_convex(const Point_2 &prev, const Point_2 &curr, const Point_2 &next) const {
        return f_orientation(prev, curr, next) == LEFT_TURN;
    }

    // Returns the point corresponding to a state (i,j).
    Point_2 get_point(int i1, int i2, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {

        return f_add(pgn1[i1], Vector_2(Point_2(ORIGIN), pgn2[i2]));
    }

    // Put the outer loop of the arrangement in 'outer_boundary'
    void get_outer_loop(Arrangement_history_2 &arr, Polygon_2 &outer_boundary) const {
        Ccb_halfedge_circulator circ_start = *(arr.unbounded_face()->holes_begin());
        Ccb_halfedge_circulator circ = circ_start;

        do {
            outer_boundary.push_back(circ->source()->point());
        } while (--circ != circ_start);
    }

    // Check whether the face is on the M-sum's border. Add it to 'holes' if it is.
    template <class OutputIterator>
    void handle_face(const Arrangement_history_2 &arr, const Face_handle face, OutputIterator holes, AABB_collision_detector_2<Kernel, Container> &collision_detector) const {

        // If the face contains holes, it can't be on the Minkowski sum's border
        if (face->holes_begin() != face->holes_end()) {
            return;
        }

        Ccb_halfedge_circulator start = face->outer_ccb();
        Ccb_halfedge_circulator circ = start;

        // The face needs to be orientable
        do {
            if (!do_original_edges_have_same_direction(arr, circ)) {
                return;
            }
        } while (++circ != start);

        // When the reversed polygon 1, translated by a point inside of this face, collides with polygon 2, this cannot be a hole
        Point_2 inner_point = get_point_in_face(face);
        if (collision_detector.check_collision(inner_point)) {
            return;
        }

        // At this point, the face is a real hole, add it to 'holes'
        Polygon_2 pgn_hole;
        circ = start;

        do {
            pgn_hole.push_back(circ->source()->point());
        } while (--circ != start);

        *holes = pgn_hole;
        ++holes;
    }

    // Check whether the convolution's original edge(s) had the same direction as the arrangement's half edge
    bool do_original_edges_have_same_direction(const Arrangement_history_2 &arr, const Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            if (f_compare_xy(segment_itr->source(), segment_itr->target()) == (Comparison_result)he->direction()) {
                return false;
            }
        }

        return true;
    }

    // Return a point in the face's interior by finding a diagonal
    Point_2 get_point_in_face(const Face_handle face) const {
        Ccb_halfedge_circulator current_edge = face->outer_ccb();
        Ccb_halfedge_circulator next_edge = current_edge;
        next_edge++;

        Point_2 a, v, b;

        // Move over the face's vertices until a convex corner is encountered:
        do {
            a = current_edge->source()->point();
            v = current_edge->target()->point();
            b = next_edge->target()->point();

            current_edge++;
            next_edge++;
        } while (!is_convex(a, v, b));

        Triangle_2 ear(a, v, b);
        FT min_distance = -1;
        Point_2 min_q;

        // Of the remaining vertices, find the one inside of the "ear" with minimal distance to v:
        while (++next_edge != current_edge) {
            Point_2 q = next_edge->target()->point();
            if (ear.has_on_bounded_side(q)) {
                FT distance = squared_distance(q, v);
                if (min_distance == -1 || distance < min_distance) {
                    min_distance = distance;
                    min_q = q;
                }
            }
        }

        // If there was no vertex inside of the ear, return it's centroid.
        // Otherwise, return a point between v and min_q.
        if (min_distance == -1) {
            return centroid(ear);
        } else {
            return midpoint(v, min_q);
        }
    }
};

} // namespace CGAL

#endif
