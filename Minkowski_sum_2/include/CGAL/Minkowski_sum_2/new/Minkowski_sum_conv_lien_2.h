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

    // Kernel functors:
    typedef typename Kernel::Construct_translated_point_2 Translate_point_2;
    typedef typename Kernel::Construct_vector_2 Construct_vector_2;
    typedef typename Kernel::Construct_direction_2 Construct_direction_2;
    typedef typename Kernel::Orientation_2 Compute_orientation_2;
    typedef typename Kernel::Compare_xy_2 Compare_xy_2;
    typedef typename Kernel::Counterclockwise_in_between_2 Ccw_in_between_2;
    typedef typename Kernel::Compare_x_2 Compare_x_2;
    typedef typename Kernel::Compute_x_2 Compute_x_2;
    typedef typename Kernel::Compute_y_2 Compute_y_2;

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
    typedef std::pair<int, int> StatePair;

    // Data members:
    Translate_point_2 f_add;
    Construct_vector_2 f_vector;
    Construct_direction_2 f_direction;
    Compute_orientation_2 f_orientation;
    Compare_xy_2 f_compare_xy;
    Ccw_in_between_2 f_ccw_in_between;
    Compute_x_2 f_compute_x;
    Compute_y_2 f_compute_y;

    typename Traits_2::Compare_y_at_x_2 f_compare_y_at_x;
    typename Traits_2::Compare_x_2 f_compare_x;

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
        _aabb_collision_detector = new AABBCollisionDetector<Kernel, Container>(pgn2, inversed_p1);

        // compute the reduced convolution
        Segments_list reduced_conv;
        buildReducedConvolutionFiberGrid(pgn1, pgn2, reduced_conv);

        // split the segments at their intersection points
        Arrangement_history_2 arr;
        CGAL::insert(arr, reduced_conv.begin(), reduced_conv.end());

        // trace outer loop
        markOutsideLoop(arr, sum_bound);

        // turn pgn1 by 180 degrees
        Polygon_2 rotated_pgn1 = transform(Aff_transformation_2<Kernel>(ROTATION, 0, -1), pgn1);

        // trace holes
        for (Face_iterator itr = arr.faces_begin(); itr != arr.faces_end(); ++itr) {
            handleFace(arr, itr, rotated_pgn1, pgn2, sum_holes);
        }

        delete _aabb_collision_detector;

        return sum_holes;
    }

private:

    void markOutsideLoop(Arrangement_history_2 &arr, Polygon_2 &out_bound) {
        Face_iterator ub_face = arr.unbounded_face();
        Hole_iterator holes_itr = ub_face->holes_begin();
        Ccb_halfedge_circulator circ_start = *holes_itr;
        Ccb_halfedge_circulator circ = circ_start;

        do {
            out_bound.push_back(circ->source()->point());
            --circ;
        } while (circ != circ_start);
    }

    template <class OutputIterator>
    void handleFace(Arrangement_history_2 &arr, Face_handle itr, const Polygon_2 &reverse_pgn1, const Polygon_2 &pgn2, OutputIterator holes) {
        if (itr->holes_begin() != itr->holes_end()) {
            return;
        }

        Ccb_halfedge_circulator start = itr->outer_ccb();
        Ccb_halfedge_circulator circ = start;

        // orientation check
        do {
            if (!checkTripNotSameDirWithSegment(arr, circ)) {
                return;
            }

            ++circ;
        } while (circ != start);

        // collision detection check.
        bool coll_detect = checkCollisionDetection(arr, start, reverse_pgn1, pgn2);

        if (coll_detect) {
            return;
        }

        // mark as hole
        circ = start;
        Polygon_2 pgn_hole;

        do {
            pgn_hole.push_back(circ->source()->point());
            --circ;
        } while (circ != start);

        *holes = pgn_hole;
        ++holes;
    }

    std::vector<Direction_2> fillPolyDirs(const Polygon_2 &pgn1) const {
        std::vector<Direction_2> outVec;
        unsigned int n = pgn1.size();

        for (int i = 0; i < n-1; ++i) {
            outVec.push_back(f_direction(f_vector(pgn1[i], pgn1[i+1])));
        }
        outVec.push_back(f_direction(f_vector(pgn1[n-1], pgn1[0])));

        return outVec;
    }

    // Gets point corresponding to a state (i,j) if exists, creates this point if asked for first time.
    Point_2 addGetPoint(int i1, int i2, boost::unordered_map<std::pair<int, int>, Point_2> &points_map, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Point_2 result;

        if (points_map.count(StatePair(i1, i2)) == 0) {
            result = f_add(pgn1[i1], Vector_2(Point_2(ORIGIN), pgn2[i2]));
            points_map[StatePair(i1, i2)] = result;
        } else {
            result = points_map[StatePair(i1, i2)];
        }

        return result;
    }

    // Builds the reduced convolution using the fiber grid approach. For each
    // starting vertex, try to add out-going next states (two states). If a
    // visited vertex is reached then do not explore. This is a BFS like
    // iteration beginning from each vertex in the first column of the fiber
    // grid.
    void buildReducedConvolutionFiberGrid(const Polygon_2 &pgn1, const Polygon_2 &pgn2, Segments_list &reduced_conv) const {
        unsigned int n1 = pgn1.size();
        unsigned int n2 = pgn2.size();

        // Init the direcions of both polygons
        std::vector<Direction_2> p1_dirs = fillPolyDirs(pgn1);
        std::vector<Direction_2> p2_dirs = fillPolyDirs(pgn2);

        boost::unordered_set<StatePair> visited_vertices_set;
        std::queue<StatePair> state_queue;
        boost::unordered_map<std::pair<int, int>, Point_2> points_map;

        // Init the queue with vertices from the first column
        for (int i = n1-1; i >= 0; --i) {
            state_queue.push(StatePair(i, 0));
        }

        while (state_queue.size() > 0) {
            StatePair curr_state = state_queue.front();
            state_queue.pop();

            int i1 = curr_state.first;
            int i2 = curr_state.second;

            if (visited_vertices_set.count(curr_state) > 0) {
                continue;
            }
            visited_vertices_set.insert(curr_state);

            // add two outgoing edges:
            int next_p1 = (i1+1) % n1;
            int next_p2 = (i2+1) % n2;
            int prev_p1 = (n1+i1-1) % n1;
            int prev_p2 = (n2+i2-1) % n2;

            StatePair next_state_p1 = StatePair(next_p1, i2);
            StatePair next_state_p2 = StatePair(i1, next_p2);

            // add geometric entites of the transition from state (i,j) to (i+1,j) and (i,j+1), if they are in the reduced convolution.

            bool is_end_coincide;
            bool is_start_coincide;

            // Add an edge from Q
            if (checkSwept(p1_dirs[prev_p1], p1_dirs[i1], p2_dirs[i2], is_start_coincide, is_end_coincide) && !is_end_coincide) {
                state_queue.push(next_state_p2);

                if (!checkReflex(pgn1[prev_p1], pgn1[i1], pgn1[next_p1])) {
                    Point_2 start_point = addGetPoint(i1, i2, points_map, pgn1, pgn2);
                    Point_2 end_point = addGetPoint(i1, next_p2, points_map, pgn1, pgn2);

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point, end_point), Segment_Data_Label(state(i1, i2), state(i1, next_p2), cres, 1));

                    reduced_conv.push_back(conv_seg);
                }
            }

            // Add an edge from P
            if (checkSwept(p2_dirs[prev_p2], p2_dirs[i2], p1_dirs[i1], is_start_coincide, is_end_coincide) && !is_start_coincide) {
                state_queue.push(next_state_p1);

                if (!checkReflex(pgn2[prev_p2], pgn2[i2], pgn2[next_p2])) {
                    Point_2 start_point = addGetPoint(i1, i2, points_map, pgn1, pgn2);
                    Point_2 end_point = addGetPoint(next_p1, i2, points_map, pgn1, pgn2);

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(typename Traits_2_A::Segment_2(start_point, end_point), Segment_Data_Label(state(i1, i2), state(next_p1, i2), cres, 0));
                    reduced_conv.push_back(conv_seg);
                }
            }
        }
    }

    AABBCollisionDetector<Kernel, Container> *_aabb_collision_detector;

    /*
    This version assumes poly1 is reflected through origin. (as called from nested loops filter)
    */
    bool checkCollisionDetection(Arrangement_history_2 &arr, Halfedge_handle &handle, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Point_2 mid_point = findInsidePoint(arr, handle);
        Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, mid_point)), pgn1);
        _aabb_collision_detector->setTranslationPoint(mid_point);
        return _aabb_collision_detector->checkCollision(t_pgn1, pgn2);
    }

    Point_2 findInsidePoint(Arrangement_history_2 &arr, Halfedge_handle &handle) const {
        Ccb_halfedge_circulator currHandle = handle->ccb();
        Ccb_halfedge_circulator nextHandle = currHandle;
        ++nextHandle;

        while (currHandle->direction() != nextHandle->direction()) {
            ++currHandle;
            ++nextHandle;

            if (checkReflex(currHandle->source()->point(), currHandle->target()->point(), nextHandle->target()->point())) {
                break;
            }
        }

        Point_2 p = currHandle->source()->point();
        Point_2 p2 = currHandle->target()->point();
        Point_2 work_point = p2;

        Ccb_halfedge_circulator best_edge = handle;
        bool has_some_point = false;
        Ccb_halfedge_circulator circ = nextHandle;
        Ccb_halfedge_circulator end = handle;

        bool shoot_upwards = (currHandle->direction() == ARR_LEFT_TO_RIGHT);

        if (nextHandle->curve().is_vertical()) {
            work_point = CGAL::midpoint(p, p2);
        }

        if (currHandle->curve().is_vertical()) {
            p = nextHandle->source()->point();
            p2 = nextHandle->target()->point();
            work_point = CGAL::midpoint(p, p2);
            ++best_edge;
            ++circ;
            ++end;
        }

        ++circ;

        while (circ != end) {
            Base_Segment_2 circ_curve = circ->curve();

            if (f_compare_x(work_point, circ_curve.min()) != f_compare_x(work_point, circ_curve.max())) {
                // we have an edge with same x range as endpoint of
                bool above_first = (f_compare_y_at_x(work_point, circ_curve) == SMALLER);

                if (has_some_point) {
                    bool under_best;
                    Base_Segment_2 best_edge_curve = best_edge->curve();

                    if (f_compare_x(best_edge_curve.min(), circ_curve.min()) != f_compare_x(best_edge_curve.max(), circ_curve.min())) {
                        under_best = f_compare_y_at_x(circ_curve.min(), best_edge_curve) == SMALLER;
                    } else {
                        under_best = f_compare_y_at_x(best_edge_curve.min(), circ_curve) != SMALLER;
                    }

                    if ((shoot_upwards && above_first && under_best) || (!shoot_upwards && !above_first && !under_best)) {
                        best_edge = circ;
                    }
                } else {
                    has_some_point = true;
                    best_edge = circ;
                }
            }

            ++circ;
        }

        if (best_edge->curve().is_vertical()) {
            Base_Segment_2 best_edge_curve = best_edge->curve();
            typename Kernel::FT x0 = f_compute_x(work_point);
            typename Kernel::FT y_point = f_compute_y(work_point);

            if (shoot_upwards) {
                typename Kernel::FT y_best = f_compute_y(best_edge_curve.min());
                typename Kernel::FT y = (y_best - y_point) / 2 + y_point;
                return Point_2(x0, y);
            } else {
                typename Kernel::FT y_best = f_compute_y(best_edge_curve.min());
                typename Kernel::FT y = (y_point - y_best) / 2 + y_best;
                return Point_2(x0, y);
            }

            return work_point;
        }

        Base_Segment_2 best_edge_curve = best_edge->curve();
        typename Kernel::FT x0 = f_compute_x(work_point);
        typename Kernel::FT x1 = f_compute_x(best_edge_curve.min());
        typename Kernel::FT x2 = f_compute_x(best_edge_curve.max());
        typename Kernel::FT alpha = (x0 - x2) / (x1 - x2);

        typename Kernel::FT y_best = alpha * f_compute_y(best_edge_curve.min()) + (1 - alpha) * f_compute_y(best_edge_curve.max());
        typename Kernel::FT y_point = f_compute_y(work_point);
        typename Kernel::FT y = (y_best - y_point) / 2 + y_point;

        return Point_2(x0, y);
    }

    bool checkTripNotSameDirWithSegment(Arrangement_history_2 &arr, Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            Segment_2 segment = *segment_itr;

            CGAL::Comparison_result c1 = segment.label()._orientation;
            CGAL::Comparison_result c2 = (CGAL::Comparison_result)he->direction();
            bool same_dir = (c1 != c2);

            if (same_dir) {
                return true;
            }
        }

        return false;
    }

    bool checkReflex(const Point_2 &prev, const Point_2 &curr, const Point_2 &next) const {
        CGAL::Orientation res_ori = f_orientation(prev, curr, next);
        return ((res_ori == RIGHT_TURN) || (res_ori == COLLINEAR));
    }

    bool checkSwept(Direction_2 &dir_start, Direction_2 &dir_end, Direction_2 &dir_new, bool &isStartConcide, bool &isEndConcide) const {
        isStartConcide = dir_new == dir_start;
        isEndConcide = dir_end == dir_new;

        return isStartConcide || f_ccw_in_between(dir_new, dir_start, dir_end) || isEndConcide;
    }
};

} // namespace internal
} // namespace CGAL

#endif
