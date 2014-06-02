#ifndef CGAL_MINKOWSKI_SUM_REDUCED_CONV_H
#define CGAL_MINKOWSKI_SUM_REDUCED_CONV_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>

#include <fstream>
#include <ostream>
#include <list>
#include <set>
#include <utility>
#include <algorithm>
#include <iterator>
#include <valarray>
#include <boost/unordered_set.hpp>
#include <queue>

#include <CGAL/Minkowski_sum_2/new/Arr_SegmentData_traits.h>
#include <CGAL/Minkowski_sum_2/new/ICollisionDetector.h>
#include <CGAL/Minkowski_sum_2/new/NaiveCollisionDetector.h>
#include <CGAL/Minkowski_sum_2/new/SweepCollisionDetection.h>
#include <CGAL/Minkowski_sum_2/new/AABB_Collision_detector.h>

#include <boost/unordered_map.hpp>
#include <boost/timer.hpp>

namespace CGAL {

struct Less_than_handle {
    template <typename Type>
    bool operator()(Type s1, Type s2) const {
        return (&(*s1) < & (*s2));
    }
};

template <class HalfedgeBase_>
class Arr_map_halfedge : public HalfedgeBase_ {
public:
    bool visited;
    bool isDegenerate;
    int loopNumber;
};

template <class Traits_,
          class HalfedgeBase_ = Arr_halfedge_base<typename Traits_::X_monotone_curve_2> > class Arr_my_extended_dcel :
    public Arr_dcel_base<Arr_vertex_base<typename Traits_::Point_2>,
    Arr_map_halfedge<HalfedgeBase_>,
    Arr_face_base> {
public:

    template<typename T>
    class rebind {
        typedef typename HalfedgeBase_::template rebind
        <typename T::X_monotone_curve_2> Rebind_halfedge;
        typedef typename Rebind_halfedge::other Halfedge_base;

    public:

        typedef Arr_my_extended_dcel<T, Halfedge_base> other;
    };
};

template <class Kernel_, class Container_>
class Minkowski_sum_by_convolution_lien_2 {
public:

    typedef Kernel_ Kernel;
    typedef CGAL::Polygon_2<Kernel, Container_> Polygon_2;

public:

    // Kernel types:
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Vector_2 Vector_2;
    typedef typename Kernel::Direction_2 Direction_2;

    // Kernel functors:
    typedef typename Kernel::Equal_2 Equal_2;
    typedef typename Kernel::Construct_translated_point_2 Translate_point_2;
    typedef typename Kernel::Construct_vector_2 Construct_vector_2;
    typedef typename Kernel::Construct_direction_2 Construct_direction_2;
    typedef typename Kernel::Construct_opposite_line_2 Opposite_line_2;
    typedef typename Kernel::Orientation_2 Compute_orientation_2;
    typedef typename Kernel::Compare_xy_2 Compare_xy_2;
    typedef typename Kernel::Counterclockwise_in_between_2 Ccw_in_between_2;
    typedef typename Kernel::Angle_2 Compute_Angle_2;
    typedef typename Kernel::Compare_x_2 Compare_x_2;
    typedef typename Kernel::Is_vertical_2 Is_vertical_2;
    typedef typename Kernel::Compute_x_2 Compute_x_2;
    typedef typename Kernel::Compute_y_2 Compute_y_2;

    // Polygon-related types:
    typedef typename Polygon_2::Vertex_circulator Vertex_circulator;
    typedef std::pair<Vertex_circulator, unsigned int> Vertex_ref;
    typedef std::pair<Vertex_ref, Vertex_ref> Anchor;
    typedef std::list<Anchor> Anchors_queue;

    // Traits-related types:
    typedef Arr_segment_traits_2<Kernel> Traits_2_A;
    typedef Arr_SegmentData_traits<Traits_2_A> Traits_2;

    typedef typename Traits_2_A::Segment_2 Base_Segment_2;
    typedef typename Traits_2::X_monotone_curve_2 Segment_2;
    typedef std::list<Segment_2> Segments_list;

    typedef CGAL::Arr_my_extended_dcel<Traits_2> Dcel;

    typedef CGAL::Arrangement_with_history_2<Traits_2, Dcel> Arrangement_history_2;
    typedef typename Arrangement_history_2::Halfedge Halfedge;
    typedef typename Arrangement_history_2::Vertex Vertex;
    typedef typename Arrangement_history_2::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement_history_2::Halfedge_iterator Halfedge_iterator;
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

    typedef std::set<Halfedge_handle, Less_than_handle> Edges_set;
    typedef std::set<Face_handle, Less_than_handle> Faces_set;

    // Data members:
    Equal_2 f_equal;
    Translate_point_2 f_add;
    Construct_vector_2 f_vector;
    Construct_direction_2 f_direction;
    Opposite_line_2 f_opp_line;
    Compute_orientation_2 f_orientation;
    Compare_xy_2 f_compare_xy;
    Ccw_in_between_2 f_ccw_in_between;
    Compute_Angle_2 f_angle;
    Is_vertical_2 f_is_vertical;
    Compute_x_2 f_compute_x;
    Compute_y_2 f_compute_y;

    typename Traits_2::Compare_endpoints_xy_2 f_compare_endpoints_xy;
    typename Traits_2::Compare_y_at_x_2 f_compare_y_at_x;
    typename Traits_2::Compare_x_2 f_compare_x;

    friend class ConvSegMapper;

    struct ConvSegment {
        Halfedge_handle _he;
        ConvSegment(Halfedge_handle &he): _he(he) {}
        ConvSegment() {}
        bool getVisited() const {
            return _he->visited;
        }

        bool getDegenerate() const {
            return _he->isDegenerate;
        }

        int getLoopNum() {
            return _he->loopNumber;
        }

        Vertex_handle getSrc() {
            return (_he->source());
        }

        Vertex_handle getDst() {
            return (_he->target());
        }

        bool operator<(const ConvSegment &rhs) const {
            return Less_than_handle()(_he, rhs._he);
        }

        bool operator==(const ConvSegment &rhs) const {
            return !Less_than_handle()(_he, rhs._he) && !Less_than_handle()(rhs._he, _he);
        }
    };

    struct ConvSegMapper {
        Arrangement_history_2 *_arr;
        Minkowski_sum_by_convolution_lien_2 *_mink;
        ConvSegMapper(Arrangement_history_2 *arr, Minkowski_sum_by_convolution_lien_2 *mink): _arr(arr), _mink(mink) {
        }

        ConvSegment getSegment(const Halfedge_handle &he) {
            return ConvSegment(_mink->getDirAgreeingHalfedge(*_arr, he));
        }

        Direction_2 getConvSegDir(const ConvSegment &seg) const {
            return _mink->getHalfedgeDir(seg._he);
        }

        void markVisited(ConvSegment &convSeg, int id) {
            _mink->setEdgeVisited(*convSeg._he, true, id);
        }

        bool isBBiggerThenAWithReagrdToC(const ConvSegment &a, const ConvSegment &b, const ConvSegment &c) const {
            Direction_2 dir_a = getConvSegDir(a);
            Direction_2 dir_b = getConvSegDir(b);
            Direction_2 dir_c = getConvSegDir(c);
            return _mink->isDirImproving(dir_a, dir_c, dir_b);
        }

        void getNeighbouringSegments(Vertex_handle v, std::list<ConvSegment> &outSegments, std::list<ConvSegment> &inSegmets) {
            std::list<Halfedge_handle> inList, outList;
            _mink->getEdgesFromVertex(*_arr, v, inList, outList);
            typename std::list<Halfedge_iterator>::const_iterator itr;

            for (itr = inList.begin(); itr != inList.end(); ++itr) {
                inSegmets.push_back(getSegment(*itr));
            }

            for (itr = outList.begin(); itr != outList.end(); ++itr) {
                outSegments.push_back(getSegment(*itr));
            }
        }

        static bool getSegVisited(const ConvSegment &seg) {
            return (seg.getVisited() == true || seg.getDegenerate());
        }

        static bool getSegNotVisited(const ConvSegment &seg) {
            return (seg.getVisited() == false && !seg.getDegenerate());
        }

        double getSignedAngle(const ConvSegment &enter, const ConvSegment &exit) {
            return _mink->getSignedAngle(enter._he, exit._he);
        }

        void filterNonVisitedSegments(const std::list<ConvSegment> &inputList, std::list<ConvSegment> &outList, std::list<ConvSegment> &visitedSegmentsList) {
            int out_size = count_if(inputList.begin(), inputList.end(), &ConvSegMapper::getSegNotVisited);
            outList.resize(out_size);
            remove_copy_if(inputList.begin(), inputList.end(), outList.begin(), &ConvSegMapper::getSegVisited);
            visitedSegmentsList.resize(inputList.size() - out_size);
            remove_copy_if(inputList.begin(), inputList.end(), visitedSegmentsList.begin(), &ConvSegMapper::getSegNotVisited);
        }

        struct SegCompare {
            ConvSegment _incomingSeg;
            ConvSegMapper *_mapperInstance;
            SegCompare(ConvSegment incomingSeg, ConvSegMapper *mapperInstance): _incomingSeg(incomingSeg), _mapperInstance(mapperInstance) {}
            bool operator()(ConvSegment a, ConvSegment b) {
                return _mapperInstance->isBBiggerThenAWithReagrdToC(a, b, _incomingSeg);
            }
        };

        ConvSegment getMaximalEdge(std::list<ConvSegment> outgoingSegs, ConvSegment incomingSeg) {
            SegCompare seg_comp(incomingSeg, this);
            return *max_element(outgoingSegs.begin(), outgoingSegs.end(), seg_comp);
        }

        bool checkLoopClosed(const Vertex_handle &v, ConvSegment &startingSeg, int id) {
            Halfedge_handle h;
            bool notFound = _mink->checkOutgoingNotVisited(*_arr, *v, h, id);

            if (h != Halfedge_handle()) {
                startingSeg = getSegment(h);
            }

            return !notFound;
        }

        template<typename T> void fillEdgesSet(T &edges_set) {
            Edge_iterator itr;

            for (itr = _arr->edges_begin(); itr != _arr->edges_end(); ++itr) {
                _mink->setEdgeVisited(*itr, false, -1);

                if (!itr->isDegenerate) {
                    edges_set.insert(getSegment(itr));
                }
            }
        }

        ConvSegment getOuterSegment() {
            Face_iterator startFace = _arr->unbounded_face();
            Halfedge_iterator perimiterFace = *(startFace -> holes_begin());
            return getSegment(perimiterFace);
        }

        void removeSegFromArr(const ConvSegment &seg) {
            _arr->remove_edge(seg._he);
        }

        void removeRangeFromArr(std::list<ConvSegment> &segsToRemove) {
            for (typename std::list<ConvSegment>::iterator itr = segsToRemove.begin(); itr != segsToRemove.end(); ++itr) {
                removeSegFromArr(*itr);
            }
        }
    };

    class TraversalManager;
    friend class ConvMovement;

    struct ConvMovement {
        ConvMovement(ConvSegMapper *mapperInstance, int id, ConvSegment &startedge, TraversalManager *managerInstance): _id(id), _mapperInstance(mapperInstance), _managerInstance(managerInstance) {
            _currEdge = startedge;
            _traversedEdges.push_back(_currEdge);
            _mapperInstance->markVisited(startedge, id);
            _traversing = true;
            _loop = NULL;
        }

        void Traverse() {
            while (_traversing) {
                Vertex_handle dst = _currEdge.getDst();
                std::list<ConvSegment> inList;
                std::list<ConvSegment> outList;
                _mapperInstance->getNeighbouringSegments(dst, outList, inList);
                std::list<ConvSegment> filteredSegments, visitedSegments;
                _mapperInstance->filterNonVisitedSegments(outList, filteredSegments, visitedSegments);

                if (checkCloseLoop()) {
                    closeLoopEvent(filteredSegments);
                }

                if (!_traversing) {
                    break;
                }

                // Check if movement may proceed.
                if (filteredSegments.size() > 0) {
                    ConvSegment c = _mapperInstance->getMaximalEdge(filteredSegments, _currEdge);
                    // Check that there is no visited outgoing edge, with diffrenet loop id (ie closed loop),
                    // which is better than c.
                    typename std::list<ConvSegment>::iterator itr = visitedSegments.begin();

                    for (; itr != visitedSegments.end(); ++itr) {
                        if (_mapperInstance->isBBiggerThenAWithReagrdToC(c, *itr, _currEdge) && itr->getLoopNum() != _id && itr->getLoopNum() != 0) {
                            closeEvent();
                            _traversing = false;
                            return;
                        }
                    }

                    double angle = _mapperInstance->getSignedAngle(_currEdge, c);
                    _angles.push_back(angle);
                    _anglesSum += angle;
                    _currEdge = c;
                    _mapperInstance->markVisited(c, _id);
                    _traversedEdges.push_back(c);
                } else {
                    closeEvent();
                    _traversing = false;
                }
            }
        }
        void closeLoopEvent(std::list<ConvSegment> &outGoingOptions) {
            _managerInstance->handleCloseLoopEvent(*_loop, outGoingOptions);
            _angles.pop_back();
        }

        void closeEvent() {
            _managerInstance->handleCloseEvent();
        }

        bool checkCloseLoop() {
            ConvSegment loopStart;
            bool found = _mapperInstance->checkLoopClosed(_currEdge.getDst(), loopStart, _id);

            if (found) {
                if (_loop != NULL) {
                    delete _loop;
                }

                _loop = new std::pair<ConvSegment, ConvSegment>(loopStart, _currEdge);
                double angle = _mapperInstance->getSignedAngle(_currEdge, loopStart);
                _angles.push_back(angle);
            }

            return found;
        }

        void stopTraverse() {
            _traversing = false;
        }

        std::list<ConvSegment> &getTraversedEdges() {
            return _traversedEdges;
        }

        std::list<double> &getAngles() {
            return _angles;
        }
    private:
        int _id;
        ConvSegment _currEdge;
        std::list<ConvSegment> _traversedEdges;
        std::list<double> _angles;
        ConvSegMapper *_mapperInstance;
        TraversalManager *_managerInstance;
        bool _traversing;
        double _anglesSum;
        std::pair<ConvSegment, ConvSegment> *_loop;
    };

    friend class TraversalManager;
    class EdgesStore;
    class LoopsTracker;
    struct TraversalManager {
        TraversalManager(ConvSegMapper *mapperInstance): _mapperInstance(mapperInstance) {
            _loopId = 0;
            _activeMovment = NULL;
        }
        void traverseLoops() {
            EdgesStore edges_db(_mapperInstance);
            _currEdgesStore = &edges_db;
            int s = edges_db.getSize();
            ConvSegment startEdge = _mapperInstance->getOuterSegment();
            traceLoop(edges_db, startEdge, true);

            while (!edges_db.isEmpty()) {
                int s = edges_db.getSize();
                ++_loopId;
                startEdge = edges_db.getEdge();
                traceLoop(edges_db, startEdge, false);
            }
        }

        void traceLoop(EdgesStore &edgesDB, ConvSegment &startSeg, bool isOuter) {
            _activeMovment = new ConvMovement(_mapperInstance, _loopId, startSeg, this);

            _activeLoopTracker = new LoopsTracker(&(_activeMovment->getTraversedEdges()), &(_activeMovment->getAngles()), isOuter);

            _activeMovment->Traverse();

            delete _activeMovment;
            delete _activeLoopTracker;
        }

        void handleCloseLoopEvent(std::pair<ConvSegment, ConvSegment> &loop, std::list<ConvSegment> &outGoingOptions) {
            _activeLoopTracker->addLoop(loop);

            // check for improvment: (TODO: maybe refactor)
            // This checks if the first segment in loop (actually second edge goes into first) is a worse choice then *itr , w.r.t second.
            // This means the loop might be improved.
            if (_activeLoopTracker->hasLoops()) {
                bool stopItrating = true;

                for (typename std::list<ConvSegment>::iterator itr = outGoingOptions.begin(); itr != outGoingOptions.end(); ++itr) {
                    if (_mapperInstance->isBBiggerThenAWithReagrdToC(loop.first, *itr, loop.second)) {
                        stopItrating = false;
                        break;
                    }
                }

                if (stopItrating) {
                    handleCloseEvent();
                    _activeMovment->stopTraverse();
                }
            }
        }

        void handleCloseEvent() {
            // remove all traversed edges from segments list

            std::list<ConvSegment> &traversedEdges = _activeMovment->getTraversedEdges();
            _currEdgesStore->removeRange(traversedEdges.begin(), traversedEdges.end());

            // remove anything but loop from arrangment.
            std::list<ConvSegment> &toRemoveFromArr = _activeLoopTracker->getNonLoopSegmentsList();

            _mapperInstance->removeRangeFromArr(toRemoveFromArr);
        }

    private:
        ConvSegMapper *_mapperInstance;
        ConvMovement *_activeMovment;
        LoopsTracker *_activeLoopTracker;
        EdgesStore *_currEdgesStore;
        int _loopId;
    };

    struct LoopsTracker {
        LoopsTracker(std::list<ConvSegment> *segmentsList, std::list<double> *anglesList, bool outerLoop): _segmentsList(segmentsList), _anglesList(anglesList), _outerLoop(outerLoop) {
            _hasLoop = false;
            _loopsCreated = false;
        }

        void addLoop(std::pair<ConvSegment, ConvSegment> &loop) {
            typename std::list<ConvSegment>::iterator loop_begin = find(_segmentsList->begin(), _segmentsList->end(), (loop.first));
            typename std::list<ConvSegment>::iterator loop_end = find(_segmentsList->begin(), _segmentsList->end(), (loop.second));
            double angleSum = sumLoopAngles(loop_begin, loop_end);

            if ((angleSum > 0 && _outerLoop) || (angleSum < 0 && !_outerLoop)) {
                // update loop
                _hasLoop = true;
                _loopBegin = loop_begin;
                _loopEnd = loop_end;
            }
        }

        double sumLoopAngles(typename std::list<ConvSegment>::iterator &loopBegin, typename std::list<ConvSegment>::iterator &loopEnd) {
            int begin = std::distance(_segmentsList->begin(), loopBegin);
            int end = std::distance(_segmentsList->begin(), loopEnd);
            std::list<double>::iterator itr_begin = _anglesList->begin();
            advance(itr_begin, begin);
            std::list<double>::iterator itr_end = _anglesList->begin();
            advance(itr_end, end + 1);
            double sum = 0;

            for (std::list<double>::iterator itr = itr_begin; itr != itr_end; ++itr) {
                sum += *itr;
            }

            return sum;
        }

        std::list<ConvSegment> &getLoopSegmentsList() {
            if (!_loopsCreated) {
                createLoopsSegmentsLists();
            }

            return _loopSegmentsList;
        }

        std::list<ConvSegment> &getNonLoopSegmentsList() {
            if (!_loopsCreated) {
                createLoopsSegmentsLists();
            }

            return _nonLoopSegmentsList;
        }

        void createLoopsSegmentsLists() {
            copy(_segmentsList->begin(), _segmentsList->end(), back_inserter(_nonLoopSegmentsList));

            if (_hasLoop) {
                typename std::list<ConvSegment>::iterator loop_begin = find(_nonLoopSegmentsList.begin(), _nonLoopSegmentsList.end(), *_loopBegin);
                typename std::list<ConvSegment>::iterator loop_end = find(_nonLoopSegmentsList.begin(), _nonLoopSegmentsList.end(), *_loopEnd);
                _loopSegmentsList.splice(_loopSegmentsList.begin(), _nonLoopSegmentsList, loop_begin, ++loop_end);
            }
        }

        bool hasLoops() {
            return _hasLoop;
        }

    private:
        typename std::list<ConvSegment>::iterator _loopBegin, _loopEnd;
        std::list<ConvSegment> *_segmentsList;
        std::list<double> *_anglesList;
        std::list<ConvSegment> _loopSegmentsList;
        std::list<ConvSegment> _nonLoopSegmentsList;
        bool _outerLoop;
        bool _hasLoop;
        bool _loopsCreated;
    };

    struct EdgesStore {
        EdgesStore(ConvSegMapper *mapper) {
            _edgesSet = new std::set<ConvSegment>();
            mapper->fillEdgesSet(*_edgesSet);
        }

        template <typename InputIterator> void removeRange(InputIterator start, InputIterator end) {
            for (InputIterator itr = start; itr != end; ++itr) {
                _edgesSet->erase(*itr);
            }
        }

        template <typename InputIterator> void removeRange_bak(InputIterator start, InputIterator end) {
            std::vector<typename InputIterator::value_type> vec(start, end);
            sort(vec.begin(), vec.end());
            std::vector<typename InputIterator::value_type> *seg_set = new std::vector<typename InputIterator::value_type>();
            set_difference(_edgesSet->begin(), _edgesSet->end(), vec.begin(), vec.end(), back_inserter(*seg_set));
            delete _edgesSet;
            _edgesSet = new std::set<ConvSegment>(seg_set->begin(), seg_set->end());
        }

        ConvSegment getEdge() {
            return *(_edgesSet->begin());
        }

        bool isEmpty() {
            return (_edgesSet->size() == 0);
        }

        int getSize() {
            return _edgesSet->size();
        }

    private:
        std::set<ConvSegment> *_edgesSet;
    };

    friend class DegenerateCassesManager;
    struct DegenerateCassesManager {
        DegenerateCassesManager(Arrangement_history_2 *arr, Minkowski_sum_by_convolution_lien_2 *mink, Polygon_2 *poly1, Polygon_2 *poly2, bool isActive): _arr(arr), _mink(mink), _poly1(poly1), _poly2(poly2), _active(isActive) {
        }

        void markDegenerateEdges() {
            Edge_iterator itr = _arr->edges_begin();

            for (; itr != _arr->edges_end(); ++itr) {
                _mink->setEdgeVisited(*itr, false, -1);

                if (_active) {
                    if (_mink->checkDegenerateEdgeOppositeSegments(*_arr, itr)) {
                        if (!_mink->checkSegmentCollisionDetection(*_arr, itr->curve(), *_poly1, *_poly2)) {
                            _mink->setEdgeDegenerate(*itr, true);
                        } else {
                            _mink->setEdgeDegenerate(*itr, false);
                        }
                    } else {
                        _mink->setEdgeDegenerate(*itr, false);
                    }
                } else {
                    _mink->setEdgeDegenerate(*itr, false);
                }
            }
        }

        void findDegenerateBorderVertices() {
            if (!_active) {
                return;
            }

            Vertex_iterator itr = _arr->vertices_begin();

            for (; itr != _arr->vertices_end(); ++itr) {
                if (_mink->checkDegenarateVertexIsIntersectionOfThreeSegments(*_arr, itr))
                {
                    Point_2 p_end = itr->point();

                    if (!_mink->checkCollisionDetection(*_arr, itr->point(), *_poly1, *_poly2)) {
                        _degenerate_points_list.push_back(itr->point());
                    }
                }
            }
        }

        void addDegenerateVerticesToArr() {
            typename std::list<Point_2>::iterator itr = _degenerate_points_list.begin();

            for (; itr != _degenerate_points_list.end(); ++itr) {
                CGAL::insert_point(*_arr, *itr);
            }
        }

    private:
        Arrangement_history_2 *_arr;
        Minkowski_sum_by_convolution_lien_2 *_mink;
        std::list<Point_2> _degenerate_points_list;
        Polygon_2 *_poly1;
        Polygon_2 *_poly2;
        bool _active;
    };

public:

    /*! Default constructor. */
    Minkowski_sum_by_convolution_lien_2() {
        // Obtain kernel functors.
        Kernel ker;

        f_equal = ker.equal_2_object();
        f_add = ker.construct_translated_point_2_object();
        f_vector = ker.construct_vector_2_object();
        f_direction = ker.construct_direction_2_object();
        f_opp_line = ker.construct_opposite_line_2_object();
        f_orientation = ker.orientation_2_object();
        f_compare_xy = ker.compare_xy_2_object();
        f_ccw_in_between = ker.counterclockwise_in_between_2_object();
        f_angle = ker.angle_2_object();
        f_compare_endpoints_xy = Traits_2().compare_endpoints_xy_2_object();
        f_is_vertical = ker.is_vertical_2_object();
        f_compare_x = Traits_2().compare_x_2_object();
        f_compare_y_at_x = Traits_2().compare_y_at_x_2_object();
        f_compute_x = ker.compute_x_2_object();
        f_compute_y = ker.compute_y_2_object();
    }

    template <class OutputIterator>
    OutputIterator operator()(const Polygon_2 &pgn1,
                              const Polygon_2 &pgn2,
                              Polygon_2 &sum_bound,
                              OutputIterator sum_holes) {
        CGAL_precondition(pgn1.is_simple());
        CGAL_precondition(pgn2.is_simple());
        CGAL_precondition(pgn1.orientation() == CGAL::COUNTERCLOCKWISE);
        CGAL_precondition(pgn2.orientation() == CGAL::COUNTERCLOCKWISE);
        Polygon_2 revP1 = revPoly(pgn1);
        Polygon_2 p2 = pgn2;
        boost::timer t_abb;
        _aabb_collision_detector = new AABBCollisionDetector<Kernel_, Container_>(p2, revP1);
        double aabb_build_time = t_abb.elapsed();
        Segments_list reduced_conv;
        boost::timer t1;
        buildReducedConvolutionFiberGrid(pgn1, pgn2, reduced_conv);
        Arrangement_history_2 arr;
        boost::timer t2;
        buildArrangementFromConv(reduced_conv, arr);
        const Minkowski_sum_by_convolution_lien_2 *ptr = this;
        boost::timer t4;
        DegenerateCassesManager degHandler(&arr, const_cast <Minkowski_sum_by_convolution_lien_2 *>(ptr), const_cast <Polygon_2 *>(&pgn1), const_cast <Polygon_2 *>(&pgn2), true);
        degHandler.findDegenerateBorderVertices();
        degHandler.markDegenerateEdges();

        boost::timer t3;
        Polygon_2 reverse_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Rotation(), 0, -1), pgn1);

        // trace outer loop
        markOutsideLoop(arr, sum_bound);

        // trace holes

        // Original code here !
        for (Face_iterator itr = arr.faces_begin(); itr != arr.faces_end(); ++itr) {
            handleFace(arr, itr, reverse_pgn1, pgn2, sum_holes);
        }

        std::list<Halfedge_handle> removeList;

        // remove all non marked edges
        for (Edge_iterator itr = arr.edges_begin(); itr != arr.edges_end(); ++itr) {
            if ((!itr->visited) && (!itr->isDegenerate)) {
                removeList.push_back(itr);
            }
        }

        for (typename std::list<Halfedge_handle>::iterator itr = removeList.begin(); itr != removeList.end(); ++itr) {
            arr.remove_edge(*itr);
        }

        degHandler.addDegenerateVerticesToArr();

        delete _aabb_collision_detector;
        return (sum_holes);
    }

    void markOutsideLoop(Arrangement_history_2 &arr) {
        Face_iterator ub_face = arr.unbounded_face();
        Hole_iterator holes_itr = ub_face->holes_begin();
        Ccb_halfedge_circulator circ_start = *holes_itr;
        Ccb_halfedge_circulator circ = circ_start;

        do {
            setEdgeVisited(*circ, true, 0);
            ++circ;
        } while (circ != circ_start);
    }

    void markOutsideLoop(Arrangement_history_2 &arr, Polygon_2 &out_bound) {
        Face_iterator ub_face = arr.unbounded_face();
        Hole_iterator holes_itr = ub_face->holes_begin();
        Ccb_halfedge_circulator circ_start = *holes_itr;
        Ccb_halfedge_circulator circ = circ_start;

        do {
            setEdgeVisited(*circ, true, 0);
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
            setEdgeVisited(*circ, true, 0);
            pgn_hole.push_back(circ->source()->point());
            --circ;
        } while (circ != start);

        *holes = pgn_hole;
        ++holes;
    }

    /*!
     * Compute the boundery Minkowski sum of two simple polygons.
     * The result is represented as
     * the outer boundary of the Minkowski sum (which is always a simple polygon).
     * \param pgn1 The first polygon.
     * \param pgn2 The second polygon.
     * \param sum_bound Output: A polygon respresenting the outer boundary
     * of the Minkowski sum.
     *
     * \pre Both input polygons are simple.
     * \return A past-the-end iterator for the holes in the sum.
     */
    template <class OutputIterator>
    OutputIterator operator()(const Polygon_2 &pgn1,
                              const Polygon_2 &pgn2,
                              Polygon_2 &sum_bound) const {
        CGAL_precondition(pgn1.is_simple());
        CGAL_precondition(pgn2.is_simple());

        Segments_list reduced_conv;
        buildReducedConvolution(pgn1, pgn2, reduced_conv);
    }

    void fillPolyDirs(const Polygon_2 &pgn1, std::vector<Direction_2> &outVec) const {
        unsigned int n1 = pgn1.size();

        for (int i = 0; i < (n1 - 1); ++i) {
            outVec[i] = f_direction(f_vector(pgn1[i], pgn1[i + 1]));
        }

        outVec[n1 - 1] = f_direction(f_vector(pgn1[n1 - 1], pgn1[0]));
    }
    void buildReducedConvolution(const Polygon_2 &pgn1, const Polygon_2 &pgn2, Segments_list &reduced_conv) const {
        unsigned int n1 = pgn1.size();
        unsigned int n2 = pgn2.size();
        Vertex_circulator vert_p1, vert_p2, prev_p1, prev_p2, next_p1, next_p2;
        vert_p1 = prev_p1 = next_p1 = pgn1.vertices_circulator();
        vert_p2 = prev_p2 = next_p2 = pgn2.vertices_circulator();

        --prev_p1;
        --prev_p2;
        ++next_p1;
        ++next_p2;
        bool is_end_coincide;
        bool is_start_coincide;

        boost::unordered_map<std::pair<int, int>, Point_2> points_map;

        for (unsigned int i1 = 0; i1 < n1; ++i1) {
            for (unsigned int i2 = 0; i2 < n2; ++i2) {
                points_map[std::pair<int, int>(i1, i2)] = f_add(*vert_p1, Vector_2(Point_2(ORIGIN), *vert_p2));
                ++vert_p2;
            }

            ++vert_p1;
        }

        std::vector<Direction_2> p1_dirs(n1);
        std::vector<Direction_2> p2_dirs(n2);

        fillPolyDirs(pgn1, p1_dirs);
        fillPolyDirs(pgn2, p2_dirs);

        vert_p1 = pgn1.vertices_circulator();
        vert_p2 = pgn2.vertices_circulator();

        for (unsigned int i1 = 0; i1 < n1; ++i1) {
            for (unsigned int i2 = 0; i2 < n2; ++i2) {

                Point_2 start_point = points_map[std::pair<int, int>(i1, i2)];
                int prev_i1 = i1 - 1;

                if (prev_i1 == -1) {
                    prev_i1 = n1 - 1;
                }

                int prev_i2 = i2 - 1;

                if (prev_i2 == -1) {
                    prev_i2 = n2 - 1;
                }

                if (!checkReflex(*prev_p1, *vert_p1, *next_p1) && checkSwept(p1_dirs[prev_i1], p1_dirs[i1], p2_dirs[i2], is_start_coincide, is_end_coincide)) {
                    int cyc_ind = i2;
                    ++cyc_ind;

                    if (cyc_ind == n2) {
                        cyc_ind = 0;
                    }

                    Point_2 end_point = points_map[std::pair<int, int>(i1, cyc_ind)];

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point, end_point), cres);

                    if (!is_end_coincide) {
                        reduced_conv.push_back(conv_seg);
                    }
                }

                if (!checkReflex(*prev_p2, *vert_p2, *next_p2) && checkSwept(p2_dirs[prev_i2], p2_dirs[i2], p1_dirs[i1], is_start_coincide, is_end_coincide)) {
                    int cyc_ind = i1;
                    ++cyc_ind;

                    if (cyc_ind == n1) {
                        cyc_ind = 0;
                    }

                    Point_2 end_point = points_map[std::pair<int, int>(cyc_ind, i2)];

                    CGAL::Comparison_result cres = f_compare_xy(start_point, end_point);
                    Segment_2 conv_seg = Segment_2(Traits_2_A::Segment_2(start_point, end_point), cres);

                    if (!is_start_coincide) {
                        reduced_conv.push_back(conv_seg);
                    }
                }

                prev_p2 = vert_p2;
                vert_p2 = next_p2;
                ++next_p2;
            }

            prev_p1 = vert_p1;
            vert_p1 = next_p1;
            ++next_p1;
        }
    }

    // Increse a cyclic integer counter with limit lim.
    int cyclicInc(int i, int lim) const {
        i = i + 1;

        if (i >= lim) {
            i = 0;
        }

        return i;
    }

    // Decrease a cyclic integer counter with limit lim.
    int cyclicDec(int i, int lim) const {
        i = i - 1;

        if (i < 0) {
            i = lim - 1;
        }

        return i;
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

    // Builds the reduced convolution using the fiber grid approach. for each starting vertex, try to add out-going next states(two states).
    // If a visited vertex is reached then do not explore. This is a BFS like iteration beggining from each vertex in the first column of the
    // fiber grid.
    void buildReducedConvolutionFiberGrid(const Polygon_2 &pgn1, const Polygon_2 &pgn2, Segments_list &reduced_conv) const {
        unsigned int n1 = pgn1.size();
        unsigned int n2 = pgn2.size();
        int i1 = 0;
        int i2 = 0;

        bool is_end_coincide;
        bool is_start_coincide;

        // Init the direcions of both polygons.
        std::vector<Direction_2> p1_dirs(n1);
        std::vector<Direction_2> p2_dirs(n2);

        fillPolyDirs(pgn1, p1_dirs);
        fillPolyDirs(pgn2, p2_dirs);

        boost::unordered_set<StatePair > visited_vertices_set;
        std::queue<StatePair > state_queue;
        boost::unordered_map<std::pair<int, int>, Point_2> points_map;

        // init the queue with vertices from the first column
        for (int i = n1 - 1; i >= 0; --i) {
            state_queue.push(StatePair(i, 0));
        }

        while (state_queue.size() > 0) {
            StatePair curr_state = state_queue.front();
            state_queue.pop();

            i1 = curr_state.first;
            i2 = curr_state.second;

            if (visited_vertices_set.count(curr_state) > 0) {
                continue;
            }

            visited_vertices_set.insert(curr_state);

            // add two outgoing edges:
            int next_p1 = cyclicInc(i1, n1);
            int next_p2 = cyclicInc(i2, n2);
            int prev_p1 = cyclicDec(i1, n1);
            int prev_p2 = cyclicDec(i2, n2);

            StatePair next_state_p1 = StatePair(next_p1, i2);
            StatePair next_state_p2 = StatePair(i1, next_p2);

            // add geometric entites of the transition from state (i,j) to (i+1,j) and (i,j+1), if they are in the reduced convolution.

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

private:
    SweepCollisionDetector<Kernel, Container_> collision_detector;
    AABBCollisionDetector<Kernel, Container_> *_aabb_collision_detector;
    /*
        Performs the stage 3 of filtering: removing nested loops which are not in the correct orientation.
    */
    void nestedLoopsFilter(Arrangement_history_2 &arr, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Face_iterator startFace = arr.unbounded_face();
        Hole_iterator hi;
        Ccb_halfedge_circulator perimiterFace;

        for (hi = startFace->holes_begin(); hi != startFace->holes_end(); ++hi) {
            perimiterFace = *hi;
        }

        // now for each hole in main face we will determine it's orientation.
        nestedLoopsFilterRec(arr, perimiterFace, false, pgn1, pgn2);
    }

    /*
    For each loop we represent it by the halfedge which is the twin of the face loop.
    ie this edge is part of the clockwise oriented halfedges loop outside the face.
    */
    void nestedLoopsFilterRec(Arrangement_history_2 &arr, Halfedge_handle &handle, bool inwards, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        std::list<Halfedge_handle> holesEdges;
        std::list<Halfedge_handle> after_removal_hole_edges;
        std::list<Face_handle> faces_list;

        Face_iterator startFace = arr.unbounded_face();
        Halfedge_handle inside_face_edge = handle->twin();
        Face_handle container_face = inside_face_edge->face();

        Face_iterator face_itr = arr.faces_begin();

        for (; face_itr != arr.faces_end(); ++face_itr) {
            if (!(face_itr->is_unbounded()) && face_itr != container_face) {
                Halfedge_handle h_e = face_itr->outer_ccb()->twin();
                holesEdges.push_back(h_e);
            }
        }

        while (holesEdges.size() > 0) { // remove loops from faces
            Halfedge_handle he = holesEdges.front();
            holesEdges.pop_front();

            if (!he->isDegenerate) {
                if (!removeAllNonConformingLoops(arr, he, inwards, holesEdges, pgn1, pgn2)) {
                    after_removal_hole_edges.push_back(he);
                }
            }
        }

        Hole_iterator hi;
        hi = startFace->holes_begin();
        container_face = (*hi)->twin()->face();
        Ccb_halfedge_circulator outside_face_itr;

        // push initial list of holes to list.
        for (hi = container_face->holes_begin(); hi != container_face->holes_end(); ++hi) { // remove degenrate cases which are false
            outside_face_itr = *hi;
            Halfedge_handle h_e = outside_face_itr; //->twin();
            holesEdges.push_back(h_e);
        }

        while (holesEdges.size() > 0) {
            Halfedge_handle he = holesEdges.front();
            holesEdges.pop_front();

            if (!he->isDegenerate) {
                if (!removeAllNonConformingLoops(arr, he, inwards, holesEdges, pgn1, pgn2)) {
                    Hole_iterator hi;
                    after_removal_hole_edges.push_back(he);
                }
            }
        }
    }

    // find all faces who has handle's face sorrounding them, but are not holes.
    void findSemiHoles(Arrangement_history_2 &arr, Halfedge_handle &handle, std::list<Halfedge_handle> &semi_holes) const {
        Faces_set faces_in_face;
        Face_handle outside_face = handle->face();
        Ccb_halfedge_circulator circ = handle->twin()->ccb();
        Ccb_halfedge_circulator curr = circ;

        do {
            // check if the twin edge is the regular outside face or a different one. if so we are in half island.
            Face_handle curr_sec_face = curr->twin()->face();

            if (curr_sec_face != outside_face) {
                typename Faces_set::iterator itr = faces_in_face.find(curr_sec_face);

                if (itr == faces_in_face.end()) {
                    faces_in_face.insert(curr_sec_face);
                }
            }
        } while (++curr != circ);

        for (typename Faces_set::iterator itr = faces_in_face.begin(); itr != faces_in_face.end(); ++itr) {
            Face_handle h = *itr;
            semi_holes.push_back((h->outer_ccb()->twin()));
        }
    }

    bool removeAllNonConformingLoops(Arrangement_history_2 &arr, Halfedge_handle &handle, bool inwards, std::list<Halfedge_handle> &holesEdges, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        // check if loop is in the right direction ((a || b) && !(a && b))
        bool a = !checkTripSameDirWithSegment(arr, handle);
        bool b = inwards;
        // nor a,b
        bool conforming_loop = !((a || b) && !(a && b));

        if (!conforming_loop) {
            removeFaceLoop(arr, handle->twin());
            return true;
        }

        // Check collision detection criterion
        bool coll_detect = checkCollisionDetection(arr, handle, pgn1, pgn2);

        if (coll_detect) {
            removeFaceLoop(arr, handle->twin());
            return true;
        }

        return false;
    }

    Polygon_2 revPoly(const Polygon_2 &input) const {
        Polygon_2 out;
        typename Polygon_2::Vertex_iterator itr = input.vertices_begin();

        for (; itr != input.vertices_end(); ++itr) {
            out.push_back(Point_2(-itr->x(), -itr->y()));
        }

        if (out.orientation() == CGAL::CLOCKWISE) {
            out.reverse_orientation();
        }

        return out;
    }

    AABBCollisionDetector<Kernel, Container_> *getColDetect() const {
        return _aabb_collision_detector;
    }

    /*
    This version assumes poly1 is reflected through origin. (as called from nested loops filter)
    */
    bool checkCollisionDetection(Arrangement_history_2 &arr, Halfedge_handle &handle, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        AABBCollisionDetector<Kernel, Container_> *collision_detector = getColDetect();
        Point_2 mid_point = findInsidePoint(arr, handle);
        Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, mid_point)), pgn1);
        collision_detector->setTranslationPoint(mid_point);
        return collision_detector->checkCollision(t_pgn1, pgn2);
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

    /*
        This version reflects poly 1.
    */
    bool checkCollisionDetection(Arrangement_history_2 &arr, Point_2 &point, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        AABBCollisionDetector<Kernel, Container_> *collision_detector = getColDetect();
        Point_2 p = point;
        Polygon_2 r_pgn1 = revPoly(pgn1);
        Polygon_2 t_pgn1 = transform(typename Kernel::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, p)), r_pgn1);
        collision_detector->setTranslationPoint(p);
        return collision_detector->checkCollision(t_pgn1, pgn2);
    }

    bool checkSegmentCollisionDetection(Arrangement_history_2 &arr, Segment_2 &seg, const Polygon_2 &pgn1, const Polygon_2 &pgn2) const {
        Point_2 mid_point = CGAL::midpoint(seg.source(), seg.target());
        return checkCollisionDetection(arr, mid_point, pgn1, pgn2);
    }

    void removeFaceLoop(Arrangement_history_2 &arr, Halfedge_handle &handle) const {
        std::list<Halfedge_handle> remove_list;
        Ccb_halfedge_circulator circ = handle->ccb();
        remove_list.push_front(handle);
        ++circ;

        while (circ != handle) {
            remove_list.push_front(circ);
            ++circ;
        }

        for (typename std::list<Halfedge_handle>::iterator itr = remove_list.begin(); itr != remove_list.end(); ++itr) {
            arr.remove_edge(*itr);
        }
    }
    void buildArrangementFromConv(const Segments_list &reduced_conv, Arrangement_history_2 &arr) const {
        CGAL_precondition(arr.is_empty());
        CGAL::insert(arr, reduced_conv.begin(), reduced_conv.end());
    }

    void constructOrientableLoops(Arrangement_history_2 &arr) const {
        const Minkowski_sum_by_convolution_lien_2 *ptr = this;
        ConvSegMapper mapper(&arr, const_cast <Minkowski_sum_by_convolution_lien_2 *>(ptr));
        TraversalManager manager(&mapper);
        manager.traverseLoops();
    }

    bool checkDegenerateEdgeOppositeSegments(Arrangement_history_2 &arr, Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);

        std::list<Direction_2> segments_dir_list;

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            Segment_2 segment = *segment_itr;
            Direction_2 seg_dir = f_direction(f_vector(segment.source(), segment.target()));
            segments_dir_list.push_back(seg_dir);
        }

        segments_dir_list.sort();
        typename std::list<Direction_2>::iterator end = unique(segments_dir_list.begin(), segments_dir_list.end());
        int i = distance(segments_dir_list.begin(), end);
        return i > 1;
    }

    bool checkDegenarateVertexIsIntersectionOfThreeSegments(Arrangement_history_2 &arr, Vertex_handle vh) {
        Halfedge_around_vertex_circulator itr = vh->incident_halfedges();
        Halfedge_around_vertex_circulator start = itr;
        int count_degree = 0;

        do {
            ++count_degree;
        } while (++itr != start);

        if (count_degree <= 2) {
            return false;
        }

        Originating_curve_iterator segment_itr;
        std::list<Segment_2 *> orig_segments_list;
        std::list<Direction_2> segments_dir_list;

        // Handle the standard case where we have two intersecting convolution segments.
        if (count_degree == 4) {
            do {
                for (segment_itr = arr.originating_curves_begin(itr); segment_itr != arr.originating_curves_end(itr); ++segment_itr) {
                    Segment_2 segment = *segment_itr;
                    orig_segments_list.push_back(&(*segment_itr));
                }
            } while (++itr != start);

            orig_segments_list.sort();
            typename std::list<Segment_2 *>::iterator end = unique(orig_segments_list.begin(), orig_segments_list.end());
            int i = distance(orig_segments_list.begin(), end);

            if (i == 2) { // this is two curves crossing case.
                return false;
            }
        }

        do {

            for (segment_itr = arr.originating_curves_begin(itr); segment_itr != arr.originating_curves_end(itr); ++segment_itr) {
                Segment_2 segment = *segment_itr;
                Direction_2 seg_dir = f_direction(f_vector(segment.source(), segment.target()));
                segments_dir_list.push_back(seg_dir);
            }
        } while (++itr != start);

        segments_dir_list.sort();
        typename std::list<Direction_2>::iterator end = unique(segments_dir_list.begin(), segments_dir_list.end());
        int i = distance(segments_dir_list.begin(), end);
        return i > 2;
    }

    // Gets the he that agrees in direction with the convolution segment.
    Halfedge_handle getDirAgreeingHalfedge(Arrangement_history_2 &arr, const Halfedge_handle &he) const {
        Halfedge_handle curr_halfedge = he;

        if (!checkTripSameDirWithSegment(arr, he)) {
            curr_halfedge = curr_halfedge->twin();
        }

        return curr_halfedge;
    }

    // Gets list of incoming and outgoing edges(as defined by directions of segments in the convolution) from the vertex.
    void getEdgesFromVertex(Arrangement_history_2 &arr, Vertex_handle v_src, std::list<Halfedge_handle> &inList, std::list<Halfedge_handle> &outList) const {
        outList.clear();
        inList.clear();
        Halfedge_around_vertex_circulator itr = v_src->incident_halfedges();
        Halfedge_around_vertex_circulator start = itr;

        do {
            Halfedge_handle curr_edge = getDirAgreeingHalfedge(arr, itr);

            if ((curr_edge->source()) == v_src) {
                outList.push_back(curr_edge);
            } else {
                inList.push_back(curr_edge);
            }
        } while (++itr != start);
    }

    // Returns the direction of a half edge
    Direction_2 getHalfedgeDir(const Halfedge_handle &he) const {
        Direction_2 dir = f_direction((f_vector(he->source()->point(), he->target()->point())));
        return dir;
    }

    double getSignedAngle(const Halfedge_handle &h_enter, const Halfedge_handle &h_exit) const {
        Direction_2 dir_enter = getHalfedgeDir(h_enter);
        Direction_2 dir_exit = getHalfedgeDir(h_exit);
        Vector_2 vec_enter = dir_enter.vector();
        Vector_2 vec_exit = dir_exit.vector();
        Point_2 org(CGAL::ORIGIN);
        Vector_2 origin_vec(org, org);
        Orientation sign_or = f_orientation(vec_enter, vec_exit);
        float sign = 0.f;

        if (sign_or == CGAL::LEFT_TURN) {
            sign = 1;
        } else if (sign_or == CGAL::RIGHT_TURN) {
            sign = -1;
        } else {
            sign = 0;
        }

        double prod = CGAL::to_double(vec_enter * vec_exit);

        if (f_equal(vec_enter, origin_vec) || f_equal(vec_exit, origin_vec)) {
            return 0;
        }

        double len1 = sqrt(CGAL::to_double(vec_enter.squared_length()));
        double len2 = sqrt(CGAL::to_double(vec_exit.squared_length()));
        double p = prod / (len1 * len2);
        p = min((double)(1), p);
        p = max((double)(-1), p);
        double ang = acos(p);
        return sign * ang;
    }

    Halfedge_handle traverseNextHalfedge() {
    }

    void traceOrientableLoops(Arrangement_history_2 &arr, Edges_set &edges_set) const {
        int loop_counter = 0;
        std::list<Halfedge_iterator> temp_segments;

        double angles_sum = 0;

        // get an edge that is surly on the outside border. we will start from the external loop.
        Face_iterator startFace = arr.unbounded_face();
        Halfedge_iterator perimiterFace = *(startFace -> holes_begin());
        Halfedge_iterator curr_halfedge = perimiterFace;//(*edges_set.begin());

        if (!checkTripSameDirWithSegment(arr, curr_halfedge)) {
            curr_halfedge = curr_halfedge->twin();
        }

        Halfedge_iterator start_halfedge = curr_halfedge;
        temp_segments.push_back(curr_halfedge);
        setEdgeVisited(*curr_halfedge, true, loop_counter);
        Halfedge_handle close_loop_handle;
        Halfedge_handle before_close_loop_handle;
        bool isOrientable = true;
        bool close_loop_found = false;
        bool tracing_holes = false;

        while (isOrientable) {
            // after first loop we are starting to trace the holes.
            if (loop_counter > 0) {
                tracing_holes = true;
            }

            /////////////////////////DEBUG
            Point_2 p_source = curr_halfedge->source()->point();
            Point_2 p_end = curr_halfedge->target()->point();
            double x1 = CGAL::to_double(p_source.x());
            double y1 = CGAL::to_double(p_source.y());
            double x2 = CGAL::to_double(p_end.x());
            double y2 = CGAL::to_double(p_end.y());
            int size_edges_set = edges_set.size();

            /////////////////////////DEBUG
            Vertex_iterator v_target = curr_halfedge->target();
            p_source = v_target->point();
            double x = CGAL::to_double(p_source.x());
            double y = CGAL::to_double(p_source.y());
            Halfedge_around_vertex_circulator itr = v_target->incident_halfedges();
            bool next_edge_found;
            bool temp_close_loop_found = false;

            Halfedge_iterator next_halfedge = getLargestExitingClockwiseEdge(arr, curr_halfedge, *v_target, next_edge_found, loop_counter, temp_close_loop_found, close_loop_handle);

            if (!next_halfedge->visited && next_edge_found) {
                setEdgeVisited(*next_halfedge, true, loop_counter);
            }

            // if any ending of loop was found remember it.
            if (temp_close_loop_found) {
                before_close_loop_handle = curr_halfedge;
                close_loop_found = true;
            }

            // Check if we are stuck or we have met a visited edge which closes the loop:
            if (!next_edge_found) {
                Halfedge_handle temp_twin = next_halfedge->twin();

                if ((next_halfedge->loopNumber == loop_counter) && (next_halfedge != curr_halfedge) && (temp_twin != curr_halfedge)) {
                    // a loop closes with part of the edges. we have to remove edges till next_halfedge from arrangment,
                    // and clear the list.
                    removeHalfEdgesArrPart(arr, edges_set, temp_segments, next_halfedge);
                    //removeHalfEdgesArrPart(arr,edges_set,temp_segments,curr_halfedge);
                    nextLoop(temp_segments, arr, edges_set, curr_halfedge, start_halfedge, isOrientable, loop_counter);
                }

                else {
                    if (close_loop_found) {
                        removeHalfEdgesArrPartEnd(arr, edges_set, temp_segments, close_loop_handle, before_close_loop_handle);
                        nextLoop(temp_segments, arr, edges_set, curr_halfedge, start_halfedge, isOrientable, loop_counter);
                        close_loop_found = false;
                    } else {
                        // we are stuck.
                        removeHalfEdgesArr(arr, edges_set, temp_segments);
                        nextLoop(temp_segments, arr, edges_set, curr_halfedge, start_halfedge, isOrientable, loop_counter);
                    }
                }
            } else {
                if (close_loop_found) {
                    // We have closed the loop, attempted to expand it and got stuck so remove everything but the loop.
                    // THIS IS WRONG ??????????????????????????????????????????????????????????????????????????????????????????????

                    removeHalfEdgesArrPartEnd(arr, edges_set, temp_segments, close_loop_handle, before_close_loop_handle);
                    nextLoop(temp_segments, arr, edges_set, curr_halfedge, start_halfedge, isOrientable, loop_counter);
                    close_loop_found = false;
                } else {
                    curr_halfedge = next_halfedge;
                    temp_segments.push_back(curr_halfedge);
                    int size_temp_segments = temp_segments.size();

                    // check if loop has closed (?)
                    if (next_halfedge == start_halfedge) {
                        removeHalfedgesSet(arr, edges_set, temp_segments);
                        nextLoop(temp_segments, arr, edges_set, curr_halfedge, start_halfedge, isOrientable, loop_counter);
                    }
                }
            }
        }
    }

    void nextLoop(std::list<Halfedge_iterator> &temp_segments, Arrangement_history_2 &arr, Edges_set &edges_set, Halfedge_iterator &curr_halfedge, Halfedge_iterator &start_halfedge, bool &isOrientable, int &loop_counter) const {
        temp_segments.clear();

        if (edges_set.size() > 0) {
            ++loop_counter;
            curr_halfedge = (*edges_set.begin());

            while (curr_halfedge->visited && edges_set.size() > 0) {
                typename Edges_set::iterator testItem = edges_set.find(curr_halfedge);

                if (testItem == edges_set.end()) {
                    testItem = edges_set.find(((curr_halfedge)->twin()));

                    if (testItem != edges_set.end()) {
                        edges_set.erase(testItem);
                    }
                } else {
                    edges_set.erase(testItem);
                }

                curr_halfedge = (*edges_set.begin());
            }

            if (!checkTripSameDirWithSegment(arr, curr_halfedge)) {
                curr_halfedge = curr_halfedge->twin();
            }

            setEdgeVisited(*curr_halfedge, true, loop_counter);
            temp_segments.push_back(curr_halfedge);
            start_halfedge = curr_halfedge;
        } else {
            isOrientable = false;
        }
    }

    void setEdgeVisited(Halfedge &he, bool value, int id) const {
        he.visited = value;
        he.twin()->visited = value;
        he.loopNumber = id;
        he.twin()->loopNumber = id;
    }

    void setEdgeDegenerate(Halfedge &he, bool value) const {
        he.isDegenerate = value;
        he.twin()->isDegenerate = value;
    }

    bool getEdgeDegenerate(Halfedge_handle &he) const {
        return he.isDegenerate;
    }

    bool getEdgeVisited(Halfedge_handle &he) const {
        return he.visited;
    }

    // Removes the edges matching to halfedges from the set
    void removeHalfedgesSet(Arrangement_history_2 &arr, Edges_set &edges_set, std::list<Halfedge_iterator> &temp_segments) const {
        removeHalfEdgesInner(arr, edges_set, temp_segments, false, false, false, *(temp_segments.begin()), *(temp_segments.begin()));
    }

    // Removes the edges matching to halfedges from the set and arrangment
    void removeHalfEdgesArr(Arrangement_history_2 &arr, Edges_set &edges_set, std::list<Halfedge_iterator> &temp_segments) const {
        removeHalfEdgesInner(arr, edges_set, temp_segments, true, false, false, *(temp_segments.begin()), *(temp_segments.begin()));
    }

    // Remove edges before loop and keep loop
    void removeHalfEdgesArrPart(Arrangement_history_2 &arr, Edges_set &edges_set, std::list<Halfedge_iterator> &temp_segments, Halfedge_handle &partition) const {
        removeHalfEdgesInner(arr, edges_set, temp_segments, false, true, false, partition, *(temp_segments.begin()));
    }

    // Remove edges before loop and keep loop and remove edges after loop.
    void removeHalfEdgesArrPartEnd(Arrangement_history_2 &arr, Edges_set &edges_set, std::list<Halfedge_iterator> &temp_segments, Halfedge_handle &partition, Halfedge_handle &end_partition) const {
        removeHalfEdgesInner(arr, edges_set, temp_segments, false, true, true, partition, end_partition);
    }

    // remove half edges from set or set and arrangment.
    void removeHalfEdgesInner(Arrangement_history_2 &arr, Edges_set &edges_set, std::list<Halfedge_iterator> &temp_segments, bool remove_arr, bool range, bool end_range, Halfedge_handle &partition_itr, Halfedge_handle &partition_end_itr) const {
        bool part_reached = false;

        typename Edges_set::iterator partItem = edges_set.find(partition_itr);

        if (partItem == edges_set.end()) {
            partItem = edges_set.find(((*partition_itr).twin()));
        }

        typename Edges_set::iterator part_item_end = edges_set.find(partition_end_itr);

        if (part_item_end == edges_set.end()) {
            part_item_end = edges_set.find(((*partition_end_itr).twin()));
        }

        Halfedge_handle h_e_start;

        if (partItem != edges_set.end()) {
            h_e_start = *partItem;
        }

        Halfedge_handle h_e_finish;

        if (part_item_end != edges_set.end()) {
            h_e_finish = *part_item_end;
        }

        bool mark_change_part = false;

        for (typename std::list<Halfedge_iterator>::iterator itr = temp_segments.begin(); itr != temp_segments.end(); ++itr) {
            if (range && !part_reached && ((*itr == h_e_start) || (*itr == h_e_start->twin()))) {
                part_reached = true;
            }

            // If we have to close the loop check if range end has arrived
            if (range && end_range && part_reached && (((*itr) == h_e_finish) || (*itr == h_e_finish->twin()))) {
                mark_change_part = true;
            }

            typename Edges_set::iterator testItem = edges_set.find(*itr);

            if (testItem == edges_set.end()) {
                testItem = edges_set.find(((*itr)->twin()));

                if (testItem != edges_set.end()) {
                    edges_set.erase(testItem);

                    if (remove_arr || (range && !part_reached)) {
                        arr.remove_edge(*itr);
                    }
                }
            } else {
                edges_set.erase(testItem);

                if (remove_arr || (range && !part_reached)) {
                    arr.remove_edge(*itr);
                }
            }

            if (mark_change_part) {
                part_reached = false;
            }
        }
    }

    Halfedge_handle getLargestExitingClockwiseEdge(Arrangement_history_2 &arr, Halfedge_iterator &curr_halfedge, Vertex &v_target, bool &next_edge_found, int loop_number, bool &close_loop_found, Halfedge_handle &handle_close_loop) const {
        Point_2 p_source = v_target.point();
        double x1 = CGAL::to_double(p_source.x());
        double y1 = CGAL::to_double(p_source.y());
        Halfedge_around_vertex_circulator itr = v_target.incident_halfedges();
        Halfedge_around_vertex_circulator start = itr;
        Halfedge_around_vertex_circulator tmp = itr;
        ++itr;

        if (itr == tmp) {
            next_edge_found = false;
            return itr;
        }

        bool found = false;
        Direction_2 entering_dir = f_direction((f_vector(curr_halfedge->source()->point(), curr_halfedge->target()->point())));

        // find first edge which is valid in direction with respect to original edges.
        //
        int count = 0;
        bool res;

        while ((!(res = checkTripSameDirWithSegment(arr, ((itr->twin())))) || (curr_halfedge == itr) || (itr->visited)) && (itr != start)) {
            ++itr;
            ++ count;
        }

        Halfedge_handle maybe_visited;

        if (count == 1) {

            close_loop_found = false;

            if (!checkOutgoingNotVisited(arr, v_target, maybe_visited, loop_number)) {
                close_loop_found = true;
                handle_close_loop = maybe_visited;
            }
        }

        if (!res || itr->visited) {
            next_edge_found = false;
            return itr->twin();
        }

        if (count != 1) {
            close_loop_found = false;

            if (!checkOutgoingNotVisited(arr, v_target, maybe_visited, loop_number)) {
                close_loop_found = true;
                handle_close_loop = maybe_visited;
            }
        }

        // Now we have the first edge, check if we can improve
        Halfedge_around_vertex_circulator next_edge_itr = itr;
        ++next_edge_itr;
        count = 0;
        Direction_2 min_edge_dir = f_direction((f_vector(itr->twin()->source()->point(), itr->twin()->target()->point())));

        while (next_edge_itr != itr) {
            if (checkTripSameDirWithSegment(arr, ((next_edge_itr->twin())))) {
                Direction_2 new_edge_dir = f_direction((f_vector(next_edge_itr->twin()->source()->point(), next_edge_itr->twin()->target()->point())));

                // If the new edge improves the old one, ie it is larger and satisfies the improvment rule, which is under consideration.
                // currently, being the anticlocwise most if we are in right halfplane of the entering direction, and then choose anti clockwise most from left plane.
                if (isDirImproving(min_edge_dir, entering_dir, new_edge_dir) &&
                        (!next_edge_itr->visited)
                   ) {
                    itr = next_edge_itr;
                    min_edge_dir = f_direction((f_vector(itr->twin()->source()->point(), itr->twin()->target()->point())));
                }
            }

            ++count;
            ++next_edge_itr;
        }

        if (close_loop_found) {
            Direction_2 new_edge_dir = f_direction((f_vector(maybe_visited->twin()->source()->point(), maybe_visited->twin()->target()->point())));

            if (isDirImproving(min_edge_dir, entering_dir, new_edge_dir)) {
                next_edge_found = true;
                return maybe_visited->twin();
            }
        }

        next_edge_found = true;
        return itr->twin();
    }

    bool isDirImproving(Direction_2 &min_edge_dir, Direction_2 &entering_dir, Direction_2 &new_edge_dir) const {
        Direction_2 opp_enter = -entering_dir;

        if ((opp_enter == min_edge_dir)) { // if minimal dir equals -entering dir
            return true;
        } else {
            return f_ccw_in_between(new_edge_dir, opp_enter, min_edge_dir);
        }
    }

    /*
        Checks that the edge leads to a vertex which an outgoing visited edge has been ie returns false if we close a loop.
        h returns the edge which begins the loop
    */
    bool checkOutgoingNotVisited(Arrangement_history_2 &arr, Vertex &v_target, Halfedge_handle &h, int loop_number) const {
        Point_2 p_source = v_target.point();

        Halfedge_around_vertex_circulator itr = v_target.incident_halfedges();
        Halfedge_around_vertex_circulator start = itr;

        do {
            if (checkTripSameDirWithSegment(arr, ((itr->twin()))) && itr->visited && itr->loopNumber == loop_number) {
                h = itr;
                return false;
            }
        } while (++itr != start);

        return true;
    }

    bool checkTripSameDirWithSegment_bak(Arrangement_history_2 &arr, Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;// = arr.originating_curves_begin ( *he);
        Direction_2 start_he_dir = f_direction((f_vector(he->source()->point(), he->target()->point())));

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            Segment_2 segment = *segment_itr;
            Direction_2 start_seg_dir = f_direction(f_vector(segment.source(), segment.target()));

            if (start_seg_dir == start_he_dir) {
                return true;
            }
        }

        return false;
    }

    bool checkTripSameDirWithSegment(Arrangement_history_2 &arr, Halfedge_handle he) const {
        Originating_curve_iterator segment_itr;

        for (segment_itr = arr.originating_curves_begin(he); segment_itr != arr.originating_curves_end(he); ++segment_itr) {
            Segment_2 segment = *segment_itr;

            CGAL::Comparison_result c1 = f_compare_endpoints_xy(segment);
            CGAL::Comparison_result c2 = (CGAL::Comparison_result)he->direction();
            bool same_dir = (c1 == c2);

            if (same_dir) {
                return true;
            }
        }

        return false;
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

    bool checkSwept(const Point_2 &prev, const Point_2 &curr, const Point_2 &next, const Point_2 &start, const Point_2 &end, bool &isStartConcide, bool &isEndConcide) const {
        Direction_2 dir_start = f_direction(f_vector(prev, curr));
        Direction_2 dir_end = f_direction(f_vector(curr, next));
        Direction_2 dir_new = f_direction(f_vector(start, end));
        isStartConcide = dir_new == dir_start;
        isEndConcide = dir_end == dir_new;

        return isStartConcide || f_ccw_in_between(dir_new, dir_start, dir_end) || isEndConcide;
    }

    bool checkSwept(Direction_2 &dir_start, Direction_2 &dir_end, Direction_2 &dir_new, bool &isStartConcide, bool &isEndConcide) const {
        isStartConcide = dir_new == dir_start;
        isEndConcide = dir_end == dir_new;

        return isStartConcide || f_ccw_in_between(dir_new, dir_start, dir_end) || isEndConcide;
    }
};

}

#endif
