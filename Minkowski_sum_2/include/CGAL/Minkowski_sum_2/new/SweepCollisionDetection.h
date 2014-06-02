#ifndef SWEEPCOLLISIONDETECTOR_HEADER
#define SWEEPCOLLISIONDETECTOR_HEADER
#include "ICollisionDetector.h"
#include <CGAL/intersections.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_default_overlay_traits_base.h>
#include <CGAL/Object.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_2_utils.h>
#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>

#include <list>
#include <algorithm>

/* A simple sweep-line visitor that determines if there are intersections
* in the interiors of the given curve set.
*/
namespace CGAL {

enum Segment_color {MY_RED, MY_BLUE};
struct Segment_Data {
    Segment_Data(Segment_color color, int min_id, int max_id): _color(color), _min_id(min_id), _max_id(max_id) {}
    Segment_color _color;
    int _min_id, _max_id;
    bool operator==(const Segment_Data &rhs) const {
        //return _color == rhs._color;
        return (this == &rhs);
    }
};

template <class Traits_> class Sweep_line_do_curves_x_visitor_ : public Sweep_line_empty_visitor<Traits_> {
    typedef Traits_ Traits_2;
    typedef Sweep_line_do_curves_x_visitor_<Traits_2> Self;

    typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits_2::Point_2 Point_2;

    typedef Sweep_line_empty_visitor<Traits_2> Base;
    typedef typename Base::Event Event;
    typedef typename Base::Subcurve Subcurve;
    typedef typename Base::Status_line_iterator Status_line_iterator;

    typedef CGAL::Sweep_line_2<Traits_2, Self> Sweep_line_2;

protected:

    // Data members:
    bool m_found_x; // Have we found an intersection so far.
    bool m_had_overlap_no_cross;
public:

    Sweep_line_do_curves_x_visitor_() :
        m_found_x(false), m_had_overlap_no_cross(false) {
    }

    template <class CurveIterator>
    void sweep(CurveIterator begin, CurveIterator end) {
        std::vector<X_monotone_curve_2> curves_vec;
        std::vector<Point_2> points_vec;

        curves_vec.reserve(std::distance(begin, end));
        make_x_monotone(begin,
                        end,
                        std::back_inserter(curves_vec),
                        std::back_inserter(points_vec),
                        this-> traits());

        // Perform the sweep.
        Sweep_line_2 *sl = reinterpret_cast<Sweep_line_2 *>(this->sweep_line());

        sl->sweep(curves_vec.begin(),
                  curves_vec.end(),
                  points_vec.begin(),
                  points_vec.end());
    }

    void update_event(Event *e ,
                      Subcurve *sc1 ,
                      Subcurve *sc2 ,
                      bool is_new) {
        //m_found_x = true;
    }

    void update_event(Event * /* e */,
                      Subcurve * /* sc1 */) {
        // m_found_x = true;
    }

    void update_event(Event * /* e */,
                      const Point_2 & /* end_point */,
                      const X_monotone_curve_2 & /* cv */,
                      Arr_curve_end /* cv_end */,
                      bool /* is_new */) {
    }

    void update_event(Event * /* e */,
                      const Point_2 & /* pt */,
                      bool /* is_new */) {
    }

    template <class XCurveIterator>
    void sweep_xcurves(XCurveIterator begin, XCurveIterator end) {
        // Perform the sweep.
        Sweep_line_2 *sl = reinterpret_cast<Sweep_line_2 *>(this->sweep_line());

        sl->sweep(begin, end);
    }

    void found_overlap(Subcurve *sc1 ,
                       Subcurve *sc2 ,
                       Subcurve *ov_sc) {

        if (sc1->last_curve().data().front()._color != sc2->last_curve().data().front()._color) {
            Sweep_line_2 *sl = reinterpret_cast<Sweep_line_2 *>(this->sweep_line());
            const Traits_2 *traits = sl->traits();

            typename Traits_2::Compare_endpoints_xy_2 t_compare_endpoints_xy_2_obj = traits->compare_endpoints_xy_2_object();
            CGAL::Comparison_result c1 = t_compare_endpoints_xy_2_obj(sc1->last_curve());
            CGAL::Comparison_result c2 = t_compare_endpoints_xy_2_obj(sc2->last_curve());
            bool same_dir = (c1 == c2);
            m_found_x = same_dir;
            m_had_overlap_no_cross = !same_dir;
        }


        //m_found_x = true;

    }

    bool after_handle_event(Event *event ,
                            Status_line_iterator iter,
                            bool flag) {
        Sweep_line_2 *sl = reinterpret_cast<Sweep_line_2 *>(this->sweep_line());

        if (m_found_x) {
            sl->stop_sweep();
            return true;
        }

        // check if there was an intersection event:
        if (((event->is_intersection() || event->is_weak_intersection() || (event->is_left_end() && event->is_right_end()))) && (event->number_of_left_curves() + event->number_of_right_curves() == 4)) {
            Sweep_line_2 *sl = reinterpret_cast<Sweep_line_2 *>(this->sweep_line());
            const Traits_2 *traits = sl->traits();
            typename Traits_2::Compare_endpoints_xy_2 t_compare_endpoints_xy_2_obj = traits->compare_endpoints_xy_2_object();

            // get all curves ordered by cyclic order. from left top counter clockwise.
            std::list<Subcurve *> ordered_list;
            ordered_list.insert(ordered_list.begin(), event->left_curves_begin(), event->left_curves_end());
            ordered_list.insert(ordered_list.begin(), event->right_curves_rbegin(), event->right_curves_rend());

            // mark for each edge whether it's incoming or outgoing.
            std::vector<bool> incoming_edges(ordered_list.size(), false);
            typename std::list<Subcurve *>::iterator itr = ordered_list.begin();
            int i = 0;

            for (; itr != ordered_list.end(); ++itr) {
                if (i < event->number_of_left_curves()) {
                    if (t_compare_endpoints_xy_2_obj((*itr)->last_curve()) == CGAL::SMALLER) {
                        incoming_edges[i] = true;
                    }
                } else {
                    if (t_compare_endpoints_xy_2_obj((*itr)->last_curve()) == CGAL::LARGER) {
                        incoming_edges[i] = true;
                    }
                }

                ++i;
            }

            if (ordered_list.size() == 4) {
                // normal intersection case
                // check for alterations
                typename std::list<Subcurve *>::iterator itr1, itr2, itr3;
                itr1 = ordered_list.begin();
                itr2 = itr1;
                ++itr2;
                itr3 = itr2;
                ++itr3;
                bool c = ((*itr1)->last_curve().data().front()._color != (*itr2)->last_curve().data().front()._color);

                if (((*itr1)->last_curve().data().front()._color != (*itr2)->last_curve().data().front()._color) &&
                        ((*itr2)->last_curve().data().front()._color != (*itr3)->last_curve().data().front()._color)) {
                    // we have alternating edges
                    m_found_x = true;
                    sl->stop_sweep();
                    return true;
                } else {
                    // either 1 and 2 are the same or 1 and for are the same colors.
                    //Subcurve* r1,r2,b1,b2;
                    if ((*itr1)->last_curve().data().front()._color == (*itr2)->last_curve().data().front()._color) {
                        // 1==2
                        /* std::list<Subcurve>::iterator itr =ordered_list.front();
                              r1 = itr; r2 = ++itr; b1 = ++itr; b2 = ++itr;*/
                        if (incoming_edges[0] || incoming_edges[2]) {
                            m_found_x = true;
                            sl->stop_sweep();
                            return true;
                        }
                    } else {// 1 == 4
                        /* std::list<Subcurve>::iterator itr =ordered_list.front();
                          r2 = itr; b1 = ++itr; b2 = ++itr; r1 = ++itr;*/
                        if (incoming_edges[1] || incoming_edges[3]) {
                            m_found_x = true;
                            sl->stop_sweep();
                            return true;
                        }
                    }


                }
            } else {
                // Maybe a bug here when we have 2 overlap not overlapping and two outgoing which cause overlap.
                // this case needs to be studied and handled here.
                //       /
                //   ----
                //       \
                //
                /* int k=7;
                      ++k;
                      */
            }

        }


        return true;
    }

    bool found_intersection() {
        return (m_found_x);
    }

    bool had_overlap_no_cross() {
        return m_had_overlap_no_cross;
    }
};

template <class Traits_> class Colored_traits : public Arr_consolidated_curve_data_traits_2<Traits_, Segment_Data> {

public:
    typedef Arr_consolidated_curve_data_traits_2<Traits_, Segment_Data> Base;
    typedef typename Base::Intersect_2 Base_intersect_2;
    //typedef CGAL::Arr_consolidated_curve_data_traits_2<Traits_,Segment_color> Data_traits_2;
    typedef typename Colored_traits::Curve_2 Colored_segment_2;
    typedef typename Colored_traits::X_monotone_curve_2 X_monotone_colored_segment_2;
    typedef typename Base::Compare_xy_2 Base_compare_xy_2;
    typedef typename Base::Construct_min_vertex_2 Base_construct_min_vertex_2;
    typedef typename Base::Construct_max_vertex_2 Base_construct_max_vertex_2;
    typedef typename Base::Point_2 Base_point_2;

    class Intersect_2 {
    protected:
        //! The base traits.
        const Arr_consolidated_curve_data_traits_2<Traits_, Segment_Data> *m_traits;

        /*! Constructor.
         * The constructor is declared protected to allow only the functor
         * obtaining function, which is a member of the nesting class,
         * constructing it.
         */
        Intersect_2(const Arr_consolidated_curve_data_traits_2<Traits_, Segment_Data> *traits) : m_traits(traits) {}

        //! Allow its functor obtaining function calling the protected constructor.
        friend class Colored_traits<Traits_>;

    public:
        template<class OutputIterator>
        OutputIterator operator()(const X_monotone_colored_segment_2 &xcv1,
                                  const X_monotone_colored_segment_2 &xcv2,
                                  OutputIterator oi) {
            // In case the curves originate from the same arrangement, they are
            // obviously interior-disjoint.
            //if (xcv1.data().front() == xcv2.data().front())
            if (xcv1.data().front()._color == xcv2.data().front()._color) {
                return (oi);
            }

            typename Base::Intersect_2 int_obj = typename Base::Intersect_2(m_traits);
            //OutputIterator oi_temp;
            OutputIterator oi_end = int_obj(xcv1, xcv2, oi);

            if (oi == oi_end) {
                return oi_end;
            }

            const std::pair<Base_point_2, unsigned int> *xp_point;

            for (; oi != oi_end; ++oi) {
                xp_point = object_cast<std::pair<Base_point_2, unsigned int> > (&(*oi));

                if (xp_point != NULL) {
                    *oi = CGAL::make_object(std::make_pair(Point_2(xp_point->first),
                                                           xp_point->second));
                }
            }

            //if (xcv1.color() == RB_OVERLAP || xcv2.color() == RB_OVERLAP)
            return (oi);
        }

    };

    /*! Obtain an Intersect_2 functor object. */
    Intersect_2 intersect_2_object() const {
        return Intersect_2(this);
    }

    class Ex_point_2 {
    public:

        typedef Base_point_2 Base_p;

    protected:

        Base_p m_base_pt; // The base point.

    public:
        int id;
        /*! Default constructor. */
        Ex_point_2() :
            m_base_pt(), id(-1)

        {}

        /*! Constructor from a base point. */
        Ex_point_2(const Base_p &pt) :
            m_base_pt(pt), id(-1)

        {}


        /*! Get the base point (const version). */
        const Base_p &base() const {
            return (m_base_pt);
        }

        /*! Get the base point (non-const version). */
        Base_p &base() {
            return (m_base_pt);
        }

        /*! Casting to a base point (const version). */
        operator const Base_p &() const {
            return (m_base_pt);
        }

        /*! Casting to a base point (non-const version). */
        operator Base_p &() {
            return (m_base_pt);
        }


    };

    typedef Ex_point_2 Point_2;


    //class Ex_point_2 : public Point_2
    //{
    //public:
    // Ex_point_2(Point_2& p)
    // {

    // }
    // int id;
    //};

    class Compare_xy_2 {
    protected:
        //! The base operator.
        Base_compare_xy_2 m_base_cmp_xy;

        /*! Constructor.
         * The constructor is declared protected to allow only the functor
         * obtaining function, which is a member of the nesting class,
         * constructing it.
         */
        Compare_xy_2(const Base_compare_xy_2 &base) :
            m_base_cmp_xy(base) {
        }

        //! Allow its functor obtaining function calling the protected constructor.
        friend class Colored_traits<Traits_>;

    public:
        Comparison_result operator()(const Point_2 &p1, const Point_2 &p2) const {
//      if (&p1 == &p2)
            //Ex_point_2 p1_e,p2_e;
            //p1_e = *dynamic_cast<const Ex_point_2 *>(&p1);
            //p2_e = *dynamic_cast<const Ex_point_2 *>(&p2);
            if ((p1.id == p2.id) && (p1.id != -1)) {
                return EQUAL;
            }

            return m_base_cmp_xy(p1, p2);
        }
    };

    /*! Obtain a Construct_min_vertex_2 functor object. */
    Compare_xy_2 compare_xy_2_object() const {
        //Base::Compare_xy_2 obj();
        return (Compare_xy_2(((Base *)this)->compare_xy_2_object()));
    }



    /*! A functor that obtains the left endpoint of an x-monotone curve. */
    class Construct_min_vertex_2 {
    protected:
        //! The base operators.
        Base_construct_min_vertex_2 m_base_min_v;
        //Base_equal_2 m_base_equal;

        /*! Constructor.
         * The constructor is declared protected to allow only the functor
         * obtaining function, which is a member of the nesting class,
         * constructing it.
         */
        Construct_min_vertex_2(const Base_construct_min_vertex_2 &base_min_v
                              ) :
            m_base_min_v(base_min_v) {
        }

        //! Allow its functor obtaining function calling the protected constructor.
        friend class Colored_traits<Traits_>;
        typedef typename Colored_traits::X_monotone_curve_2 X_monotone_curve_2;

    public:
        Point_2 operator()(const X_monotone_curve_2 &xcv) {
            Point_2 min_p = m_base_min_v(xcv);

            if (xcv.data().size() > 1) {
                min_p.id = -1;
            } else {
                min_p.id = xcv.data().front()._min_id;
            }

            return (min_p);
        }
    };

    /*! Obtain a Construct_min_vertex_2 functor object. */
    Construct_min_vertex_2 construct_min_vertex_2_object() const {
        return
            (Construct_min_vertex_2(((Base *)this)->construct_min_vertex_2_object()));
    }

    /*! A functor that obtains the right endpoint of an x-monotone curve. */
    class Construct_max_vertex_2 {
    protected:
        //! The base operators.
        Base_construct_max_vertex_2 m_base_max_v;


        /*! Constructor.
         * The constructor is declared protected to allow only the functor
         * obtaining function, which is a member of the nesting class,
         * constructing it.
         */
        Construct_max_vertex_2(const Base_construct_max_vertex_2 &base_max_v) :
            m_base_max_v(base_max_v)

        {}

        //! Allow its functor obtaining function calling the protected constructor.
        friend class Colored_traits<Traits_>;
        typedef typename Colored_traits::X_monotone_curve_2 X_monotone_curve_2;

    public:
        Point_2 operator()(const X_monotone_curve_2 &xcv) const {
            Point_2 max_p = m_base_max_v(xcv);

            if (xcv.data().size() > 1) {
                max_p.id = -1;
            } else {
                max_p.id = xcv.data().front()._max_id;
            }

            return (max_p);
        }
    };

    /*! Obtain a Construct_min_vertex_2 functor object. */
    Construct_max_vertex_2 construct_max_vertex_2_object() const {
        return
            (Construct_max_vertex_2(((Base *)this)->construct_max_vertex_2_object()));
    }



};


template <class Kernel_, class Container_> class SweepCollisionDetector : public ICollisionDetector< Kernel_, Container_> {



public:


    SweepCollisionDetector() {}
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel_> Traits_2;
    typedef Colored_traits<Traits_2> Data_traits_2;
    //typedef CGAL::Arr_consolidated_curve_data_traits_2<Traits_2,Segment_color> Data_traits_2;

    typedef typename Data_traits_2::Curve_2 Colored_segment_2;
    typedef Arrangement_2<Data_traits_2> Colored_arr_2;
    typedef typename CGAL::Polygon_2<Kernel_>::Edge_const_iterator Edge_iterator ;
    typedef typename CGAL::Polygon_2<Kernel_>::Traits::Segment_2 Segment_2 ;
    //typedef typename Polygon_2::Vertex_circulator Vertex_circulator;
    //typedef typename

protected:
    Traits_2 m_traits;


public:

    virtual bool checkCollision(const typename CGAL::Polygon_2<Kernel_> &p, const typename CGAL::Polygon_2<Kernel_> &q) {
        typename Traits_2::Compare_endpoints_xy_2 cmp_obj = m_traits.compare_endpoints_xy_2_object();

        if (p.has_on_bounded_side(*(q.vertices_begin())) || q.has_on_bounded_side(*(p.vertices_begin()))) {
            return true;
        }

        //std::list<Segment_2> edges(p.edges_begin(),p.edges_end());
        std::list<Colored_segment_2> edges;
        Edge_iterator itr = p.edges_begin();
        int i = 0;
        int n = p.size();

        for (; itr != p.edges_end(); ++itr, ++i) {
            //Traits_2::

            //CGAL::Comparison_result res = CGAL::compare((itr)->source(),(itr)->target());
            CGAL::Comparison_result res = cmp_obj(*itr);

            if (res != CGAL::SMALLER) {
                if (i == n - 1) {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_RED, 0, i)));
                } else {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_RED, i + 1, i)));
                }
            } else {
                if (i == n - 1) {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_RED, i, 0)));
                } else {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_RED, i, i + 1)));
                }
            }

        }

        itr = q.edges_begin();
        n = q.size();
        int j = 0;
        int q_first = i;

        for (; itr != q.edges_end(); ++itr, ++i, ++j) {
            //CGAL::Comparison_result res = CGAL::compare((itr)->source(),(itr)->target());
            CGAL::Comparison_result res = cmp_obj(*itr);

            if (res != CGAL::SMALLER) {
                if (j == n - 1) {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_BLUE, q_first, i)));
                } else {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_BLUE, i + 1, i)));
                }
            } else {
                if (j == n - 1) {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_BLUE, i, q_first)));
                } else {
                    edges.push_back(Colored_segment_2(*itr, Segment_Data(MY_BLUE, i, i + 1)));
                }
            }


            //edges.push_back(Colored_segment_2(*itr,Segment_Data(MY_BLUE,2*i,2*i+1)));
        }

        //edges.insert(edges.begin(),q.edges_begin(),q.edges_end());
        CGAL::Sweep_line_do_curves_x_visitor_<Data_traits_2> visitor;
        Sweep_line_2<Data_traits_2, Sweep_line_do_curves_x_visitor_<Data_traits_2> > sweep_line(&visitor);
        visitor.attach((void *)&sweep_line);
        visitor.sweep(edges.begin(), edges.end());

        if (visitor.found_intersection()) {
            return true;
        }

        if (visitor.had_overlap_no_cross()) {
            return false;
        }

        //else
        //{
        //
        //  //if (CGAL::do_intersect(p,q))
        //  //  return true;
        //  //else return false;
        //
        //}
        return false;
        /*bool intersect = false;
        if (CGAL::do_intersect(p,q))
            return true;
        else return false;*/

        /*
        for (Edge_iterator itr1 = p.edges_begin();itr1!=p.edges_end();++itr1)
        {
            for (Edge_iterator itr2 = q.edges_begin();itr2!=q.edges_end();++itr2)
            {
                //intersect = intersect CGAL::Do_intersect(*itr1,*itr2);
                 CGAL::Object result = CGAL::intersection(*itr1,*itr2);
                 if (const CGAL::Point_2<Kernel> *ipoint = CGAL::object_cast<CGAL::Point_2<Kernel> >(&result)) {
                    // handle the point intersection case with *ipoint.
                    // check if the intersection point is the vertex of one of the edges
                     if (((*ipoint)==itr1->source()) || ((*ipoint)==itr1->target())
                         ||((*ipoint)==itr2->source()) || ((*ipoint)==itr2->target()))
                        intersect = false;
                     else
                        intersect = true;
                 } else
                if (const CGAL::Segment_2<Kernel> *iseg = CGAL::object_cast<CGAL::Segment_2<Kernel> >(&result)) {
                    // handle the segment intersection case with *iseg. this is the case where the edges are touching
                    intersect = false;
                } else {
                // handle the no intersection case.
                }

                if (intersect)
                    return intersect;
            }

        }

        if (p.bounded_side(*q.vertices_begin())== CGAL::ON_BOUNDED_SIDE || q.bounded_side(*p.vertices_begin())== CGAL::ON_BOUNDED_SIDE)
            return true;
        return intersect;
        */

    }






};

}
#endif
