#ifndef ARR_SEGMENTDATA_HEADER
#define ARR_SEGMENTDATA_HEADER

namespace CGAL {
namespace internal {

typedef std::pair<int, int> state;

struct Segment_Data_Label {
    Segment_Data_Label(state min_id, state max_id, CGAL::Comparison_result orientation, int origin): _min_id(min_id), _max_id(max_id), _orientation(orientation), _origin(origin) {}
    Segment_Data_Label() {
        _min_id = state(-1, -1);
        _max_id = state(-1, -1);
        _orientation = CGAL::SMALLER;
        _origin = -1;
    }

    state _min_id, _max_id;
    CGAL::Comparison_result _orientation;
    int _origin;

    bool operator==(const Segment_Data_Label &rhs) const {
        return (this == &rhs);
    }
};

template <class Traits_> class Arr_SegmentData_traits : public Traits_ {

public:

    typedef Traits_ Base;
    typedef typename Base::Intersect_2 Base_intersect_2;
    typedef typename Base::Split_2 Base_split_2;
    typedef typename Base::Compare_xy_2 Base_compare_xy_2;
    typedef typename Base::Construct_min_vertex_2 Base_construct_min_vertex_2;
    typedef typename Base::Construct_max_vertex_2 Base_construct_max_vertex_2;
    typedef typename Base::Point_2 Base_point_2;

    typedef typename Base::Curve_2 Base_Curve_2;
    typedef typename Base::X_monotone_curve_2 Base_x_monotone_curve_2;
    typedef typename Base::Compare_y_at_x_2 Compare_y_at_x_2;
    typedef typename Base::Compare_x_2 Compare_x_2;
    typedef typename Base::Compare_endpoints_xy_2 Compare_endpoints_xy_2;

    typedef typename Base::Has_left_category Has_left_category;
    typedef Tag_false Has_merge_category;

    class Ex_point_2 {

    public:

        typedef Base_point_2 Base_p;

    protected:

        Base_p m_base_pt; // The base point.

    public:
        state id;
        /*! Default constructor. */
        Ex_point_2() :
            m_base_pt(), id(state(-1, -1))

        {}

        /*! Constructor from a base point. */
        Ex_point_2(const Base_p &pt) :
            m_base_pt(pt), id(state(-1, -1))

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

    class X_monotone_curve_2 : public Base_x_monotone_curve_2 {

    private:

        Segment_Data_Label _label;

    public:

        /*! Default constructor. */
        X_monotone_curve_2() {
        }

        /*! Constructor from a base x-monotone curve. */
        X_monotone_curve_2(const Base_x_monotone_curve_2 &p) :
            Base_x_monotone_curve_2(p),
            _label() {
        }

        /*! Constructor from an x-monotone curve an a label. */
        X_monotone_curve_2(const Base_x_monotone_curve_2 &p,
                           const Segment_Data_Label &label) :
            Base_x_monotone_curve_2(p),
            _label(label) {
        }

        /*! Get the label (const version). */
        const Segment_Data_Label &label() const {
            return (_label);
        }

        const Segment_Data_Label &data() const {
            return (_label);
        }

        /*! Get the label (non-const version). */
        Segment_Data_Label &label() {
            return (_label);
        }

        Segment_Data_Label &data() {
            return (_label);
        }

        /*! Set the label. */
        void set_label(const Segment_Data_Label &label) {
            _label = label;
            return;
        }
    };

    typedef X_monotone_curve_2 Curve_2;

    class Make_x_monotone_2 {

    public:

        /*!
         * Cut the given curve into x-monotone subcurves and insert them into the
         * given output iterator. As segments are always x_monotone, only one
         * object will be contained in the iterator.
         * \param cv The curve.
         * \param oi The output iterator, whose value-type is Object.
         * \return The past-the-end iterator.
         */
        template<class OutputIterator>
        OutputIterator operator()(const Curve_2 &cv, OutputIterator oi) const {
            // Wrap the segment with an object.
            *oi = make_object(cv);
            ++oi;
            return (oi);
        }
    };

    /*! Get a Make_x_monotone_2 functor object. */
    Make_x_monotone_2 make_x_monotone_2_object() const {
        return Make_x_monotone_2();
    }

    class Split_2 {

    private:

        const Traits_ *m_traits;

    public:

        /*! Constructor. */
        Split_2(const Base *_base) :
            m_traits(_base) {
        }

        /*!
         * Split a given x-monotone curve at a given point into two sub-curves.
         */
        void operator()(const X_monotone_curve_2 &cv, const Point_2 &p,
                        X_monotone_curve_2 &c1, X_monotone_curve_2 &c2) const {
            // Split the base curve into two.
            m_traits->split_2_object()(cv, p, c1, c2);

            if (cv.data()._orientation == CGAL::SMALLER) {
                c1.data()._min_id = cv.data()._min_id;
                c2.data()._max_id = cv.data()._max_id;
            } else {
                c1.data()._min_id = cv.data()._max_id;
                c2.data()._max_id = cv.data()._min_id;
            }

            c1.data()._orientation = cv.data()._orientation;
            c2.data()._orientation = cv.data()._orientation;

            return;
        }
    };

    /*! Get a Split_2 functor object. */
    Split_2 split_2_object() const {
        return (Split_2(this));
    }

    class Intersect_2 {

    protected:
        //! The base traits.
        const Traits_ *m_traits;

        /*! Constructor.
         * The constructor is declared protected to allow only the functor
         * obtaining function, which is a member of the nesting class,
         * constructing it.
         */
        Intersect_2(const Traits_* traits) : m_traits(traits) {}

        //! Allow its functor obtaining function calling the protected constructor.
        friend class Arr_SegmentData_traits<Traits_>;

    public:

        template<class OutputIterator>
        OutputIterator operator()(const X_monotone_curve_2 &xcv1,
                                  const X_monotone_curve_2 &xcv2,
                                  OutputIterator oi) {
            state no_state = state(-1, -1);

            // In case the curves originate from the same arrangement, they are
            // obviously interior-disjoint.
            //if (xcv1.data().front() == xcv2.data().front())
            //  if (xcv1.data().front()._color == xcv2.data().front()._color)
            //   return (oi);
            //Arr_SegmentData_traits<Traits_>::count++;

            // Ignore edges which originates from the same polygon and are added to the same vertex.
            if (xcv1.data()._origin == xcv2.data()._origin) {
                if ((xcv1.data()._origin == 0 && ((xcv1.data()._min_id.second == xcv2.data()._min_id.second)) ||
                        (xcv1.data()._origin == 1 && ((xcv1.data()._min_id.first == xcv2.data()._min_id.first))))) {
                    return oi;
                }
            }

            Construct_min_vertex_2 min_vertex_obj = m_traits->construct_min_vertex_2_object();
            Construct_max_vertex_2 max_vertex_obj = m_traits->construct_max_vertex_2_object();

            int endp_coll = 0;

            if ((xcv1.data()._min_id == xcv2.data()._min_id) && (xcv2.data()._min_id != no_state)) {
                return oi;

                if (xcv1.data()._orientation == CGAL::SMALLER) {
                    *oi = CGAL::make_object(std::make_pair(min_vertex_obj(xcv1), 1));
                } else {
                    *oi = CGAL::make_object(std::make_pair(max_vertex_obj(xcv1), 1));
                }

                ++oi;
                return oi;
            }

            if ((xcv1.data()._min_id == xcv2.data()._max_id) && (xcv2.data()._max_id != no_state)) {
                return oi;

                if (xcv1.data()._orientation == CGAL::SMALLER) {
                    *oi = CGAL::make_object(std::make_pair(min_vertex_obj(xcv1), 1));
                } else {
                    *oi = CGAL::make_object(std::make_pair(max_vertex_obj(xcv1), 1));
                }

                ++oi;
                return oi;
            }

            if ((xcv1.data()._max_id == xcv2.data()._min_id) && (xcv2.data()._min_id != no_state)) {
                return oi;

                if (xcv1.data()._orientation == CGAL::SMALLER) {
                    *oi = CGAL::make_object(std::make_pair(max_vertex_obj(xcv1), 1));
                } else {
                    *oi = CGAL::make_object(std::make_pair(min_vertex_obj(xcv1), 1));
                }

                ++oi;
                return oi;
            }

            if ((xcv1.data()._max_id == xcv2.data()._max_id) && (xcv2.data()._max_id != no_state)) {
                return oi;

                if (xcv1.data()._orientation == CGAL::SMALLER) {
                    *oi = CGAL::make_object(std::make_pair(max_vertex_obj(xcv1), 1));
                } else {
                    *oi = CGAL::make_object(std::make_pair(min_vertex_obj(xcv1), 1));
                }

                ++oi;
                return oi;
            }

            if (endp_coll == 1) {
                return oi;
            }

            // code ONLY FOR SEGMENTS !!!!!!!!

            std::list<CGAL::Object> base_objs;

            m_traits->intersect_2_object()(xcv1, xcv2, std::back_inserter(base_objs));

            if (base_objs.empty()) {
                return (oi);
            }

            // Attach labels to the intersection objects.
            std::list<CGAL::Object>::iterator obj_it;
            const std::pair<Base_point_2, unsigned int> *base_pt;
            const Base_x_monotone_curve_2 *base_xcv;

            for (obj_it = base_objs.begin(); obj_it != base_objs.end(); ++obj_it) {
                base_pt =
                    object_cast<std::pair<Base_point_2, unsigned int> > (&(*obj_it));

                if (base_pt != NULL) {
                    // Attach an invalid label to an itersection point.
                    *oi = CGAL::make_object
                          (std::make_pair(Point_2(base_pt->first), base_pt->second));
                    ++oi;
                } else {
                    base_xcv = object_cast<Base_x_monotone_curve_2> (&(*obj_it));
                    CGAL_assertion(base_xcv != NULL);

                    // Attach a merged label to the overlapping curve.
                    *oi = CGAL::make_object(X_monotone_curve_2(*base_xcv, Segment_Data_Label(no_state, no_state, xcv1.data()._orientation, -1)));
                    ++oi;
                }
            }

            return (oi);
        }
    };

    /*! Obtain an Intersect_2 functor object. */
    Intersect_2 intersect_2_object() const {
        return Intersect_2(this);
    }

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
        friend class Arr_SegmentData_traits<Traits_>;

    public:

        Comparison_result operator()(const Point_2 &p1, const Point_2 &p2) const {
            if ((p1.id == p2.id) && (p1.id != state(-1, -1))) {
                return EQUAL;
            }

            return m_base_cmp_xy(p1, p2);
        }
    };

    /*! Obtain a Construct_min_vertex_2 functor object. */
    Compare_xy_2 compare_xy_2_object() const {
        return (Compare_xy_2(((Base *)this)->compare_xy_2_object()));
    }

    /*! A functor that obtains the left endpoint of an x-monotone curve. */
    class Construct_min_vertex_2 {

    protected:

        //! The base operators.
        Base_construct_min_vertex_2 m_base_min_v;

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
        friend class Arr_SegmentData_traits<Traits_>;

    public:

        Point_2 operator()(const X_monotone_curve_2 &xcv) {
            Point_2 min_p = m_base_min_v(xcv);

            if (xcv.label()._orientation == CGAL::SMALLER) {
                min_p.id = xcv.data()._min_id;
            } else {
                min_p.id = xcv.data()._max_id;
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
        friend class Arr_SegmentData_traits<Traits_>;

    public:

        Point_2 operator()(const X_monotone_curve_2 &xcv) const {
            Point_2 max_p = m_base_max_v(xcv);

            if (xcv.label()._orientation == CGAL::SMALLER) {
                max_p.id = xcv.data()._max_id;
            } else {
                max_p.id = xcv.data()._min_id;
            }

            return (max_p);
        }
    };

    /*! Obtain a Construct_min_vertex_2 functor object. */
    Construct_max_vertex_2 construct_max_vertex_2_object() const {
        return (Construct_max_vertex_2(((Base *)this)->construct_max_vertex_2_object()));
    }

    Compare_y_at_x_2 compare_y_at_x_2_object() const {
        return (Compare_y_at_x_2(((Base *)this)->compare_y_at_x_2_object()));
    }

    Compare_x_2 compare_x_2_object() const {
        return (Compare_x_2(((Base *)this)->compare_x_2_object()));
    }

    Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const {
        return (Compare_endpoints_xy_2(((Base *)this)->compare_endpoints_xy_2_object()));
    }
};

} // namespace internal
} // namespace CGAL

#endif
