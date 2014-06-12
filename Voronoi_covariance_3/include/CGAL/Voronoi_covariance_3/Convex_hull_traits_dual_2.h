#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_2_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_2_H

namespace CGAL
{
  namespace Voronoi_covariance_3
  {
      // Traits classes for Convex_hull_2
      // XY
    template <class R>
        class Traits_xy_dual {
            public:
                typedef typename R::Plane_3 Point_2;

                struct Equal_2 {
                    typedef typename R::RT        RT;

                    bool operator() (Point_2 p, Point_2 q) {
                        RT diffa = p.a() * q.d() - q.a() * p.d();
                        RT diffb = p.b() * q.d() - q.b() * p.d();

                        return CGAL::is_zero(diffa) && CGAL::is_zero(diffb);
                    }
                };

                struct Left_turn_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_xy_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_yx_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        // TODO
                        return true;
                    }
                };

                struct Orientation_2 {
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return CGAL::COLLINEAR;
                    }
                };


                Equal_2 equal_2_object () const {
                    typedef bool result_type;

                    return Equal_2();
                }

                Left_turn_2 left_turn_2_object () const {
                    return Left_turn_2();
                }

                Less_xy_2 less_xy_2_object () const {
                    return Less_xy_2();
                }

                Less_yx_2 less_yx_2_object () const {
                    return Less_yx_2();
                }

                Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object () const {
                    return Less_signed_distance_to_line_2();
                }

                Orientation_2 orientation_2_object () const {
                    return Orientation_2();
                }
        };

      // XZ
    template <class R>
        class Traits_xz_dual {
            public:
                typedef typename R::Plane_3 Point_2;

                struct Equal_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) {
                        RT diffa = p.a() * q.d() - q.a() * p.d();
                        RT diffb = p.b() * q.d() - q.b() * p.d();

                        return CGAL::is_zero(diffa) && CGAL::is_zero(diffb);
                    }
                };

                struct Left_turn_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_xy_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_yx_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        // TODO
                        return true;
                    }
                };

                struct Orientation_2 {
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return CGAL::COLLINEAR;
                    }
                };


                Equal_2 equal_2_object () const {
                    return Equal_2();
                }

                Left_turn_2 left_turn_2_object () const {
                    return Left_turn_2();
                }

                Less_xy_2 less_xy_2_object () const {
                    return Less_xy_2();
                }

                Less_yx_2 less_yx_2_object () const {
                    return Less_yx_2();
                }

                Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object () const {
                    return Less_signed_distance_to_line_2();
                }

                Orientation_2 orientation_2_object () const {
                    return Orientation_2();
                }
        };

      // YZ
    template <class R>
        class Traits_yz_dual {
            public:
                typedef typename R::Plane_3 Point_2;

                struct Equal_2 {
                    typedef typename R::RT        RT;

                    bool operator() (Point_2 p, Point_2 q) {
                        RT diffa = p.a() * q.d() - q.a() * p.d();
                        RT diffb = p.b() * q.d() - q.b() * p.d();

                        return CGAL::is_zero(diffa) && CGAL::is_zero(diffb);
                    }
                };

                struct Left_turn_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_xy_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_yx_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q) const {
                        // TODO
                        return true;
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef bool result_type;

                    bool operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        // TODO
                        return true;
                    }
                };

                struct Orientation_2 {
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        // TODO
                        return CGAL::COLLINEAR;
                    }
                };

                Equal_2 equal_2_object () const {
                    return Equal_2();
                }

                Left_turn_2 left_turn_2_object () const {
                    return Left_turn_2();
                }

                Less_xy_2 less_xy_2_object () const {
                    return Less_xy_2();
                }

                Less_yx_2 less_yx_2_object () const {
                    return Less_yx_2();
                }

                Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object () const {
                    return Less_signed_distance_to_line_2();
                }

                Orientation_2 orientation_2_object () const {
                    return Orientation_2();
                }
        };
  } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_4_H

