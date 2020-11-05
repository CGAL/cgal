// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jocelyn Meyron
//


#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_2_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_2_H

#include <CGAL/license/Convex_hull_3.h>


// Traits classes used during the computation of the dual
// convex hull when the dual points are coplanar
// The traits are used by convex_hull_2

// Each traits class is based on the concept ConvexHullTraits_2

namespace CGAL
{
  namespace Convex_hull_3
  {
      // XY
    template <class R>
        class Traits_xy_dual {
            public:
                typedef typename R::Plane_3 Point_2;

                struct Equal_2 {
                    typedef bool result_type;
                    typedef typename R::RT        RT;

                    result_type operator() (Point_2 p, Point_2 q) {
                        RT diffa = p.a() * q.d() - q.a() * p.d();
                        RT diffb = p.b() * q.d() - q.b() * p.d();

                        return CGAL::is_zero(diffa) && CGAL::is_zero(diffb);
                    }
                };

                struct Left_turn_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffapq = p.a() * q.d() - q.a() * p.d();
                        RT diffbpr = p.b() * r.d() - r.b() * p.d();
                        RT diffbpq = p.b() * q.d() - q.b() * p.d();
                        RT diffapr = p.a() * r.d() - r.a() * p.d();

                        RT prod = diffapq * diffbpr - diffbpq * diffapr;

                        if (CGAL::is_positive(q.d() * r.d())) {
                            return CGAL::is_negative(prod);
                        } else {
                            return CGAL::is_positive(prod);
                        }
                    }
                };

                struct Less_xy_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffa = q.d() * p.a() - p.d() * q.a();
                        RT diffb = q.d() * p.b() - p.d() * q.b();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffa)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffa) && CGAL::is_positive(diffb));
                            }
                        }

                        if (CGAL::is_negative(diffa)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffa) && CGAL::is_negative(diffb));
                    }
                };

                struct Less_yx_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffa = q.d() * p.a() - p.d() * q.a();
                        RT diffb = q.d() * p.b() - p.d() * q.b();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffb)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffb) && CGAL::is_positive(diffa));
                            }
                        }

                        if (CGAL::is_negative(diffb)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffb) && CGAL::is_negative(diffa));
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        RT A = p.a() * q.d() - q.a() * p.d();
                        RT B = p.b() * q.d() - q.b() * p.d();

                        RT dist = s.d() * (B * r.a() - A * r.b()) - r.d() * (B * s.a() - A * s.b());

                        if (CGAL::is_positive(r.d() * s.d())) {
                            return CGAL::is_positive(dist);
                        }

                        return CGAL::is_negative(dist);
                    }
                };

                struct Orientation_2 {
                    typedef typename R::RT        RT;
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffapq = p.a() * q.d() - q.a() * p.d();
                        RT diffbpr = p.b() * r.d() - r.b() * p.d();
                        RT diffbpq = p.b() * q.d() - q.b() * p.d();
                        RT diffapr = p.a() * r.d() - r.a() * p.d();

                        RT prod = diffapq * diffbpr - diffbpq * diffapr;

                        if (CGAL::is_zero(prod)) {
                            return CGAL::COLLINEAR;
                        }

                        if (CGAL::is_positive(q.d() * r.d())) {
                            if (CGAL::is_positive(prod)) {
                                return CGAL::RIGHT_TURN;
                            } else {
                                return CGAL::LEFT_TURN;
                            }
                        }

                        if (CGAL::is_positive(prod)) {
                                return CGAL::LEFT_TURN;
                        }

                        return CGAL::RIGHT_TURN;
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
                        RT diffc = p.c() * q.d() - q.c() * p.d();

                        return CGAL::is_zero(diffa) && CGAL::is_zero(diffc);
                    }
                };

                struct Left_turn_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffapq = p.a() * q.d() - q.a() * p.d();
                        RT diffcpr = p.c() * r.d() - r.c() * p.d();
                        RT diffcpq = p.c() * q.d() - q.c() * p.d();
                        RT diffapr = p.a() * r.d() - r.a() * p.d();

                        RT prod = diffapq * diffcpr - diffcpq * diffapr;

                        if (CGAL::is_positive(q.d() * r.d())) {
                            return CGAL::is_negative(prod);
                        } else {
                            return CGAL::is_positive(prod);
                        }
                    }
                };

                struct Less_xy_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffa = q.d() * p.a() - p.d() * q.a();
                        RT diffc = q.d() * p.c() - p.d() * q.c();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffa)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffa) && CGAL::is_positive(diffc));
                            }
                        }

                        if (CGAL::is_negative(diffa)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffa) && CGAL::is_negative(diffc));
                    }
                };

                struct Less_yx_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffa = q.d() * p.a() - p.d() * q.a();
                        RT diffc = q.d() * p.c() - p.d() * q.c();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffc)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffc) && CGAL::is_positive(diffa));
                            }
                        }

                        if (CGAL::is_negative(diffc)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffc) && CGAL::is_negative(diffa));
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        RT A = p.a() * q.d() - q.a() * p.d();
                        RT B = p.c() * q.d() - q.c() * p.d();

                        RT dist = s.d() * (B * r.a() - A * r.c()) - r.d() * (B * s.a() - A * s.c());

                        if (CGAL::is_positive(r.d() * s.d())) {
                            return CGAL::is_positive(dist);
                        }

                        return CGAL::is_negative(dist);
                    }
                };

                struct Orientation_2 {
                    typedef typename R::RT        RT;
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffapq = p.a() * q.d() - q.a() * p.d();
                        RT diffcpr = p.c() * r.d() - r.c() * p.d();
                        RT diffcpq = p.c() * q.d() - q.c() * p.d();
                        RT diffapr = p.a() * r.d() - r.a() * p.d();

                        RT prod = diffapq * diffcpr - diffcpq * diffapr;

                        if (CGAL::is_zero(prod)) {
                            return CGAL::COLLINEAR;
                        }

                        if (CGAL::is_positive(q.d() * r.d())) {
                            if (CGAL::is_positive(prod)) {
                                return CGAL::RIGHT_TURN;
                            } else {
                                return CGAL::LEFT_TURN;
                            }
                        }

                        if (CGAL::is_positive(prod)) {
                                return CGAL::LEFT_TURN;
                        }

                        return CGAL::RIGHT_TURN;
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
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) {
                        RT diffb = p.b() * q.d() - q.b() * p.d();
                        RT diffc = p.c() * q.d() - q.c() * p.d();

                        return CGAL::is_zero(diffb) && CGAL::is_zero(diffc);
                    }
                };

                struct Left_turn_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffbpq = p.b() * q.d() - q.b() * p.d();
                        RT diffcpr = p.c() * r.d() - r.c() * p.d();
                        RT diffcpq = p.c() * q.d() - q.c() * p.d();
                        RT diffbpr = p.b() * r.d() - r.b() * p.d();

                        RT prod = diffbpq * diffcpr - diffcpq * diffbpr;

                        if (CGAL::is_positive(q.d() * r.d())) {
                            return CGAL::is_negative(prod);
                        } else {
                            return CGAL::is_positive(prod);
                        }
                    }
                };

                struct Less_xy_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffb = q.d() * p.b() - p.d() * q.b();
                        RT diffc = q.d() * p.c() - p.d() * q.c();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffb)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffb) && CGAL::is_positive(diffc));
                            }
                        }

                        if (CGAL::is_negative(diffb)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffb) && CGAL::is_negative(diffc));
                    }
                };

                struct Less_yx_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q) const {
                        RT diffb = q.d() * p.b() - p.d() * q.b();
                        RT diffc = q.d() * p.c() - p.d() * q.c();

                        if (CGAL::is_positive(p.d() * q.d())) {
                            if (CGAL::is_positive(diffc)) {
                                return true;
                            } else {
                                return (CGAL::is_zero(diffc) && CGAL::is_positive(diffb));
                            }
                        }

                        if (CGAL::is_negative(diffc)) {
                            return true;
                        }

                        return (CGAL::is_zero(diffc) && CGAL::is_negative(diffb));
                    }
                };

                struct Less_signed_distance_to_line_2 {
                    typedef typename R::RT        RT;
                    typedef bool result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r,Point_2 s) const {
                        RT A = p.b() * q.d() - q.b() * p.d();
                        RT B = p.c() * q.d() - q.c() * p.d();

                        RT dist = s.d() * (B * r.b() - A * r.c()) - r.d() * (B * s.b() - A * s.c());

                        if (CGAL::is_positive(r.d() * s.d())) {
                            return CGAL::is_positive(dist);
                        }

                        return CGAL::is_negative(dist);
                    }
                };

                struct Orientation_2 {
                    typedef typename R::RT        RT;
                    typedef CGAL::Orientation result_type;

                    result_type operator() (Point_2 p, Point_2 q, Point_2 r) const {
                        RT diffbpq = p.b() * q.d() - q.b() * p.d();
                        RT diffcpr = p.c() * r.d() - r.c() * p.d();
                        RT diffcpq = p.c() * q.d() - q.c() * p.d();
                        RT diffbpr = p.b() * r.d() - r.b() * p.d();

                        RT prod = diffbpq * diffcpr - diffcpq * diffbpr;

                        if (CGAL::is_zero(prod)) {
                            return CGAL::COLLINEAR;
                        }

                        if (CGAL::is_positive(q.d() * r.d())) {
                            if (CGAL::is_positive(prod)) {
                                return CGAL::RIGHT_TURN;
                            } else {
                                return CGAL::LEFT_TURN;
                            }
                        }

                        if (CGAL::is_positive(prod)) {
                                return CGAL::LEFT_TURN;
                        }

                        return CGAL::RIGHT_TURN;
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
  } // namespace Convex_hull_3
} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_2_H

