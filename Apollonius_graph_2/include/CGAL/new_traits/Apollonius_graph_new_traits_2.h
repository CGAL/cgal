// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.2-branch/Apollonius_graph_2/include/CGAL/Apollonius_graph_traits_2.h $
// $Id: Apollonius_graph_traits_2.h 28567 2006-02-16 14:30:13Z lsaboret $
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_NEW_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_NEW_TRAITS_2_H


// FIXME: We include the old traits class file for now to get the functors.
#include <CGAL/Apollonius_graph_traits_2.h>

CGAL_BEGIN_NAMESPACE

// Would be *much* more convenient
// static inline
// Sign operator- (const Sign &s)
// {
//     return opposite (s);
// }
// 
// static inline
// Sign operator* (const Sign &s1, const Sign &s2)
// {
//     return static_cast<Sign> (static_cast<int> (s1) * static_cast<int> (s2));
// }

//-----------------------------------------------------------------------
//                        Conflict Base
//-----------------------------------------------------------------------
template < class K, class Method_tag >
class Conflict_2
{
public:
    typedef CGAL::Inverted_weighted_point<K>  Inverted_weighted_point;
    typedef typename K::FT                    FT;	
    typedef Sign                              result_type;
    struct Arity {};

protected:

    Sign orientation(const Inverted_weighted_point &p1, 
                     const Inverted_weighted_point &p2,
                     const Inverted_weighted_point &p3) const
    {
        return CGAL::sign(det3x3_by_formula<FT>(
                    p1.p(), p2.p(), p3.p(),
                    p1.x(), p2.x(), p3.x(),
                    p1.y(), p2.y(), p3.y()));
    }

    Sign radical_intersection(const Inverted_weighted_point &p1,
                              const Inverted_weighted_point &p2,
                              const Inverted_weighted_point &p3, int i) const
    {
        Sign s =  CGAL::sign(
                square(det3x3_by_formula<FT>(
                        p1.p(), p1.weight(), p1.y(),
                        p2.p(), p2.weight(), p2.y(),
                        p3.p(), p3.weight(), p3.y()))
                + square(det3x3_by_formula<FT>(
                        p1.p(), p1.x(), p1.weight(),
                        p2.p(), p2.x(), p2.weight(),
                        p3.p(), p3.x(), p3.weight()))
                - square(det3x3_by_formula<FT>(
                        p1.p(), p1.x(), p1.y(),
                        p2.p(), p2.x(), p2.y(),
                        p3.p(), p3.x(), p3.y())));

        if (s != ZERO || i < 0) return s;

        // perturbation
        switch (i) {
            case (1) : return radical_side(p3, p2, p1, 3);
            case (2) : return radical_side(p1, p3, p2, 3);
            case (3) : return radical_side(p1, p2, p3, 3);
            default :
                CGAL_assertion_msg(false, "Case should not hapen");
                return ZERO;
        }
    }	

    Sign radical_side(const Inverted_weighted_point &p1,
                      const Inverted_weighted_point &p2,
                      const Inverted_weighted_point &p3, int i) const
    {
        CGAL_assertion(i == -1 || i == 1 || i == 2 || i == 3);

        Sign s =  CGAL::sign(-1 * 
                ( det2x2_by_formula<FT>(
                      p1.p(), p1.x(), 
                      p2.p(), p2.x()) 
                  * det3x3_by_formula<FT>(
                      p1.p(), p1.weight(), p1.x(),
                      p2.p(), p2.weight(), p2.x(),
                      p3.p(), p3.weight(), p3.x())
                  + det2x2_by_formula<FT>(
                      p1.p(), p1.y(),
                      p2.p(), p2.y())
                  * det3x3_by_formula<FT>(
                      p1.p(), p1.weight(), p1.y(),
                      p2.p(), p2.weight(), p2.y(),
                      p3.p(), p3.weight(), p3.y())));

        if (s != ZERO || i < 1) return s; 

        // perturbation
        return POSITIVE;
    }

    Sign power_test(const Inverted_weighted_point &p1,
            const Inverted_weighted_point &p2,
            const Inverted_weighted_point &p3, int i) const
    {
        CGAL_assertion(i == -1 || i == 1 || i == 2 || i == 3);

        Sign s;
        FT xTest = det2x2_by_formula<FT>(
                p1.p(), p1.x(),
                p2.p(), p2.x());
        if (xTest != 0) {
            s = CGAL::sign(-1 * 
                    xTest 
                    * det3x3_by_formula<FT>(
                        p1.p(), p1.x(), p1.weight(),
                        p2.p(), p2.x(), p2.weight(),
                        p3.p(), p3.x(), p3.weight()));
        } else {
            s = CGAL::sign(-1 * 
                    det2x2_by_formula<FT>(
                        p1.p(), p1.y(),
                        p2.p(), p2.y())
                    * det3x3_by_formula<FT>(
                        p1.p(), p1.y(), p1.weight(),
                        p2.p(), p2.y(), p2.weight(),
                        p3.p(), p3.y(), p3.weight()));
        }
        if (s != ZERO || i < 1) return s;

        switch (i) {
            case (1) : return ordered_on_line (p1, p2, p3) ? NEGATIVE : POSITIVE;
            case (2) : return ordered_on_line (p2, p1, p3) ? NEGATIVE : POSITIVE;
            case (3) : return NEGATIVE;
            default : CGAL_assertion_msg (false, "this should not happen.");
        }

        // perturbation
        bool ool;
        if (i == 1) {
            ool = ordered_on_line(p2, p1, p3);
        } else if (i == 2) {
            ool = ordered_on_line(p1, p2, p3);
        } else if (i == 3) {
            ool = ordered_on_line(p1, p3, p2); 
        } else 	{		
            CGAL_assertion_msg(false, "this does not happen.");
            ool = false;
        }

        if (ool) return NEGATIVE;
            
        return POSITIVE;
    }

    Sign ordered_on_line_test(const Inverted_weighted_point &p1,
                              const Inverted_weighted_point &p2) const
    {
        FT detVal = det2x2_by_formula<FT>(
                p1.p(), p1.x(), 
                p2.p(), p2.x());
        if (detVal != 0) return CGAL::sign(detVal);
            
        return CGAL::sign(det2x2_by_formula<FT>(
                    p1.p(), p1.y(),
                    p2.p(), p2.y()));
    }

    bool ordered_on_line(const Inverted_weighted_point &p1,
                         const Inverted_weighted_point &p2,
                         const Inverted_weighted_point &p3) const
    {
        if (ordered_on_line_test(p1, p2) == POSITIVE)
           return ordered_on_line_test(p2, p3) == POSITIVE;
        return ordered_on_line_test(p3, p2) == POSITIVE;
    }	
};


//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag > 
class Vertex_conflict_new_2 
{
public:
    typedef typename K::Site_2                Site_2;
    typedef typename K::FT                    FT;	
    typedef Sign                              result_type;
    struct Arity {};

private:
    
    inline
    bool is_less (const Site_2 &p0, const Site_2 &p1) const
    {
        if (p0.weight() < p1.weight()) return true;
        if (p0.weight() > p1.weight()) return false;

        if (p0.x() < p1.x()) return true;
        if (p0.x() > p1.x()) return false;

        return p0.y() < p1.y();
    }

    inline
    int max_radius(const Site_2	&p0, const Site_2 &p1,
            const Site_2 &p2, const Site_2 &p3) const
    {
        int i = 0;
        const Site_2 *p = &p0;

        if (is_less (*p, p1)) { i = 1; p = &p1; }
        if (is_less (*p, p2)) { i = 2; p = &p2; }
        if (is_less (*p, p3)) { i = 3; }

        return i;
    }

    inline
    Sign predicate (const Site_2 &p1, const Site_2 &p2,
            const Site_2 &p3, const Site_2 &q, bool perturb) const
    {	
        FT xq = q.x() - p1.x(), yq = q.y() - p1.y(), wq = q.weight() - p1.weight();
        FT aq = xq * xq + yq * yq - wq * wq;
        
        // q is hiding p1
        if (sign (aq) != POSITIVE){
            // I BELIEVE MENELAOS RETURNS -1 in this case even when degernate 
            //if (sign (aq) == ZERO && ! perturb) return ZERO;
            return NEGATIVE;
        }

        FT x2 = p2.x() - p1.x(), y2 = p2.y() - p1.y(), w2 = p2.weight() - p1.weight();
        FT a2 = x2 * x2 + y2 * y2 - w2 * w2;

        CGAL_assertion (a2 > 0);

        FT x3 = p3.x() - p1.x(), y3 = p3.y() - p1.y(), w3 = p3.weight() - p1.weight();
        FT a3 = x3 * x3 + y3 * y3 - w3 * w3;

        CGAL_assertion (a3 > 0);

        FT ax3q = a3 * xq - x3 * aq; 
        FT ax2q = a2 * xq - x2 * aq;
        FT ax23 = a2 * x3 - x2 * a3;

        FT ay23 = a2 * y3 - y2 * a3;
        FT ay2q = a2 * yq - y2 * aq;
        FT ay3q = a3 * yq - y3 * aq;

        FT axw23q = ax23 * wq - ax2q * w3 + ax3q * w2;
        FT ayw23q = ay23 * wq - ay2q * w3 + ay3q * w2;

        FT axy23q = y2 * ax3q - y3 * ax2q + yq * ax23;

        // orientation
        Sign orient = CGAL::sign(axy23q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign orient1 = sign (ax23);

            Sign power_test = (orient1 == ZERO ?
                    Sign (sign (ay23) * sign (ayw23q)) :
                    Sign (orient1 * sign (axw23q)));

            if (power_test != ZERO || ! perturb) return Sign (- power_test);

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) return NEGATIVE;

            Sign o23, o2q, o3q;

            if (orient1 == ZERO) {
                o23 = sign (ay23);
                o2q = sign (ay2q);
                o3q = sign (ay3q);
            } else {
                o23 = sign (ax23);
                o2q = sign (ax2q);
                o3q = sign (ax3q);
            }

            if (o23 != o2q) return i == 2 ? NEGATIVE : POSITIVE;

            if (o23 == o3q) return i == 1 ? NEGATIVE : POSITIVE;

            return i == 0 ? NEGATIVE : POSITIVE;
        }
 
        // radical side 
        FT rs23q = ax23 * axw23q + ay23 * ayw23q;
        Sign radSide = CGAL::sign (rs23q);

        if (radSide == ZERO || radSide != orient) return orient;
       
        // radical intersection
        Sign radInt = enum_cast<Sign> (CGAL::compare (axw23q * axw23q + ayw23q * ayw23q, axy23q * axy23q));

        // radical intersection degenerate
        if (radInt == ZERO) {
            Sign radSideQ = sign (ax23 * axw23q + ay23 * ayw23q);
            
            CGAL_assertion (radSideQ != ZERO);

            if (! perturb) return (radSideQ == orient) ? ZERO : orient;

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) { 
                radInt = radSideQ;
            } else if (i == 2) {
                radInt = Sign (- sign (ax2q * axw23q + ay2q * ayw23q));
                if (radInt == ZERO) return NEGATIVE;
            } else if (i == 1) {
                radInt = sign (ax3q * axw23q + ay3q * ayw23q);
                if (radInt == ZERO) return NEGATIVE;
            } else {
                CGAL_assertion (i == 0);
                Sign radSide1 = Sign (- sign (ax2q * axw23q + ay2q * ayw23q));      
                if (radSide1 == ZERO) return NEGATIVE;

                Sign radSide2 = sign (ax3q * axw23q + ay3q * ayw23q);	
                if (radSide2 == ZERO) return NEGATIVE;

                radInt = Sign (- (radSideQ + radSide1 + radSide2));
            }
        }
        
        CGAL_assertion (! perturb || radInt != ZERO);

        if (radInt == NEGATIVE) return orient;
        
        return Sign (- radSide);
    }
    

    inline
    Sign predicate(const Site_2 &p1, const Site_2 &p2, 
            const Site_2 &q, bool perturb) const
    {
        // NOTE:***************************************
        // * the perturb boolean variable is not used 
        // * for consistancy with Menelaos
        // NOTE:***************************************
        FT x2 = p2.x() - p1.x(), y2 = p2.y() - p1.y(), w2 = p2.weight() - p1.weight();
        FT xq =  q.x() - p1.x(), yq =  q.y() - p1.y(), wq =  q.weight() - p1.weight();

        FT xw2q = x2 * wq - xq * w2;
        FT yw2q = y2 * wq - yq * w2;
        FT xy2q = x2 * yq - xq * y2;
        
        // orientation
        Sign orient = CGAL::sign(xy2q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign o12 = sign (x2);
            Sign o1q, o2q;

            Sign power_test;
            if (o12 != ZERO) {
                power_test = Sign (o12 * sign (xw2q));
                 
                // this results is consistant with Menelaos
                if (power_test != ZERO) return Sign(- power_test);

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return Sign(- power_test);

                o1q = sign (xq);
                o2q = sign (q.x() - p2.x());
            } else {
                o12 = sign (y2);
                power_test = Sign (o12 * sign (yw2q));

                // this results is consistant with Menelaos
                if (power_test != ZERO) return Sign(- power_test);

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return Sign(- power_test);

                o1q = sign (yq);
                o2q = sign (q.y() - p2.y());
            }

            if (o1q != o12) return POSITIVE;
            if (o2q == o12) return POSITIVE;

            return NEGATIVE;
        }

        // radical side 
        FT rs12q = x2 * xw2q + y2 * yw2q;
        Sign radSide = CGAL::sign (rs12q);

        if (radSide == ZERO || radSide == orient) {
            return Sign(- orient);
        }

        // radical intersection
        Sign radInt = enum_cast<Sign> (CGAL::compare (xw2q * xw2q + yw2q * yw2q, xy2q * xy2q));

        // radical intersection degerate
        if (radInt == ZERO) {
            CGAL_assertion (radSide != ZERO);
            
            // this result is consistant with the perturb on off idea
            //if (! perturb) return (radSide == orient) ? ZERO : orient;

            FT rs2q1 = (p2.x() - q.x()) * xw2q + (p2.y() - q.y()) * yw2q;
            Sign radSide1 = sign (rs2q1);
            if (radSide1 == ZERO) return NEGATIVE;
            
            FT rsq12 = xq * xw2q + yq * yw2q;
            Sign radSide2 = sign (rsq12);
            if (radSide2 == ZERO) return NEGATIVE;
 
            return Sign (- (radSide1 * radSide2));
        }

        CGAL_assertion (! perturb || radInt != ZERO);

        if (radInt == POSITIVE) return orient;
        return radSide;
    }

public:
    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &p3, const Site_2 &q, bool perturb = true) const
    {	
        Sign newPred = predicate (p1, p2, p3, q, perturb);	
        CGAL_assertion (! perturb || newPred != ZERO);
        return newPred;
    }

    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &q, bool perturb = true) const
    {
        Sign newPred = predicate (p1, p2, q, perturb);	
        CGAL_assertion (! perturb || newPred != ZERO);
        return newPred;
    }
};

//-----------------------------------------------------------------------
//                     Edge Conflict Base
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Edge_conflict_2 : public Conflict_2<K, Method_tag>
{
public:
    typedef CGAL::Inverted_weighted_point<K>  Inverted_weighted_point;
    typedef bool                              result_type;
    struct Arity {};

protected:

    bool edge_conflict_test(const Inverted_weighted_point &p2,
                            const Inverted_weighted_point &p3,
                            const Inverted_weighted_point &p4,
                            const Inverted_weighted_point &q,
                            bool b, int i23Q, int i24Q) const
    {

        // orientations
        Sign orient23Q = this->orientation(p2, p3, q);
        Sign orient42Q = this->orientation(p4, p2, q);
        Sign orient234 = this->orientation(p2, p3, p4);

        // radical intersections
        Sign radInt23Q = radical_intersection(p2, p3, q, -1);
        Sign radInt24Q = radical_intersection(p2, p4, q, -1);

        // radical side
        Sign radSid2Q3 = radical_side(p2, q, p3, -1);
        Sign radSid2Q4 = radical_side(p2, q, p4, -1);

        // order of a line
        bool oolQ24 = ordered_on_line(q, p2, p4);
        bool oolQ23 = ordered_on_line(q, p2, p3); 

        if(b)
        {
            if (sign (q.p()) != POSITIVE ) return true;
            // degenerate case
            if(orient234 == ZERO && orient23Q == ZERO && orient42Q == ZERO)
                return (oolQ23 || oolQ24);	
            // non degenerate case
            else if (! ((radInt23Q != NEGATIVE && radSid2Q3 == NEGATIVE) && 
                        (radInt24Q != NEGATIVE && radSid2Q4 == NEGATIVE)))
                return true;
            else if (orient234 != NEGATIVE)
                return orient23Q != POSITIVE && orient42Q != POSITIVE;
            else
                return orient23Q != POSITIVE || orient42Q != POSITIVE;
        }
        else
        {
            CGAL_assertion (sign (q.p()) == POSITIVE);
            // degenerate case
            if(orient234 == ZERO && orient23Q == ZERO && orient42Q == ZERO)
                return (oolQ23 && oolQ24);		
            // non degenerate case	
            else if (! ((radInt23Q != NEGATIVE && radSid2Q3 == NEGATIVE) && 
                        (radInt24Q != NEGATIVE && radSid2Q4 == NEGATIVE)))
                return false;
            else if (orient234 != NEGATIVE) 
                return orient23Q != POSITIVE || orient42Q != POSITIVE;
            else 
                return orient23Q != POSITIVE && orient42Q != POSITIVE;
        }
    }
};


//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Finite_edge_interior_conflict_new_2 : public Edge_conflict_2<K, Method_tag>
{
public:
    typedef CGAL::Inverted_weighted_point<K> Inverted_weighted_point;
    typedef CGAL::Weighted_point_inverter<K> Weighted_point_inverter;
    typedef typename K::Site_2               Site_2;
    typedef typename K::Point_2              Point_2;
    typedef bool                             result_type;
    struct Arity {};

    inline
    bool operator()(const Site_2& p1, const Site_2& p2, 
            const Site_2& q, bool b) const
    {
        Weighted_point_inverter inverter(p1);
        return edge_conflict_test(inverter(p2), 
                Inverted_weighted_point (Site_2(Point_2(0,0),0),1),
                Inverted_weighted_point (Site_2(Point_2(0,0),0),1),
                inverter(q), b, 2, 2);
    }

    inline
    bool operator()(const Site_2& p1, const Site_2& p2, 
            const Site_2& p3, const Site_2& q, bool b) const
    {
        Weighted_point_inverter inverter(p2);
        return edge_conflict_test(inverter(p1), 
                Inverted_weighted_point(Site_2(Point_2(0,0),0),1),
                inverter(p3), inverter(q), b, 2, 1);
    }

    inline
    bool operator()(const Site_2& p1, const Site_2& p2,
            const Site_2& p3, const Site_2& p4, const Site_2& q, bool b) const
    {
        Weighted_point_inverter inverter(p2);
        return edge_conflict_test(inverter(p1), inverter(p4), inverter(p3), 
                inverter(q), b, 1, 1);
    }
};

//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Infinite_edge_interior_conflict_new_2 : public Edge_conflict_2<K, Method_tag>
{
public:
    typedef CGAL::Weighted_point_inverter<K> Weighted_point_inverter;
    typedef CGAL::Inverted_weighted_point<K> Inverted_weighted_point;
    typedef typename K::Site_2               Site_2;
    typedef typename K::Point_2              Point_2;
    typedef bool                             result_type;
    typedef Arity_tag<5>                     Arity;

    inline
    bool operator()(const Site_2& p2, const Site_2& p3, 
            const Site_2& p4, const Site_2& q, bool b) const
    {
        Weighted_point_inverter inverter(p2);
        return edge_conflict_test(
                Inverted_weighted_point(Site_2(Point_2(0,0),0), 1),
                inverter(p4), inverter(p3), inverter(q), b, 1, 1);
    }
};


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
template < class Rep, class MTag = Ring_tag >
class Apollonius_graph_new_traits_2
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
private:  
  typedef Apollonius_graph_new_traits_2<Rep,MTag>       Self;
  typedef Apollonius_graph_kernel_wrapper_2<Rep>        Kernel;

public:
  typedef Rep                                           R;
  typedef MTag                                          Method_tag;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Site_2                       Site_2;

  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Rep::Segment_2                       Segment_2;

  typedef typename Kernel::Object_2                     Object_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::RT                           RT;


public:
  // OBJECT CONSTRUCTION & ASSIGNMENT
  //---------------------------------
  typedef typename Kernel::Construct_object_2     Construct_object_2;
  typedef typename Kernel::Assign_2               Assign_2;

  // CONSTRUCTIONS
  //--------------
  // vertex and dual site
  typedef CGAL::Construct_Apollonius_vertex_2<Kernel>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_site_2<Kernel>
  /*                                        */ Construct_Apollonius_site_2;


  // PREDICATES
  //-----------
  typedef CGAL::Ag2_compare_x_2<Kernel>                 Compare_x_2;
  typedef CGAL::Ag2_compare_y_2<Kernel>                 Compare_y_2;
  typedef CGAL::Ag2_compare_weight_2<Kernel>            Compare_weight_2;
  typedef CGAL::AG2_Orientation_test_2<Kernel,MTag>     Orientation_2;
  typedef CGAL::Is_hidden_2<Kernel,MTag>                Is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<Kernel,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_new_2<Kernel,MTag>          Vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_new_2<Kernel,MTag>
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_new_2<Kernel,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<Kernel,MTag>       Is_degenerate_edge_2;


public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  // OBJECT CONSTRUCTION & ASSIGNMENT
  Assign_2
  assign_2_object() const {
    return Assign_2();
  }

  Construct_object_2
  construct_object_2_object() const { 
    return Construct_object_2();
  }


  // CONSTRUCTIONS
  //--------------
  Construct_Apollonius_vertex_2
  construct_Apollonius_vertex_2_object() const { 
    return Construct_Apollonius_vertex_2();
  }

  Construct_Apollonius_site_2
  construct_Apollonius_site_2_object() const {
    return Construct_Apollonius_site_2();
  }


  // PREDICATES
  //-----------
  Compare_x_2
  compare_x_2_object() const {
    return Compare_x_2();
  }

  Compare_y_2
  compare_y_2_object() const {
    return Compare_y_2();
  }

  Compare_weight_2
  compare_weight_2_object() const {
    return Compare_weight_2();
  }

  Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
  }

  Oriented_side_of_bisector_2
  oriented_side_of_bisector_2_object() const {
    return Oriented_side_of_bisector_2();
  }

  Vertex_conflict_2
  vertex_conflict_2_object() const {
    return Vertex_conflict_2();
  }

  Finite_edge_interior_conflict_2
  finite_edge_interior_conflict_2_object() const {
    return Finite_edge_interior_conflict_2();
  }

  Infinite_edge_interior_conflict_2
  infinite_edge_interior_conflict_2_object() const {
    return Infinite_edge_interior_conflict_2();
  }

  Is_degenerate_edge_2
  is_degenerate_edge_2_object() const {
    return Is_degenerate_edge_2();
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_NEW_TRAITS_2_H
