// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Kernel_checker.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================

#ifndef CGAL_KERNEL_CHECKER_H
#define CGAL_KERNEL_CHECKER_H

#include <CGAL/basic.h>
#include <CGAL/multirep.h>
#include <utility>
#include <iostream>

CGAL_BEGIN_NAMESPACE


// Class used by Kernel_checker.
template <class O1, class O2, 
          class Conv>
class Special_construction
{
    O1 o1;
    O2 o2;

public:

    Special_construction(const O1 &oo1 = O1(), const O2 &oo2 = O2())
	: o1(oo1), o2(oo2) {}

    // we have another result type now ...
    typedef MultiRep<typename O1::result_type, typename O2::result_type, Conv> result_type;
    
    // the arity can stay ...
    typedef typename O1::Arity       Arity;
    
    // --------------------------------------------------------------------------------        
    // attention: we have to handle MultiRep arguments in the templated
    // operators, but also other stuff (for instance ints in Construct_vertex)
    // problem: does this work on VC++ ?
    
    template<class FR, class SR, class Converter>
    const FR& get_first_rep(const MultiRep<FR,SR,Converter>& obj) const
    { return obj.get_first_rep(); }
    
    template<class T>
    const T& get_first_rep(const T& other_obj) const
    { return other_obj; }
    
    template<class FR, class SR, class Converter>
    const SR& get_second_rep(const MultiRep<FR,SR,Converter>& obj) const
    { return obj.get_second_rep(); }
    
    template<class T>
    const T& get_second_rep(const T& other_obj) const
    { return other_obj; }
    
    // --------------------------------------------------------------------------------       
    

    template <class A1>
    result_type
    operator()(const A1 &a1) const
    { 
	
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1));
	typename O2::result_type res2 = o2(get_second_rep(a1));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);
    }

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {	
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), get_first_rep(a3));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), get_second_rep(a3));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), 
	                                   get_first_rep(a3), get_first_rep(a4));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), 
	                                   get_second_rep(a3), get_second_rep(a4));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), 
	                                   get_first_rep(a3), get_first_rep(a4), get_first_rep(a5));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), 
	                                   get_second_rep(a3), get_second_rep(a4), get_second_rep(a5));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2); 
    }
    
    template <class A1, class A2, class A3, class A4, class A5, class A6>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), 
	                                   get_first_rep(a3), get_first_rep(a4),
					   get_first_rep(a5), get_first_rep(a6));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), 
	                                   get_second_rep(a3), get_second_rep(a4),
					   get_second_rep(a5), get_second_rep(a6));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);  
    }
    
    template <class A1, class A2, class A3, class A4, class A5,class A6,class A7>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), 
	                                   get_first_rep(a3), get_first_rep(a4),
					   get_first_rep(a5), get_first_rep(a6),
					   get_first_rep(a7));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), 
	                                   get_second_rep(a3), get_second_rep(a4),
					   get_second_rep(a5), get_second_rep(a6),
					   get_second_rep(a7));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);   
    }        

    template <class A1, class A2, class A3, class A4, 
              class A5, class A6, class A7, class A8>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
    {
	// construct the two representations ...
	typename O1::result_type res1 = o1(get_first_rep(a1), get_first_rep(a2), 
	                                   get_first_rep(a3), get_first_rep(a4),
					   get_first_rep(a5), get_first_rep(a6),
					   get_first_rep(a7), get_first_rep(a8));
	typename O2::result_type res2 = o2(get_second_rep(a1), get_second_rep(a2), 
	                                   get_second_rep(a3), get_second_rep(a4),
					   get_second_rep(a5), get_second_rep(a6),
					   get_second_rep(a7), get_second_rep(a8));
	
	// should we add some kind of check here as well ?
	
	// we use the special constructor here ...
        return result_type(res1, res2);    
    }        

    // Same thing with more arguments...
};



// Class used by Kernel_checker.
template <class O1, class O2>
class Predicate_checker
{
    O1 o1;
    O2 o2;

public:

    Predicate_checker(const O1 &oo1 = O1(), const O2 &oo2 = O2())
	: o1(oo1), o2(oo2) {}

    // the result type and arity of predicates can stay ...
    
    typedef typename O1::result_type result_type;
    typedef typename O1::Arity       Arity;

    template <class A1>
    result_type
    operator()(const A1 &a1) const
    {     
#if defined(CGAL_KCH_DEBUG)    
        std::cout << "Predicate_checker(1):";
	a1.output(std::cout);
	std::cout << "\n";
#endif         
	typename O1::result_type res1 = o1(a1);
	typename O2::result_type res2 = o2(a1.get_second_rep());
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {
#if defined(CGAL_KCH_DEBUG) 
        std::cout << "Predicate_checker(2):";   
        a1.output(std::cout); 
	a2.output(std::cout);
	std::cout << "\n";
#endif     
	typename O1::result_type res1 = o1(a1, a2);
	typename O2::result_type res2 = o2(a1.get_second_rep(), a2.get_second_rep());
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
#if defined(CGAL_KCH_DEBUG)  
        std::cout << "Predicate_checker(3):";    
        a1.output(std::cout); a2.output(std::cout); a3.output(std::cout);
	std::cout << "\n";
#endif     
	typename O1::result_type res1 = o1(a1, a2, a3);
	typename O2::result_type res2 = o2(a1.get_second_rep(), a2.get_second_rep(),
	                                   a3.get_second_rep());
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
#if defined(CGAL_KCH_DEBUG)   
        std::cout << "Predicate_checker(4):";  
        a1.output(std::cout); a2.output(std::cout); 
	a3.output(std::cout); a4.output(std::cout);
#endif    
	typename O1::result_type res1 = o1(a1, a2, a3, a4);
	typename O2::result_type res2 = o2(a1.get_second_rep(), a2.get_second_rep(),
	                                   a3.get_second_rep(), a4.get_second_rep());
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4, a5);
	typename O2::result_type res2 = o2(a1.get_second_rep(), a2.get_second_rep(),
	                                   a3.get_second_rep(), a4.get_second_rep(), 
					   a5.get_second_rep());	
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << ", " << a5 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }

    template <class A1, class A2, class A3, class A4, class A5, class A6>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6) const
    {
	typename O1::result_type res1 = o1(a1, a2, a3, a4, a5, a6);	
	typename O2::result_type res2 = o2(a1.get_second_rep(), a2.get_second_rep(),
	                                   a3.get_second_rep(), a4.get_second_rep(), 
					   a5.get_second_rep(), a6.get_second_rep());	
	
	if (res1 != res2)
	{
	    std::cerr << "Kernel_checker error : " << res1 << " != " << res2
		      << " for the inputs : " << std::endl;
	    std::cerr << a1 << ", " << a2 << ", " << a3 << ", " << a4
		      << ", " << a5 << std::endl;
#ifdef __GNUG__
	    std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
	    CGAL_kernel_assertion(false);
	}
	return res1;
    }
    // Same thing with more arguments...
};



// offer the ability to use different converters for different objects ... 
// (might be important for some traits classes with vector type = direction type ?) 

template<class Converter>
struct Default_converter_traits {
  typedef Converter   Converter_object_2;
  typedef Converter   Converter_object_3;          

  typedef Converter   Converter_point_2;
  typedef Converter   Converter_vector_2;
  typedef Converter   Converter_direction_2;
  typedef Converter   Converter_line_2;
  typedef Converter   Converter_ray_2;
  typedef Converter   Converter_segment_2;
  typedef Converter   Converter_triangle_2;
  typedef Converter   Converter_circle_2;
  typedef Converter   Converter_iso_rectangle_2;

  typedef Converter   Converter_point_3;
  typedef Converter   Converter_vector_3;
  typedef Converter   Converter_direction_3;
  typedef Converter   Converter_line_3;
  typedef Converter   Converter_plane_3;
  typedef Converter   Converter_ray_3;
  typedef Converter   Converter_segment_3;
  typedef Converter   Converter_triangle_3;
  typedef Converter   Converter_tetrahedron_3;
  typedef Converter   Converter_sphere_3;
  typedef Converter   Converter_iso_cuboid_3;  
};

// For now, we inherit all geometric objects and constructions from K1, and
// just overload the predicates.


template <class K1, class K2, 
          class Conv>
class Kernel_checker
  : public K1
{
    typedef K1     Kernel1;
    typedef K2     Kernel2;

    Kernel2 k2;

    typedef Conv   c;
    
    // the FT/RT stuff is still a problem ...
    // we derive it at the moment from K1

    // get the conversion functors from Conv ...
    
public:

    typedef typename Conv::Converter_object_2    Converter_object_2;
    typedef typename Conv::Converter_object_3    Converter_object_3;          

    typedef typename Conv::Converter_point_2     Converter_point_2;
    typedef typename Conv::Converter_vector_2    Converter_vector_2;
    typedef typename Conv::Converter_direction_2 Converter_direction_2;
    typedef typename Conv::Converter_line_2      Converter_line_2;
    typedef typename Conv::Converter_ray_2       Converter_ray_2;
    typedef typename Conv::Converter_segment_2   Converter_segment_2;
    typedef typename Conv::Converter_triangle_2  Converter_triangle_2;
    typedef typename Conv::Converter_circle_2    Converter_circle_2;
    typedef typename Conv::Converter_iso_rectangle_2   Converter_iso_rectangle_2;

    typedef typename Conv::Converter_point_3     Converter_point_3;
    typedef typename Conv::Converter_vector_3    Converter_vector_3;
    typedef typename Conv::Converter_direction_3 Converter_direction_3;
    typedef typename Conv::Converter_line_3      Converter_line_3;
    typedef typename Conv::Converter_plane_3     Converter_plane_3;
    typedef typename Conv::Converter_ray_3       Converter_ray_3;
    typedef typename Conv::Converter_segment_3   Converter_segment_3;
    typedef typename Conv::Converter_triangle_3  Converter_triangle_3;
    typedef typename Conv::Converter_tetrahedron_3   Converter_tetrahedron_3;
    typedef typename Conv::Converter_sphere_3    Converter_sphere_3;
    typedef typename Conv::Converter_iso_cuboid_3   Converter_iso_cuboid_3;

    
    // now the typedefs; we use the MultiRep class template ...
    typedef CGAL::MultiRep<typename K1::Object_2,typename K2::Object_2, Converter_object_2>          Object_2;
    typedef CGAL::MultiRep<typename K1::Object_3,typename K2::Object_3, Converter_object_3>          Object_3;          

    typedef CGAL::MultiRep<typename K1::Point_2,typename K2::Point_2, Converter_point_2>             Point_2;
    typedef CGAL::MultiRep<typename K1::Vector_2,typename K2::Vector_2, Converter_vector_2>          Vector_2;
    typedef CGAL::MultiRep<typename K1::Direction_2,typename K2::Direction_2, Converter_direction_2> Direction_2;
    typedef CGAL::MultiRep<typename K1::Line_2,typename K2::Line_2, Converter_line_2>                Line_2;
    typedef CGAL::MultiRep<typename K1::Ray_2,typename K2::Ray_2, Converter_ray_2>                   Ray_2;
    typedef CGAL::MultiRep<typename K1::Segment_2,typename K2::Segment_2, Converter_segment_2>       Segment_2;
    typedef CGAL::MultiRep<typename K1::Triangle_2,typename K2::Triangle_2, Converter_triangle_2>    Triangle_2;
    typedef CGAL::MultiRep<typename K1::Circle_2,typename K2::Circle_2, Converter_circle_2>          Circle_2;
    typedef CGAL::MultiRep<typename K1::Iso_rectangle_2,typename K2::Iso_rectangle_2, Converter_iso_rectangle_2>   Iso_rectangle_2;

    typedef CGAL::MultiRep<typename K1::Point_3,typename K2::Point_3, Converter_point_3>             Point_3;
    typedef CGAL::MultiRep<typename K1::Vector_3,typename K2::Vector_3, Converter_vector_3>          Vector_3;
    typedef CGAL::MultiRep<typename K1::Direction_3,typename K2::Direction_3, Converter_direction_3> Direction_3;
    typedef CGAL::MultiRep<typename K1::Line_3,typename K2::Line_3, Converter_line_3>                Line_3;
    typedef CGAL::MultiRep<typename K1::Plane_3,typename K2::Plane_3, Converter_plane_3>             Plane_3;
    typedef CGAL::MultiRep<typename K1::Ray_3,typename K2::Ray_3, Converter_ray_3>                   Ray_3;
    typedef CGAL::MultiRep<typename K1::Segment_3,typename K2::Segment_3, Converter_segment_3>       Segment_3;
    typedef CGAL::MultiRep<typename K1::Triangle_3,typename K2::Triangle_3, Converter_triangle_3>    Triangle_3;
    typedef CGAL::MultiRep<typename K1::Tetrahedron_3,typename K2::Tetrahedron_3, Converter_tetrahedron_3>       Tetrahedron_3;
    typedef CGAL::MultiRep<typename K1::Sphere_3,typename K2::Sphere_3, Converter_sphere_3>          Sphere_3;
    typedef CGAL::MultiRep<typename K1::Iso_cuboid_3,typename K2::Iso_cuboid_3, Converter_iso_cuboid_3> Iso_cuboid_3;
    

    // typedef std::pair<K1::Point_2, K2::Point_2>  Point_2;
    // ...  Same thing for all objects.

#define CGAL_check_pred(X, Y) \
    typedef Predicate_checker<typename K1::X, typename K2::X> X; \
    X Y() const { return X(K1::Y(), k2.Y()); }
    
#define CGAL_special_cons(X, Y) \
    typedef Special_construction<typename K1::X, typename K2::X, Conv> X; \
    X Y() const { return X(K1::Y(), k2.Y()); }    

#define CGAL_Kernel_pred(Y,Z) CGAL_check_pred(Y, Z)
#define CGAL_Kernel_cons(Y,Z) CGAL_special_cons(Y,Z)

public:

#include <CGAL/Kernel/interface_macros.h>
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_CHECKER_H
