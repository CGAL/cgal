#ifndef VISIBILITY_COMPLEX_2_H
#define VISIBILITY_COMPLEX_2_H

#include <CEP/Visibility_complex/Visibility_complex_antichain.h>
#include <CEP/Visibility_complex/Visibility_complex_iterators.h>
#include <CEP/Visibility_complex/Visibility_complex_sweep_iterator.h>

#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class _Tr , 
	   class It   = Visibility_complex_items ,
	   class Flip = Visibility_complex_flip_traits >
class Visibility_complex_2 
{
public:
    // -------------------------------------------------------------------------
    typedef _Tr                                          Gt;
    typedef Visibility_complex_2<_Tr,It,Flip>                 Self;
    typedef Visibility_complex_antichain<_Tr,It,Flip>         Antichain;
    typedef typename Antichain::Disk_handle                    Disk_handle;
    // -------------------------------------------------------------------------
    typedef typename Antichain::Vertex                            Vertex;
    typedef typename Antichain::Vertex_handle                     Vertex_handle;
    typedef typename Antichain::Edge                              Edge;
    typedef typename Antichain::Edge_handle                       Edge_handle;
    typedef typename Antichain::Face                              Face;
    typedef typename Antichain::Face_handle                       Face_handle;
    // -------------------------------------------------------------------------
    typedef typename Antichain::Flip_traits                       Flip_traits;
    typedef typename Antichain::Ccw_traits                        Ccw_traits;
    typedef typename Antichain::Cw_traits                         Cw_traits;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_vertex_iterator<Self ,Vertex, 
					       Vertex&,Vertex_handle> 
						         Vertex_iterator;
    typedef Visibility_complex_vertex_iterator<Self ,Vertex,
					     const Vertex&,const Vertex_handle>
							 Vertex_const_iterator;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_edge_iterator<Self ,Edge, 
					     Edge&,Edge_handle> 
						         Edge_iterator;
    typedef Visibility_complex_edge_iterator<Self ,Edge,
					     const Edge&,const Edge_handle>
							 Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_face_iterator<Self ,Face, 
					     Face&,Face_handle> 
						         Face_iterator;
    typedef Visibility_complex_face_iterator<Self ,Face,
					     const Face&,const Face_handle>
							 Face_const_iterator;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_linear_sweep_iterator<Antichain,Vertex,
						     Vertex&,Vertex_handle,
						typename Gt::Is_upward_directed>
							 Linear_sweep_iterator;
    typedef Visibility_complex_linear_sweep_iterator<Antichain,Vertex,
					      const Vertex&,const Vertex_handle,
						typename Gt::Is_upward_directed>
						    Linear_sweep_const_iterator;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_sweep_iterator<Antichain,Vertex,
					      Vertex&,Vertex_handle>
							 Sweep_iterator;
    typedef Visibility_complex_sweep_iterator<Antichain,Vertex,
					      const Vertex&,const Vertex_handle>
							 Sweep_const_iterator;
    // -------------------------------------------------------------------------

private:
    Antichain* _antichain;
    std::list<Vertex_handle> vertices;
public:
    // -------------------------------------------------------------------------
    Visibility_complex_2() : _antichain(0) { }
    template <class InputIterator , class ConstraintIt>
    Visibility_complex_2(InputIterator first, InputIterator last,
			 ConstraintIt  firstc,ConstraintIt  lastc);
    template <class InputIterator>
    Visibility_complex_2(InputIterator first, InputIterator last) 
    { 
	std::list<Vertex> L;
	*this = Visibility_complex_2(first,last,L.begin(),L.end());
    }
    // -------------------------------------------------------------------------
    // Functions from Antichain
    Antichain* antichain() const      { return _antichain; }
    void purge(Vertex_handle v)       { }
    void sweep(Vertex_handle v)       { antichain()->sweep(v); }
    Vertex_handle pop_minimal(bool b) { return antichain()->pop_minimal(b); }
    bool is_valid()        const      { return antichain()->is_valid(); }
    int size()             const      { return vertices.size(); }
    bool is_on_convex_hull(Vertex_handle v) const 
    { return antichain()->is_on_convex_hull(v); }
    // -------------------------------------------------------------------------
    Vertex_iterator       vertices_begin()       
    { return Vertex_iterator(vertices.begin());  }
    Vertex_const_iterator vertices_begin() const 
    { return Vertex_iterator(vertices.begin());  }
    Vertex_iterator       vertices_end()         
    { return Vertex_iterator(vertices.end());    }
    Vertex_const_iterator vertices_end()   const 
    { return Vertex_iterator(vertices.end());    }
    // -------------------------------------------------------------------------
    Edge_iterator edges_begin() { return Edge_iterator(vertices.begin()); }
    Edge_const_iterator edges_begin() const 
	{ return Edge_iterator(vertices.begin()); }
    Edge_iterator edges_end()   { return Edge_iterator(vertices.end());   }
    Edge_const_iterator edges_end() const  
	{ return Edge_iterator(vertices.end());   }
    // -------------------------------------------------------------------------
    Face_iterator faces_begin() { return Face_iterator(vertices.begin()); }
    Face_const_iterator faces_begin() const
	{ return Face_iterator(vertices.begin()); }
    Face_iterator faces_end()   { return Face_iterator(vertices.end());   }
    Face_const_iterator faces_end() const
	{ return Face_iterator(vertices.end());   }
    // -------------------------------------------------------------------------
    Vertex_handle find(const Vertex& v);
    // -------------------------------------------------------------------------
    // Views  FIXME
    Disk_handle     left_back_view(Vertex_handle b);
    Disk_handle     left_front_view(Vertex_handle b);
    Disk_handle     right_back_view(Vertex_handle b);
    Disk_handle     right_front_view(Vertex_handle b);
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class InputIterator , class ConstraintIt >
Visibility_complex_2<_Gtr,It,Flip>::
Visibility_complex_2(InputIterator first, InputIterator last ,
		     ConstraintIt  firstc,ConstraintIt  lastc)
{
    // -------------------------------------------------------------------------
    // The antichain used to sweep the complex.
    _antichain = new Antichain(first,last,firstc,lastc);
    _antichain->set_linear_space(false);
    Sweep_iterator v(_antichain) , vend(_antichain,0);
    for ( ; v != vend ; ++v ) vertices.push_back(&(*v));
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
// search for v in the visibility complex. Takes time O(n).
template< class _Tr , class It , class Flip >
typename Visibility_complex_2<_Tr,It,Flip>::Vertex_handle
Visibility_complex_2<_Tr,It,Flip>::find(const Vertex& v)
{
    // -------------------------------------------------------------------------
    // Find an edge on the signed disk v.source_object() 
    typename Antichain::Edge_iterator e = _antichain->edges_begin(); 
    while (e != _antichain->edges_end() && 
	   (v.source_object() != e->object() || v.is_left_xx() != e->sign()))
	e++;
    if (e == _antichain->edges_end()) return 0;
    // -------------------------------------------------------------------------
    // iterate over the vertices tangent to v.source_object()
    typename Ccw_traits::Cc cc;
    Vertex_handle w = e->sup();
    while (*w != v) {
	w = cc(w,e->object());
	if (w == e->sup()) return 0;
    }
    return w;
    // -------------------------------------------------------------------------
}

// -------------------------------------------------------------------------
// -------------------------- Views ----------------------------------------
// -------------------------------------------------------------------------

template< class _Tr , class It , class Flip >
typename Visibility_complex_2<_Tr,It,Flip>::Disk_handle
Visibility_complex_2<_Tr,It,Flip>::right_front_view(Vertex_handle b)
{
    return (b->is_xx_right()) ? b->target_object() :
				b->ccw_target_edge()->dr()->front_object();
}

template< class _Tr , class It , class Flip >
typename Visibility_complex_2<_Tr,It,Flip>::Disk_handle
Visibility_complex_2<_Tr,It,Flip>::right_back_view(Vertex_handle b)
{
    return (b->is_right_xx()) ? b->source_object() :
				b->cw_source_edge()->dr()->back_object();
}

template< class _Tr , class It , class Flip >
typename Visibility_complex_2<_Tr,It,Flip>::Disk_handle
Visibility_complex_2<_Tr,It,Flip>::left_front_view(Vertex_handle b)
{
    return (b->is_xx_left()) ? b->target_object() :
			       b->cw_target_edge()->ur()->front_object();
}

template< class _Tr , class It , class Flip >
typename Visibility_complex_2<_Tr,It,Flip>::Disk_handle
Visibility_complex_2<_Tr,It,Flip>::left_back_view(Vertex_handle b)
{
    return (b->is_left_xx()) ? b->source_object() :
			       b->ccw_source_edge()->ur()->back_object();
}

// -----------------------------------------------------------------------------

template <class InputIterator, class OutputIterator,
	  class Traits, class Items>
OutputIterator
visibility_complex_2(InputIterator first,
		     InputIterator last, 
		     OutputIterator result, 
		     Visibility_complex_antichain<Traits,Items>& a)
{
    typedef Visibility_complex_antichain<Traits,Items> Antichain;
    typename Antichain::Vertex_handle min = a.pop_minimal(true);
    while (min != 0) {
	*result++ = *min;
	a.sweep(min);
	min = a.pop_minimal(true);
    }
    return result;
}

// -----------------------------------------------------------------------------

#ifdef VISIBILITY_COMPLEX_POLYGON_TRAITS_H
template <class InputIterator, class OutputIterator,
	  class _R, class Items>
OutputIterator
__visibility_complex_2(InputIterator first, InputIterator last, 
		       OutputIterator result, 
		       Polygon_2<Polygon_traits_2<_R>,std::list<Point_2<_R> > >, 
		       Items)
{
    typedef Visibility_complex_polygon_traits<_R> Traits;
    Visibility_complex_antichain<Traits,Items> a(first,last);
    visibility_complex_2(first,last,result,a);
    return result;
}
#endif

#ifdef VISIBILITY_COMPLEX_CIRCLE_TRAITS_H
template <class InputIterator, class OutputIterator,
	  class _R, class Items>
OutputIterator
__visibility_complex_2(InputIterator first, InputIterator last, 
		       OutputIterator result, Circle_2<_R>, Items)
{
    typedef Visibility_complex_circle_traits<_R> Traits;
    Visibility_complex_antichain<Traits,Items> a(first,last);
    visibility_complex_2(first,last,result,a);
    return result;
}
#endif

template <class InputIterator, class OutputIterator, class Items>
OutputIterator
__visibility_complex_2(InputIterator first, InputIterator last, 
		       OutputIterator result, Items t)
{
    typedef typename 
      std::iterator_traits<InputIterator>::value_type object_type;

    __visibility_complex_2(first,last,result,object_type(),t);
    return result;
}

template <class InputIterator, class OutputIterator>
OutputIterator
visibility_complex_2(InputIterator first, InputIterator last, 
		     OutputIterator result)
{
    __visibility_complex_2(first,last,result,Visibility_complex_items());
    return result;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
