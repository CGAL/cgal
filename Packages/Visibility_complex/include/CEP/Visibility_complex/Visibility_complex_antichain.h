#ifndef VISIBILITY_COMPLEX_ANTICHAIN_H
#define VISIBILITY_COMPLEX_ANTICHAIN_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_FUNCTION_OBJECTS_H
#include <CGAL/function_objects.h>
#endif

#ifndef CGAL_IN_PLACE_LIST_H
#include <CGAL/In_place_list.h>
#endif

#ifndef VISIBILITY_COMPLEX_ITEMS_H
#include <CEP/Visibility_complex/Visibility_complex_items.h>
#endif 

#ifndef VISIBILITY_COMPLEX_FLIP_TRAITS_H
#include <CEP/Visibility_complex/Visibility_complex_flip_traits.h>
#endif 

#ifndef VISIBILITY_COMPLEX_ANTICHAIN_ITERATORS_H
#include <CEP/Visibility_complex/Visibility_complex_antichain_iterators.h>
#endif 

#ifndef VISIBILITY_COMPLEX_CCW_CW_TRAITS_H
#include <CEP/Visibility_complex/Visibility_complex_ccw_cw_traits.h>
#endif 

#ifndef VISIBILITY_COMPLEX_FUNCTION_OBJECTS
#include <CEP/Visibility_complex/Visibility_complex_function_objects.h>
#endif

#ifndef VISIBILITY_COMPLEX_SWEEP_ITERATOR_H
#include <CEP/Visibility_complex/Visibility_complex_sweep_iterator.h>
#endif

#include <queue>
#include <list>
#include <set>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class Vertex_base>
class Visibility_complex_vertex
    : public Vertex_base, 
      public In_place_list_base< Visibility_complex_vertex<Vertex_base> > 
{
public:
    typedef Visibility_complex_vertex< Vertex_base> Self;
    typedef typename Vertex_base::Vertex_handle     Vertex_handle;
    typedef typename Vertex_base::Edge_handle       Edge_handle;
    typedef typename Vertex_base::Disk_handle    Disk_handle;
    typedef typename Vertex_base::Bitangent_2       Bitangent_2;
    typedef typename Vertex_base::Type              Type;

    Visibility_complex_vertex() {}
    Visibility_complex_vertex(Type t , Disk_handle start , 
				       Disk_handle finish) 
	: Vertex_base(t,start,finish) {}
    Visibility_complex_vertex(Edge_handle start , Edge_handle finish)
	: Vertex_base(start,finish)   {}
    Visibility_complex_vertex(const Bitangent_2& b) 
	: Vertex_base(b) {}
    Visibility_complex_vertex( const Vertex_base& v)   // down cast
        : Vertex_base(v) {}
    Self& operator=( const Self& v) {
        this->Vertex_base::operator=(v);
        return *this;
    }
    void set_color(bool b) { _color = b; }
    bool color() const { return _color; }
private:
    bool _color;
};

// -----------------------------------------------------------------------------

template <class Edge_base>
class Visibility_complex_edge
    : public Edge_base, 
      public In_place_list_base< Visibility_complex_edge<Edge_base> > 
{
public:
    typedef Visibility_complex_edge<Edge_base> Self;
    typedef typename Edge_base::Vertex_handle  Vertex_handle;
    typedef typename Edge_base::Disk_handle Disk_handle;

    Visibility_complex_edge() {}                   
    Visibility_complex_edge(bool s, Disk_handle p)
	: Edge_base(s,p) {}
    Visibility_complex_edge(Vertex_handle v0 , Vertex_handle v1 , 
			    Disk_handle p)
	: Edge_base(v0,v1,p) {}
    Visibility_complex_edge( const Edge_base& h)
        : Edge_base(h) {}
    Self& operator=( const Self& e) {
        this->Edge_base::operator=(e);
        return *this;
    }
};

// -----------------------------------------------------------------------------

template < class Face_base>
class Visibility_complex_face
    : public Face_base, public In_place_list_base<
                            Visibility_complex_face< Face_base> > {
public:
    typedef Visibility_complex_face< Face_base>  Self;
    Visibility_complex_face() {}
    Visibility_complex_face(const Face_base& f) : Face_base(f) {}
    Self& operator=( const Self& f) {
        this->Face_base::operator=(f);
        return *this;
    }
};

// -----------------------------------------------------------------------------

template < class _Gtr , 
	   class It   = Visibility_complex_items ,
	   class Flip = Visibility_complex_flip_traits >
class Visibility_complex_antichain 
 : public In_place_list< Visibility_complex_edge<typename It::template Edge_wrapper< Visibility_complex_antichain<_Gtr, It, Flip> >::Edge>, false>
{
public:
    // -------------------------------------------------------------------------
    typedef _Gtr                                          Gt;
    typedef typename _Gtr::Disk                           Disk; 
    typedef typename _Gtr::Point_2                        Point_2;
    typedef typename _Gtr::Bitangent_2                    Bitangent_2;
    typedef Bitangent_2                                   BT;
    typedef typename BT::Disk_handle                      Disk_handle;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_antichain<_Gtr,It,Flip>    Self;
    typedef Visibility_complex_antichain<_Gtr,It,Flip>    Antichain;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_left_ccw_traits<Self>      Left_ccw_traits;
    typedef Visibility_complex_right_ccw_traits<Self>     Right_ccw_traits;
    typedef Visibility_complex_left_cw_traits<Self>       Left_cw_traits;
    typedef Visibility_complex_right_cw_traits<Self>      Right_cw_traits;
    typedef Right_ccw_traits                              Ccw_traits;
    typedef Right_cw_traits                               Cw_traits;
    // -------------------------------------------------------------------------
    /*
    typedef typename Visibility_complex_backward_flip_traits:: Flip_wrapper<Self>
							  Flip_wrapper;
    */
    typedef typename  Flip::template Flip_wrapper<Self>             Flip_wrapper;
  typedef  Flip_wrapper FW;
    typedef typename FW::Flip_traits            Flip_traits;
    // -------------------------------------------------------------------------
    typedef typename It::template Vertex_wrapper<Self>             Vertex_wrapper;
    typedef typename Vertex_wrapper::Vertex               Vertex_base;
    typedef Visibility_complex_vertex< Vertex_base>       Vertex;
    typedef Vertex*                                       Vertex_handle;
    typedef typename In_place_list<Vertex,false>::iterator         Minimals_iterator;
    // -------------------------------------------------------------------------
    typedef typename It::template Edge_wrapper<Self>               Edge_wrapper;
    typedef typename Edge_wrapper::Edge                   Edge_base;
    typedef Visibility_complex_edge< Edge_base>           Edge;
    typedef Edge*                           		  Edge_handle;
    typedef typename In_place_list<Edge,false>::iterator           Edge_iterator;
    typedef typename In_place_list<Edge,false>::reverse_iterator   Edge_reverse_iterator;
    typedef typename In_place_list<Edge,false>::const_iterator     Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef typename It::template Face_wrapper<Self>               Face_wrapper;
    typedef typename Face_wrapper::Face                   Face_base;
    typedef Visibility_complex_face< Face_base>           Face;
    typedef Face*                           		  Face_handle;
    typedef Vc_antichain_face_iterator<Self,Face,
				       Face&,Face_handle> Face_iterator;
    typedef Vc_antichain_face_iterator<Self,Face,
				       const Face&,const Face_handle> 
							  Face_const_iterator;
    // -------------------------------------------------------------------------
    typedef Vc_antichain_vertex_iterator<Self,Vertex, Vertex&,
					 Vertex_handle, typename Ccw_traits::Sup> 
						       Vertex_iterator;
    typedef Vc_antichain_vertex_iterator<Self,Vertex, const Vertex&,
				        const Vertex_handle, 
                                        typename Ccw_traits::Sup> 
						       Vertex_const_iterator;
    typedef Vc_antichain_vertex_iterator<Self,Vertex, Vertex&,
                                         Vertex_handle, typename Cw_traits::Sup> 
						       Vertex_cw_iterator;
    typedef Vc_antichain_vertex_iterator<Self,Vertex, const Vertex&,
				         const Vertex_handle,
                                         typename Cw_traits::Sup> 
						       Vertex_cw_const_iterator;
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
    typedef typename BT::Type_util               Type_util;
    // -------------------------------------------------------------------------
private:
    bool                        _valid;
    bool                        _straight_sweep;
    bool                        _linear_space;
    Face_handle                 _infinite_face;
    In_place_list<Vertex,false> _minimals_ccw;
    In_place_list<Vertex,false> _minimals_cw;

#ifdef DEBUG
    map<long,int>               Key;
#endif

public :

    // -------------------------------------------------------------------------
    Visibility_complex_antichain() 
	: _valid(false) , _straight_sweep(false) , _linear_space(true)
    { }
    template < class InputIterator ,class ConstraintIt >
    Visibility_complex_antichain(InputIterator first, InputIterator last,
				 ConstraintIt  firstc,ConstraintIt lastc); 
    ~Visibility_complex_antichain() { destroy(); }
    // -------------------------------------------------------------------------
    // Options when sweeping
    bool is_valid()               const { return _valid;          }
    bool is_straight()            const { return _straight_sweep; }
    void make_straight()                { _straight_sweep = true; }
    bool linear_space()           const { return _linear_space;   }
    void set_linear_space(bool b)       { _linear_space = b;      }
    // -------------------------------------------------------------------------
    // Identifying convex-hull vertices
    Face_handle infinite_face() const   { return _infinite_face;  }
    void create_infinite_face() { _infinite_face = new Face;  }
    bool is_on_convex_hull(Vertex_handle v) const;
    // -------------------------------------------------------------------------
    // Iterator pairs for trversing the sink of the faces of the antichain
    // These vertices form the Greedy pseudo-triangulation
    Vertex_iterator       vertices_begin()
	{ return Vertex_iterator(this,begin()); } 
    Vertex_const_iterator vertices_begin() const
	{ return Vertex_const_iterator(this,begin()); } 
    Vertex_iterator   vertices_end()
	{ return Vertex_iterator(this,end()); } 
    Vertex_const_iterator   vertices_end() const
	{ return Vertex_const_iterator(this,end()); } 
    // -------------------------------------------------------------------------
    // Iterator pairs for trversing the sources of the faces of the antichain
    // These vertices form the dual Greedy pseudo-triangulation
    Vertex_cw_iterator       cw_vertices_begin()
	{ return Vertex_cw_iterator(this,begin()); } 
    Vertex_cw_const_iterator cw_vertices_begin() const
	{ return Vertex_cw_const_iterator(this,begin()); } 
    Vertex_cw_iterator   cw_vertices_end()
	{ return Vertex_cw_iterator(this,end()); } 
    Vertex_cw_const_iterator   cw_vertices_end() const
	{ return Vertex_cw_const_iterator(this,end()); } 
    // -------------------------------------------------------------------------
    Edge_iterator       edges_begin()       { return begin(); } 
    Edge_const_iterator edges_begin() const { return begin(); } 
    Edge_iterator       edges_end()         { return end();   } 
    Edge_const_iterator edges_end()   const { return end();   } 
    Edge_reverse_iterator edges_rbegin()    { return rbegin(); } 
    Edge_reverse_iterator edges_rend()      { return rend(); } 
    // -------------------------------------------------------------------------
    Face_iterator       faces_begin()     {return Face_iterator(this,begin());}
    Face_const_iterator faces_begin()const{return Face_iterator(this,begin());}
    Face_iterator       faces_end()       { return Face_iterator(this,end()); }
    Face_const_iterator faces_end()  const{ return Face_iterator(this,end()); }
    // -------------------------------------------------------------------------
    Minimals_iterator minimals_begin()    { return _minimals_ccw.begin(); }
    Minimals_iterator minimals_end()      { return _minimals_ccw.end(); }
    Minimals_iterator cw_minimals_begin() { return _minimals_cw.begin(); }
    Minimals_iterator cw_minimals_end()   { return _minimals_cw.end();   }
    Minimals_iterator minimals_begin(Ccw_traits) { return minimals_begin();    }
    Minimals_iterator minimals_end  (Ccw_traits) { return minimals_end();      }
    Minimals_iterator minimals_begin(Cw_traits)  { return cw_minimals_begin(); }
    Minimals_iterator minimals_end  (Cw_traits)  { return cw_minimals_end();   }
    // -------------------------------------------------------------------------
    // Iterator pair for linear sweep
    Linear_sweep_iterator       sweep_begin() 
	{ return Linear_sweep_iterator(this); }
    Linear_sweep_const_iterator sweep_begin() const 
	{ return Linear_sweep_const_iterator(this); }
    Linear_sweep_iterator       sweep_end()
	{ return Linear_sweep_iterator(this,0); }
    Linear_sweep_const_iterator sweep_end() const 
	{ return Linear_sweep_const_iterator(this,0); }
    // -------------------------------------------------------------------------
    // Testing minimality
    template < class _Tr > 
    bool is_minimal(const Vertex_handle& v, _Tr tr) const;
    bool is_minimal(const Vertex_handle& v) const
    { return is_minimal(v,Ccw_traits()); }
    template < class _Tr > 
    bool is_xx_minimal(const Vertex_handle& v,_Tr tr) const;
    bool is_right_minimal(const Vertex_handle& v) const
    { return is_xx_minimal(v,Right_ccw_traits()); }
    bool is_left_minimal(const Vertex_handle& v)  const
    { return is_xx_minimal(v,Left_ccw_traits()); }
    // -------------------------------------------------------------------------
    Vertex_handle pop_minimal(bool finite = false);
    template < class _Tr > 
    bool is_swept_regular(const Vertex_handle& v , _Tr tr) const;
    template < class _Tr > 
    bool is_swept_constraint(const Vertex_handle& v , _Tr tr) const;
    template < class _Tr > 
    bool is_swept(const Vertex_handle& v , _Tr tr) const;
    bool is_swept(const Vertex_handle& v) const 
    { return is_swept(v,Ccw_traits()); }
    // -------------------------------------------------------------------------
    // Adding a minimal
    void push_back_minimal(Vertex_handle v ,Ccw_traits);
    void push_back_minimal(Vertex_handle v ,Cw_traits);
    void push_back_minimal(Vertex_handle v) {push_back_minimal(v,Ccw_traits());}
    // -------------------------------------------------------------------------
    // Removing a minimal
    void erase_minimal(Vertex_handle v, Ccw_traits) { _minimals_ccw.erase(v); }
    void erase_minimal(Vertex_handle v, Cw_traits)  { _minimals_cw.erase(v);  }
    void erase_minimal(Vertex_handle v) { erase_minimal(v, Ccw_traits()); }
    // Depreciated - for backward compatibility
    template < class _Tr > 
    void erase_minimal(Minimals_iterator v,_Tr tr) { erase_minimal(&(*v),tr); }
    void erase_minimal(Minimals_iterator v) { erase_minimal(v,Ccw_traits()); }
    // -------------------------------------------------------------------------
    // Removing all minimals
    template < class _Tr > void clear_minimals(_Tr tr);
    void clear_minimals() { clear_minimals(Ccw_traits()); }
    // -------------------------------------------------------------------------
    // Sweeping - flipping operations
    template < class _Tr > void sweep(Vertex_handle v, _Tr tr);
    template < class _Tr > void sweep_good(Vertex_handle v , _Tr tr);
    template < class _Tr > void sweep_all_minimals(_Tr tr);
    template < class _Tr > void sweep_good_all_minimals(_Tr tr);
    void sweep(Vertex_handle v) { sweep(v,Ccw_traits()); }
    void sweep_all_minimals()   { sweep_all_minimals(Ccw_traits()); }
    // -------------------------------------------------------------------------
    void set_constraint(Vertex_handle v);
    void unset_constraint(Vertex_handle v);
    void add_constraint(Vertex_handle v);
    void remove_constraint(Vertex_handle v);
    // -------------------------------------------------------------------------
protected:
    // -------------------------------------------------------------------------
    // Compute the flipped bitangent phi(v)
    template < class _Tr >
    Vertex_handle compute_phi(Vertex_handle v , _Tr)  /*const*/;
    // -------------------------------------------------------------------------
    // Update the antichain while adding a constraint
    template < class _Tr > void sweep_constraint(const Vertex_handle& v, const _Tr&) const;
    template < class _Tr > void sweep_regular   (const Vertex_handle& v, const _Tr&) const;
    // -------------------------------------------------------------------------
    // Initialisation methods
    template<class InputIterator , class ConstraintIt> 
    void compute_graph(InputIterator first, InputIterator last,
		       ConstraintIt  firstc,ConstraintIt  lastc);
    template < class _Tr > void compute_minimals(_Tr tr);
    template < class _Tr > void compute_vertices(_Tr tr);
    template < class _Tr > void fix_extreme_edges(_Tr tr) const;
    template < class _Tr > void glue_ccw_cw(_Tr tr) const;
    template < class _Tr > void initialize_convex_hull(_Tr tr) const;
    // -------------------------------------------------------------------------
    // Method used during intialization, to compute candidates during the
    // Bentley-Ottmann rotational sweep.
  //   template<class _Tr> bool is_candidate(const Face_handle& fa) const;
template < class _Tr >
bool
is_candidate(const Face_handle& f, const _Tr&) const
{
    typename _Tr::Dl dl; typename _Tr::Ur ur; typename _Tr::Sup sup;
    if (f == 0 || f->top_edge() == 0 || f->bottom_edge() == 0) return false;
    if (sup(f) != 0 && sup(f)->is_constraint()) return is_minimal(sup(f),_Tr());
    if (f->top_edge()->object() == 0 || f->bottom_edge()->object() == 0)
	return false;
    return (ur(f->bottom_edge()) == f && dl(f->top_edge()) == f);
}
    // -------------------------------------------------------------------------
#ifdef DEBUG
public:
    void print(Vertex_handle v) {
	//cout << v << " " << flush;
	if (v == 0) return; 
	cout << *v << " , (";
	if (v->is_left_xx()) cout << "+"; else cout << "-";
	cout << Key[long(v->source_object())] << ",";
	if (v->is_xx_left()) cout << "+"; else cout << "-";
	cout << Key[long(v->target_object())] << ")" ;
	/*
	cout << " minLR(" << flush 
	     << is_xx_minimal(v,Left_ccw_traits()) << "," << flush 
	     << is_xx_minimal(v,Right_ccw_traits()) << ")";
	cout << " cw_minLR(" << flush 
	     << is_xx_minimal(v,Left_cw_traits()) << "," << flush 
	     << is_xx_minimal(v,Right_cw_traits()) << ")";
	if (is_swept(v)) cout << " swept";
	else cout << " not swept";
	*/
    }
    void print(Edge_handle e, bool coord = true) {
	cout << e << " " << flush;
	if (e == 0) return; 
	if (e->object() == 0) {
	    if (e == e->sup()->target_cusp_edge())
		 cout << "-" << Key[long(e->sup()->target_object())] << "c";
	    else cout << "+" << Key[long(e->sup()->source_object())] << "c";
	}
	else {
	    if (e->sign()) cout << "+"; else cout << "-";
	    cout << Key[long(e->object())] ;
	}
	cout << " {";
	if (e->sign()) 
	     cout << e->dl() << "," 
		  << e->ur() << "," 
		  << e->ul();
	else cout << e->dl() << "," 
		  << e->dr() << "," 
		  << e->ul();
	cout << "}";
	if (coord) {
	    cout << " [";
	    cout << *e->begin() << "," << *--e->end() << "]";
	}
    }
    void print(Face_handle f) {
	cout << f << " " << flush;
	if (f == 0) return; 
	cout << "[" ; 
	print(f->bottom_edge()); cout << ",";
	print(f->top_edge());    cout << "]";
	cout << " inf = " << f->inf() << " , sup = " << f->sup();
    }
    template < class _Tr >
    void print_left_minimal(Vertex_handle v ,_Tr tr);
    template < class _Tr >
    void print_right_minimal(Vertex_handle v ,_Tr tr);
#endif
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------
// Computing the Antichain using a Bentley-Ottman sweep

template < class _Gtr , class It , class Flip >
template < class InputIterator , class ConstraintIt >
Visibility_complex_antichain<_Gtr,It,Flip>::
Visibility_complex_antichain(InputIterator first, InputIterator last ,
			     ConstraintIt  firstc,ConstraintIt  lastc)
{
    // -------------------------------------------------------------------------
    // Default values : linear space topological sweep
    _straight_sweep = false;
    _linear_space = true;
    //_linear_space = false;
    // -------------------------------------------------------------------------
    if (first == last) { _valid = false; return; }
    compute_graph(first,last,firstc,lastc);
    if (!is_valid()) return;
    // -------------------------------------------------------------------------
    compute_vertices(Ccw_traits());
    compute_vertices(Cw_traits());
    // -------------------------------------------------------------------------
    fix_extreme_edges(Ccw_traits());
    fix_extreme_edges(Cw_traits());
    // -------------------------------------------------------------------------
    glue_ccw_cw(Ccw_traits());
    glue_ccw_cw(Cw_traits());
    // -------------------------------------------------------------------------
    initialize_convex_hull(Ccw_traits());
    initialize_convex_hull(Cw_traits());
    // -------------------------------------------------------------------------
    compute_minimals(Ccw_traits());
    compute_minimals(Cw_traits());
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class InputIterator , class ConstraintIt >
void
Visibility_complex_antichain<_Gtr,It,Flip>::
compute_graph(InputIterator first, InputIterator last,
	      ConstraintIt  firstc,ConstraintIt  lastc)
{
    // -------------------------------------------------------------------------
    // The pseudo-triangulation is valid if no pair of objects intersect
    CGAL_expensive_precondition_code( Do_intersect<Self> do_intersect; );
    _valid = true;
    // -------------------------------------------------------------------------
    // X-structure containing the nodes of the antichain == edges of 
    // Visibility graph
    typedef std::priority_queue<Edge_handle,
    				std::vector<Edge_handle>,
				Less_edge_handle<_Gtr> >      Xstructure;
    Xstructure XE;
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Map to keep track of the edge with opposite sign
    std::map<long,Edge_handle> opposite;
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Creating the edges and pushing them in the priority queue
    typename Ccw_traits::Set_adjacent_faces_one_to_one set_adjacent_faces;
    // Regular edges - two per disk
    for (InputIterator it = first; it != last ; ++it ) {
	Edge_handle neg = new Edge(false,&(*it)); 
	Edge_handle pos = new Edge(true,&(*it));
	push_back(*neg); push_back(*pos); 
	XE.push(neg);  XE.push(pos);
	set_adjacent_faces(neg,new Face,new Face,0);
	set_adjacent_faces(pos,new Face,0,0);
	neg->dl()->set_front_view(neg); 
	neg->dr()->set_back_view(pos);
	opposite[long(neg)] = pos;
	opposite[long(pos)] = neg;
	neg->set_opposite(pos);
    }
    // Constraint edges - two per constraint
    // We make sure the constraint is upward directed
    typename _Gtr::Is_upward_directed is_upward_directed;
    for (ConstraintIt c = firstc; c != lastc ; ++c ) {
	CGAL_precondition(c->sup() == 0);
	if (!c->is_constraint())       set_constraint(&(*c));
	if (!c->pi()->is_constraint()) set_constraint(c->pi());
	Vertex_handle d = (is_upward_directed(*c)) ? &(*c) : c->pi();
	push_back(*d->target_cusp_edge());
	push_back(*d->source_cusp_edge());
	XE.push(d->target_cusp_edge());  
	XE.push(d->source_cusp_edge());
    }
    if (XE.size() == 2) { _valid = false; return; }
    // -------------------------------------------------------------------------
#ifdef DEBUG
    int index = 1;
    for (Edge_iterator e = edges_begin(); e != edges_end(); ++e) 
    {
	if (!e->sign()) {
	    Key[long(e->object())] = index; 
	    Key[long(&(*e))] = -index;
	}
	else {
	    Key[long(&(*e))] = +index; 
	    ++index;
	}
    }
#endif
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Y-structure containing at each moment the ordered list of object
    // intersected by the sweeping line.
    // The faces are ordered by their front view
    // An object is represented by a pair of faces (f1,f2) such that
    // f1->front_object() == f2->back_object()
    typedef std::set<Face_handle,Less_face_handle<Self> >  Ystructure;
    Ystructure YE;
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Splice function object used when inserting a constraint
    typename Ccw_traits::Splice splice_ccw;
    typename Cw_traits::Splice  splice_cw;
    // -------------------------------------------------------------------------
    // The infinite face helps us to identify the convex-hull vertices
    _infinite_face = new Face;
    _infinite_face->set_bottom_edge(XE.top());
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Start of the Bentley-Ottman sweep
    Face_handle infinite = NULL;
    while ( !XE.empty() && _valid ) {
	Edge_handle pe = XE.top(); XE.pop();
	// ---------------------------------------------------------------------
	// Treating a node with in_degree 2 and out_degree 1
	if (!pe->sign()) {
	    // -----------------------------------------------------------------
	    // Inserting a Regular Edge
	    if (pe->object() != 0) { // regular Edge
		// -------------------------------------------------------------
		// Insert two new faces in YE
		Face_handle upleft = new Face; upleft->set_front_view(pe);
		typename Ystructure::iterator upleft_it = YE.upper_bound(upleft);
		delete upleft;
		if ( upleft_it == YE.end() ) upleft = infinite;
		else { upleft = *upleft_it; YE.erase(upleft_it); }
		set_adjacent_faces(pe,pe->dl(),pe->dr(),upleft);
		if (upleft != 0) {
		    // ---------------------------------------------------------
		    CGAL_expensive_precondition_code(
			_valid = !do_intersect(pe->ul()->back_view(),pe);
			if (_valid == false) return;
			_valid = !do_intersect(pe->ul()->front_view(),pe);
			if (_valid == false) return;
		    );
		    // ---------------------------------------------------------
		    pe->dr()->set_front_view(pe->ul()->front_view());
		    pe->dl()->set_back_view (pe->ul()->back_view());
		}
		YE.insert(pe->dl()); YE.insert(pe->dr());
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	    // Inserting a cusp Edge from an xx-left constraint
	    else if (pe->sup()->is_xx_left()) { 
		// -------------------------------------------------------------
		// Find the place where to insert the constraint
		Face_handle right = new Face; right->set_front_view(pe);
		typename Ystructure::iterator right_it = YE.upper_bound(right);
		CGAL_precondition( right_it != YE.end() );
		delete right; right = *right_it;
		// -------------------------------------------------------------
		CGAL_expensive_precondition_code(
		    _valid = !do_intersect(right->front_view(),pe);
		    if (_valid == false) return;
		);
		// -------------------------------------------------------------
		// Splitting the boundary of the target object of pe->sup()
		CGAL_precondition(right->back_view()->object() != 0);
		splice_ccw(right->back_view(),pe->sup());
		// -------------------------------------------------------------
		pe->dl()->set_back_view (pe->sup()->cw_target_edge());
		pe->dl()->set_front_view(pe);
		right->set_back_view(pe);
		YE.insert(pe->dl());
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	    // Inserting a cusp Edge from an xx-right constraint
	    else {
		// -------------------------------------------------------------
		// Find the place where to insert the constraint
		Face_handle left = new Face; left->set_front_view(pe);
		typename Ystructure::iterator left_it = YE.upper_bound(left);
		CGAL_precondition( left_it != YE.end() );
		delete left; left = *left_it;
		YE.erase(left_it);
		// -------------------------------------------------------------
		CGAL_expensive_precondition_code(
		    _valid = !do_intersect(left->back_view(),pe);
		    if (_valid == false) return;
		);
		// -------------------------------------------------------------
		CGAL_precondition(left->front_view()->object() != 0);
		splice_ccw(left->front_view(),pe->sup());
		// -------------------------------------------------------------
		pe->dr()->set_back_view (pe);
		pe->dr()->set_front_view(pe->sup()->ccw_target_edge());
		left->set_front_view(pe);
		YE.insert(left); YE.insert(pe->dr());
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	}
	// ---------------------------------------------------------------------
	// Treating a node with in_degree 1 and out_degree 2
	else {
	    // -----------------------------------------------------------------
	    // Inserting a regular Edge
	    if (pe->object() != 0) {
		// -------------------------------------------------------------
		// Find the two adjacent faces in YE that connect to pe and
		// erase them from YE.
		Face_handle tmpf = new Face; tmpf->set_front_view(pe);
		typename Ystructure::iterator it1 = YE.find(tmpf);
		typename Ystructure::iterator it2 = it1; ++it2;
		Face_handle upleft = *it1; Face_handle upright = *it2;
		delete tmpf; YE.erase(it1); YE.erase(it2);
		// -------------------------------------------------------------
		// Update existing faces
		if ( XE.empty() ) {
		    delete pe->dl();
		    set_adjacent_faces(pe,0,upright,upleft);
		}
		else { 
		    set_adjacent_faces(pe,pe->dl(),upright,upleft);
		    pe->dl()->set_front_view(upright->front_view());
		    pe->dl()->set_back_view (upleft->back_view());
		    // ---------------------------------------------------------
		    CGAL_expensive_precondition_code(
			_valid = !do_intersect(pe->dl()->back_view(),
					       pe->dl()->front_view());
			if (_valid == false) return;
		    );
		    // ---------------------------------------------------------
		    YE.insert(pe->dl());
		}
		// -------------------------------------------------------------
		infinite = pe->dl();
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	    // Inserting a tail cusp Edge from an left-xx constraint
	    else if (pe->sup()->is_left_xx()) {
		// -------------------------------------------------------------
		// Find the face in YE with pe->sup() as a front view.
		Face_handle upleft = new Face; 
		upleft->set_front_view(pe->sup()->target_cusp_edge());
		typename Ystructure::iterator upleft_it = YE.find(upleft);
		typename Ystructure::iterator right = upleft_it; ++right;
		delete upleft; upleft = *upleft_it;
		YE.erase(upleft_it);
		// -------------------------------------------------------------
		// Update degenerate face emanating from edge.
		upleft->set_sup(pe->sup());
		upleft->set_inf(pe->sup()->pi());
		CGAL_precondition(upleft->back_view()->object() != 0);
		splice_ccw(upleft->back_view(),pe->sup());
		// -------------------------------------------------------------
		/* delete pe->ul(); */ set_adjacent_faces(pe,0,0,upleft);
		CGAL_precondition(pe->dl() == 0 && pe->ur() == 0);
		// -------------------------------------------------------------
		// Update the back view of right.
		(*right)->set_back_view(pe->sup()->cw_source_edge());
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	    // Inserting a tail cusp Edge from a right-xx constraint
	    else {
		// -------------------------------------------------------------
		// Find the face in YE with pe->sup() as a front view.
		Face_handle left = new Face; 
		left->set_front_view(pe->sup()->target_cusp_edge());
		typename Ystructure::iterator left_it = YE.find(left);
		typename Ystructure::iterator upright = left_it; ++upright;
		delete left; left = *left_it;
		// -------------------------------------------------------------
		// Set upright as the degenerate face of the edge pe.
		(*upright)->set_sup(pe->sup());
		(*upright)->set_inf(pe->sup()->pi());
		CGAL_precondition((*upright)->front_view()->object() != 0);
		splice_ccw((*upright)->front_view(),pe->sup());
		// -------------------------------------------------------------
		/* delete pe->ur(); */ set_adjacent_faces(pe,0,*upright,0);
		CGAL_precondition(pe->dl() == 0 && pe->ul() == 0);
		// -------------------------------------------------------------
		// Update the back view of left and insert the face.
		YE.erase(left_it); YE.erase(upright); 
		left->set_front_view(pe->sup()->ccw_source_edge());
		YE.insert(left); 
		// -------------------------------------------------------------
	    }
	    // -----------------------------------------------------------------
	}
	if (XE.empty()) _infinite_face->set_top_edge(pe);
    } // end while ( !XE.empty() )
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Splicing the boundary of objects by inserting the pi of the constraints
    // The fact that the edges come in pairs allows us to recover the edge with
    // opposite sign on the same disk.
    for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e) {
	if (e->object() != 0 && e->sup() != 0) {
	    Edge_handle eo = opposite[long(&(*e))];
	    CGAL_precondition(eo == e->opposite());

	    Vertex_handle v = e->sup();
	    splice_cw(eo,v->pi());
	    while (v->ccw_edge(e->object())->sup() != 0) {
		splice_cw(v->pi()->ccw_edge(e->object()),
			  v->ccw_edge(e->object())->sup()->pi());
		v = v->ccw_edge(e->object())->sup();
	    }
	}
    }
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
/*
// put inline 
template < class _Gtr , class It , class Flip >
template < class _Tr >
bool
Visibility_complex_antichain<_Gtr,It,Flip>::
is_candidate(const Face_handle& f) const
{
    typename _Tr::Dl dl; typename _Tr::Ur ur; typename _Tr::Sup sup;
    if (f == 0 || f->top_edge() == 0 || f->bottom_edge() == 0) return false;
    if (sup(f) != 0 && sup(f)->is_constraint()) return is_minimal(sup(f),_Tr());
    if (f->top_edge()->object() == 0 || f->bottom_edge()->object() == 0)
	return false;
    return (ur(f->bottom_edge()) == f && dl(f->top_edge()) == f);
}
*/
// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
void
Visibility_complex_antichain<_Gtr,It,Flip>::compute_vertices(_Tr tr)
{
    // -------------------------------------------------------------------------
    // All the operators from this methods are taken from the algorithm traits 
    // class _Tr. Two traits classes are given : Ccw_traits and Cw_traits to
    // compute G and G_* respectively
    typename _Tr::Sup     sup; typename _Tr::Inf     inf;
    typename _Tr::Set_sup set_sup;// typename _Tr::Set_inf set_inf;
    typename _Tr::Vertex_creator vc; 
    typename _Tr::Set_adjacent_faces set_adjacent_old_faces;
    typename _Tr::Set_adjacent_faces_one_to_one  set_adjacent_faces;
    typename _Tr::Dr      dr;  typename _Tr::Dl      dl;
    typename _Tr::Ur      ur;  typename _Tr::Ul      ul;
    typename _Tr::Cw_target_edge cw_target_edge;
    typename _Tr::Cw_source_edge cw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Splice splice;
    typedef _Tr TR;
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Saving the Edge --> Face pointers because we will loose them during the
    // rotational sweep below
    Self a;
    for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e) {
	Edge_handle tmp = new Edge; tmp->set_sign(e->sign());
	a.push_back(*tmp);
	if (e->sign()) 
	     set_adjacent_old_faces(&a.back(),dl(&(*e)),ur(&(*e)),ul(&(*e)));
	else set_adjacent_old_faces(&a.back(),dl(&(*e)),dr(&(*e)),ul(&(*e)));
    }
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // death is the event queue. 
    // It contains the vertices which are about to appear.
    // Pushing first candidates to initialize death vertices
    std::set<Vertex_handle,typename _Tr::Less_vertex_handle> death;

    for (Face_iterator fit = faces_begin() ; fit != faces_end() ; ++fit) {
	Face_handle f = &(*fit);
	if (is_candidate(f, _Tr())) { // af: gcc needs it, for bcc it must not be there
	    if (sup(f) == 0) set_sup(f,vc(f->bottom_edge(),f->top_edge()));
	    death.insert(sup(f));
	}
    }
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Rotational sweep € la Bentley-Otman
    while ( death.size() != 0 ) {
	Vertex_handle vmin = *death.begin();
	death.erase(death.begin());
	// ---------------------------------------------------------------------
	Edge_handle top = (vmin->is_constraint()) ?  cw_target_edge(vmin) : 
						     inf(vmin)->top_edge();
	Edge_handle bot = (vmin->is_constraint()) ?  cw_source_edge(vmin) : 
						     inf(vmin)->bottom_edge();
	// ---------------------------------------------------------------------
	// A new Vertex ( == bitangent) has been found 
	// splice the two edges defining it
	if (vmin != sup(top)) splice(top,vmin);
	if (vmin != sup(bot)) splice(bot,vmin);
	// ---------------------------------------------------------------------
	// The two possible adjacent candidates must be erased
	Face_handle ef[6];
	if (bot->sign()) { ef[0] = dl(bot); ef[1] = ur(bot); ef[2] = ul(bot); }
	else             { ef[0] = dl(bot); ef[1] = dr(bot); ef[2] = ul(bot); }
	if (top->sign()) { ef[3] = dl(top); ef[4] = ur(top); ef[5] = ul(top); }
	else             { ef[3] = dl(top); ef[4] = dr(top); ef[5] = ul(top); }
	for (int i = 0 ; i < 6 ; i++) 
	  if (is_candidate(ef[i], _Tr()) && sup(ef[i]) != vmin) {
		death.erase(sup(ef[i]));
		if (!sup(ef[i])->is_constraint()) set_sup(ef[i],0);
	    }
	// ---------------------------------------------------------------------
	// Update the antichain
	CGAL_precondition(sup(vmin) == 0);
	sweep_constraint(vmin, _Tr());
	set_adjacent_old_faces(top,0,0,0);
	set_adjacent_old_faces(bot,0,0,0);
	// ---------------------------------------------------------------------
	// Computing new candidates by looking at the adjacent faces of
	// bot and top.
	bot = ccw_source_edge(vmin);
	top = ccw_target_edge(vmin);
	if (sup(bot) != 0 && sup(bot)->is_constraint() && 
	    is_minimal(sup(bot),tr)) 
	    death.insert(sup(bot));
	if (sup(top) != 0 && sup(top)->is_constraint() && 
	    is_minimal(sup(top),tr))
	    death.insert(sup(top));
	Face_handle af[] = { dl(bot) , ur(bot) , dl(top) , ur(top) };
	for (int i = 0 ; i < 4 ; i++) 
	  if (is_candidate(af[i], _Tr()) &&
		sup(af[i]) != vmin && sup(af[i]) != vmin->pi()) {
		if (sup(af[i]) == 0) set_sup(af[i],vc(af[i]->bottom_edge(),
						      af[i]->top_edge()));
		death.insert(sup(af[i]));
	    }
	// ---------------------------------------------------------------------
    } // end while (death.size() != 0) 

    // -------------------------------------------------------------------------
    // Recovering the initial Edge --> Face pointers with a
    Edge_iterator ea = a.edges_begin();
    for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e,++ea) {
	if (e->sign()) 
	     set_adjacent_faces(&(*e),dl(&(*ea)),ur(&(*ea)),ul(&(*ea)));
	else set_adjacent_faces(&(*e),dl(&(*ea)),dr(&(*ea)),ul(&(*ea)));
    }
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
void
Visibility_complex_antichain<_Gtr,It,Flip>::compute_minimals(_Tr tr)
{
    // -------------------------------------------------------------------------
    // Operators used in this method.
    typename _Tr::Sup sup; 
    typedef typename _Tr::Right_traits Right_traits;
    typename Right_traits::Target_cusp_edge target_cusp_edge;
    typename Right_traits::Top_edge         top_edge;
    // -------------------------------------------------------------------------
    // Computing the set of minimal bitangents
    Vertex_handle ll = 0, rr = 0;
    for (Face_iterator fi = faces_begin(); fi != faces_end() ; ++fi) {
	Face_handle f = &(*fi);
	if (sup(f) != 0 && is_minimal(sup(f),tr)) {
	    if (!sup(f)->is_constraint() || 
		top_edge(f) == target_cusp_edge(sup(f)))  {
		Vertex_handle v = sup(f);
		if (v != ll && v != rr) push_back_minimal(sup(f),tr);
		if (ll == 0 && v->is_left_left() && 
			       is_on_convex_hull(v))      ll = v;
		else if (rr == 0 && v->is_right_right() && 
				    is_on_convex_hull(v)) rr = v;
	    }
	}
    }
    // -------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
// Find topmost and bottommost edge.
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::
initialize_convex_hull(_Tr /*tr*/) const
{
    typename _Tr::Set_adjacent_faces            set_adjacent_faces;
    typename _Tr::Dr   dr;  typename _Tr::Dl  dl;
    typename _Tr::Ur   ur;  typename _Tr::Ul  ul;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Sup sup;

    Edge_handle e = infinite_face()->bottom_edge();
    do {
	e = ccw_source_edge(sup(e));
	set_adjacent_faces(e,dl(e),dr(e),infinite_face());
    } while (e != infinite_face()->bottom_edge()); 

    e = infinite_face()->top_edge();
    do {
	e = ccw_target_edge(sup(e));
	set_adjacent_faces(e,infinite_face(), ur(e),ul(e));
    } while (e != infinite_face()->top_edge()); 
}

// ----------------------------------------------------------------------------- 
// The antichain has 3n-1 faces whereas the pseudo-triangulation has only 3n-3
// bitangents. There are two faces which do not have their sup pointer set at
// this point. These are ul(bottommost) and dr(topmost). This function sets the
// sink of the two supplementary faces.
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::fix_extreme_edges(_Tr /*tr*/) const // warning tr is never used 
{
    //--------------------------------------------------------------------------
    // Operators used in this method.
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    typename _Tr::Set_sup set_sup; typename _Tr::Set_inf set_inf; 
    typename _Tr::Dr   dr;  typename _Tr::Dl  dl;
    typename _Tr::Ur   ur;  typename _Tr::Ul  ul;
    typename _Tr::CcL  ccL; typename _Tr::CcR ccR;
    typename _Tr::Set_adjacent_faces set_adjacent_faces;
    typename _Tr::Ccw_edge ccw_edge;
    typename _Tr::Splice splice;
    // ------------------------------------------------------------------------- 
    Edge_handle bot = infinite_face()->top_edge();    // Bottommost edge
    Edge_handle top = infinite_face()->bottom_edge(); // Topmost edge
    set_adjacent_faces(bot,infinite_face(),ur(bot),ul(bot));
    set_adjacent_faces(top,dl(top),dr(top),infinite_face());
    // ------------------------------------------------------------------------- 
    // A face may appear twice in the antichain due to our identification 
    // v = v->pi()->pi(). If this is the case, delete one of the two and update
    // pointers. 
    if (ur(bot) == dr(top) && ul(bot) != dl(top)) {
	Edge_handle t = ul(bot)->top_edge();
	set_inf(dl(top),inf(ul(bot)));
	if (t->sign()) set_adjacent_faces(t,dl(top),ur(t),ul(t));
	else           set_adjacent_faces(t,dl(top),dr(t),ul(t));
	dl(top)->set_top_edge(t);
	set_adjacent_faces(bot,dr(bot),ur(bot),dl(top));
	// delete ul(bot);
    }
    else if (ul(bot) == dl(top) && ur(bot) != dr(top)) {
	Edge_handle t = dr(top)->bottom_edge();
	set_inf(ur(bot),inf(dr(top)));
	if (!t->sign()) set_adjacent_faces(t,dl(t),dr(t),ur(bot));
	else            set_adjacent_faces(t,dl(t),ur(bot),ul(t));
	ur(bot)->set_bottom_edge(t);
	set_adjacent_faces(top,dl(top),ur(bot),ul(top));
	// delete dr(top);
    }
    // ------------------------------------------------------------------------- 
    // Set the sink of the two faces on the convex-hull.
    if (sup(ur(bot)) != 0 && sup(ul(bot)) == 0) {
	Vertex_handle v = sup(ur(bot));
	Vertex_handle w = inf(dr(top));
	if (*v == *w->pi()) set_sup(ul(bot),w);
	else {
	    do { v = ccR(v); } while (v->is_constraint());
	    CGAL_precondition(v != 0 && inf(v) != 0);
	    set_sup(ul(bot),inf(inf(v)));
	}
	sup(ur(bot))->set_pi(sup(ul(bot)));
    }
    if (sup(dl(top)) != 0 && sup(dr(top)) == 0) {
	Vertex_handle v = sup(dl(top));
	Vertex_handle w = inf(ul(bot));
	if (*v == *w->pi()) set_sup(dr(top),w);
	else {
	    do { v = ccL(v); } while (v->is_constraint());
	    CGAL_precondition(v != 0 && inf(v) != 0);
	    set_sup(dr(top),inf(inf(v)));
	}
	sup(dl(top))->set_pi(sup(dr(top)));
    }
    // ------------------------------------------------------------------------- 
    // Due to the problem of the face appearing twice in the antichain, it is
    // possible that sup(bot) or sup(top) is 0. We fix this.
    //CGAL_precondition(sup(bot) != 0 || sup(top) != 0);
    if (sup(bot) == 0) {
	CGAL_precondition(sup(top) != 0);
	CGAL_precondition(inf(bot) != 0);
	splice(ccw_edge(inf(bot),bot->object()), sup(top)->pi());
    }
    else if (sup(top) == 0) {
	CGAL_precondition(sup(bot) != 0);
	CGAL_precondition(inf(top) != 0);
	splice(ccw_edge(inf(top),top->object()), sup(bot)->pi());
    }
    // ------------------------------------------------------------------------- 
    CGAL_precondition(sup(bot) != 0);
    CGAL_precondition(sup(top) != 0);
}

// ----------------------------------------------------------------------------- 
// This function glues the last arc of G to the first arc of G_*.
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::glue_ccw_cw(_Tr /*tr*/) const // warning tr is never used
{
    //--------------------------------------------------------------------------
    // Operators used in this method.
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    typename _Tr::Ccw_edge ccw_edge; typename _Tr::Cc   cc;
    typename _Tr::Splice splice;
    // -------------------------------------------------------------------------
    // Let eo be the edge on the same object as e but with opposite sign.
    // Let vmax           be the maximal bitangent leaving the disk of e.
    // Let vmin = sup(eo) be the minimal bitangent leaving the disk of eo.
    // We set the pointer sup(ccw_edge(vmax,e->object())) = pi(vmin)
    for (Edge_const_iterator e = edges_begin(); e != edges_end() ; ++e) {
	if (e->object() != 0) {
	    Edge_const_iterator eo = e; if (e->sign()) --eo; else ++eo;
	    Vertex_handle vmin = inf(&(*eo));
	    Vertex_handle vmax = sup(&(*e));
	    while (cc(vmax,e->object()) != 0 && *vmax != *vmin->pi())
		vmax = cc(vmax,e->object());
	    if (*vmax == *vmin->pi()) {
		if (vmax != vmin->pi()) delete vmin->pi();
	       	vmin->set_pi(vmax);
	    }
	    else splice(ccw_edge(vmax,e->object()),vmin->pi());
	    Edge_handle f = ccw_edge(vmax,e->object());  // Fixed 050701
	    if (f != &(*e)) f->set_adjacent_faces(0,0,0);// Fixed 050701
	    //ccw_edge(vmax,e->object())->set_adjacent_faces(0,0,0);
	    CGAL_precondition(inf(ccw_edge(vmax,e->object())) == vmax);
	    splice(ccw_edge(vmin->pi(),e->object()),sup(&(*eo))->pi());
	}
    }
    // ------------------------------------------------------------------------- 
}
    
// -----------------------------------------------------------------------------
// Identifying vertices of the convex-hull

template < class _Gtr , class It , class Flip >
inline bool 
Visibility_complex_antichain<_Gtr,It,Flip>::
is_on_convex_hull(Vertex_handle v) const
{
    return ((v->is_left_left()   &&  v->cw_source_edge() != 0 && 
	     v->cw_source_edge()->dr() == infinite_face()) ||
	    (v->is_right_right() &&  v->cw_target_edge() != 0 &&
	     v->cw_target_edge()->ul() == infinite_face()));
}

// -----------------------------------------------------------------------------
// Minimality Testing

template < class _Gtr , class It , class Flip >
template < class _Tr > 
inline bool 
Visibility_complex_antichain<_Gtr,It,Flip>::
is_minimal(const Vertex_handle& v, _Tr /*tr*/)  const
{
  typename  _Tr::Left_traits lt;
  typename _Tr::Right_traits rt;
  return (is_xx_minimal(v,lt) && 
	    is_xx_minimal(v,rt));
}

template < class _Gtr , class It , class Flip >
template < class _Tr > 
inline bool 
Visibility_complex_antichain<_Gtr,It,Flip>::
is_xx_minimal(const Vertex_handle& v, _Tr /*tr*/ ) const
{ 
    typename _Tr::Cw_source_edge cw_source_edge;
    typename _Tr::Bottom_edge    bottom_edge;
    typename _Tr::Top_edge       top_edge;
    typename _Tr::Is_left_xx     is_left_xx;
    typename _Tr::Left_cw_traits left_cw_traits;
    typename _Tr::Right_cw_traits right_cw_traits;
    if (v == 0 || cw_source_edge(v) == 0) return false;

    if (is_on_convex_hull(v)) {
	typename _Tr::CcL ccL; typename _Tr::CwR cwR;
	if (is_left_xx(v)) 
	     return cw_source_edge(v) == top_edge(infinite_face());
	else if (ccL(cwR(v)) == v) 
	     return is_xx_minimal(cwR(v), left_cw_traits); 
	else return is_xx_minimal(cwR(v), right_cw_traits);
    }

    if (v->is_constraint() || v->pi()->is_constraint()) {
	typename _Tr::Ul ul; typename _Tr::Dr dr;
	typename _Tr::Ur ur; typename _Tr::Dl dl; 
	if (is_left_xx(v)) {
	    Face_handle f = ur(cw_source_edge(v));
	    if (f == 0) f = ul(cw_source_edge(v));
	    return (f != 0 && bottom_edge(f) == cw_source_edge(v));
	}
	else {
	    Face_handle f = dl(cw_source_edge(v));
	    if (f == 0) f = dr(cw_source_edge(v));
	    return (f != 0 && top_edge(f)    == cw_source_edge(v));
	}
    }

    typename _Tr::Sup sup; typename _Tr::Inf inf; 
    Face_handle f = inf(v);
    return (f != 0 && bottom_edge(f) != 0 && v == sup(bottom_edge(f)));
}

#ifdef DEBUG
template < class _Gtr , class It , class Flip >
template < class _Tr > 
inline void
Visibility_complex_antichain<_Gtr,It,Flip>::
print_left_minimal(Vertex_handle v, _Tr tr) 
{ 
    cerr << "Not implemented" << endl;
}

template < class _Gtr , class It , class Flip >
template < class _Tr > 
inline void
Visibility_complex_antichain<_Gtr,It,Flip>::
print_right_minimal(Vertex_handle v, _Tr tr) 
{ 
    cerr << "Not implemented" << endl;
}
#endif

// -----------------------------------------------------------------------------
// Methods to manage the list of minimals

template < class _Gtr , class It , class Flip >
template < class _Tr > 
void 
Visibility_complex_antichain<_Gtr,It,Flip>::clear_minimals(_Tr tr) {
    Minimals_iterator first = minimals_begin(tr);
    Minimals_iterator last  = minimals_end(tr);
    while (first != last) erase_minimal(first++,tr);
}

template < class _Gtr , class It , class Flip >
void          
Visibility_complex_antichain<_Gtr,It,Flip>::
push_back_minimal(Vertex_handle v,Ccw_traits)
{
    Less_bitangent<_Gtr> lv;
    Minimals_iterator w = (is_straight()) ? std::lower_bound(minimals_begin(),
							minimals_end(), 
							*v,lv) :
					    minimals_end();
    _minimals_ccw.insert(w,*v);
}

template < class _Gtr , class It , class Flip >
void          
Visibility_complex_antichain<_Gtr,It,Flip>::
push_back_minimal(Vertex_handle v,Cw_traits)
{
    Greater_bitangent<_Gtr> lv;
    Minimals_iterator w = (is_straight()) ? std::lower_bound(cw_minimals_begin(),
							cw_minimals_end(), 
							*v,lv) :
					    cw_minimals_end();
    _minimals_cw.insert(w,*v);
}

template < class _Gtr , class It , class Flip >
typename Visibility_complex_antichain<_Gtr,It,Flip>::Vertex_handle
Visibility_complex_antichain<_Gtr,It,Flip>::pop_minimal(bool finite/* = false*/)
{
    if (minimals_begin() == minimals_end()) return 0;
    if (finite == false) {
	Vertex_handle min = &(*minimals_begin());
	return min;
    }
    else {
	Vertex_handle min = pop_minimal(false);
	if (min != 0 && is_on_convex_hull(min)) return min;
	while (minimals_begin() != minimals_end() && 
	       min->sup() != 0 && !is_on_convex_hull(min)) {
	    erase_minimal(min);
	    min = pop_minimal(false);
	}
	if (min != 0 && is_on_convex_hull(min)) return min;
	return (min == 0 || min->sup() != 0) ? 0 : min;
    }
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::sweep(Vertex_handle v , _Tr tr)
{
    //--------------------------------------------------------------------------
    // Operators used in this method.
    typename _Tr::Inf inf; typename _Tr::Sup sup; 
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Ccw_edge  ccw_edge; typename _Tr::Cw_edge   cw_edge; 
    typename _Tr::Cc cc;   
    typename _Tr::CcR ccR; typename _Tr::CcL ccL;
    typename _Tr::CwR cwR; typename _Tr::CwL cwL;
    typename _Tr::Splice splice;
    typename _Tr::Source_object source_object;
    typename _Tr::Target_object target_object;
    typename _Tr::Merge  merge(this);
    //--------------------------------------------------------------------------
    CGAL_precondition(is_minimal(v,tr));
    CGAL_precondition(cwL(v) != 0 && cwR(v) != 0);
    //--------------------------------------------------------------------------
    // Erasing the edges that are about to be swept
    erase(cw_source_edge(v)); 
    erase(cw_target_edge(v));
    // -------------------------------------------------------------------------
    // Removing cwR(v) and cwL(v) from the minimal list corresponding to the
    // opposite orientation. Adding v.
    typename _Tr::Cw_traits cw_traits;
    erase_minimal(v,tr); push_back_minimal(v,cw_traits);
    if (is_minimal(cwL(v),cw_traits)) erase_minimal(cwL(v),cw_traits);
    if (cwR(v) != cwL(v) && is_minimal(cwR(v),cw_traits))
	erase_minimal(cwR(v),cw_traits);
    // -------------------------------------------------------------------------
    // Fixing operator pi on the source object of v
    if ( cc(cwR(v)->pi(),source_object(v)) != 0 &&
	*cc(cwR(v)->pi(),source_object(v)) == *v->pi()) 
	v->set_pi(cc(cwR(v)->pi(),source_object(v)));
    if (cwL(v->pi()) != 0 && *cwL(v->pi()) == *cwR(v)->pi()) 
	cwR(v)->set_pi(cwL(v->pi()));
    // -------------------------------------------------------------------------
    // Arranging Edge <---> Vertex pointer on the source object of v
    if (ccw_edge(cwR(v)->pi(),source_object(v)) == 0) 
	splice(cw_target_edge(v->pi()),cwR(v)->pi());
    splice(ccw_edge(cwR(v)->pi(),source_object(v)),v->pi());
    if (linear_space() && !cwR(v)->pi()->is_constraint() && 
			  inf(cwR(v)->pi()) == 0) // merge if pi(cwR(v)) notin G
			  //!is_on_convex_hull(cwR(v)->pi()))
	merge(cw_edge (cwR(v)->pi(),source_object(v)),
	      ccw_edge(cwR(v)->pi(),source_object(v)));
    // -------------------------------------------------------------------------
    // Fixing operator pi on the target object of v
    if ( cc(cwL(v)->pi(),target_object(v)) != 0 &&
        *cc(cwL(v)->pi(),target_object(v)) == *v->pi()) 
	v->set_pi(cc(cwL(v)->pi(),target_object(v)));
    if (cwR(v->pi()) != 0 && *cwR(v->pi()) == *cwL(v)->pi()) 
	cwL(v)->set_pi(cwR(v->pi()));
    // -------------------------------------------------------------------------
    // Arranging Edge <---> Vertex pointer on the target object of v
    if (ccw_edge(cwL(v)->pi(),target_object(v)) == 0) 
	splice(cw_source_edge(v->pi()),cwL(v)->pi());
    splice(ccw_edge(cwL(v)->pi(),target_object(v)),v->pi());
    if (linear_space() && !cwL(v)->pi()->is_constraint() && 
			  inf(cwL(v)->pi()) == 0) // merge if pi(cwL(v)) notin G
			  //!is_on_convex_hull(cwR(v)->pi()))
	merge(cw_edge (cwL(v)->pi(),target_object(v)),
	      ccw_edge(cwL(v)->pi(),target_object(v)));
    // -------------------------------------------------------------------------
    // Erasing unneeded element to keep the storage linear
    if (linear_space() && !is_on_convex_hull(v) && !v->is_constraint()) 
    { delete sup(v->pi()); delete inf(v->pi()); }
    //--------------------------------------------------------------------------
    // Flipping v and updating the antichain
    if (v->is_constraint() && v->pi()->is_constraint()) {
	// ---------------------------------------------------------------------
	//erase(target_cusp_edge(v));
	//erase(source_cusp_edge(v));
	sweep_constraint(v, _Tr());
	//push_back(*target_cusp_edge(v->pi()));
	//push_back(*source_cusp_edge(v->pi()));
	// ---------------------------------------------------------------------
    }
    else {
	// ---------------------------------------------------------------------
	unset_constraint(v);
	unset_constraint(v->pi());
	compute_phi(v,tr);       // Flip bitangent in pseudo-triangulation
	sweep_regular(v, _Tr());   // Modify the antichain
	// ---------------------------------------------------------------------
    }
    //--------------------------------------------------------------------------
    // Pushing the two new created edges 
    push_back(*ccw_source_edge(v)); 
    push_back(*ccw_target_edge(v));
    //--------------------------------------------------------------------------
    // Updating minimals
    if (is_minimal(ccL(v),tr))                     push_back_minimal(ccL(v),tr);
    if (ccR(v) != ccL(v) && is_minimal(ccR(v),tr)) push_back_minimal(ccR(v),tr);
    //--------------------------------------------------------------------------
    //Removing phis(v) and merging the edges incident to phis
    if (linear_space() && !is_on_convex_hull(v) && 
	inf(v) != 0 && inf(inf(v)) != 0) {
	Vertex_handle phis = inf(inf(v));
	unset_constraint(phis);
	unset_constraint(phis->pi());
	merge(cw_source_edge(phis),ccw_source_edge(phis));
	merge(cw_target_edge(phis),ccw_target_edge(phis));
	delete phis; delete inf(v);
    }
}

template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::sweep_all_minimals(_Tr tr)
{
    std::list<Vertex_handle> mins;
    for (Minimals_iterator m = minimals_begin(tr); m != minimals_end(tr); ++m)
	mins.push_back(&(*m));
    typename std::list<Vertex_handle>::iterator clic = mins.begin();
    for ( ; clic != mins.end() ; ++clic) 
	if (is_minimal(*clic,tr)) sweep(*clic,tr);
}

// -----------------------------------------------------------------------------

/*
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::sweep_good(Vertex_handle v,_Tr tr)
{
    typename _Tr::CcR ccR; typename _Tr::CcL ccL; typename _Tr::Sup sup;
//    typename _Tr::CwR cwR; typename _Tr::CwL cwL; 
    // -------------------------------------------------------------------------
    sweep(v,tr);
    // -------------------------------------------------------------------------
    // If ccR(v) and/or ccL(v) are geometrically equal to v and if ccR(v) and/or
    // ccL(v) are minimal, we sweep them.
    bool source_is_point = _Gtr().is_point(*v->source_object());
    bool target_is_point = _Gtr().is_point(*v->target_object());
    if (!source_is_point && !target_is_point) return;
    CGAL_precondition(ccR(v) != 0 && ccL(v) != 0 && sup(sup(v)) != 0);
    Vertex_handle w[] = { ccR(v) , ccL(v) , sup(sup(v)) };
    for (int i = 0; i < 3 ; i++)
	if (v->source_object() == w[i]->source_object() &&
	    v->target_object() == w[i]->target_object() &&
	    (source_is_point || v->is_left_xx() == w[i]->is_left_xx()) &&
	    (target_is_point || v->is_xx_left() == w[i]->is_xx_left()) &&
	    is_minimal(w[i],tr))
	    sweep(w[i],tr);
    // -------------------------------------------------------------------------
}
*/
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::sweep_good(Vertex_handle v,_Tr tr)
{
    typename _Tr::CcR ccR; typename _Tr::CcL ccL; typename _Tr::Sup sup;
//    typename _Tr::CwR cwR; typename _Tr::CwL cwL; 
    // -------------------------------------------------------------------------
    sweep(v,tr);
    // -------------------------------------------------------------------------
    // If ccR(v) and/or ccL(v) are geometrically equal to v and if ccR(v) and/or
    // ccL(v) are minimal, we sweep them.
    CGAL_precondition(ccR(v) != 0 && ccL(v) != 0 && sup(sup(v)) != 0);
    Vertex_handle w[] = { ccR(v) , ccL(v) , sup(sup(v)) };
    typename _Gtr::Equal_as_segments equal_as_segments;
    for (int i = 0; i < 3 ; i++)
	if (equal_as_segments(*v,*w[i]) && is_minimal(w[i],tr))
	    sweep(w[i],tr);
    // -------------------------------------------------------------------------
}
template < class _Gtr , class It , class Flip >
template < class _Tr >
void
Visibility_complex_antichain<_Gtr,It,Flip>::sweep_good_all_minimals(_Tr tr)
{
    std::list<Vertex_handle> mins;
    for (Minimals_iterator m = minimals_begin(tr); m != minimals_end(tr); ++m)
	mins.push_back(&(*m));
    typename std::list<Vertex_handle>::iterator clic = mins.begin();
    for ( ; clic != mins.end() ; ++clic) 
	if (is_minimal(*clic,tr)) sweep_good(*clic,tr);
}

// ----------------------------------------------------------------------------- 

template < class _Gtr , class It , class Flip >
template < class _Tr >
typename Visibility_complex_antichain<_Gtr,It,Flip>::Vertex_handle
Visibility_complex_antichain<_Gtr,It,Flip>::
compute_phi(Vertex_handle v, _Tr tr) //const
{
    // -------------------------------------------------------------------------
    typename _Tr::Sup sup; 
    typename _Tr::Set_sup set_sup; typename _Tr::Set_inf set_inf;
    typename _Tr::Splice splice;
    typename _Tr::Cw_source_edge   cw_source_edge;
    typename _Tr::Cw_target_edge   cw_target_edge;
    typename _Tr::CcL ccL; typename _Tr::CcR ccR; 
    typename _Tr::Is_left_xx is_left_xx;
    // -------------------------------------------------------------------------
    Vertex_handle phiv; // The new Vertex that we must compute
    // -------------------------------------------------------------------------
    // Vertex has already been swept.
    if (sup(v) != 0) {
	phiv  = sup(sup(v)); 
	if (!is_on_convex_hull(v)) { 
	    splice(cw_source_edge(phiv),phiv); 
	    splice(cw_target_edge(phiv),phiv); 
	}
    }
    // -------------------------------------------------------------------------
    // The pi of the Vertex has already been swept. We use the formula:
    // phi(v)->pi() = phi(v->pi()).
    else if (sup(v->pi()) != 0) {
	if (is_on_convex_hull(v)) {
	    if (is_left_xx(v)  && ccL(v->pi()) != 0)
		ccR(v)->set_pi(ccL(v->pi()));
	    if (!is_left_xx(v) && ccR(v->pi()) != 0) 
		ccL(v)->set_pi(ccR(v->pi()));
	}
	phiv = sup(sup(v->pi()))->pi();
	if (!is_on_convex_hull(v)) { 
	    splice(cw_source_edge(phiv),phiv); 
	    splice(cw_target_edge(phiv),phiv); 
	}
    }
    // -------------------------------------------------------------------------
    // We flip v in the current pseudo-triangulation. 
    // To compute phi(v) we walk on the incident pseudo-triangles. 
    //else phiv = Chi2_strategy(this)(v,tr);
    else phiv = Flip_traits(this)(v,tr);
    // -------------------------------------------------------------------------
    // We splice the arcs left and right if b is not on the convex hull
    //if (!is_on_convex_hull(v)) { splice(right,phiv); splice(left,phiv); }
    // -------------------------------------------------------------------------
    // Creating a new face with source v and sink phi(v)
    if (sup(v) == 0) set_sup(v,new Face);
    set_inf(sup(v),v);
    set_sup(sup(v),phiv);
    // -------------------------------------------------------------------------
    // Creating a new face with source pi(v) and sink pi(phi(v))
    if (sup(v->pi()) == 0)  set_sup(v->pi(),new Face);
    set_inf(sup(v->pi()),v->pi());
    set_sup(sup(v->pi()),phiv->pi());
    // -------------------------------------------------------------------------
    CGAL_precondition(sup(v) != 0 && sup(sup(v)) != 0);
    return phiv;      
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
void
Visibility_complex_antichain<_Gtr,It,Flip>::set_constraint(Vertex_handle v) 
{
    typedef Ccw_traits _Tr;
    typename _Tr::Set_sup set_sup; typename _Tr::Set_inf set_inf;
    typename _Tr::Set_adjacent_faces_one_to_one set_adjacent_faces;
    typename _Tr::Set_target_cusp_edge set_target_cusp_edge;
    typename _Tr::Set_source_cusp_edge set_source_cusp_edge;
    typename _Tr::Target_cusp_edge target_cusp_edge;
    typename _Tr::Source_cusp_edge source_cusp_edge;
    typename _Tr::Is_left_xx is_left_xx;
    typename _Tr::Is_xx_left is_xx_left;
    // -------------------------------------------------------------------------
    typename _Tr::Sup sup; typename _Tr::Inf inf;
    Face_handle fs = sup(v);         Face_handle fi = inf(v);
    Face_handle pifs = sup(v->pi()); Face_handle pifi = inf(v->pi());
    // -------------------------------------------------------------------------
    // New edge for rays emanating from b->target() where b = this
    // with angle \theta such that \theta(b) - \pi <  \theta < \theta(b)
    // New edge for rays emanating from b->source() where b = this
    // with angle \theta such that \theta(b) - \pi <  \theta < \theta(b)
    set_target_cusp_edge(v,new Edge);
    set_source_cusp_edge(v,new Edge); 
    target_cusp_edge(v)->set_sign(false);
    set_sup(target_cusp_edge(v),v);
    set_inf(target_cusp_edge(v),v->pi());
    source_cusp_edge(v)->set_sign(true);
    set_sup(source_cusp_edge(v),v);
    set_inf(source_cusp_edge(v),v->pi());
    // -------------------------------------------------------------------------
    // New face directed from b->target_object() to b, the sink is b
    Face_handle f0 = new Face; 
    set_sup(f0,v);
    set_inf(f0,v->pi());
    f0->set_top_edge(target_cusp_edge(v)); 
    if (is_xx_left(v)) set_adjacent_faces(target_cusp_edge(v),f0,0,0);
    else               set_adjacent_faces(target_cusp_edge(v),0,f0,0);
    // New face directed from b->source_object() to b, the sink is b
    Face_handle f1 = new Face; 
    set_sup(f1,v);
    set_inf(f1,v->pi());
    f1->set_bottom_edge(source_cusp_edge(v)); 
    if (is_left_xx(v)) set_adjacent_faces(source_cusp_edge(v),0,0,f1);
    else               set_adjacent_faces(source_cusp_edge(v),0,f1,0);
    // -------------------------------------------------------------------------
    set_sup(v,fs);         set_inf(v,fi);
    set_sup(v->pi(),pifs); set_inf(v->pi(),pifi);
    // -------------------------------------------------------------------------
    // Do the same with b->pi()
    //if (!v->pi()->is_constraint()) set_constraint(v->pi());
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
void
Visibility_complex_antichain<_Gtr,It,Flip>::remove_constraint(Vertex_handle v) 
{
    if (!v->is_constraint()) return;
    // -------------------------------------------------------------------------
  //  erase(v->target_cusp_edge());
//    erase(v->source_cusp_edge());
    unset_constraint(v);
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
void
Visibility_complex_antichain<_Gtr,It,Flip>::add_constraint(Vertex_handle v) 
{
    if (v->is_constraint()) return;
    // -------------------------------------------------------------------------
    set_constraint(v);
    set_constraint(v->pi());
    //push_back(*v->target_cusp_edge());
    //push_back(*v->source_cusp_edge());
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
void
Visibility_complex_antichain<_Gtr,It,Flip>::unset_constraint(Vertex_handle v) 
{
    // -------------------------------------------------------------------------
    delete v->target_cusp_face(); delete v->source_cusp_face();
    delete v->target_cusp_edge(); delete v->source_cusp_edge();
//    if (v->pi()->is_constraint()) unset_constraint(v->pi());
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
void
Visibility_complex_antichain<_Gtr,It,Flip>::
sweep_regular(const Vertex_handle& v, const _Tr&) const
{
    // -------------------------------------------------------------------------
    // The operators used by this method
    typename _Tr::Sup sup; typename _Tr::Inf inf; 
    typename _Tr::Set_inf set_inf; 
    typename _Tr::Set_adjacent_faces_one_to_one set_adjacent_faces;
    typename _Tr::Dl dl; typename _Tr::Dr dr;
    typename _Tr::Ul ul; typename _Tr::Ur ur;
    typename _Tr::Is_left_xx is_left_xx;
    typename _Tr::Cw_source_edge cw_source_edge;
      typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Cw_target_edge cw_target_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    // -------------------------------------------------------------------------
    // The four edges adjacent to v
    Edge_handle e  = cw_source_edge(v);
    Edge_handle ep = ccw_source_edge(v);
    Edge_handle f  = cw_target_edge(v);
    Edge_handle fp = ccw_target_edge(v);
    // -------------------------------------------------------------------------
    // The six faces adjacent to v.
    Face_handle f0,f1,f2,f3;
    if (e->sign()) { f0 = ul(e); f1 = dl(e); }
    else { 
	f0 = dl(e); if (f0 == 0 || sup(f0) == v) f0 = dl(f); 
	f1 = dr(e); 
    }
    if (f->sign()) { 
	f2 = ur(f); if (f2 == 0 || sup(f2) == v) f2 = ur(e); 
	f3 = ul(f); 
    }
    else           { f2 = dr(f); f3 = ul(f); }
    Face_handle f4 = sup(v); set_inf(f4,v);
    Face_handle f5 = inf(v);
    // -------------------------------------------------------------------------
    // Due to our identification pi^2(v) == v, a face can appear twice in the
    // antichain. As a consequence a face can appear twice in the antichain.
    // This implies that the pointers Edge <---> Face are not necessarily
    // reversible. 
    // Our solution: make the Edge --> Face pointers correct and use the
    // infinite face and the two Edges below.
    Edge_handle botf4 = f4->bottom_edge();
    Edge_handle topf4 = f4->top_edge();
    // -------------------------------------------------------------------------
    // Updating the Edge <--> Face pointers for the two new edges of the
    // antichain.
    if (linear_space() && !is_on_convex_hull(v)) {
	set_adjacent_faces(e,0,0,0);
	set_adjacent_faces(f,0,0,0);
    }
    if (is_left_xx(v)) set_adjacent_faces(ep,f4,f3,f0);
    else               set_adjacent_faces(ep,f0,f4,f3);
    set_adjacent_faces(fp,f1,f2,f4);
    // -------------------------------------------------------------------------
    // Fix the pointers due to the problem mentionned above. This only occurs
    // for faces whose sink are on the convex-hull.
    CGAL_precondition(f4 != 0);
    CGAL_precondition(f5 != 0);
    if (is_on_convex_hull(v)) {
	if (f2 == f4) { f4->set_top_edge(ep); f4->set_bottom_edge(fp); }
	else {
	    Face_handle w = (is_left_xx(v)) ? f2 : f0;
	    if (botf4 != 0 && w->bottom_edge() != botf4) 
		f4->set_bottom_edge(botf4);
	    if (topf4 != 0 && w->top_edge()    != topf4) 
		f4->set_top_edge(topf4);
	}
	if (f5->bottom_edge() != e) {
	    f5->set_top_edge(infinite_face()->bottom_edge());
	    return;
	}
	else if (f5->top_edge() != f) {
	    f5->set_bottom_edge(infinite_face()->top_edge());
	    return;
	}
    }
    // -------------------------------------------------------------------------
    // The face inf(v) is no longer swept.
    f5->set_bottom_edge(0);
    f5->set_top_edge   (0);
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
void
Visibility_complex_antichain<_Gtr,It,Flip>::
sweep_constraint(const Vertex_handle& v, const _Tr&) const
{
    // -------------------------------------------------------------------------
    // The operators used by this method
    typename _Tr::Sup sup; 
    typename _Tr::Set_adjacent_faces set_adjacent_old_faces;
    typename _Tr::Set_adjacent_faces_one_to_one set_adjacent_faces;
    typename _Tr::Dl dl; typename _Tr::Dr dr;
    typename _Tr::Ul ul; typename _Tr::Ur ur;
    typename _Tr::Target_cusp_face target_cusp_face;
    typename _Tr::Source_cusp_face source_cusp_face;
    typename _Tr::Is_left_xx is_left_xx; typename _Tr::Is_xx_left is_xx_left;
    typename _Tr::Cw_source_edge cw_source_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Cw_target_edge cw_target_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    // -------------------------------------------------------------------------
    // The four edges adjacent to v
    Edge_handle e  = cw_source_edge(v);
    Edge_handle ep = ccw_source_edge(v);
    Edge_handle f  = cw_target_edge(v);
    Edge_handle fp = ccw_target_edge(v);
    // -------------------------------------------------------------------------
    // The four regular faces adjacent to v
    Face_handle f0,f1,f2,f3;
    if (e->sign()) { f0 = ul(e);                               f1 = dl(e); }
    else           { f0 = dl(e); if (sup(f0) == v) f0 = dl(f); f1 = dr(e); }
    if (f->sign()) { f2 = ur(f); if (sup(f2) == v) f2 = ur(e); f3 = ul(f); }
    else           { f2 = dr(f);                               f3 = ul(f); }
    // -------------------------------------------------------------------------
    // The four degenerate faces adjacent to v
    Face_handle a  = source_cusp_face(v);
    Face_handle ap = source_cusp_face(v->pi());
    Face_handle b  = target_cusp_face(v);
    Face_handle bp = target_cusp_face(v->pi());
    // -------------------------------------------------------------------------
    Edge_handle fix = 0;
    if (is_on_convex_hull(v)) {
	if (e->sign() && e != f2->bottom_edge())    fix = f2->bottom_edge();
	else if (!e->sign() && f != f0->top_edge()) fix = f0->top_edge();
    }
    // -------------------------------------------------------------------------
    // Updating the Edge <--> Face pointers for the two new edges of the
    // antichain.
    typename Vertex::Type_util type;
    switch (type(is_left_xx(v),is_xx_left(v))) {
	case Vertex::LL:
	    set_adjacent_old_faces(e,f1,f2,f0);
	    set_adjacent_old_faces(f,a,b,f3);
	    set_adjacent_faces(ep,ap,f3,bp);
	    set_adjacent_faces(fp,f1,f2,f0);
	break;
	case Vertex::RR:
	    set_adjacent_old_faces(e,a,f1,b);
	    set_adjacent_old_faces(f,f0,f2,f3);
	    set_adjacent_faces(ep,f0,f2,f3);
	    set_adjacent_faces(fp,f1,ap,bp);
	break;
	case Vertex::LR:
	    set_adjacent_old_faces(e,f1,b,f0);
	    set_adjacent_old_faces(f,a,f2,f3);
	    set_adjacent_faces(ep,f2,f3,bp);
	    set_adjacent_faces(fp,f1,ap,f0);
	break;
	case Vertex::RL:
	    set_adjacent_old_faces(e,a,f1,f2);
	    set_adjacent_old_faces(f,f0,b,f3);
	    set_adjacent_faces(ep,f0,ap,f3);
	    set_adjacent_faces(fp,f1,f2,bp);
	break;
    }
    // -------------------------------------------------------------------------
    if (is_on_convex_hull(v) && fix != 0) {
	if (fix->sign()) set_adjacent_faces(fix,dl(fix),ur(fix),ul(fix));
	else             set_adjacent_faces(fix,dl(fix),dr(fix),ul(fix));
    }
    // -------------------------------------------------------------------------
    if (a != 0) a->set_top_edge(0);
    if (b != 0) b->set_bottom_edge(0);
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
bool
Visibility_complex_antichain<_Gtr,It,Flip>::
is_swept_regular(const Vertex_handle& v, _Tr /*tr*/) const
{
    // -------------------------------------------------------------------------
    CGAL_precondition(!v->is_constraint());
    // -------------------------------------------------------------------------
    // The operators used by this method
    typename _Tr::Sup sup; 
    typename _Tr::Dl dl; typename _Tr::Dr dr;
    typename _Tr::Ul ul; typename _Tr::Ur ur;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Is_left_xx is_left_xx; typename _Tr::Is_xx_left is_xx_left;
    // -------------------------------------------------------------------------
    // A vertex v has been swept iff. the two following conditions are
    // satisfied:
    // (1) sup(v) is not 0
    // (2) the faces adjacent to ccw_source_edge(v) and ccw_target_edge(v) are
    // correct. We check this with the help of the faces adjacent to
    // cw_source_edge(v) and cw_target_edge(v) which are assumed to be correct.
    // -------------------------------------------------------------------------
    // Checking (1).
    Face_handle f4 = sup(v); 
    if (f4 == 0) return false;
    // -------------------------------------------------------------------------
    // Checking (2).
    // -------------------------------------------------------------------------
    // The four edges adjacent to v
    Edge_handle e  = cw_source_edge(v);
    Edge_handle ep = ccw_source_edge(v);
    Edge_handle f  = cw_target_edge(v);
    Edge_handle fp = ccw_target_edge(v);
    if (e == 0 || ep == 0 || f == 0 || fp == 0) return false;
    // -------------------------------------------------------------------------
    // The remaining faces adjacent to v different from inf(v).
    Face_handle f0,f1,f2,f3;
    if (e->sign()) { f0 = ul(e); f1 = dl(e); }
    else { 
	f0 = dl(e); if (f0 == 0 || sup(f0) == v) f0 = dl(f); 
	f1 = dr(e); 
    }
    if (f->sign()) { 
	f2 = ur(f); if (f2 == 0 || sup(f2) == v) f2 = ur(e); 
	f3 = ul(f); 
    }
    else           { f2 = dr(f); f3 = ul(f); }
    // -------------------------------------------------------------------------
    if (is_left_xx(v)) {
	if (f4 != dl(ep) || f3 != ur(ep) || f0 != ul(ep)) return false;
    }
    else {
	if (f0 != dl(ep) || f4 != dr(ep) || f3 != ul(ep)) return false;
    }
    if (is_xx_left(v)) {
	if (f1 != dl(fp) || f2 != ur(fp) || f4 != ul(fp)) return false;
    }
    else {
	if (f1 != dl(fp) || f2 != dr(fp) || f4 != ul(fp)) return false;
    }
    // -------------------------------------------------------------------------
    return true;
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
bool
Visibility_complex_antichain<_Gtr,It,Flip>::
is_swept_constraint(const Vertex_handle& v, _Tr /*tr*/) const
{
    // -------------------------------------------------------------------------
    CGAL_precondition(v->is_constraint());
    // -------------------------------------------------------------------------
    // The operators used by this method
    typename _Tr::Sup sup; 
    typename _Tr::Dl dl; typename _Tr::Dr dr;
    typename _Tr::Ul ul; typename _Tr::Ur ur;
    typename _Tr::Cw_source_edge  cw_source_edge;
    typename _Tr::Cw_target_edge  cw_target_edge;
    typename _Tr::Ccw_source_edge ccw_source_edge;
    typename _Tr::Ccw_target_edge ccw_target_edge;
    typename _Tr::Target_cusp_face target_cusp_face;
    typename _Tr::Source_cusp_face source_cusp_face;
    typename _Tr::Is_left_xx is_left_xx; typename _Tr::Is_xx_left is_xx_left;
    // -------------------------------------------------------------------------
    // See is_swept_regular for the criterion. As sup(v) is not defined we check
    // only point (2).
    // -------------------------------------------------------------------------
    // The four edges adjacent to v
    Edge_handle e  = cw_source_edge(v);
    Edge_handle ep = ccw_source_edge(v);
    Edge_handle f  = cw_target_edge(v);
    Edge_handle fp = ccw_target_edge(v);
    if (e == 0 || ep == 0 || f == 0 || fp == 0) return false;
    // -------------------------------------------------------------------------
    // The four regular faces adjacent to v
    Face_handle f0,f1,f2,f3;
    if (e->sign()) { 
	if (ul(e) == 0 || ur(e) == 0) return false;
	f0 = ul(e); f1 = dl(e); 
    }
    else { 
	if (dl(e) == 0 || dr(e) == 0) return false;
	f0 = dl(e); if (sup(f0) == v) f0 = dl(f); f1 = dr(e); 
    }
    if (f->sign()) { 
	if (ul(f) == 0 || ur(f) == 0) return false;
	f2 = ur(f); if (sup(f2) == v) f2 = ur(e); f3 = ul(f); 
    }
    else { 
	if (dl(f) == 0 || dr(f) == 0) return false;
	f2 = dr(f); f3 = ul(f); 
    }
    // -------------------------------------------------------------------------
    // The four degenerate faces adjacent to v
    Face_handle a  = source_cusp_face(v);
    Face_handle b  = target_cusp_face(v);
    Face_handle ap = source_cusp_face(v->pi());
    Face_handle bp = target_cusp_face(v->pi());
    if (a == 0 || b == 0 || ap == 0 || bp == 0) return false;
    // -------------------------------------------------------------------------
    // Updating the Edge <--> Face pointers for the two new edges of the
    // antichain.
    typename Vertex::Type_util type;
    switch (type(is_left_xx(v),is_xx_left(v))) {
	case Vertex::LL:
	    if (ap != dl(ep) || f3 != ur(ep) || bp != ul(ep)) return false;
	    if (f1 != dl(fp) || f2 != ur(fp) || f0 != ul(fp)) return false;
	break;
	case Vertex::RR:
	    if (f0 != dl(ep) || f2 != dr(ep) || f3 != ul(ep)) return false;
	    if (f1 != dl(fp) || ap != dr(fp) || bp != ul(fp)) return false;
	break;
	case Vertex::LR:
	    if (f2 != dl(ep) || f3 != ur(ep) || bp != ul(ep)) return false;
	    if (f1 != dl(fp) || ap != dr(fp) || f0 != ul(fp)) return false;
	break;
	case Vertex::RL:
	    if (f0 != dl(ep) || ap != dr(ep) || f3 != ul(ep)) return false;
	    if (f1 != dl(fp) || f2 != ur(fp) || bp != ul(fp)) return false;
	break;
    }
    return true;
}

// -----------------------------------------------------------------------------

template < class _Gtr , class It , class Flip >
template < class _Tr>
bool
Visibility_complex_antichain<_Gtr,It,Flip>::
is_swept(const Vertex_handle& v, _Tr tr) const
{
    if (!v->is_constraint()) return is_swept_regular(v,tr);
    return (is_swept_constraint(v,tr) && is_swept_constraint(v->pi(),tr)); 
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
