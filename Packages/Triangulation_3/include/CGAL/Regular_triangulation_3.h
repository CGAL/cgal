#ifndef CGAL_REGULAR_TRIANGULATION_3_H
#define CGAL_REGULAR_TRIANGULATION_3_H

#include <set.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Oriented_side_power_test_3.h>

#define DEBUG_INSERT 0
#define DEBUG_STAR   0
#define DEBUG_SIDE_OF_SPHERE   0

template < class Gt, class Tds>
class CGAL_Regular_triangulation_3 : public CGAL_Triangulation_3<Gt,Tds>
{
public:
  
  typedef typename Gt::Weighted_point Weighted_point;

  typedef typename CGAL_Triangulation_3<Gt,Tds>::Cell_handle Cell_handle;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Vertex_handle Vertex_handle;

  typedef typename CGAL_Triangulation_3<Gt,Tds>::Cell Cell;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Vertex Vertex;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Facet Facet;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Edge Edge;

  typedef typename CGAL_Triangulation_3<Gt,Tds>::Cell_circulator Cell_circulator;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Cell_iterator Cell_iterator;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Facet_iterator Facet_iterator;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Edge_iterator Edge_iterator;
  typedef typename CGAL_Triangulation_3<Gt,Tds>::Vertex_iterator Vertex_iterator;

  typedef typename CGAL_Triangulation_3<Gt,Tds>::Locate_type Locate_type;

  CGAL_Regular_triangulation_3()
    : CGAL_Triangulation_3<Gt,Tds>() {}
  
  CGAL_Regular_triangulation_3(const Gt & gt)
  : CGAL_Triangulation_3<Gt,Tds>(gt) {}
  
  // copy constructor duplicates vertices and cells
  CGAL_Regular_triangulation_3(const CGAL_Regular_triangulation_3<Gt,Tds> 
			       & rt) : CGAL_Triangulation_3<Gt,Tds>(rt)
    { 
      CGAL_triangulation_postcondition( is_valid(true) );  
    }
  
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#else
#if defined(LIST_H) || defined(__SGI_STL_LIST_H)
  int
  insert(list<Point>::const_iterator first,
         list<Point>::const_iterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // LIST_H
#if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
  int
  insert(vector<Point>::const_iterator first,
         vector<Point>::const_iterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // VECTOR_H
#ifdef ITERATOR_H
  int
  insert(istream_iterator<Point, ptrdiff_t> first,
         istream_iterator<Point, ptrdiff_t> last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // ITERATOR_H
  
  int insert(Point* first,
	     Point* last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif 

// CGAL_TEMPLATE_MEMBER_FUNCTIONS

  /// +++  Insertion  +++ //

  // -- 1 -- //
  Vertex_handle insert( const Weighted_point &p ) 
    {
      Cell_handle start;
      if ( dimension() >= 1 ) {
	// there is at least one finite "cell" (or facet or edge)
	start = infinite_vertex()->cell()
	  ->neighbor(infinite_vertex()->cell()->index(infinite_vertex()));
      }
      else {
	start = NULL;
      }
      return insert( p, start );
    }  

  // -- 2 -- //
  Vertex_handle insert( const Weighted_point & p, Cell_handle start ) {

    Vertex_handle v;
    
    if (DEBUG_INSERT){
      cout << "dimension dela triangulation avant d'inserer le nouveau point : ";
    }
    
    switch (dimension()) {
    case 3: 
      {
	if (DEBUG_INSERT){
	  cout << "3" << endl;
	}
	
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate( p, start, lt, li, lj);
	
	switch (lt) {
	case OUTSIDE_CONVEX_HULL:
	case CELL:
	case FACET:
	case EDGE:
	  {
	    if (DEBUG_INSERT){
	      cout << "j'essaye d'inserer un nouveau point en dimension 3\n";
	      cout << "-- dim = 3 -- " << p << endl;
	    }
	    v = new Vertex(p);
	    if (DEBUG_INSERT){
	      cout << "vertex : " << v->point() << endl;
	    }
	    // a voir : comment compter les sommets redondants ? (***)
	    // set_number_of_vertices(number_of_vertices()+1);
	    set<void*, less<void*> > conflicts;
	    set<void*, less<void*> > deleted_points;
	    Cell_handle aconflict;
	    int ineighbor;
	    if (in_conflict_3(p, c)) {
	      find_conflicts_3(conflicts, p, c, aconflict, ineighbor);
	      deleted_points = 
		star_region_delete_points(conflicts,&(*v),&(*aconflict),
					  ineighbor);
	      // debug
	      if (DEBUG_INSERT){
		set<void*, less<void*> >::const_iterator it;
		cout << "##### Points supprimes #####\n";
		for( it = deleted_points.begin(); it != deleted_points.end(); ++it) {
		  Point *p_tmp = (Point *) *it;
		  cout << p_tmp->point() << endl;
		}
		cout << "############################\n";
	      }

	      // a voir : comment compter les sommets redondants ? (***)
	      set_number_of_vertices(number_of_vertices()+1);
	    }
	    // else : traiter le cas des points redondants a stocker dans
	    // la face associee pour le cas d'une future suppression
	    return v;
	  }	 
	case VERTEX:
	  return c->vertex(li);
	default :
	  CGAL_triangulation_assertion(false);  // impossible
	}
	break;
      }
    case 2:
      {
	if (DEBUG_INSERT){
	  cout << "2" << endl;
	}
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate( p, start, lt, li, lj);
	switch (lt) {
	case OUTSIDE_CONVEX_HULL:
	case CELL:
	case FACET:
	case EDGE:
	  {
	    // 2D : particular case
	    // not yet implemented
	  }
	case VERTEX:
	  return c->vertex(li);
	case OUTSIDE_AFFINE_HULL:
	  {
	    // if the 2d triangulation is Regular, the 3d
	    // triangulation will be Regular 
	    if (DEBUG_INSERT){
	      cout << "j'insere un point en dehors de l'enveloppe affine\n";
	      cout << "-- dim = 2 -- " << p << endl;
	    }
	    return
	      CGAL_Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p); 
	  }
	}
      }
    default : // dim <= 1
      if (DEBUG_INSERT){
	if (dimension() == 1) cout << "1" << endl;
	else cout << "0" << endl;
	cout << "j'insere un point dans la triangulation de base\n";
	cout << "-- dim <= 1 -- " << p << endl;
      }
      return CGAL_Triangulation_3<Gt,Tds>::insert(p);
    }
    return CGAL_Triangulation_3<Gt,Tds>::insert(p); // to avoid warning with egcs
  }
  

  set<void*, less<void*> > 
  star_region_delete_points( set<void*, less<void*> > & region, Vertex* v,
			     Cell* c, int li)
    // region is a set of connected cells
    // c belongs to region and has facet i on the boundary of region 
    // replaces the cells in region  
    // by linking v to the boundary of region
    // deleted weighted points that are not in the triangulation
    // anymore, pts will be the list of deleted points
    {
      set<void*, less<void*> > vert;
      Cell *c_tmp;
      Vertex *v_tmp;
      Vertex_handle vh;
      set<void*, less<void*> > pts;
      Point *p;
      int i=0;
 
      if (DEBUG_STAR)
	cout << "@ star_region_delete_points\n";
      // for each cell to be deleted, keep vertices
      set<void*, less<void*> >::const_iterator it;
      for( it = region.begin(); it != region.end(); ++it) {
	c_tmp = (Cell *) *it;
	if (DEBUG_STAR) {
	  cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
	  cout << "Region a detruire : cellule " << i++ << endl;
	}
	// store vertices
	for (int i=0; i<4 ; i++){
	  vh = c_tmp->vertex(i);
	  if (DEBUG_STAR)
	    cout << *vh << endl;
	  if ( (vert.find((void*) &(*(vh))) == vert.end()) &&
	       (! is_infinite(v_tmp)) ) {
	    if (DEBUG_STAR)
	      cout << "stocke le sommet : " << *vh << endl;
	    vert.insert( (void*) &( *(vh) ) );
	  }
	}
      }
      if (DEBUG_STAR)
	cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";


      if (DEBUG_STAR)
	cout << "j'appelle star_region\n";

      // Create the new faces and delete old ones
      _tds.star_region( region, v, c, li );
      
      // get the vertices incident to v
      set<Vertex*, less<Vertex*> > inc_vert;
      incident_vertices(v, inc_vert);

      // for each vertex, check if it is a vertex incident to v
      // if not, delete it
      if (DEBUG_STAR)
	cout << "je passe en revue les sommets\n";
      for( it = vert.begin(); it != vert.end(); ++it) {
	v_tmp = (Vertex *) *it;
	if (DEBUG_STAR)
	  cout << *v_tmp << endl;
	if ( (inc_vert.find( v_tmp )) == inc_vert.end() ) {
	  // vertex has to be deleted and point to be stored
	  p = new Point( v_tmp->point() );
	  if (DEBUG_STAR) {
	    cout << "-- Suppression du point -- " << *p << endl;
	  }
	  pts.insert( p );
	  set_number_of_vertices(number_of_vertices()-1);
	  delete(*it);
	}
      }
      
      // returns list of deleted points
      return pts;
    }
    


private:
  bool in_conflict_3(const Weighted_point & p, Cell_handle c)
    {
      return ( side_of_sphere( c, p ) ==  CGAL_ON_BOUNDED_SIDE );
    }

  void find_conflicts_3(set<void*, less<void*> > &conflicts, 
			const Weighted_point & p,
			Cell_handle c, Cell_handle & ac, int & i)
    {
      if ( ( conflicts.find( (void *) &(*c) ) ) != conflicts.end() )
	{
	  return;   // c was already found
	}
      
      (void) conflicts.insert( (void *) &(*c) );
      
      for ( int j=0; j<4; j++ ) {
	if ( in_conflict_3( p, c->neighbor(j) ) ) {
	  find_conflicts_3(conflicts, p, c->neighbor(j), ac, i);
	}
	else {
	  ac = c;
	  i = j;
	}
      }      
    }
 
  
  void find_conflicts_2(set<void*, less<void*> > &conflicts, 
			const Weighted_point & p,
			Cell_handle c, Cell_handle & ac, int & i)
    {
      // find_conflicts_2 : not yet implemented
      ;
    }

public:

  CGAL_Bounded_side
  side_of_sphere( Cell_handle c, const Weighted_point &p) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    int i3;
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {  
      // Cas d'une face de dimension finie
      CGAL_Oriented_side o = geom_traits().power_test_3(c->vertex(0)->point(),
					       c->vertex(1)->point(),
					       c->vertex(2)->point(),
					       c->vertex(3)->point(),p);
      return ( (o == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
	       (o == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
	       CGAL_ON_BOUNDARY );
    }
    else {
      // cas d'une face infinie
      int i0,i1,i2;
      if ( (i3%2) == 1 ) {
	i0 = (i3+1)&3;
	i1 = (i3+2)&3;
	i2 = (i3+3)&3;
      }
      else {
	i0 = (i3+2)&3;
	i1 = (i3+1)&3;
	i2 = (i3+3)&3;
      }

      bool coll = false;
      CGAL_Oriented_side side;
      // on verifie que p n'est pas collineaire a 2 autres points
      // dans ce cas, on traite un cas en dimension 1
      if (geom_traits().collinear(c->vertex(i0)->point(),
				  c->vertex(i1)->point(), p))
	{
	  if (DEBUG_SIDE_OF_SPHERE) {
	    cout << "le point est collineaire a i0(";
	    cout << c->vertex(i0)->point() << ") et i1(";
	    cout << c->vertex(i1)->point() << ")\n";
	  }
	  coll = true;
	  side = geom_traits().power_test_3 
	    ( c->vertex(i0)->point(), c->vertex(i1)->point(), p );
	}
      else if (geom_traits().collinear(c->vertex(i0)->point(),
				  c->vertex(i2)->point(), p))
	{
	  if (DEBUG_SIDE_OF_SPHERE) {
	    cout << "le point est collineaire a i0(";
	    cout << c->vertex(i0)->point() << ") et i2(";
	    cout << c->vertex(i2)->point() << ")\n";
	  }
	  coll = true;
	  side = geom_traits().power_test_3 
	    ( c->vertex(i0)->point(), c->vertex(i2)->point(), p );
	}
      else if (geom_traits().collinear(c->vertex(i1)->point(),
				       c->vertex(i2)->point(), p))
	{
	  if (DEBUG_SIDE_OF_SPHERE) {
	    cout << "le point est collineaire a i1(";
	    cout << c->vertex(i1)->point() << ") et i2(";
	    cout << c->vertex(i2)->point() << ")\n";
	  }
	  coll = true;
	  side = geom_traits().power_test_3 
	    ( c->vertex(i1)->point(), c->vertex(i2)->point(), p );
	}
	 
      if (coll) 
	{
	  return ( (side == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
		   (side == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
		   CGAL_ON_BOUNDARY );	   
	}

      // teste de quel cote du plan defini par (i0, i1, i2) p repose
      CGAL_Orientation
	o = geom_traits().orientation(c->vertex(i0)->point(),
				      c->vertex(i1)->point(),
				      c->vertex(i2)->point(),
				      p);
      
      switch (o) {
      case CGAL_POSITIVE: 
	{
	  return CGAL_ON_BOUNDED_SIDE;
	}
      case CGAL_NEGATIVE:
	{
	  return CGAL_ON_UNBOUNDED_SIDE;
	}
      case CGAL_ZERO:
	{
	  // i0, i1, i2, p - coplanar
	  //cout << "les points sont coplanaires\n";
	  CGAL_Oriented_side s = 
	    geom_traits().power_test_3
	    ( c->vertex(i0)->point(),
	      c->vertex(i1)->point(),
	      c->vertex(i2)->point(),
	      p );
	  
	  return ( (s == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
		   (s == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
		   CGAL_ON_BOUNDARY );
	}
      } 
    }
    return CGAL_ON_UNBOUNDED_SIDE;// to avoid warning with egcs
  }// end side of sphere




  CGAL_Bounded_side side_of_circle( const Facet &f, 
				    const Weighted_point &p) const;
  CGAL_Bounded_side side_of_circle( Cell_handle c, int i, 
				    const Weighted_point &p) const;

  bool is_valid(bool verbose = false, int level = 0) const
    {
      if ( ! tds().is_valid(verbose,level) ) {
	if (verbose) { cerr << "invalid data structure" << endl; }
	CGAL_triangulation_assertion(false); return false;
      }
      if ( &(*infinite_vertex()) == NULL ) {
	if (verbose) { cerr << "no infinite vertex" << endl; }
	CGAL_triangulation_assertion(false); return false;
      }

      int i;

      switch ( dimension() ) {
      case 3:
	{
	  Cell_iterator it;
	  for ( it = finite_cells_begin(); it != cells_end(); ++it ) {
	    is_valid_finite((*it).handle());
	    for ( i=0; i<4; i++ ) {
	      if ( side_of_sphere
		   ( (*it).handle(), 
		     it->vertex( (it->neighbor(i))->index((*it).handle() ) )
		     ->point() )
		   == CGAL_ON_BOUNDED_SIDE ) {
		if (verbose) { 
		  cerr << "non-empty sphere " << endl;
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	  }
	  break;
	}
	//      case 2:
	//	{
	//	  Facet_iterator it;
	//	  for ( it = finite_facets_begin(); it != facets_end(); ++it ) {
	//	    is_valid_finite((*it).first);
	//	    for ( i=0; i<2; i++ ) {
	//	      if ( side_of_circle
	//		   ( (*it).first, 3,
	//		     (*it).first
	//		     ->vertex( (((*it).first)->neighbor(i))->index((*it).first) )->point() )
	//		   == CGAL_ON_BOUNDED_SIDE ) {
	//		if (verbose) { 
	//		  cerr << "non-empty circle " << endl;
	//		}
	//		CGAL_triangulation_assertion(false); return false;
	//	      }
	//	    }
	//	  }
	//	  break;
	//	}
      case 1:
	{
	  Edge_iterator it;
	  for ( it = finite_edges_begin(); it != edges_end(); ++it ) {
	    is_valid_finite((*it).first);
	  }
	  break;
	}
      }
      if (verbose) { cerr << "Regular valid triangulation" << endl;}
      return true;
    }
};

///////////////////////////  IMPLEMENTATION //////////////////////////

// +++  Insertion  +++ //

// +++  Determination de la region en conflit  +++ //

#endif CGAL_REGULAR_TRIANGULATION_3_H
