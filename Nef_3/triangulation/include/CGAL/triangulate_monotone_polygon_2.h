#ifndef TRIANGULATE_MONOTONE_POLYGON_2_H
#define TRIANGULATE_MONOTONE_POLYGON_2_H

#include <CGAL/Indirect_not_less_yx_2.h>
#include <CGAL/Unique_hash_map.h>
#include <deque>
#include <vector>

#undef _DEBUG 
#define _DEBUG 5
#include <CGAL/Nef_2/debug.h>

#ifdef _DEBUG_WINDOW
extern CGAL::Window_stream W;
#endif

// Implementaion of the algorithm from pp 55--58 of "Computational Geometry 
// Algorithms and  Applications" by de Berg, van Kreveld, Overmars, and 
// Schwarzkopf for triangulating a y-monotone polygon.

namespace CGAL {

enum ChainId { LEFTCHAIN, RIGHTCHAIN, EXTREMEPOINT};

#ifdef _DEBUG_WINDOW_DISABLED
template <class PointCirculator, class ChainMap>
inline void color_chains_vertices( PointCirculator polygon, 
				   const ChainMap& C) {
  Color left = PURPLE, right = VIOLET, extreme = RED;
  PointCirculator c(polygon), done(c);
  CGAL_For_all( c, done)
    W << (C[c] == LEFTCHAIN ? left: C[c] == RIGHTCHAIN ? right: extreme) 
      << **c;
}
#endif

template <class PointIterator, class Traits>
bool is_vertex_visible( PointIterator uk, PointIterator uj, PointIterator ukp,
			ChainId chain, const Traits& traits) {
  // check if the vertex uk is visible from uj, in a ccw oriented polygon.
  // ukp is the adjacent vertex to uk, nearest to uj on the same chain
  CGAL_assertion( chain == LEFTCHAIN || chain == RIGHTCHAIN);
  Orientation proper = (chain == LEFTCHAIN ? RIGHT_TURN: LEFT_TURN);
  return( traits.orientation_2_object()( *uk, *uj, *ukp) == proper);
}


template <class InputCirculator, class OutputIterator, class Traits>
void triangulate_monotone_polygon_2( InputCirculator polygon,
				     OutputIterator diagonals,
				     const Traits& traits) {

  typedef typename Traits::Diagonal Diagonal;
  typedef std::vector<InputCirculator> Circulator_vector;
  typedef typename Circulator_vector::iterator Circulator_iterator;

#ifdef _DEBUG_WINDOW
  typedef typename InputCirculator::value_type Point_2;
  typedef typename Point_2::R K;
  typedef typename K::Segment_2 Segment_2;
  Point_2 pause;
#endif

  // put the input iterators in a vector, and sort it lexicographically
  // decreasing on the yx coordinates

  Circulator_vector U; // TODO: especify vector size
  InputCirculator c(polygon), done(c);
  CGAL_NEF_TRACEN("y-monotone polygon:");
  CGAL_For_all( c, done) {
    CGAL_NEF_TRACEN(*c);
    U.push_back(c);
  }
  std::sort( U.begin(), U.end(), Indirect_not_less_yx_2<Traits>(traits));

  // identify top and bottom vertices

  InputCirculator top = *(U.begin()), bottom = *(--U.end()); // TO VERIFY

  // map every ignput iterator to a chain identifier (left, right or extreme) 

  Unique_hash_map< InputCirculator, ChainId> Chain(RIGHTCHAIN);
  for( c = top; c != bottom; ++c)
    Chain[c] = LEFTCHAIN;
  Chain[top] = Chain[bottom] = EXTREMEPOINT;

  // create the diagonals that triangulate the monotone polygon

  std::deque<Circulator_iterator> S;
  Circulator_iterator uj(U.begin()), ujp, un(--U.end());
  S.push_front(uj); // push u1 on S
  uj++;
  S.push_front(uj); // push u2 on S
  ujp = uj;
  uj++;             // now ujp = u2, uj = u3

  while( uj != un) { // for each uj, j=3..n-1

#ifdef _DEBUG_WINDOW_DISABLED
    color_chains_vertices( polygon, Chain);
    W << YELLOW << **uj << ORANGE << **S.front(); 
    W >> pause;
#endif

    CGAL_assertion( !S.empty());
    if( Chain[*uj] != Chain[*S.front()]) {
      while( S.size() > 1) { 

	// Connect uj to all vertices on stack but the bottom,
	// which is already conneted to uj by an edge of P

	Circulator_iterator uk = S.front();
	S.pop_front();
	CGAL_NEF_TRACEN( "Diagonal { " << **uk << ", " << **uj << " } (diff)");
#ifdef _DEBUG_WINDOW
	W << RED << Segment_2( **uk, **uj); 
	W >> pause;
#endif
	*diagonals++ = Diagonal( *uk, *uj);
      }
      S.pop_front();
      CGAL_assertion( S.empty());

      // uj and ujp bound the untriangulated part of P and so must be kept
      // on the stack

      S.push_front(ujp);
      S.push_front(uj);
      CGAL_assertion( !traits.less_yx_2_object()( **ujp, **uj));
    }
    else {

      // Connect uk with all the visible vertices on the stack, but the one
      // on the top since they are already connected by an edge of P.
      // Visible vertices from uk, if any, are placed consecutively
      // on the top of the stack

      Circulator_iterator ukp = S.front();
      S.pop_front(); 
      CGAL_assertion( !S.empty());
      Circulator_iterator uk = S.front();

#ifdef _DEBUG_WINDOW_DISABLED
      W << RED << **uk << GREEN << **uj << BLUE << **ukp;
      W >> pause;
      color_chains_vertices( polygon, Chain);
#endif
      while( !S.empty() && 
	     is_vertex_visible( *uk, *uj, *ukp, Chain[*uj], traits)) {
	S.pop_front();
	CGAL_NEF_TRACEN( "Diagonal { " << **uk << ", " << **uj << " } (same)");
#ifdef _DEBUG_WINDOW
	W << RED << Segment_2( **uk, **uj); 
	W >> pause;
#endif
#ifdef _DEBUG_WINDOW_DISABLED
	if( !S.empty()) {
	  W << RED << **uk << GREEN << **uj << BLUE << **ukp;
	  W >> pause;
	  color_chains_vertices( polygon, Chain);
	}
#endif
	*diagonals++ = Diagonal( *uk, *uj);
	ukp = uk;
	uk = S.front();
      }

      // uk and ukp now bound the untriangulated part of P and so, must
      // be kept on the stack

      S.push_front(ukp);
      S.push_front(uj);
      CGAL_assertion( !traits.less_yx_2_object()( **ukp, **uj));
    }
    ujp = uj;
    uj++;
  }

  // un is already connected to the bottom and top vertices on the stack
  // by edges of P

  CGAL_assertion( S.size() >= 2);

  // Connect un to all vertices on the stack but the first and last one

  S.pop_front();
  while( S.size() > 1) {
    Circulator_iterator uk = S.front();
    S.pop_front();
    CGAL_NEF_TRACEN( "Diagonal { " << **uk << ", " << **un << " } (final)");
#ifdef _DEBUG_WINDOW
    W << RED << Segment_2( **uk, **un); 
    W >> pause;
#endif
    *diagonals++ = Diagonal( *uk, *uj);
  }
  S.pop_front();
  CGAL_assertion( S.empty());

  return;
}

}

#endif // TRIANGULATE_MONOTONE_POLYGON_2_H
