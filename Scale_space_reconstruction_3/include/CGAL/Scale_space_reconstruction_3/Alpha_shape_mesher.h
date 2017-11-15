// Copyright (C) 2013 INRIA - Sophia Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s):      Thijs van Lankveld, Simon Giraudot


#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_ALPHA_SHAPE_MESHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_ALPHA_SHAPE_MESHER_H

#include <CGAL/license/Scale_space_reconstruction_3.h>

#include <CGAL/Scale_space_reconstruction_3/Shape_construction_3.h>

#include <CGAL/Union_find.h>

namespace CGAL
{

namespace Scale_space_reconstruction_3
{
  
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *
 *  Surface mesher for scale space reconstruction based on
 *  `CGAL::Alpha_shape_3`.
 * 
 *  The surface can be constructed either for a fixed neighborhood
 *  radius, or for a dynamic radius. When constructing the surface for
 *  exactly one neighborhood radius, it is faster to set
 *  `FixedSurface` to `Tag_true`. If the correct neighborhood radius
 *  should be changed or estimated multiple times, it is faster to set
 *  `FixedSurface` to `Tag_false`.
 *
 *  It is undefined whether a surface with fixed radius may have its radius
 *  changed, but if so, this will likely require more computation time than
 *  changing the radius of a dynamic surface. In either case, it is possible to
 *  change the point set while maintaining the same radius.
 *
 *  The surface can be stored either as an unordered collection of triangles, 
 *  or as a collection ordered by shells. A shell is a maximally connected
 *  component of the surface where connected facets are locally oriented
 *  towards the same side of the surface.
 *
 *  \cgalModels CGAL::Scale_space_reconstruction_3::Mesher
 *
 *  \tparam Geom_traits is the geometric traits class. It must be a
 *  model of `DelaunayTriangulationTraits_3`. It must have a
 *  `RealEmbeddable` field number type. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.
 *  \tparam FixedSurface determines whether the surface is expected to
 *  be constructed for a fixed neighborhood radius. It must be a
 *  `Boolean_tag` type. The default value is `Tag_true`. Note that the
 *  value of this parameter does not change the result but only has an
 *  impact on the run-time.
 */
template <typename Geom_traits, typename FixedSurface = Tag_true>
class Alpha_shape_mesher
{
public:
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3                        Point;          ///< defines the point type.
  
  typedef CGAL::cpp11::array< unsigned int, 3 >       Facet;                 ///< defines a triple of point indices indicating a triangle of the surface.
private:
  typedef std::list< Facet >                         Facetset;              ///< defines a collection of triples.
  // Note that this is a list for two reasons: iterator validity for the shell iterators, and memory requirements for the expected huge collections.

public:
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                            Facet_iterator;        ///< defines an iterator over the triples.
  typedef const unspecified_type                      Facet_const_iterator;  ///< defines a constant iterator over the triples.
#else // DOXYGEN_RUNNING
  typedef Facetset::iterator                         Facet_iterator;
  typedef Facetset::const_iterator                   Facet_const_iterator;
#endif // DOXYGEN_RUNNING

private:
  typedef std::vector< Facet_iterator >              FacetIterSet;

private:
  
  // Constructing the surface.
  typedef CGAL::Shape_construction_3<Geom_traits, FixedSurface> Shape_construction_3;

  typedef typename Shape_construction_3::Shape           Shape;
  typedef typename Shape_construction_3::Triangulation   Triangulation;

  typedef typename Shape::Vertex_handle               Vertex_handle;
  typedef typename Shape::Cell_handle                 Cell_handle;
  typedef typename Shape::Facet                       SFacet;
  typedef typename Shape::Edge                        Edge;
  typedef std::pair<Vertex_handle, Vertex_handle>     VEdge;
    
  typedef typename Shape::Vertex_iterator             Vertex_iterator;
  typedef typename Shape::Cell_iterator               Cell_iterator;
  typedef typename Shape::Facet_iterator              SFacet_iterator;
  typedef typename Shape::Edge_iterator               Edge_iterator;

  typedef typename Shape::Finite_cells_iterator       Finite_cells_iterator;
  typedef typename Shape::Finite_facets_iterator      Finite_facets_iterator;
  typedef typename Shape::Finite_edges_iterator       Finite_edges_iterator;
  typedef typename Shape::Finite_vertices_iterator    Finite_vertices_iterator;

  typedef typename Shape::Facet_circulator            SFacet_circulator;
  
  typedef typename Shape::All_cells_iterator          All_cells_iterator;

  typedef typename Shape::Classification_type         Classification_type;
  
  typedef std::map<SFacet, unsigned int> Map_facet_to_shell;
  typedef typename CGAL::cpp11::array<std::set<SFacet>, 2 >   Bubble;

  bool _separate_shells;
  bool _force_manifold;
  FT _border_angle;

  // The shape must be a pointer, because the alpha of a Fixed_alpha_shape_3
  // can only be set at construction and its assignment operator is private.
  // We want to be able to set the alpha after constructing the scale-space
  // reconstructer object.
  Shape*          _shape;

  // The surface. If the surface is collected per shell, the triples of the
  // same shell are stored consecutively.
  Facetset       _surface;

  // The shells can be accessed through iterators to the surface.
  FacetIterSet   _shells;

  // If the surface is forced to be manifold, removed facets are stored
  Facetset _garbage;

  // Map TDS facets to shells
  Map_facet_to_shell _map_f2s;
  unsigned int _index;

  std::vector<Bubble> _bubbles;
  std::map<SFacet, std::size_t> _map_f2b;

  FT _squared_radius;
  
public:

  /**
   *  Constructs an alpha shape mesher.
   *
   *  \param squared_radius \f$\alpha\f$ parameter of the alpha shape algorithm.
   *  \param separate_shells determines whether to collect the surface per shell. 
   *  \param force_manifold determines if the surface is forced to be 2-manifold.
   *  \param border_angle sets the maximal angle between two facets
   *  such that the edge is seen as a border.
   *
   *  If the output is forced to be 2-manifold, some almost flat
   *  volume bubbles are detected. To do so, border edges must be
   *  estimated.
   *
   *  An edge adjacent to 2 regular facets is considered as a border
   *  if it is also adjacent to a singular facet or if the angle
   *  between the two regular facets is lower than this parameter
   *  (set to 45° by default).
   *
   *  \note `border_angle` is not used if `force_manifold` is set to false.
   */
  Alpha_shape_mesher (FT squared_radius,
                      bool separate_shells = false,
                      bool force_manifold = false,
                      FT border_angle = 45.)
    : _separate_shells (separate_shells),
      _force_manifold (force_manifold),
      _border_angle (border_angle),
      _shape (NULL),
      _squared_radius (squared_radius)
  {

  }

  /// \cond SKIP_IN_MANUAL
  template <typename InputIterator, typename OutputIterator>
  void operator() (InputIterator begin, InputIterator end, OutputIterator output)
  {
    clear_surface();

    _shape = Shape_construction_3().construct (begin, end, _squared_radius);
    
    // If shells are not separated and no manifold constraint is given,
    // then the quick collection of facets can be applied
    if (!_separate_shells && !_force_manifold)
    {
      collect_facets_quick ();
      _shells.push_back (_surface.begin ());
    }
    else
    {
      // Init shell index
      _index = 0;

      // Collect all surface meshes from the alpha-shape in a fashion similar to ball-pivoting.
      // Reset the facet handled markers.
      for( All_cells_iterator cit = _shape->all_cells_begin(); cit != _shape->all_cells_end(); ++cit )
        cit->info() = 0;

      if (_force_manifold)
      {
        // If manifold surface is wanted, small volumes (bubbles) are first detected in
        // order to know which layer to ignore when collecting facets
        detect_bubbles();
      }

      collect_facets ();

      if (_force_manifold)
      {
        // Even when taking into account facets, some nonmanifold features might remain
        fix_nonmanifold_edges();
        fix_nonmanifold_vertices();
      }
    }

    for (Facet_iterator it = _surface.begin(); it != _surface.end(); ++ it)
    {
      cpp11::array<std::size_t, 3> f = {{ std::size_t((*it)[0]), std::size_t((*it)[1]), std::size_t((*it)[2]) }};
      *(output ++) = f;
    }
  }
  /// \endcond

  /// gives the number of triangles of the surface.
  std::size_t number_of_triangles() const { return _surface.size(); }

  /// gives an iterator to the first triple in the surface.
  Facet_const_iterator surface_begin() const { return _surface.begin(); }
  /// gives an iterator to the first triple in the surface.
  /** \warning Changes to the surface may change its topology.
   */
  Facet_iterator surface_begin() { return _surface.begin(); }
    
  /// gives a past-the-end iterator of the triples in the surface.
  Facet_const_iterator surface_end() const { return _surface.end(); }
  /// gives a past-the-end iterator of the triples in the surface.
  /** \warning Changes to the surface may change its topology.
   */
  Facet_iterator surface_end() { return _surface.end(); }


  /// gives the number of shells of the surface.
  std::size_t number_of_shells() const {
    return _shells.size();
  }

  /// gives an iterator to the first triple in a given shell.
  /** \param shell is the index of the shell to access.
   *
   *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
   */
  Facet_const_iterator shell_begin( std::size_t shell ) const
  {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
  }
  /// gives an iterator to the first triple in a given shell.
  /** \param shell is the index of the shell to access.
   *
   *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
   *
   *  \warning Changes to a shell may invalidate the topology of the surface.
   */
  Facet_iterator shell_begin( std::size_t shell )
  {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
  }

  /// gives a past-the-end iterator of the triples in a given shell.
  /** \param shell is the index of the shell to access.
   *
   *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
   */
  Facet_const_iterator shell_end( std::size_t shell ) const
  {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
      return _surface.end();
    return _shells[ shell+1 ];
  }

  /// gives a past-the-end iterator of the triples in a given shell.
  /** \param shell is the index of the shell to access.
   *
   *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
   *
   *  \warning Changes to a shell may invalidate the topology of the surface.
   */
  Facet_iterator shell_end( std::size_t shell )
  {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
        return _surface.end();
    return _shells[ shell+1 ];
  }
  
  /// gives an iterator to the first triple of the garbage facets
  /// that may be discarded if 2-manifold output is required.
  Facet_const_iterator garbage_begin() const { return _garbage.begin(); }
  /// gives an iterator to the first triple of the garbage facets
  /// that may be discarded if 2-manifold output is required.
  Facet_iterator garbage_begin() { return _garbage.begin(); }
    
  /// gives a past-the-end iterator of the triples of the garbage facets
  /// that may be discarded if 2-manifold output is required.
  Facet_const_iterator garbage_end() const { return _garbage.end(); }
  /// gives a past-the-end iterator of the triples of the garbage facets
  /// that may be discarded if 2-manifold output is required.
  Facet_iterator garbage_end() { return _garbage.end(); }

private:
  
  void deinit_shape()
  {
    if (_shape != NULL)
    {
      delete _shape;
      _shape = NULL;
    }
  }
  
  void clear_surface()
  {
    _shells.clear();
    _surface.clear();
    _garbage.clear();
    deinit_shape();
  }
  
  void collect_facets ()
  {
    // We check each of the facets: if it is not handled and either regular or singular,
    // we start collecting the next surface from here.
    for (Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit)
    {
      switch(_shape->classify (*fit))
      {
        case Shape::REGULAR:
          // Build a surface from the outer cell.
          if (_shape->classify(fit->first) == Shape::EXTERIOR)
            collect_shell (*fit);
          else
            collect_shell (_shape->mirror_facet (*fit));
          break;
        case Shape::SINGULAR:
          // Build a surface from both incident cells.
          collect_shell (*fit);
          collect_shell (_shape->mirror_facet (*fit));
          break;
        default:
          break;
      }
    }
  }

  void collect_facets_quick ()
  {
    // Collect all facets from the alpha-shape in an unordered fashion.
    for (Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit)
    {
      switch (_shape->classify(*fit))
      {
        case Shape::REGULAR:
          // Collect the outer cell.
          if (_shape->classify(fit->first) == Shape::EXTERIOR)
            _surface.push_back (ordered_facet_indices (*fit));
          else
            _surface.push_back (ordered_facet_indices (_shape->mirror_facet(*fit)));
          break;
        case Shape::SINGULAR:
          // Collect both incident cells.
          _surface.push_back (ordered_facet_indices (*fit));
          _surface.push_back (ordered_facet_indices (_shape->mirror_facet (*fit)));
          break;
        default:
          break;
      }
    }
  }
  
  inline bool is_handled( Cell_handle c, unsigned int li ) const
  {
    switch( li ) {
      case 0: return ( c->info()&1 ) != 0;
      case 1: return ( c->info()&2 ) != 0;
      case 2: return ( c->info()&4 ) != 0;
      case 3: return ( c->info()&8 ) != 0;
    }
    return false;
  }
  inline bool is_handled( const SFacet& f ) const { return is_handled( f.first, f.second ); }

  inline void mark_handled( Cell_handle c, unsigned int li )
  {
    switch( li ) {
      case 0: c->info() |= 1; return;
      case 1: c->info() |= 2; return;
      case 2: c->info() |= 4; return;
      case 3: c->info() |= 8; return;
    }
  }
  inline void mark_handled( SFacet f ) { mark_handled( f.first, f.second ); }
 
  inline void mark_opposite_handled( SFacet f )
  {

    Classification_type cl = _shape->classify (f);

    // If cell is singular, simply mark mirror facet as handled
    if (cl == Shape::SINGULAR)
    {
      SFacet mirror = _shape->mirror_facet (f);
      mark_handled (mirror);
    }
    // If cell is regular, get corresponding bubble and mark
    // facets of the other layer of the bubble as handled
    else if (cl == Shape::REGULAR)
    {
      SFacet fac = (_shape->classify (f.first) == Shape::EXTERIOR)
	? f
	: _shape->mirror_facet (f);

      typename std::map<SFacet, std::size_t>::iterator
	search = _map_f2b.find (fac);

      if (search == _map_f2b.end ())
	return;
      
      unsigned int layer = (_bubbles[search->second][0].find (fac) == _bubbles[search->second][0].end ())
	? 0 : 1;

      typename std::set<SFacet>::iterator it = _bubbles[search->second][layer].begin ();

      // If bubble has already been handled, no need to do it again
      if (is_handled (*it))
	return;
      
      for (;it != _bubbles[search->second][layer].end (); ++ it)
      {
        _garbage.push_back (ordered_facet_indices (*it));
        mark_handled (*it);
      }
      
    }


  }


  inline Facet ordered_facet_indices( const SFacet& f ) const
  {
    if( (f.second&1) == 0 )
      return make_array<unsigned int>( f.first->vertex( (f.second+2)&3 )->info(),
                                       f.first->vertex( (f.second+1)&3 )->info(),
                                       f.first->vertex( (f.second+3)&3 )->info() );
    else
      return make_array<unsigned int>( f.first->vertex( (f.second+1)&3 )->info(),
                                       f.first->vertex( (f.second+2)&3 )->info(),
                                       f.first->vertex( (f.second+3)&3 )->info() );
  }

  void collect_shell( const SFacet& f)
  {
    collect_shell (f.first, f.second);
  }
  void collect_shell( Cell_handle c, unsigned int li)
  {
    // Collect one surface mesh from the alpha-shape in a fashion similar to ball-pivoting.
    // Invariant: the facet is regular or singular.


    // To stop stack overflows: use own stack.
    std::stack<SFacet> stack;
    stack.push( SFacet(c, li) );

    SFacet f;
    Cell_handle n, p;
    int ni, pi;
    Vertex_handle a;
    Classification_type cl;
    bool processed = false;

    while( !stack.empty() ) {
      f = stack.top();
      stack.pop();

      // Check if the cell was already handled.
      // Note that this is an extra check that in many cases is not necessary.
      if( is_handled(f) )
        continue;

      // The facet is part of the surface.
      CGAL_assertion( !_shape->is_infinite(f) );
      mark_handled(f);
      // Output the facet as a triple of indices.
      _surface.push_back( ordered_facet_indices(f) );
      if( !processed ) {
        if (_separate_shells || _shells.size () == 0)
        {
          _shells.push_back( --_surface.end() );
          _index ++;
        }
        processed = true;
      }

      // Save in which shell the facet is stored
      _map_f2s[f] = _index-1;

      // If the surface is required to be manifold,
      // the opposite layer should be ignored
      if (_force_manifold)
        mark_opposite_handled (f);
		
      // Pivot over each of the facet's edges and continue the surface at the next regular or singular facet.
      for( int i = 0; i < 4; ++i ) {
        // Skip the current facet.
        if( i == f.second || is_handled(f.first, i) )
          continue;

        // Rotate around the edge (starting from the shared facet in the current cell) until a regular or singular facet is found.
        n = f.first;
        ni = i;
        a = f.first->vertex( f.second );
        cl = _shape->classify( SFacet(n, ni) );
	    
        while( cl != Shape::REGULAR && cl != Shape::SINGULAR ) {
          p = n;
          n = n->neighbor(ni);
          ni = n->index(a);
          pi = n->index(p);
          a = n->vertex(pi);
          cl = _shape->classify( SFacet(n, ni) );
        }

        // Continue the surface at the next regular or singular facet.
        stack.push( SFacet(n, ni) );
      }

    }

  }

  void detect_bubbles ()
  {
    std::set<Cell_handle> done;

    unsigned int nb_facets_removed = 0;
  
    unsigned int nb_skipped = 0;
    for (Cell_iterator cit = _shape->cells_begin (); cit != _shape->cells_end (); ++ cit)
    {
      if (_shape->is_infinite (cit))
	continue;
      if (_shape->classify (cit) != Shape::INTERIOR)
	continue;
      if (done.find (cit) != done.end ())
	continue;

      std::set<VEdge> borders;
      std::vector<Cell_handle> cells;
      std::stack<Cell_handle> todo;
      todo.push (cit);

      // Get all cells of volume and all borders
      while (!(todo.empty ()))
      {
        Cell_handle c = todo.top ();
        todo.pop ();

        if (!(done.insert (c).second))
          continue;

        cells.push_back (c);

        for (unsigned int i = 0; i < 4; ++ i)
        {
          if (_shape->classify (c->neighbor (i)) == Shape::INTERIOR)
            todo.push (c->neighbor (i));
          else
          {
            // Test if edge is proper border
            for (unsigned int j = 0; j < 3; ++ j)
            {
              unsigned int i0 = (i + j + 1)%4;
              unsigned int i1 = (i + (j+1)%3 + 1)%4;
              CGAL_assertion (i0 != i && i1 != i);
              Edge edge (c, i0, i1);

              if (_shape->classify (edge) != Shape::REGULAR)
                continue;

              VEdge vedge = (c->vertex (i0) < c->vertex (i1))
                ? std::make_pair (c->vertex (i0), c->vertex (i1))
                : std::make_pair (c->vertex (i1), c->vertex (i0));

              if (borders.find (vedge) != borders.end ())
                continue;

              SFacet_circulator start = _shape->incident_facets (edge);
              SFacet_circulator circ = start;
              unsigned int cnt = 0;
              do
              {
                if (_shape->classify (*circ) == Shape::SINGULAR
                    || _shape->classify (*circ) == Shape::REGULAR)
                  ++ cnt;
                ++ circ;
              }
              while (circ != start);

              // If edge is non-manifold, use as border
              if (cnt > 2)
              {
                borders.insert (vedge);
                continue;
              }

              // Else, if facets in cell are regular and angle is
              // under _border_angle limit, use as border
              SFacet f0 (c, i);
              SFacet f1 (c, (i + (j+2)%3 + 1)%4);

              if (_shape->classify (f0) != Shape::REGULAR
                  || _shape->classify (f1) != Shape::REGULAR)
                continue;
		      
              double angle = Geom_traits().compute_approximate_dihedral_angle_3_object()(vedge.first->point (),
                                                                                         vedge.second->point (),
                                                                                         c->vertex (i)->point (),
                                                                                         c->vertex ((i + (j+2)%3 + 1)%4)->point ());

              if (-_border_angle < angle && angle < _border_angle)
              {
                borders.insert (vedge);
              }
            }

          }
        }
      }

      int layer = -1;

      // Try to generate bubble from the volume found
      _bubbles.push_back (Bubble());
      std::set<SFacet> done;
      for (unsigned int c = 0; c < cells.size (); ++ c)
      {
        for (unsigned int ii = 0; ii < 4; ++ ii)
        {
          SFacet start = _shape->mirror_facet (SFacet (cells[c], ii));

          if (_shape->classify (start) != Shape::REGULAR)
            continue;

          if (done.find (start) != done.end ())
            continue;

          ++ layer;
	      
          std::stack<SFacet> stack;
          stack.push (start);

          SFacet f;
          Cell_handle n, p;
          int ni, pi;
          Vertex_handle a;
          Classification_type cl;

          // A bubble is well formed is the border contains one loop that
          // separates two layers.
          // If the number of layers is different than 2, the volume is completely ignored.
          while( !stack.empty() )
          {
            f = stack.top();
            stack.pop();

            if (!(done.insert (f).second))
              continue;

            if (_shape->classify (f.first) == Shape::EXTERIOR)
            {
              if (layer < 2)
              {
                _bubbles.back ()[layer].insert (f);
                _map_f2b[f] = _bubbles.size () - 1;
              }
              else
              {
                nb_facets_removed ++;
                mark_handled (f);
                _garbage.push_back (ordered_facet_indices (f));
              }
            }
            else
            {
              if (layer < 2)
              {
                _bubbles.back ()[layer].insert (_shape->mirror_facet (f));
                _map_f2b[_shape->mirror_facet (f)] = _bubbles.size () - 1;
              }
              else
              {
                nb_facets_removed ++;
                mark_handled (_shape->mirror_facet (f));
                _garbage.push_back (ordered_facet_indices (_shape->mirror_facet (f)));
              }
            }


            for( int i = 0; i < 4; ++i )
            {
              // Skip the current facet.
              if( i == f.second)
                continue;

              n = f.first;
              ni = i;
              a = f.first->vertex( f.second );
              cl = _shape->classify( SFacet(n, ni) );

              int n0 = -1, n1 = -1;
              bool n0found = false;
              for (int j = 0; j < 4; ++ j)
              {
                if (j != ni && j != f.second)
                {
                  if (n0found)
                  {
                    n1 = j;
                    break;
                  }
                  else
                  {
                    n0 = j;
                    n0found = true;
                  }
                }
              }

              VEdge vedge = (n->vertex (n0) < n->vertex (n1))
                ? std::make_pair (n->vertex (n0), n->vertex (n1))
                : std::make_pair (n->vertex (n1), n->vertex (n0));

              // If the edge is a border, propagation stops in this direction.
              if (borders.find (vedge) != borders.end ())
                continue;
		      
              while( cl != Shape::REGULAR && cl != Shape::SINGULAR ) {
                p = n;
                n = n->neighbor(ni);
                ni = n->index(a);
                pi = n->index(p);
                a = n->vertex(pi);
                cl = _shape->classify( SFacet(n, ni) );
              }
		      
              stack.push (SFacet (n, ni));

            }

          }
        }

      }
            
      // If number of layers is != 2, ignore volume and discard bubble
      if (layer != 1)
      {
        nb_skipped ++;
        for (unsigned int i = 0; i < 2; ++ i)
          for (typename std::set<SFacet>::iterator fit = _bubbles.back()[i].begin ();
               fit != _bubbles.back()[i].end (); ++ fit)
          {
            mark_handled (*fit);
            _map_f2b.erase (*fit);
            _garbage.push_back (ordered_facet_indices (*fit));
            nb_facets_removed ++;
          }
        _bubbles.pop_back ();
      }
      
    }
  }
  

  void fix_nonmanifold_edges()
  {

    typedef std::map<std::pair<VEdge, unsigned int>, std::set<Facet> > Edge_shell_map_triples;
    typedef typename Edge_shell_map_triples::iterator Edge_shell_map_triples_iterator;

    unsigned int nb_facets_removed = 0;

    unsigned int nb_nm_edges = 0;

    // Store for each pair edge/shell the incident facets
    Edge_shell_map_triples eshell_triples;
    std::map<Facet, SFacet> map_t2f;
  
    for (typename Map_facet_to_shell::iterator fit = _map_f2s.begin ();
         fit != _map_f2s.end (); ++ fit)
    {
      SFacet f = fit->first;
      Facet t = ordered_facet_indices (f);
      map_t2f[t] = f;
	      
      for (unsigned int k = 0; k < 3; ++ k)
      {
        Vertex_handle v0 = f.first->vertex ((f.second + k + 1)%4);
        Vertex_handle v1 = f.first->vertex ((f.second + (k+1)%3 + 1)%4);
        VEdge vedge = (v0 < v1) ? std::make_pair (v0, v1) : std::make_pair (v1, v0);
		  
        std::pair<Edge_shell_map_triples_iterator, bool>
          search = eshell_triples.insert (std::make_pair (std::make_pair (vedge, fit->second),
                                                          std::set<Facet>()));

        search.first->second.insert (t);
      }
    }

    for (Edge_shell_map_triples_iterator eit = eshell_triples.begin ();
         eit != eshell_triples.end (); ++ eit)
    {
      // If an edge has more than 2 incident facets for one shell, it is non-manifold
      if (eit->second.size () < 3)
	continue;

      ++ nb_nm_edges;

      Facet_iterator tit = _shells[eit->first.second];
      Facet_iterator end = (eit->first.second == _shells.size () - 1)
	? _surface.end () : _shells[eit->first.second + 1];

      // Remove facets until the edge is manifold in this shell
      while (tit != end && eit->second.size () > 2)
      {
        Facet_iterator current = tit ++;

        typename std::set<Facet>::iterator search = eit->second.find (*current);

        if (search != eit->second.end ())
        {
          if (current == _shells[eit->first.second])
            _shells[eit->first.second] = tit;

          _garbage.push_back (*current);
          _map_f2s.erase (map_t2f[*current]);
          _surface.erase (current);

          ++ nb_facets_removed;
          eit->second.erase (search);
        }

      }
	  
    }
  }

  template <typename T>
  struct operator_less
  {
    bool operator() (const T& a, const T& b) const
    {
      return &*a < &*b;
    }
  };


  void find_two_other_vertices(const SFacet& f, Vertex_handle v,
                               Vertex_handle& v1, Vertex_handle& v2)
  {
    Vertex_handle vother = f.first->vertex (f.second);
    bool v1found = false;
  
    for (unsigned int i = 0; i < 4; ++ i)
    {
      Vertex_handle vi = f.first->vertex (i);
      if (vi != v && vi != vother)
      {
        if (v1found)
        {
          v2 = vi;
          return;
        }
        else
        {
          v1 = vi;
          v1found = true;
        }
      }
    }
  }
  
  void fix_nonmanifold_vertices()
  {

    typedef ::CGAL::Union_find<SFacet> UF;
    typedef typename UF::handle UF_handle;


    typedef std::map<std::pair<Vertex_handle, unsigned int>, std::vector<SFacet> > Vertex_shell_map_facets;
    typedef typename Vertex_shell_map_facets::iterator Vertex_shell_map_facet_iterator;

    // For faster facet removal, we sort the triples of each shell as a preprocessing
    for (unsigned int i = 0; i < _shells.size (); ++ i)
    {
      Facet_iterator begin = _shells[i];
      Facet_iterator end = (i+1 == _shells.size ()) ? _surface.end () : _shells[i+1];
      
      Facetset tmp;
      tmp.splice (tmp.end(), _surface, begin, end);
      
      tmp.sort();
      _shells[i] = tmp.begin ();      
      _surface.splice(end, tmp, tmp.begin(), tmp.end());
    }

    unsigned int nb_facets_removed = 0;
    unsigned int nb_nm_vertices = 0;
    // Removing facets to fix non-manifold vertices might make some other vertices
    // become non-manifold, therefore we iterate until no facet needs to be removed.
    do
    {
      nb_nm_vertices = 0;
      nb_facets_removed = 0;

      // Store for each pair vertex/shell the incident facets
      Vertex_shell_map_facets vshell_facets;

      for (typename Map_facet_to_shell::iterator fit = _map_f2s.begin ();
	   fit != _map_f2s.end (); ++ fit)
      {
        SFacet f = fit->first;
	  
        for (unsigned int k = 0; k < 3; ++ k)
        {
          Vertex_handle v = f.first->vertex ((f.second+k+1)%4);

          std::pair<Vertex_shell_map_facet_iterator, bool>
            search = vshell_facets.insert (std::make_pair (std::make_pair (v, fit->second),
                                                           std::vector<SFacet>()));
          search.first->second.push_back (f);

        }
	  
      }

      for (Vertex_shell_map_facet_iterator fit = vshell_facets.begin ();
	   fit != vshell_facets.end (); ++ fit)
      {
        if (fit->second.size () < 2)
          continue;

        Vertex_handle vit = fit->first.first;
        unsigned int shell = fit->first.second;

        UF uf;
        std::map<SFacet, UF_handle> map_f2h;
	  
        for (unsigned int i = 0; i < fit->second.size (); ++ i)
          map_f2h.insert (std::make_pair (fit->second[i], uf.make_set (fit->second[i])));

        std::map<Vertex_handle, SFacet> map_v2f;
	    
        for (unsigned int i = 0; i < fit->second.size (); ++ i)
        {
          Vertex_handle v1, v2;
          find_two_other_vertices (fit->second[i], vit, v1, v2);
          std::pair<typename std::map<Vertex_handle, SFacet>::iterator, bool>
            insertion1 = map_v2f.insert (std::make_pair (v1, fit->second[i]));
          if (!(insertion1.second))
            uf.unify_sets (map_f2h[fit->second[i]], map_f2h[insertion1.first->second]);
          std::pair<typename std::map<Vertex_handle, SFacet>::iterator, bool>
            insertion2 = map_v2f.insert (std::make_pair (v2, fit->second[i]));
          if (!(insertion2.second))
            uf.unify_sets (map_f2h[fit->second[i]], map_f2h[insertion2.first->second]);
        }

        if (uf.number_of_sets () > 1)
        {
          ++ nb_nm_vertices;

          typedef std::map<UF_handle, std::vector<SFacet>, operator_less<UF_handle> > Map_uf_sets;
          Map_uf_sets map_h2f;
          for (unsigned int i = 0; i < fit->second.size (); ++ i)
          {
            UF_handle handle = uf.find (map_f2h[fit->second[i]]);

            std::pair<typename Map_uf_sets::iterator, bool>
              insertion = map_h2f.insert (std::make_pair (handle, std::vector<SFacet>()));

            insertion.first->second.push_back (fit->second[i]);
          }

          typename Map_uf_sets::iterator largest = map_h2f.end ();
          std::size_t nb_largest = 0;
          for (typename Map_uf_sets::iterator ufit = map_h2f.begin (); ufit != map_h2f.end (); ++ ufit)
          {
            std::size_t size = ufit->second.size ();
            if (size > nb_largest)
            {
              nb_largest = size;
              largest = ufit;
            }
          }

          std::vector<Facet> triples;

          for (typename Map_uf_sets::iterator ufit = map_h2f.begin (); ufit != map_h2f.end (); ++ ufit)
          {
            if (ufit == largest)
              continue;
            for (unsigned int i = 0; i < ufit->second.size (); ++ i)
            {
              _map_f2s.erase (ufit->second[i]);
              triples.push_back (ordered_facet_indices (ufit->second[i]));
            }
          }
          std::sort (triples.begin (), triples.end ());

          Facet_iterator tit = _shells[shell];
          Facet_iterator end = (shell == _shells.size () - 1)
            ? _surface.end () : _shells[shell + 1];

          unsigned int tindex = 0;
	      
          while (tit != end && tindex < triples.size ())
          {
            Facet_iterator current = tit ++;

            if (*current == triples[tindex])
            {
              if (current == _shells[shell])
                _shells[shell] = tit;

              _garbage.push_back (*current);
              _surface.erase (current);

              ++ nb_facets_removed;
              ++ tindex;
            }

          }
        }

      }

    }
    while (nb_nm_vertices != 0);

  }

};

  
} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_ALPHA_SHAPE_MESHER_H
