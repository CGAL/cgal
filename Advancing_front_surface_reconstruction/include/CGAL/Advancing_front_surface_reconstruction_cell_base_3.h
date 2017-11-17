// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_CELL_BASE_3_H
#define CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_CELL_BASE_3_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


#include <CGAL/Triangulation_cell_base_3.h>

namespace CGAL {

  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstruction

  The class `Advancing_front_surface_reconstruction_cell_base_3` is the default
  cell type for the class  `Advancing_front_surface_reconstruction`.

  \tparam Traits has to be a model of `DelaunayTriangulationTraits_3`.

  \tparam Cb has to be a model of `TriangulationCellBase_3`.
  */
  template < typename Traits, typename Cb = Triangulation_cell_base_3<Traits> >
  class Advancing_front_surface_reconstruction_cell_base_3 : public Cb
  {

  public:
    template < typename TDS2 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
      typedef Advancing_front_surface_reconstruction_cell_base_3<Traits,Cb2>                    Other;
    };

    typedef typename Cb::Vertex_handle Vertex_handle;
    typedef typename Cb::Cell_handle Cell_handle;

  private:

#ifdef AFSR_FACET_NUMBER
    int _facet_number[4];
#endif
    typedef double coord_type;
#ifdef AFSR_LAZY
    typedef typename CGAL::Simple_cartesian<coord_type>::Point_3 D_Point;
#endif
    //-------------------- DATA MEMBERS ---------------------------------

    coord_type* _smallest_radius_facet_tab;
    unsigned char selected_facet;
#ifdef AFSR_LAZY
    D_Point* _circumcenter;
    coord_type* _squared_radius;
#endif

    //-------------------- CONSTRUCTORS ----------------------------------

  public:

    Advancing_front_surface_reconstruction_cell_base_3()
      : Cb(),
        _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef AFSR_FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }

    Advancing_front_surface_reconstruction_cell_base_3(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
      : Cb( v0, v1, v2, v3),
        _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }

    Advancing_front_surface_reconstruction_cell_base_3(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
                                                       Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3)
      : Cb(v0,  v1,  v2, v3,
           n0,  n1,  n2, n3),
        _smallest_radius_facet_tab(NULL), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(NULL), _squared_radius(NULL)
#endif
    {
#ifdef AFSR_FACET_NUMBER
      for(int i = 0; i < 4; i++){
	_facet_number[i] = -1;
      }
#endif
    }

    //-------------------- DESTRUCTOR -----------------------------------

    inline ~Advancing_front_surface_reconstruction_cell_base_3()
    {
      if (_smallest_radius_facet_tab != NULL)
        delete[] _smallest_radius_facet_tab;
#ifdef AFSR_LAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      if (_squared_radius != NULL)
	delete _squared_radius;
#endif
    }

    //--------------------  MEMBER FUNCTIONS ----------------------------
  public:

    inline void clear()
    {
      if (_smallest_radius_facet_tab != NULL)
	delete[] _smallest_radius_facet_tab;
      _smallest_radius_facet_tab = NULL;
      selected_facet = 0;
#ifdef AFSR_LAZY
      if (_circumcenter != NULL)
	delete _circumcenter;
      _circumcenter = NULL;
      if (_squared_radius != NULL)
	delete _squared_radius;
      _squared_radius = NULL;
#endif
    }

    //-------------------------------------------------------------------
    inline coord_type smallest_radius(const int& i)
    {
      if (_smallest_radius_facet_tab == NULL)
	return -1;
      return _smallest_radius_facet_tab[i];
    }

    inline void set_smallest_radius(const int& i, const coord_type& c)
    {
      if (_smallest_radius_facet_tab == NULL)
	{
	  _smallest_radius_facet_tab = new coord_type[4];
	  for(int i = 0; i < 4; i++)
	    _smallest_radius_facet_tab[i] = -1;
	}
      _smallest_radius_facet_tab[i] = c;
    }

    // pour un controle de l'allocation memoire... utile???
    inline bool alloc_smallest_radius_tab(coord_type* ptr)
    {
      if (_smallest_radius_facet_tab==NULL)
	{
	  _smallest_radius_facet_tab = ptr;
	  return true;
	}
      return false;
    }


    //-------------------------------------------------------------------

#ifdef FACET_NUMBER
    void set_facet_number(int i, int n){}
    {
      _facet_number[i] = n;
    }

    int facet_number(int i)
    {
      return _facet_number[i];
    }
#else
    void set_facet_number(int, int){}
    int facet_number(int){return 0;}
#endif


    //-------------------------------------------------------------------

    inline void select_facet(const int& i)
    {
      selected_facet |= (1 << i);
    }

    inline void unselect_facet(const int& i)
    {
      selected_facet &= (15 - (1 << i));
    }

    inline bool is_selected_facet(const int& i)
    {
      return (selected_facet & (1 << i)) != 0;
    }

    inline bool has_facet_on_surface(const int& i)
    {
      return (selected_facet & (1 << i)) != 0;
    }

#ifdef AFSR_LAZY

    //-------------------------------------------------------------------

    inline D_Point* lazy_circumcenter()
    {
      return _circumcenter;
    }

    inline void set_lazy_circumcenter(const D_Point& p)
    {
      _circumcenter = new D_Point(p);
    }

    //-------------------------------------------------------------------

    inline coord_type* lazy_squared_radius()
    {
      return _squared_radius;
    }

    inline void set_lazy_squared_radius(const coord_type& sr)
    {
      _squared_radius = new coord_type(sr);
    }

#endif //AFSR_LAZY

  };

} // namespace CGAL

#endif // CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_CELL_BASE_3_H
