// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_CELL_BASE_3_H
#define CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_CELL_BASE_3_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>

#include <CGAL/Delaunay_triangulation_cell_base_3.h>

namespace CGAL {

  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstructionRef

  The class `Advancing_front_surface_reconstruction_cell_base_3` is the default
  cell type for the class  `Advancing_front_surface_reconstruction`.

  \tparam Traits has to be a model of `DelaunayTriangulationTraits_3`.

  \tparam Cb has to be a model of `DelaunayTriangulationCellBase_3`.
  */
  template < typename Traits, typename Cb = Delaunay_triangulation_cell_base_3<Traits> >
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
        _smallest_radius_facet_tab(nullptr), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(nullptr), _squared_radius(nullptr)
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
        _smallest_radius_facet_tab(nullptr), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(nullptr), _squared_radius(nullptr)
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
        _smallest_radius_facet_tab(nullptr), selected_facet(0)
#ifdef AFSR_LAZY
      , _circumcenter(nullptr), _squared_radius(nullptr)
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
      if (_smallest_radius_facet_tab != nullptr)
        delete[] _smallest_radius_facet_tab;
#ifdef AFSR_LAZY
      if (_circumcenter != nullptr)
        delete _circumcenter;
      if (_squared_radius != nullptr)
        delete _squared_radius;
#endif
    }

    //--------------------  MEMBER FUNCTIONS ----------------------------
  public:

    inline void clear()
    {
      if (_smallest_radius_facet_tab != nullptr)
        delete[] _smallest_radius_facet_tab;
      _smallest_radius_facet_tab = nullptr;
      selected_facet = 0;
#ifdef AFSR_LAZY
      if (_circumcenter != nullptr)
        delete _circumcenter;
      _circumcenter = nullptr;
      if (_squared_radius != nullptr)
        delete _squared_radius;
      _squared_radius = nullptr;
#endif
    }

    //-------------------------------------------------------------------
    inline coord_type smallest_radius(const int& i)
    {
      if (_smallest_radius_facet_tab == nullptr)
        return -1;
      return _smallest_radius_facet_tab[i];
    }

    inline void set_smallest_radius(const int& i, const coord_type& c)
    {
      if (_smallest_radius_facet_tab == nullptr)
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
      if (_smallest_radius_facet_tab==nullptr)
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
