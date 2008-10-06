// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
// File          : include/SoX/GAPS/Z_stack.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef SoX_GAPS_Z_STACK_H
#define SoX_GAPS_Z_STACK_H 1

/*! \file SoX/GAPS/Z_stack.h
 * \brief definition of Z_stack class template
 */

#include <CGAL/config.h>

#include <map>

#include <boost/optional.hpp>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_accessor.h>
#include <CGAL/Arrangement_2l/Z_stack_helpers.h>
#include <CGAL/Arrangement_2l/Adjacencies_3.h>

namespace SoX {

namespace Intern {

template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_cell_rep {
public:
  
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! this instance's second template parameter
    typedef DcelData Dcel_data;
    
    //! the class itself
    typedef Z_cell_rep< Surface_z_at_xy_isolator_traits, Dcel_data > Self;

    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of z_stack
    typedef Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data > Z_stack;
    
    //! type of Surface less
    typedef typename Surface_3::Surface_less_than Surface_less_than;
    
    //! type of sheet map
    typedef std::map< Surface_3, int, Surface_less_than > Sheet_map;

    //!  type of Boundary
    typedef typename Z_at_xy_isolator::Boundary Boundary;
    
    Z_cell_rep(const Z_stack& z_stack, int z) :
        _m_z_stack(z_stack),
        _m_z(z),
        _m_min_length(0) {
    }
    
private:    
    //! the z_stack
    Z_stack _m_z_stack;
    
    //! z-position in stacke
    int _m_z;

    //! sheet of surface
    Sheet_map _m_sheets;

    //! current min_length
    Boundary _m_min_length;
    
    //! surface that induces min_length
    Surface_3 _m_min_length_surface;

    //! surface sheet that induces min_length
    int _m_min_length_surface_sheet;

#if 0
    //!
    std::map< CGAL::Object, int, Object_less_than > _m_adjacent_items;
#endif

    friend class Z_cell< Surface_z_at_xy_isolator_traits, Dcel_data >;
    
};

template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_cell : public 
CGAL::Handle_with_policy< Z_cell_rep <SurfaceZAtXyIsolatorTraits, DcelData > > 
{
public:
    
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! this instance's second template parameter
    typedef DcelData Dcel_data;
    
    //! type of rep
    typedef Z_cell_rep< Surface_z_at_xy_isolator_traits, Dcel_data > Rep;
    
    //! type of base
    typedef CGAL::Handle_with_policy< Rep > Base;

    //! the class itself
    typedef Z_cell< Surface_z_at_xy_isolator_traits, Dcel_data > Self;

    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

    //! type of sheet map
    typedef typename Rep::Sheet_map Sheet_map;
    

    //! type of sheet map iterator
    typedef typename Sheet_map::const_iterator Sheet_const_iterator;

    //! type of Boundary
    typedef typename Rep::Boundary Boundary;
    
    //! type of z_stack
    typedef Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data > Z_stack;
    
    //!\name Constructor 
    //!@{

    //! standard constructor for given \c z-stack  and \c int 
    Z_cell(const Z_stack& z_stack, int z) :
        Base(Rep(z_stack, z)) {
    }
    
    //!@}
    
    //!\name Surfaces, sheets, and levels
    //!@{
    
    //! returns z-entry of cell
    int z() {
        return this->ptr()->_m_z;
    }
    
    //! returns number of involved sheets
    int number_of_sheets() const {
        return static_cast< int >(this->ptr()->_m_sheets.size());
    }
    
    //! returns sheet number of surface in given cell
    int sheet_number(const Surface_3& surface) const {
        typename Sheet_map::const_iterator it = 
            this->ptr()->_m_sheets.find(surface);
        if (it != this->ptr()->_m_sheets.end()) {
            CGAL_postcondition(it->second >= 0);
            return it->second;
        }
        return -1;
    }

    //! returns beginning of sheets
    Sheet_const_iterator sheets_begin() const {
        return this->ptr()->_m_sheets.begin();
    }

    //! returns past-the-end iterator for sheets
    Sheet_const_iterator sheets_end() const {
        return this->ptr()->_m_sheets.end();
    }
    
    // TODO random surface

    //! return surface whose z-interval is minimal in cell
    const Surface_3& minimal_length_surface() const {
        CGAL_assertion(!this->ptr()->_m_sheets.empty());
        return this->ptr()->_m_min_length_surface;
    }

    //! return sheet of surface whose z-interval is minimal in cell
    int minimal_length_surface_sheet() const {
        CGAL_assertion(!this->ptr()->_m_sheets.empty());
        return this->ptr()->_m_min_length_surface_sheet;
    }
    
    //!@}
    
    //! Approximation
    //!@{

    //! returns double approximation of z-cell
    std::pair< double, double > to_double() const {
        CGAL_assertion(!this->ptr()->_m_sheets.empty());
        Surface_3 surface = this->ptr()->_m_min_length_surface;
        int sheet = this->ptr()->_m_min_length_surface_sheet;
        return std::make_pair(
                CGAL::to_double(
                        this->ptr()->_m_z_stack._isolator(surface).
                        left_boundary(sheet)
                ),
                CGAL::to_double(
                        this->ptr()->_m_z_stack._isolator(surface).
                        right_boundary(sheet)
                )
        );
    }
    
    //!@}

private:    
    
    //!\name Modifying
    //!@{
    
    //! add sheet \c sheet of \c surface to cell
    void _add_sheet(const Surface_3& surface, int sheet, 
                    Boundary interval_length) {
        CGAL_precondition(sheet >= 0);
        this->ptr()->_m_sheets.insert(std::make_pair(surface, sheet));
        if (interval_length < this->ptr()->_m_min_length ||
            this->ptr()->_m_min_length == Boundary(0)) {
            this->ptr()->_m_min_length = interval_length;
            this->ptr()->_m_min_length_surface = surface;
            this->ptr()->_m_min_length_surface_sheet = sheet;
        }
            
    }

    
    //! makes cell unique for \c z's cell of \c z_stack
    void _make_unique_for(const Z_stack& z_stack, int z) {
        if (z_stack.id() != this->ptr()->_m_z_stack.id()) {
            this->copy_on_write();
            this->ptr()->_m_z_stack = z_stack;
        }
        this->ptr()->_m_z = z;
    }

    //!@}
    
public:
    //!\name IO
    //!@{
    
    /*!\brief
     * prints pretty-formated version of z_stack
     */
    void pretty_print(std::ostream& os) const {
        os << "(z=" << this->ptr()->_m_z << ",";
        std::pair< double, double > dp = to_double();
        os << "[" << dp.first << "," << dp.second << "],";
        for (typename Sheet_map::const_iterator iit = 
                 this->ptr()->_m_sheets.begin();
             iit != this->ptr()->_m_sheets.end(); iit++) {
            os << "<" << iit->first.id() << "," << iit->second << ">";
            typename Sheet_map::const_iterator tmp = iit;
            tmp++;
            if (tmp != this->ptr()->_m_sheets.end()) {
                os << ",";
            }
        }
        os << ")";
    }
    
    // friends
    friend class SoX::Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data >;

}; // Z_cell



template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_stack_rep {

public:
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    //! this instance's second template parameter
    typedef DcelData Dcel_data;
    
    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );

    //! the class itself
    typedef Z_stack_rep< Surface_z_at_xy_isolator_traits, Dcel_data > 
    Self;
    
    //! type of Surface less
    typedef typename Surface_3::Surface_less_than Surface_less_than;
    
    //! type of z-isolator
    typedef std::map< Surface_3, Z_at_xy_isolator, Surface_less_than > 
    Z_at_xy_isolators;

    //! type of z-cell
    typedef Z_cell< Surface_z_at_xy_isolator_traits, Dcel_data > Z_cell;

    //! type of z-cell container
    typedef std::list< Z_cell > Z_cell_container;

    //! type of const_iterator
    typedef typename Z_cell_container::const_iterator Z_cell_const_iterator;
    
    //!\name Constructors
    //!@{


    //! default constructor
    Z_stack_rep() {}
    
    /*!\brief
     * standard constructor with the help of a point \c xy
     */
    Z_stack_rep(SoX::Dcel_feature feature,
                const Dcel_data* data, const Point_2& pt) :
        _m_dcel_feature(feature),
        _m_dcel_data(data),
        _m_point(pt) {
        _m_z_at_xy_isolators.clear();
        _m_z_at_xy_empty_isolators.clear();
        _m_z_cells.clear();
    }

    //!@}

private:
    // data members
    //! type of dcel_feature
    mutable SoX::Dcel_feature _m_dcel_feature;

    //! storage for isolators
    mutable Z_at_xy_isolators _m_z_at_xy_isolators;

    //! storage for empty isolators
    mutable Z_at_xy_isolators _m_z_at_xy_empty_isolators;
    
    //! z-cells
    // taken a list since we need to update the cells
    // actually we also require RAM with the help of operator[]
    mutable Z_cell_container _m_z_cells;

    //! dcel data
    const Dcel_data* _m_dcel_data;
    
    //! stored point
    mutable Point_2 _m_point;

    // friends
    friend class Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data >;
};

} // namespace Intern

/*!\brief
 * Class to describe the order of surfaces of a planar (refinable) point.
 */
template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_stack : public 
::CGAL::Handle_with_policy< 
Intern::Z_stack_rep< SurfaceZAtXyIsolatorTraits, DcelData > > {
public:
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;

    //! this instance's second template parameter
    typedef DcelData Dcel_data;
    
    SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of rep
    typedef Intern::Z_stack_rep< Surface_z_at_xy_isolator_traits, Dcel_data > 
    Rep; 
    
    //! base type
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //! the class itself
    typedef Z_stack< Surface_z_at_xy_isolator_traits, Dcel_data > Self;

    //! type of isolator map
    typedef typename Rep::Z_at_xy_isolators Z_at_xy_isolators;
    
    //! type of Z_stack_cell
    typedef typename Rep::Z_cell Z_cell;
    
    //! type of Z_cell_container
    typedef typename Rep::Z_cell_container Z_cell_container;

    //! type of const_iterator
    typedef typename Rep::Z_cell_const_iterator Z_cell_const_iterator;
    
    //!\name Constructors
    //!@{
    
    //! defualt constructor
    Z_stack() {};

    /*!\brief
     * standard constructor from a refinable point
     */
    Z_stack(SoX::Dcel_feature feature, 
            const Dcel_data* data, const Point_2& pt) :
        Base(Rep(feature,data,pt)) {
    };

    /*!\brief
     * copies z-stack \c self, but replaces feature
     */
    Z_stack(const Self& self, 
            SoX::Dcel_feature feature, const Dcel_data* data) :
        Base(*self.ptr()) {
        this->copy_on_write();
        this->ptr()->_m_dcel_feature = feature;
        this->ptr()->_m_dcel_data = data;
    };
    
    //!@}

private:
    //!\name Query isolators
    //!@{

    /*!\brief 
     * returns \c true if z_stack already knows an isolator
     * at stored point for \c surface
     */
    bool _knows_isolator(const Surface_3& surface, bool& empty) const {
        empty = true;
        typename Z_at_xy_isolators::const_iterator it =
            this->ptr()->_m_z_at_xy_isolators.find(surface);
        if (it != this->ptr()->_m_z_at_xy_isolators.end()) {
            return true;
        } 
        // else
        empty = false;
        return (this->ptr()->_m_z_at_xy_empty_isolators.find(surface) !=
                this->ptr()->_m_z_at_xy_empty_isolators.end());
    }

    /*!\brief returns isolator for surface
     *
     * \pre _knows_isolator(surface)
     */
    Z_at_xy_isolator _isolator(const Surface_3& surface) const {
        CGAL_precondition_code(bool empty;);
        CGAL_precondition(_knows_isolator(surface, empty));
        typename Z_at_xy_isolators::const_iterator it =
            this->ptr()->_m_z_at_xy_isolators.find(surface);
        if (it != this->ptr()->_m_z_at_xy_isolators.end()) {
            return it->second;
        } 
        // else
        return this->ptr()->_m_z_at_xy_empty_isolators.find(surface)->second;
    }
    
public:
    
     /*!\brief 
     * returns \c true if z_stack already knows an isolator
     * at stored point for \c surface
     */
    inline
    bool knows_isolator(const Surface_3& surface, bool& empty) const {
        return _knows_isolator(surface, empty);
    }

    /*!\brief returns isolator for surface
     *
     * \pre knows_isolator(surface)
     */
    inline
    Z_at_xy_isolator isolator(const Surface_3& surface) const {
        return _isolator(surface);
    }
    
    //!@}

private:
    //!\name Adding surfaces
    //!@{
    
    /*!\brief
     * add surface to z_stack and maintain sequence of cells
     */
    void _add_surface(const Surface_3& surface, 
                      const Z_at_xy_isolator& isolator) const {
        
        // check whether surface influences z_stack
        if (isolator.number_of_real_roots() == 0) {
            typename Z_at_xy_isolators::iterator it =
                this->ptr()->_m_z_at_xy_empty_isolators.find(surface);
            if (it == this->ptr()->_m_z_at_xy_empty_isolators.end()) {
                this->ptr()->_m_z_at_xy_empty_isolators.insert(
                        it, std::make_pair(surface,isolator)
                );
            }
            // surface has no cut with line parallel to z-axis 
            // -> nothing else to do
            return;
        } 
        // else
        
        // first check whether z_stack aldready knows something about the curve
        typename Z_at_xy_isolators::iterator it =
            this->ptr()->_m_z_at_xy_isolators.find(surface);
        if (it == this->ptr()->_m_z_at_xy_isolators.end()) {
            this->ptr()->_m_z_at_xy_isolators.insert(
                    it, std::make_pair(surface,isolator)
            );
        } else {
            // surface already exists in z_stack -> no further changes on stack
            return;
        }
        
        // maintenance
        _update_cells(surface, isolator);
    } 
    
private:
    
    /*!\brief
     * updates the cells, i.e., find intersections!
     */
    void _update_cells(const Surface_3& surface2, 
                       const Z_at_xy_isolator& isolator2) const {

        int roots2 = isolator2.number_of_real_roots();
        
        if (this->ptr()->_m_z_cells.empty()) {
            // simple case, fast exit
            for (int i = roots2 - 1; i >= 0; i--) {
                Z_cell cell(*this, i);
                cell._add_sheet(surface2, i, isolator2.length(i));
                this->ptr()->_m_z_cells.push_front(cell);
            }
            // exit function
            return;
        }
        
        // else
        // for each 0 <= i < roots2: 
        // find position of (isol2.left_boundary(i), isol2.right_boundary(i)) 
        // in sequence of cells
        
        // we try to apply filters
        // FACE: No intersection can take place
        // EDGE/VERTEX: No projected intersection -> no equal_z-test
        // EDGE: If mult == 1 => only one intersection
        // VERTEX: not isolated -> adjacency

        // there are situation where filters do not work
        // EDGE: If mult > 1 => equal-z
        // VERTEX: isolated
        
        SoX::Dcel_feature feature = this->ptr()->_m_dcel_feature;
        
        typedef typename Surface_3::Surface_less_than Surface_less_than;
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits >
            Restricted_cad_3;

        Restricted_cad_3 cad = 
            Restricted_cad_3::cad_cache()(
                    this->ptr()->_m_dcel_data->_rs_id()
            );
        
        
        SoX::Intern::Refineable_interval_helper< Z_at_xy_isolator > iv_helper;
        
        std::set< Surface_3, Surface_less_than > projected_sil1s;
        std::map< Surface_3, std::pair< int, 
                                        std::list< std::pair < int, int > >
            >, Surface_less_than > projected_cut1s;
        
        bool projected_sil2 = false;
        bool projected_cut2 = false;

        typedef std::pair< Surface_3, int > Sheet;
        
        //! type of Less of sheet
        typedef CGAL::Pair_lexicographical_less_than< Surface_3, int, 
            Surface_less_than, std::less< int > > Sheet_less;
        
        typedef std::map< Sheet, int, Sheet_less > Adj_x_map;

        Adj_x_map adj_x_map;
        
        if (feature != FACE) {

            // check existing projected intersections and boundaries
            projected_sil2 = 
                this->ptr()->_m_dcel_data->_has_silhouette(surface2);
            
            for (typename Z_at_xy_isolators::iterator it =
                     this->ptr()->_m_z_at_xy_isolators.begin();
                 it != this->ptr()->_m_z_at_xy_isolators.end(); it++) {
                if (this->ptr()->_m_dcel_data->_has_silhouette(it->first)) {
                    projected_sil1s.insert(it->first);
                }
                if (this->ptr()->_m_dcel_data->_has_cut(it->first, surface2)) {
                    // compute overlaps
                    std::list< std::pair< int, int > > overlaps;
                    iv_helper.refined_overlaps(it->second, isolator2,
                                               std::back_inserter(overlaps));
                    
                    // and cache them
                    projected_cut1s.insert(
                            std::make_pair(
                                    it->first,
                                    std::make_pair(
                                            // make sense for EDGE 
                                            // in the test *mult == 1
                                            (feature == SoX::VERTEX ? 
                                             -1 : 
                                             *this->ptr()->
                                             _m_dcel_data->multiplicity_of_cut(
                                                     it->first, surface2
                                             )
                                            ),
                                            overlaps
                                    )
                            )
                    );
                    projected_cut2 = true;
                }
            }
            
            // for VERTEX:
            if (feature == SoX::VERTEX && projected_cut2) { 
                // if there is a projected cut over a vertex
                
                
                // Create adjacency information at z_stack-location from 
                // neighboring cells and for one non-isolated
                // surfaces of each cell
                
                // search for neighboring cells that form an intersection
                // of surface2 with at least one other surface
                typename Restricted_cad_3::Vertex_const_handle vh;

                CGAL_assertion_code(bool check = )
                    CGAL::assign(vh,this->ptr()->_m_dcel_data->_dcel_handle());
                CGAL_assertion(check);
                
                if (!vh->is_isolated()) {
                 
                    std::pair< Self, SoX::Dcel_feature > z_stack2to = 
                        vh->data()->_z_stack_of_surface(surface2);
                    
                    boost::optional< 
                        std::pair< SoX::Dcel_feature,CGAL::Object >
                        > opt2to = vh->data()->_dcel(surface2);
                    CGAL_assertion(opt2to);
                    CGAL_assertion(
                            opt2to->first == z_stack2to.second
                    );
                    
                    typename Restricted_cad_3::
                        Halfedge_around_vertex_const_circulator circ, start =
                        vh->incident_halfedges();
                    circ = start;
                    
                    do {
                        ++circ;
                        Self z_stack = 
                            circ->data()->_z_stack_for_halfedge_handle(circ);
                        
                        for (typename Z_cell_container::iterator it =
                                 z_stack.z_cells_begin();
                             it != z_stack.z_cells_end(); it++) {
                            // going over all cells
                            
                            // check whether surface2 and at least some
                            // other curve is involved
                            int sheet2 = it->sheet_number(surface2);
                            
                            if (it->number_of_sheets() > 1 && 
                                sheet2 >= 0) {
                                
                                // compute adjacencies of surface2 towards goal
                                std::pair< Self, SoX::Dcel_feature > 
                                    z_stack2from =
                                    circ->data()->_z_stack_of_surface(
                                            surface2
                                    );
                                boost::optional< 
                                    std::pair< SoX::Dcel_feature,CGAL::Object >
                                    > opt2from = circ->data()->_dcel(surface2);
                                CGAL_assertion(opt2from);
                                CGAL_assertion(
                                        opt2from->first == z_stack2from.second
                                );
                                
                                std::pair< int, int > adj2 = 
                                    z_stack2from.first.adjacency(
                                            surface2, sheet2,
                                            opt2from->first,
                                            opt2from->second,
                                            z_stack2to.first,
                                            opt2to->first,
                                            opt2to->second
                                    );
                                
                                // TODO check adj [1,0]
                                CGAL_assertion(adj2.first >= 0 && 
                                               adj2.first <= adj2.second);
                                
                                // and adjacancies of other surfaces 
                                // towards goal
                                for (typename Z_cell::Sheet_const_iterator 
                                         sit = it->sheets_begin();
                                     sit != it->sheets_end(); sit++) {
                                    
                                    Surface_3 surface1 = sit->first;
                                    int sheet1 = sit->second;
                                    
                                    if (surface1.id() == surface2.id()) {
                                        continue;
                                    }
                                    
                                    std::pair< Self, SoX::Dcel_feature > 
                                        z_stack1from = 
                                        circ->data()->_z_stack_of_surface(
                                                surface1
                                        );
                                    boost::optional< 
                                        std::pair< SoX::Dcel_feature,
                                        CGAL::Object >
                                        > opt1from = 
                                        circ->data()->_dcel(surface1);
                                    CGAL_assertion(opt1from);
                                    CGAL_assertion(
                                            opt1from->first == 
                                            z_stack1from.second
                                    );

                                    std::pair< Self, SoX::Dcel_feature > 
                                        z_stack1to = 
                                        vh->data()->_z_stack_of_surface(
                                                surface1
                                        );
                                    boost::optional< 
                                        std::pair< SoX::Dcel_feature,
                                        CGAL::Object >
                                        > opt1to = 
                                        vh->data()->_dcel(surface1);
                                    CGAL_assertion(opt1to);
                                    CGAL_assertion(
                                            opt1to->first == 
                                            z_stack1to.second
                                    );
                                    
                                    std::pair< int, int > adj1 = 
                                        z_stack1from.first.adjacency(
                                                surface1, sheet1, 
                                                opt1from->first,
                                                opt1from->second,
                                                z_stack1to.first,
                                                opt1to->first,
                                                opt1to->second);
                                                
                                    // TODO check adj [1,0]
                                    CGAL_assertion(adj1.first >= 0 && 
                                                   adj1.first <= adj1.second);
                                    
                                    // and store other surface-sheets as key
#if 0
                                    std::cout << "s1: " << surface1.id()
                                              << "t1: " << adj1.first
                                              << std::endl;
                                    std::cout << "s2: " << surface2.id()
                                              << "t2: " << adj2.first
                                              << std::endl;
#endif
                                    adj_x_map.insert(
                                            std::make_pair(
                                                    std::make_pair(
                                                            surface1,
                                                            adj1.first), 
                                                    adj2.first)
                                    );
                                    
                                    // whenever we later see the stored 
                                    // surface-sheet-key, we've actually 
                                    // detected an intersection
                                }
                            }
                        }
                    } while (circ != start);
                }
            }
        } 
        // else: if z_stack is for FACE, we will not detect any intersection, 
        // as we require surface to be "coprime"

        // setup for filters finished
        // actual merge starts now
        int z = 0;
        int root2 = 0;
        
        typename Z_cell_container::iterator it =
            this->ptr()->_m_z_cells.begin();
       
        bool still_possible_equal = true;
        
        while (root2 < roots2 && it != this->ptr()->_m_z_cells.end()) {
            // read isolator and correspondent root
            // TASK make selection of surface1 random!
            const Surface_3& surface1 = it->minimal_length_surface();
            const int root1 = it->minimal_length_surface_sheet();
            
            const Z_at_xy_isolator& isolator1 =
                this->ptr()->_m_z_at_xy_isolators.find(surface1)->second;
            
            switch (iv_helper.overlap_or_order(isolator1, root1, 
                                               isolator2, root2)) {
            case CGAL::SMALLER: {
                //std::cout << "SMALLER root1: " << root1 << " root2: " 
                //          << root2 << std::endl;
                it->_make_unique_for(*this, z);
                
                // go to next existing cell in z_stack
                it++;
                z++;

                // the next combination can be equal again
                still_possible_equal = true;
                break;
            }
            case CGAL::EQUAL: {
                //std::cout << "EQUAL root1: " << root1 << " root2: " << root2
                //          << std::endl;
                
                bool surely_equal = false;
                bool surely_unequal = false;
            
                bool projected_sil1 =  
                                (projected_sil1s.find(surface1) != 
                                 projected_sil1s.end());
    
                // still_possible_equal can already be *false*, i.e.,
                // if in last iterator CGAL::EQAUL was detected and
                // refinement was triggered. Then CGAL:EQUAL is still
                // a wrong overlap ..so directly refine!
                
                // General ////////////////////////////////////////////////////

                // filter: intersections can only happen for non-FACES
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    still_possible_equal = (feature != SoX::FACE);
                    //std::cout << "filterA" << std::endl;
                    if (feature == SoX::FACE) {
                        //std::cout << "filterA2" << std::endl;
                        surely_unequal = true;
                    }
                }
                
                // filter: is surface2 involved in any projected cut?
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    still_possible_equal = projected_cut2;
                    //std::cout << "filterB" << std::endl;
                    if (!projected_cut2) {
                        //std::cout << "filterB2" << std::endl;
                        surely_equal = false;
                        surely_unequal = true;
                    }
                }
                
                typename 
                    std::map< Surface_3, 
                              std::pair< int, 
                                         std::list< std::pair< int, int > >
                    >, Surface_less_than >::iterator sit;
                
                int mult = -1;

                // filter: is surface1 involved in any projected cut?
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    //std::cout << "filterC" << std::endl;
                    sit = projected_cut1s.find(surface1);
                    
                    if (sit == projected_cut1s.end()) {
                        //std::cout << "filterC2" << std::endl;
                        still_possible_equal = false;
                        surely_equal = false;
                        surely_unequal = true;
                    } else {
                        //std::cout << "filterC3" << std::endl;
                        if (sit->second.second.empty()) {
                            // no overlapping found
                            still_possible_equal = false;
                            surely_equal = false;
                            surely_unequal = true;
                        } else {
                            mult = sit->second.first;
                        }
                    }
                }
                
                // Edge ///////////////////////////////////////////////////////
                
                // filter: EDGE: at most one intersection?
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    //std::cout << "filterD" << std::endl;
                    if (feature == SoX::EDGE) {
                        //std::cout << "filterD2" << "mult=" << mult << std::endl;
                        CGAL_assertion(mult > 0);
                        if (mult == 1) {
                            //std::cout << "filterD3" << std::endl;
                            std::pair< int, int > intersection =
                                iv_helper.unique_overlap(
                                        isolator1, isolator2,
                                        sit->second.second.begin(),
                                        sit->second.second.end());
                            if (intersection.first == root1 &&
                                intersection.second == root2) {
                                //std::cout << "filterD4" << std::endl;
                                surely_equal = true;
                                surely_unequal = false;
                            } else {
                                //std::cout << "filterD5" << std::endl;
                                surely_equal = false;
                                surely_unequal = true;
                            }
                        } // else complex intersections can create 
                          // the projection -> equal_z
                    }
                }
                
                // the other case (mult > 1) must be handled by equal_z
                
                // Vertex /////////////////////////////////////////////////////
                // we now turn to vertices only
                
                // filter: VERTEX using adjacencies
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    //std::cout << "filterE" << std::endl;
                    if (feature == SoX::VERTEX) {
                        //std::cout << "filterE2" << std::endl;
                        if (!adj_x_map.empty()) {
                            //std::cout << "filterE3" << std::endl;
                            // it is possible that we can derive the 
                            // intersection by adjacencies, 
                            // so let's give it a try
                            typename Adj_x_map::const_iterator
                                iit = adj_x_map.find(
                                        std::make_pair(surface1, root1)
                                );
                            if (iit != adj_x_map.end() &&
                                iit->second == root2) {
                                //std::cout << "filterE4" << std::endl;
                                // adj-intersection detected
                                surely_equal = true;
                                surely_unequal = false;
                            }
                        } 
                    }
                }
                
                // Remark: When reching this point, two z-intervals are 
                //         overlapping and the adjacency does not found 
                //         an intersection. If no silhouette of 
                //         surface1 and surface2 exists, there will be 
                //         no further intersection
                //         The other way around: There can only be an 
                //         intersection if an isolated object meets another one
                //         of "lies within a sheet that locally looks like a
                //         plane". 

                // Filter: VERTEX further intersections require silhouette
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    if (feature == SoX::VERTEX) {
                        if (!projected_sil2) {
                            
                            if (! projected_sil1) {
                                surely_equal = false;
                                surely_unequal = true;
                            }
                        }
                    }
                }
                
                // Remark: We cannot claim that the existence of a single
                //         overlap induces a real intersection,
                //         this is only true oif the point is non-singular,
                //         i.e., mult = 1 (but we do not have this information)
                //         otherwise complex intersection can create the
                //         projected intersection. Therefore it is not possible
                //         to use the following filter
                /*
                // filter: only one overlap?
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    // find possible pairings first: unique_overlap!
                }
                */
                
                //         So, here can assume that the point is singular 
                //         and at least one silhouette exists
                // Filter:
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    if (feature == SoX::VERTEX) {
                      CGAL_assertion(projected_sil1 || projected_sil2);
                    }
                }
                
                // If *at least one of) the current events is isolated,
                // we can also do not deduce whether the overlapping
                // intervals will seperate or not. So the following
                // filter is also not possible -> need answer from equal_z
                /*
                // Filter: Surface-"sheet" at vertex isolated and
                //         NOT involved in overlaps -> no intersection
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    
                }
                */
                
                // TODO move filter before adjacenceny?
                // Filter: If vertex is not a genuine vertex in cut-curve arr
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    if (feature == SoX::VERTEX) {
                        boost::optional< 
                            std::pair< SoX::Dcel_feature, CGAL::Object > >
                            dcel_pair = 
                            this->ptr()->_m_dcel_data->_dcel(
                                    surface1, surface2
                            );
                        CGAL_assertion(dcel_pair);
                        if (dcel_pair->first == SoX::EDGE) {
                            typename Restricted_cad_3::Halfedge_const_handle 
                                hh;
                            CGAL_assertion_code(bool check = )
                                CGAL::assign(hh, dcel_pair->second);
                            CGAL_assertion(check);
                            boost::optional< int > m = 
                                hh->data()->multiplicity_of_cut(
                                        surface1, surface2
                                );
                            CGAL_assertion(m);
                            if (*m == 1) {
                                std::pair< int, int > intersection =
                                    iv_helper.unique_overlap(
                                            isolator1, isolator2,
                                            sit->second.second.begin(),
                                            sit->second.second.end());
                                if (intersection.first == root1 &&
                                    intersection.second == root2) {
                                    //std::cout << "filterD4" << std::endl;
                                    surely_equal = true;
                                    surely_unequal = false;
                                } else {
                                    //std::cout << "filterD5" << std::endl;
                                    surely_equal = false;
                                    surely_unequal = true;
                                }
                            }
                        }
                    }
                }
                

#if 0
                // empty filter
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    
                }
#endif
                
                // final check (costly one!)
                if (!surely_equal && !surely_unequal && still_possible_equal) {
                    // reaching here, no filter applied

                    // what are the possibilities:
                    // (1) vertex and at least one silhouette
                    CGAL_assertion_code(
                            bool poss1 = (
                                    feature == VERTEX && 
                                    (projected_sil1 || projected_sil2)
                            );
                            
                    );
                    // (2) EDGE with mult > 1
                    CGAL_assertion_code(
                            bool poss2 = (
                                    feature == EDGE && 
                                    (mult > 1)
                            );
                    );
                    // TODO (2a) an isolated segment? Something special?
                    CGAL_assertion_code(
                            bool poss2a = false;
                    );
                    CGAL_assertion(poss1 || poss2 || poss2a);
                    
                    // Remark: It is CRUCIAL that we can decide (somehow)
                    //         whether two roots are equal
                    typename Surface_z_at_xy_isolator_traits::Equal_z equal_z
                        // TODO (object isolator1.equal_z_object())
                        ;
                    
                    // we have to ask the traits class using equal_z
                    // remove mult as upper bound of possible intersections
                    surely_equal = equal_z(
                            surface1, isolator1, root1, projected_sil1,
                            surface2, isolator2, root2, projected_sil2);
                    surely_unequal = !surely_equal;
                }
                
                // final result                
                if (surely_equal) {
                    CGAL_assertion(!surely_unequal);
                    // reach here, it's not only possible equal -> it is equal!
                    
                    // update z-id
                    it->_make_unique_for(*this, z);
                    
                    // insert it to current cell
                    it->_add_sheet(surface2, root2, isolator2.length(root2));
                    
                    // proceed to next root and next cell!
                    it++;
                    root2++;
                    z++;

                    // the next combination can be equal again
                    still_possible_equal = true;
                    
                } else {
                    CGAL_assertion(!still_possible_equal || surely_unequal);
                    CGAL_assertion(!surely_equal);
                    
                    // otherwise, we refine interval
                    isolator1.refine_interval(root1);
                    isolator2.refine_interval(root2);

                    // the next combination is equal to this one,
                    // which is cannot be equal!
                    still_possible_equal = false;
                }
                break;
            }
            case CGAL::LARGER: {
                //std::cout << "LARGER root1: " << root1 << " root2: " << root2
                //          << std::endl;
                // create new Z_cell
                Z_cell cell(*this, z);
                cell._add_sheet(surface2, root2, isolator2.length(root2));

                // and insert it before current position
                this->ptr()->_m_z_cells.insert(it, cell);
                
                // go to next root
                root2++;
                z++;

                // the next combination can be equal again
                still_possible_equal = true;
                
                break;
            }
            }
        }
        while (it != this->ptr()->_m_z_cells.end()) {
            it->_make_unique_for(*this, z);
            
            // go to next existing cell in z_stack
            it++;
            z++;
        }
        while (root2 < roots2) {
            // append cells of isolator2 at the end
            Z_cell cell(*this,z);
            cell._add_sheet(surface2, root2, isolator2.length(root2));

            // and insert it at the end
            this->ptr()->_m_z_cells.push_back(cell);
            
            // go to next root
            root2++;
            z++;
        }
        
    }
    
    //!@}
    
private:
    //!\name Merging
    //!@{
    
    /*!\brief
     * merges \c *this and \c z_stack to a new z_stack
     */
    Self _merge(const Self& z_stack, const Dcel_data* data) {
        Self tmp(z_stack);
        tmp.copy_on_write();
        tmp.ptr()->_m_dcel_data = data;
        
        for (typename Z_at_xy_isolators::iterator it =
                 this->ptr()->_m_z_at_xy_isolators.begin();
             it != this->ptr()->_m_z_at_xy_isolators.end(); it++) {
            tmp._add_surface(it->first, it->second);
        }

        for (typename Z_at_xy_isolators::iterator it =
                 this->ptr()->_m_z_at_xy_empty_isolators.begin();
             it != this->ptr()->_m_z_at_xy_empty_isolators.end(); it++) {
            tmp._add_surface(it->first, it->second);
        }
        
        return tmp;
    }
    
    //!@}
    
public:
    //!\name Information
    //!@{

    /*\brief
     * returns one of stored projected point for z_stack
     */
    Point_2 point() const {
        return this->ptr()->_m_point;
    }

    /*!\brief
     * returns the number of z_stack cells
     */
    bool is_empty() const {
        return this->ptr()->_m_z_cells.empty();
    }
    
    /*!\brief
     * returns the number of z_stack cells
     */
    int number_of_z_cells() const {
        return static_cast< int >(this->ptr()->_m_z_cells.size());
    }

    //!@}
    
    //!\name Accessors
    //!@{

    /*!\brief
     * returns \c i-th z_stack cell
     */
    Z_cell z_cell(int i) const {
        CGAL_precondition(i >= 0 && i < this->number_of_z_cells());
        typename Z_cell_container::const_iterator it = 
            this->ptr()->_m_z_cells.begin();
        std::advance(it, i);
        return *it;
    }
    
    /*\brief
     * return level of \c surface in \c i-th z-cell
     */
    int level_of_surface_in_z_cell(const Surface_3& surface, int i) const {
        CGAL_precondition(i >= 0);
        CGAL_precondition(i < this->number_of_z_cells());
        const Z_cell& cell = this->z_cell(i);
        return cell.sheet_number(surface);
    }

    /*!\brief
     * beginning of z-cells
     */
    typename Z_cell_container::iterator z_cells_begin() {
        return this->ptr()->_m_z_cells.begin();
    }
    
    /*!\brief
     * past-the-end of z-cells
     */
    typename Z_cell_container::iterator z_cells_end() {
        return this->ptr()->_m_z_cells.end();
    }

    /*!\brief
     * beginning of z-cells
     */
    Z_cell_const_iterator z_cells_begin() const {
        return this->ptr()->_m_z_cells.begin();
    }
    
    /*!\brief
     * past-the-end of z-cells
     */
    Z_cell_const_iterator z_cells_end() const {
        return this->ptr()->_m_z_cells.end();
    }
    
    /*\brief
     * returns all levels of \c surface
     */
    template < class OutputIterator >
    OutputIterator levels_of(
            const Surface_3& surface,
            OutputIterator oi
    ) const {
        int level = 0;
        for (typename Z_cell_container::const_iterator it =
                 this->z_cells_begin(); it != this->z_cells_end(); 
             it++, level++) {
            if (it->sheet_number(surface) != -1) {
                *oi++ = level;
            }
        }
        return oi;
    }
    
    
    /*!\brief
     * returns z-level for a given surface sheet
     */
    int z_level_of_sheet(const Surface_3& surface, int sheet) const {
        CGAL_precondition(sheet < this->number_of_z_cells());
        // it suffices to start at "sheet"
        for (int i = 0; i < this->number_of_z_cells(); i++) {
            if (this->level_of_surface_in_z_cell(surface, i) == sheet) {
                return i;
            }
        }
        CGAL_error_msg("Not allowed to reach here");
        /* Should not be reached */ return -1;
    }
    
#if 0
    // TASK implement other intersections predicates
    /*\brief
     * returns all pairs of intersection \c surface1 and \c surface2
     * are involved
     */
    template < class OutputIterator >
    OutputIterator intersections(
            const Surface_3& surface1, const Surface_3& surface2,
            OutputIterator oi
    ) const {
        typedef Z_cell::value_type Surface_id;
        typedef std::pair< Surface_id, Surface_id > Surface_id_pair;

        for (typename Z_cell_container::const_iterator it =
                 z_cells_begin(); it != z_cells.end(); it++) {
            Surface_id si1 = it->find(surface1);
            Surface_id si2 = it->find(surface2);
            if (si1 != it->end() && si2 != it->end()) {
                *oi++ = std::make_pair(si1, si2);
            }
        }
        return oi;
    }
    
    /*\brief
     * returns all pairs of intersection \c surface1 and \c surface2
     * are involved
     */
    template < class OutputIterator >
    OutputIterator intersections(
            const Surface_3& surface1,
            OutputIterator oi
    ) const {
        typedef Z_cell::value_type Surface_id;
        typedef std::pair< Surface_id, Surface_id > Surface_id_pair;

        for (typename Z_cell_container::const_iterator it =
                 z_cells_begin(); it != z_cells.end(); it++) {
            Surface_id si1 = it->find(surface1);
            Surface_id si2 = it->find(surface2);
            if (si1 != it->end()) {
                for (typename Z_cell::const_iterator sit =
                         it->begin(); sit != it->end(); sit++) {
                    if (sit != si1) {
                        *oi++ = std::make_pair(si1, *sit);
                    }
                }
            }
        }
        return oi;
    }

#endif
    
    /*\brief
     * returns all levels of \c surface that is involved in some intersection
     */
    template < class OutputIterator >
    OutputIterator intersection_levels_of(
            const Surface_3& surface,
            OutputIterator oi
    ) const {
        if (this->ptr()->_m_dcel_data->_has_cut()) {
            int level = 0;
            for (typename Z_cell_container::const_iterator it =
                     this->z_cells_begin(); it != this->z_cells_end(); 
                 it++, level++) {
                if (it->number_of_sheets() > 1) {
                    int sheet = it->sheet_number(surface);
                    if (sheet != -1) {
                        *oi++ = sheet;
                    }
                }
            }
        }
        return oi;
    }
    

    //!@}
    
public:
    //!\name Adjacency
    //!@{

    /*!\brief
     * Compute adjacency information when going from isolator \c *this and 
     * sheet defined by \c surface and \c sheet to \c to. Return the final
     * sheet number at \c to.
     */
    std::pair< int, int > adjacency(const Surface_3& surface, int sheet, 
                                    SoX::Dcel_feature feat_from,
                                    CGAL::Object dcel_handle_from,
                                    Self to,
                                    SoX::Dcel_feature feat_to,
                                    CGAL::Object dcel_handle_to
    ) const {
        typename Z_at_xy_isolators::const_iterator thisit =
            this->ptr()->_m_z_at_xy_isolators.find(surface);
        CGAL_assertion(thisit != this->ptr()->_m_z_at_xy_isolators.end());
        
        typename Z_at_xy_isolators::const_iterator thatit =
            to.ptr()->_m_z_at_xy_isolators.find(surface);
        if (thatit == to.ptr()->_m_z_at_xy_isolators.end()) {
            // if surface is not involved in "that" dcel feature
            // -> either not connected or infinity
            thatit = to.ptr()->_m_z_at_xy_empty_isolators.find(surface);
            CGAL_assertion(
                    thatit != to.ptr()->_m_z_at_xy_empty_isolators.end()
            );
        }
        
        typedef typename Surface_z_at_xy_isolator_traits::Adjacency Adjacency;
        
        typedef typename SoX::Adjacencies_3::Adjacency_interval 
            Adjacency_interval;
        
        if (this->id() == to.id()) {
            return std::make_pair(sheet, sheet);
        }
        // else
        
        typedef Restricted_cad_3< Surface_z_at_xy_isolator_traits >
            Restricted_cad_3;
        
        CGAL_assertion_code(
                Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(surface);
        );
        
        bool has_vertical_line_from = false;
        if (feat_from == SoX::VERTEX) {
            typename Restricted_cad_3::Vertex_const_handle vh;
            CGAL_assertion_code(bool check = )
                CGAL::assign(vh,dcel_handle_from);
            CGAL_assertion(check);
            has_vertical_line_from = 
                vh->data()->_supports_vertical_line(surface);

            CGAL_assertion(vh->data()->_rs_id() == cad.id());
        }

        CGAL_assertion_code((
        {
            if (feat_from == SoX::EDGE) {
                typename Restricted_cad_3::Halfedge_const_handle hh1, hh2;
                CGAL_assertion_code(CGAL::assign(hh1, dcel_handle_from));
                CGAL_assertion(hh1->data()->_rs_id() == cad.id());
                CGAL::Object obj = 
                    cad.locate(thisit->second.traits().point());
                CGAL_assertion(CGAL::assign(hh2, obj));
                CGAL_assertion(hh1 == hh2);
            }
        })
        );
        CGAL_assertion_code((
        {
            if (feat_from == SoX::FACE) {
                typename Restricted_cad_3::Face_const_handle fh1, fh2;
                CGAL_assertion(CGAL::assign(fh1, dcel_handle_from));
                CGAL_assertion(fh1->data()->_rs_id() == cad.id());
                CGAL::Object obj = 
                    cad.locate(thisit->second.traits().point());
                CGAL_assertion(CGAL::assign(fh2, obj));
                CGAL_assertion(fh1 == fh2);
            }
        })
        );
        
        bool has_vertical_line_to = false;
        if (feat_to == SoX::VERTEX) {
            typename Restricted_cad_3::Vertex_const_handle vh;
            CGAL_assertion_code(bool check = )
                CGAL::assign(vh,dcel_handle_to);
            CGAL_assertion(check);
            has_vertical_line_to = 
                vh->data()->_supports_vertical_line(surface);
            
            CGAL_assertion(vh->data()->_rs_id() == cad.id());
        }

        CGAL_assertion_code((
        {
            if (feat_to == SoX::EDGE) {
                typename Restricted_cad_3::Halfedge_const_handle hh1, hh2;
                CGAL_assertion_code(CGAL::assign(hh1, dcel_handle_to));
                CGAL_assertion(hh1->data()->_rs_id() == cad.id());
                CGAL::Object obj = 
                    cad.locate(thatit->second.traits().point());
                CGAL_assertion(CGAL::assign(hh2, obj));
                CGAL_assertion(hh1 == hh2 || hh1->twin() == hh2);
            }
        })
        );
        CGAL_assertion_code((
        {
            if (feat_to == SoX::FACE) {
                typename Restricted_cad_3::Face_const_handle fh1, fh2;
                CGAL_assertion(CGAL::assign(fh1, dcel_handle_to));
                CGAL_assertion(fh1->data()->_rs_id() == cad.id());
                CGAL::Object obj = 
                    cad.locate(thatit->second.traits().point());
                CGAL_assertion(CGAL::assign(fh2, obj));
                CGAL_assertion(fh1 == fh2);
            }
        })
        );
       
        
        Adjacency compute_adj; // TODO single instance
        SoX::Adjacencies_3 adj = 
            compute_adj(surface, 
                        thisit->second, feat_from, dcel_handle_from,
                        has_vertical_line_from,
                        thatit->second, feat_to, dcel_handle_to,
                        has_vertical_line_to);
        Adjacency_interval intv = adj.interval(sheet);
        return intv;
    }

    //!@}

    //!\name IO
    //!@{
    
    /*!\brief
     * prints pretty-formated version of z_stack
     */
    void pretty_print(std::ostream& os) const {
        os << "#ZS @ ";
        os << point() << ":" << std::endl;
        os << "[" << number_of_z_cells() << ",";
        for (typename Z_cell_container::const_iterator sit =
                 z_cells_begin(); sit != z_cells_end(); sit++) {
            sit->pretty_print(os);
        }
        os << " id=" << this->id() << ">";
        os << "]";
            
    }
    
    //!@}

    // friend
    friend class SoX::P_dcel_info< Surface_z_at_xy_isolator_traits >;
    // for _isolator
    friend class 
    SoX::Intern::Z_cell< Surface_z_at_xy_isolator_traits, Dcel_data >;
    friend class SoX::Restricted_cad_3< Surface_z_at_xy_isolator_traits >;
};


/*!\relates Z_stack
 * \brief
 * outputs Z_stack object to stream \c os
 */
template < class SurfaceZAtXyIsolatorTraits, class DcelData >
std::ostream& operator<<(
        std::ostream& os, 
        const Z_stack< SurfaceZAtXyIsolatorTraits, DcelData >& z_stack
) {
    z_stack.pretty_print(os);
    return os;
}
    

} // namespace SoX

#endif // SoX_GAPS_Z_STACK_H
// EOF
