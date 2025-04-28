//
//  hdvf_duality.hpp
//  HDVF
//
//  Created by umenohana on 27/08/2024.
//

#ifndef SUB_CHAIN_COMPLEX_H
#define SUB_CHAIN_COMPLEX_H

#include <vector>
#include <cassert>
#include <iostream>
#include "CGAL/OSM/OSM.hpp"
#include "CGAL/OSM/Bitboard.hpp"

/**
 * \class SubChainComplex
 * \brief Implementation of the sub chain complex relation (for duality) - subcomplex and cocomplex
 *
 * The SubChainComplex class contains all functions to build sub chain complex (for HDVF).
 *
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 */

template <typename CoefficientType, typename ComplexType> // !!! change everywhere -> CoefficientType
class SubChainComplex
{
public:
    int _dim ;
    std::vector<OSM::Bitboard> _sub ;
    std::vector<int> _nb_cells ;
    bool _full ;
private:
    // Load a set of cells and close the enumeration down (faces)
    void down_closure (const std::vector<std::vector<int> > &cells, const ComplexType& K)
    {
        std::vector<std::vector<int> > faces(cells.size()) ;
        bool rec_needed = false ;
        for (int q = 0; q < cells.size(); ++q)
        {
            for (int i : cells.at(q))
            {
                if (!_sub.at(q).isOn(i))
                {
                    _sub.at(q).setOn(i) ;
                    ++_nb_cells.at(q) ;
                    if (q>0) // Then faces must be considered
                    {
                        // Add all its faces to faces.at(q-1)
                        OSM::Chain<CoefficientType, OSM::COLUMN> bnd(K.d(i,q)) ;
                        for (typename OSM::Chain<CoefficientType, OSM::COLUMN>::const_iterator it = bnd.cbegin(); it != bnd.cend(); ++it)
                        {
                            const int c(it->first) ;
                            if (!_sub.at(q-1).isOn(c))
                            {
                                faces.at(q-1).push_back(c) ;
                                rec_needed = true ;
                            }
                        }
                    }
                }
            }
        }
        // Call recursively
        if (rec_needed)
            down_closure(faces, K);
    }
    
    // Load a set of cells (without closure)
    void set_cells (const std::vector<std::vector<int> > &cells, const ComplexType& K)
    {
        for (int q = 0; q < cells.size(); ++q)
        {
            for (int i : cells.at(q))
            {
                    _sub.at(q).setOn(i) ;
                    ++_nb_cells.at(q) ;
            }
        }
    }
    
public:
    /** \brief Constructor of a full sub CC from a complex. */
    SubChainComplex(const ComplexType& K, bool full=true)
    {
        _dim = K.dim() ;
        _sub.resize(_dim+1) ;
        _nb_cells.resize(_dim+1) ;
        // Create Bitboards
        for (int q=0; q<=K.dim(); ++q)
        {
            if (full)
            {
                _sub.at(q) = OSM::Bitboard(K.nb_cells(q),false) ; // all cells to 1
                _full = true ;
            }
            else
                _sub.at(q) = OSM::Bitboard(K.nb_cells(q),true) ; // all cells to 0
            _nb_cells.at(q) = K.nb_cells(q) ;
        }
    }
    
    
    /** \brief Constructor from an enumeration of cellId of the sub chain complex (close the complex down if necessary). */
    SubChainComplex(const ComplexType& K, const std::vector<std::vector<int> > &cells, bool close = true)
    {
        _dim = K.dim() ;
        _sub.resize(_dim+1) ;
        _nb_cells.resize(_dim+1) ;
        _full = true ;
        // Create Bitboards
        for (int q=0; q<=K.dim(); ++q)
        {
            _sub.at(q) = OSM::Bitboard(K.nb_cells(q)) ;
        }
        
        // Set bits
        if (close)
            down_closure(cells, K) ;
        else
            set_cells(cells, K) ;
        
        for (int q=0; q<=K.dim(); ++q)
        {
            _full = (_full && (_nb_cells.at(q) == _sub.at(q).size())) ;
        }
    }
    
    /** \brief Copy constructor. */
    SubChainComplex(const SubChainComplex& _otherToCopy)
    {
        _dim = _otherToCopy._dim ;
        _nb_cells = _otherToCopy._nb_cells ;
        _full = _otherToCopy._full ;
        _sub = _otherToCopy._sub ;
    }
    
    /** \brief Complement of a sub chain complex. */
    SubChainComplex complement()
    {
        SubChainComplex cSub(*this) ;
        // Complement bitboards
        for (int q=0; q<=_dim; ++q)
        {
            cSub._sub.at(q).bit_not() ;
            cSub._nb_cells.at(q) = _sub.at(q).size() - _nb_cells.at(q) ;
            cSub._full = ~_full ;
        }
        return cSub ;
    }
    
    /** \brief Check for a bit (bit i in dimension q). */
    inline bool get_Bit(int q, int i) const
    {
        return _sub.at(q).isOn(i) ;
    }
    
    /** \brief Check for a bit (bit i in dimension q). */
    inline void set_BitOn(int q, int i)
    {
        return _sub.at(q).setOn(i) ;
    }
    
    /** \brief Check for a bit (bit i in dimension q). */
    inline void set_BitOff(int q, int i)
    {
        return _sub.at(q).setOff(i) ;
    }
    
    /** \brief Getter for the bitboards of the sub chain complex. */
    inline const std::vector<OSM::Bitboard>& get_BitBoard() const
    {
        return _sub ;
    }
    
    /** \brief Getter for the bitboard of the sub chain complex in dimension q. */
    inline const OSM::Bitboard& get_BitBoard(int q) const
    {
        return _sub.at(q) ;
    }
    
//    template <typename _CoefficientType, typename _ComplexType>
//    friend ostream & operator << (ostream & out, const SubChainComplex<_CoefficientType,_ComplexType> & sub)
//    {
//        for (int q = 0; q <= sub._dim; ++q)
//            out << "dim " << q << " : " << sub._sub.at(q) << endl ;
//        return out ;
//    }
    
    friend std::ostream & operator << (std::ostream & out, const SubChainComplex & sub)
    {
        for (int q = 0; q <= sub._dim; ++q)
            out << "dim " << q << " : " << sub._sub.at(q) << std::endl ;
        return out ;
    }
};

#endif /* SUB_CHAIN_COMPLEX_H */
