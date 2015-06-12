// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
//
// Author(s)     : Fernando de Goes, Pierre Alliez

#ifndef COST_H_
#define COST_H_

template <class FT>
class Cost
{
private:
    FT m_norm;
    FT m_tang;
    FT m_max_norm;
    FT m_max_tang;
    
public:
    Cost()
    {
        m_norm = 0.0;
        m_tang = 0.0; 
        m_max_norm = 0.0;
        m_max_tang = 0.0;
    }
    
    Cost(const FT norm, const FT tang)
    {
        m_norm = norm;
        m_tang = tang;
        m_max_norm = norm;
        m_max_tang = tang;
    }
    
    ~Cost() { }
    
    Cost& operator = (const Cost& cost)
    {
        m_norm = cost.norm();
        m_tang = cost.tang();
        m_max_norm = cost.max_norm();
        m_max_tang = cost.max_tang();
        return *this;
    }
    
    const FT norm() const { return m_norm; }
    
    const FT tang() const { return m_tang; }

    const FT max_norm() const { return m_max_norm; }

    const FT max_tang() const { return m_max_tang; }

    FT finalize(const FT alpha = 0.5) const 
    { 
        return 2.0 * (alpha * m_norm + (1.0 - alpha) * m_tang);
    }
    
    void divide(const FT ratio)
    {
        assert(ratio != 0.0);
        m_norm /= ratio;
        m_tang /= ratio;
    }
    
    void add(const Cost& cost, const FT mass = 1.0)
    {
        m_norm += mass * cost.norm();
        m_tang += mass * cost.tang();
    }
    
    void reset_max()
    {
        m_max_norm = 0.0;
        m_max_tang = 0.0;
    }

    void update_max(const Cost& cost)
    {
        m_max_norm = (std::max)(m_max_norm, cost.max_norm());
        m_max_tang = (std::max)(m_max_tang, cost.max_tang());
    }

    void compute_max(const FT norm, const FT tang)
    {
        m_max_norm = (std::max)(m_max_norm, norm);
        m_max_tang = (std::max)(m_max_tang, tang);
    }
};

#endif
