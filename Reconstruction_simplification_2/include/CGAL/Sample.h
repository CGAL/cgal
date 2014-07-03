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

#ifndef SAMPLE_H
#define SAMPLE_H


/// \cond SKIP_IN_MANUAL

template <class Kernel>
class Sample
{
public:
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point;
    
private:
    Point m_point;
    FT m_mass;
    
    FT m_dist2_to_edge;
    FT m_coordinate;
    
    FT m_backup_dist2;
    FT m_backup_coord;
    
public:
    Sample(const Point& point, 
			const FT mass = 1.0)
    {
        m_mass  = mass;
        m_point = point;
        
        m_dist2_to_edge = 0.0;
        m_coordinate    = 0.0;
        
        m_backup_dist2 = 0.0;
        m_backup_coord = 0.0;
    }

    Sample(const Sample& sample)
    {
        m_mass  = sample.mass();
        m_point = sample.point();
        
        m_dist2_to_edge = 0.0;
        m_coordinate    = 0.0;
        
        m_backup_dist2 = 0.0;
        m_backup_coord = 0.0;
    }
    
    ~Sample() { }
    
    const Point& point() const { return m_point; } 
    Point& point() { return m_point; } 
    
    const FT& mass() const { return m_mass; }
    FT& mass() { return m_mass; }
    
    const FT& distance2() const { return m_dist2_to_edge; }
    FT& distance2() { return m_dist2_to_edge; }
    
    const FT& coordinate() const { return m_coordinate; }
    FT& coordinate() { return m_coordinate; }
    
    void backup()
    {
        m_backup_dist2 = m_dist2_to_edge;
        m_backup_coord = m_coordinate;
    }
    
    void restore()
    {
        m_dist2_to_edge = m_backup_dist2;
        m_coordinate = m_backup_coord;
    }
};

template <class	Sample>
class Sample_with_priority
{
public:
    typedef typename Sample::FT FT;
    
private:
    Sample* m_sample;
    FT m_priority;
    
public:
    Sample_with_priority(Sample* sample, const FT priority = 0.0)
    {
        m_sample   = sample;
        m_priority = priority;
    }
    
    Sample_with_priority(const Sample_with_priority& psample) 
    {
        m_sample   = psample.sample();
        m_priority = psample.priority();
    }
    
    ~Sample_with_priority() { }
    
    Sample_with_priority& operator = (const Sample_with_priority& psample)
    {
        m_sample   = psample.sample();
        m_priority = psample.priority();
        return *this;
    }
    
    Sample* sample() const { return m_sample; }
    
    const FT priority() const { return m_priority; }
};


//TODO: IV: What is the correct place for this struct?
template <class T>
struct greater_priority
{
    bool operator() (const T& a, const T& b) const
    {
        return ( a.priority() > b.priority() );
    }
};

#endif // SAMPLE_H
