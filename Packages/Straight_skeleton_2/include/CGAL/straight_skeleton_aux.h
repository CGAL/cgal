// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/straight_skeleton_aux.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_STRAIGHT_SKELETON_AUX_H
#define CGAL_STRAIGHT_SKELETON_AUX_H 1

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW

namespace SS_IO_AUX
{
  class ScopedDrawing
  {
    public :
    
      virtual ~ScopedDrawing() 
      {
        if ( mID != -1 ) 
          Straight_skeleton_external_undraw_object(mID) ; 
      }
      
      void Release() { mID = -1 ; }

    protected :
    
      ScopedDrawing ( int aID ) : mID(aID) {}
            
    private :
    
      int mID ;
  } ;
  
  class ScopedPointDrawing : public ScopedDrawing
  {
    public :
    
    template<class Point_2>
    ScopedPointDrawing( Point_2 const& aP, int aColor, char const* aLayer )
      :
      ScopedDrawing
      (
        Straight_skeleton_external_draw_point(  to_double( aP.x() )
                                               ,to_double( aP.y() )
                                               ,aColor
                                               ,aLayer
                                             )
      )                                         
    {}
  } ;
  
  class ScopedSegmentDrawing : public ScopedDrawing
  {
    public :
    
    template<class Segment_2>
    ScopedSegmentDrawing( Segment_2 const& aSeg, int aColor, char const* aLayer )
      :
      ScopedDrawing
      (
        Straight_skeleton_external_draw_segment(  to_double( aSeg.source().x() )
                                                 ,to_double( aSeg.source().y() )
                                                 ,to_double( aSeg.target().x() )
                                                 ,to_double( aSeg.target().y() )
                                                 ,aColor
                                                 ,aLayer
                                               )
      )                                         
    {}
  } ;
     
}
#endif

class Ref_counted_base
{
private:

  mutable long mCount ;

  Ref_counted_base( Ref_counted_base const &);
  Ref_counted_base& operator=( Ref_counted_base const &);

protected:

  Ref_counted_base(): mCount(0) {}

  virtual ~Ref_counted_base() {}

public:

#if !defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)

    inline friend void intrusive_ptr_add_ref(Ref_counted_base const* p) { ++p->mCount; }
    
    inline friend void intrusive_ptr_release(Ref_counted_base const* p)
      {
        if ( --p->mCount == 0 )
          delete p;
      }

#else

    void AddRef() const { ++mCount; }

    void Release() const
      {
        if( --mCount == 0 )
          delete this;
      }

#endif
};

CGAL_END_NAMESPACE

#if defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)

namespace boost
{

inline void intrusive_ptr_add_ref( CGAL_NS Ref_counted_base const* p ) { p->AddRef(); }
inline void intrusive_ptr_release( CGAL_NS Ref_counted_base const* p ) { p->Release(); }

}

} // namespace boost

#endif

#endif // CGAL_STRAIGHT_SKELETON_AUX_H //
// EOF //
 
