#ifndef CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE


template <class Gt>
class Apollonius_graph_face_base_2
  :  public Triangulation_face_base_2<Gt>
{
  //public:
private:
  typedef Gt                                  Geom_traits;
  typedef Triangulation_face_base_2<Gt>       Fbase;
  typedef Apollonius_graph_face_base_2<Gt>  Face_base;
  //  typedef typename Gt::Weighted_point   Weighted_point;
protected:
  void* next_face_in_list[3];
  int   next_indx_in_list[3];
  void* prev_face_in_list[3];
  int   prev_indx_in_list[3];

  inline
  void init()
  {
    for (int i = 0; i < 3; i++) {
      next_face_in_list[i] = NULL;
      prev_face_in_list[i] = NULL;
    }
  }
public:
  Apollonius_graph_face_base_2() : Fbase()
  { init(); }

  Apollonius_graph_face_base_2(void* v0, void* v1, void* v2)
    : Fbase(v0,v1,v2)
  { init(); }

  Apollonius_graph_face_base_2(void* v0, void* v1, void* v2,
			       void* n0, void* n1, void* n2)
    : Fbase(v0,v1,v2,n0,n1,n2)
  { init(); }

public:
  // methods for handling the in-place queue
  inline bool is_in_list(int i) const
  {
    CGAL_triangulation_assertion( i >= 0 && i <= 2 );
    return (next_face_in_list[i] != NULL ||
	    prev_face_in_list[i] != NULL);
  }

  inline void set_next(int i, const pair<void*,int>& next)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    CGAL_triangulation_precondition( next.first == NULL ||
				     (next.second >= 0 &&
				      next.second <= 2) );
    next_face_in_list[i] = next.first;
    next_indx_in_list[i] = next.second;
  }

  inline void set_previous(int i, const pair<void*,int>& prev)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    CGAL_triangulation_precondition( prev.first == NULL ||
				     (prev.second >= 0 &&
				      prev.second <= 2) );
    prev_face_in_list[i] = prev.first;
    prev_indx_in_list[i] = prev.second;
  }

#if 0
  inline void remove(int i)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    next_face_in_list[i] = NULL;
    prev_face_in_list[i] = NULL;
  }
#endif

  inline std::pair<void*, int> next(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    return std::pair<void*, int>(next_face_in_list[i],
				 next_indx_in_list[i]);
  }

  inline std::pair<void*, int> previous(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 2 );
    return std::pair<void*, int>(prev_face_in_list[i],
				 prev_indx_in_list[i]);
  }
};

CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_FACE_BASE_2_H
