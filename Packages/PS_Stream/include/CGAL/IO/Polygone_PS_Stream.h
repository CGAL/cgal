#ifdef CGAL_POLYGON_2_H

CGAL_BEGIN_NAMESPACE

template <class _Traits, class _Container>
PS_Stream &operator<<(PS_Stream &ps, const Polygon_2<_Traits,_Container>& p)
{
  typedef Polygon_2<_Traits,_Container>::Vertex_const_iterator VI;
  
  VI i= p.vertices_begin();
  VI end= p.vertices_end();
  ps << point_style(PS_Stream::NONE);
  ps.os() <<"/poly {newpath " <<endl;
  ps << *i;
  ps.os() <<"mt" <<endl;

  do{ 
     ps << *i;
     ps.os() << "lt" << endl;
     i++;
    }while ( i != end);
 
 ps.os() <<"closepath" << endl;
 ps.os() <<"} def" <<endl;
 ps.os() <<"gsave" <<endl;
 ps.os() << " poly fill" <<endl;
 ps.os() << "0 setgray " <<endl;
 ps.os() << " poly st" <<endl;
 ps.os() << "grestore" <<endl;
return ps;
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_2_H
