void
any_button(CGAL::Window_stream &W)
{
    double x, y;
    std::cerr << "Press any button to continue" << std::endl;
    W.read_mouse(x,y);
}

template < class TRIANGULATION > 
Vertex_handle  closest_vertex(const TRIANGULATION &T,
			      Face_handle f,
			      const Point& p)
{
  Vertex_handle v ;
  typename Gt::Compare_distance_2 cmp =
    T.geom_traits().compare_distance_2_object();

  if( T.is_infinite(f)){
    int i = f->index(T.infinite_vertex());
    Point pcwi = f->vertex(f->cw(i))->point();
    Point pccwi = f->vertex(f->ccw(i))->point();
    v =  cmp(p, pcwi, pccwi) == CGAL::SMALLER ? f->vertex(f->cw(i)) :
                                                f->vertex(f->ccw(i));
  }
  else{ 
    v = f->vertex(0);
    if (cmp(p, f->vertex(1)->point(), v->point()) == CGAL::SMALLER) 
      v = f->vertex(1);
    if (cmp(p, f->vertex(2)->point(), v->point()) == CGAL::SMALLER) 
      v = f->vertex(2);
  }
  return v;
}



template <class TRIANGULATION>
void file_input(TRIANGULATION  &T,
		Window_stream& W,
		const Options& opt)
{
    if(! opt.file_input){
        return;
    }

    std::ifstream is(opt.fname);
    CGAL::set_ascii_mode(is);

    int n, count = 0;
    is >> n;
    std::cerr << "Reading " << n << " points" << std::endl;

    if( (! opt.check) && (! opt.draw)){
	std::cerr << "Reading from iterator" << std::endl;
        std::istream_iterator<Point> begin(is);
        std::istream_iterator<Point> end;
        T.insert(begin, end);
    }else{
	std::cerr << "Reading from point to point" << std::endl;
        Point mp;

        for(; n > 0; n--){
            is >> mp;
            T.insert(mp);
            if(opt.check){
		std::cerr << "Checking validity" << std::endl;
                T.is_valid();
            }
            if(opt.draw){
                W.clear();
                W << CGAL::BLUE << T << CGAL::RED;
            }

            if(++count == 100){
                std::cerr << ".";
                count = 0;
            }
        }
    }
    std::cerr << "Done with file input" << std::endl;
}

template <class TRIANGULATION>
void container_input(TRIANGULATION &T,
		     Window_stream &W)
{
    std::list<Point> L;
    L.push_front(Point(0,0));
    L.push_front(Point(1,0));
    L.push_front(Point(1,1));

    int n = T.insert(L.begin(), L.end());
    std::cerr << n << " points inserted from a list." << std::endl;

    std::vector<Point> V(3);
    V[0] = Point(0, 0);
    V[1] = Point(0.4, 0.4);
    V[2] = Point(0.3, 0.3);

    n = T.insert(V.begin(), V.end());
    std::cerr << n << " points inserted from a vector." << std::endl;

    W.clear();
    W << T;
}

template <class TRIANGULATION>
void draw_face(Face_handle fh, TRIANGULATION &T, Window_stream &W) 
{
  if(! T.is_infinite( fh))  W << T.triangle( fh );
}

template <class TRIANGULATION>
void draw_edge(Edge e, TRIANGULATION &T,Window_stream &W) 
{
  if(! T.is_infinite( e ))  W << T.segment( e );
}
