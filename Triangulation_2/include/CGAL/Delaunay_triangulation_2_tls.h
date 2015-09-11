

    CGAL_THREAD_LOCAL_DECLARE_POD(int, maxd,30); 
    if (CGAL_THREAD_LOCAL_IS_UNINITIALIZED_POD(maxd)) {
      CGAL_THREAD_LOCAL_SET_POD(int, maxd,30);
    }
    CGAL_THREAD_LOCAL_GET_POD(int,maxd);


#   ifdef CGAL_HAS_THREADS
        CGAL_THREAD_LOCAL_DECLARE(std::vector<Face_handle> , f);
        CGAL_THREAD_LOCAL_DECLARE(std::vector<int> , i);
        CGAL_THREAD_LOCAL_DECLARE(std::vector<Vertex_handle> , w);
        if (CGAL_THREAD_LOCAL_IS_UNINITIALIZED(f)) {
          CGAL_THREAD_LOCAL_SET(f, new std::vector<Face_handle>(maxd));
          CGAL_THREAD_LOCAL_SET(i, new std::vector<int>(maxd));
          CGAL_THREAD_LOCAL_SET(w, new std::vector<Vertex_handle>(maxd));
        }
        CGAL_THREAD_LOCAL_GET(std::vector<Face_handle>, f);
        CGAL_THREAD_LOCAL_GET(std::vector<int>, i);
        CGAL_THREAD_LOCAL_GET(std::vector<Vertex_handle>, w);
#   else
      static std::vector<Face_handle> f(maxd);
      static std::vector<int> i(maxd);
      static std::vector<Vertex_handle> w(maxd);
#   endif
