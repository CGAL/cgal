
void myhandler(int x, int y, int but, CGAL::GL_win *W) {
  std::vector<double> v(3);
  v[0]=x*2; v[1]=y*2; v[2]=0;
  W->add_point_to_object(1,1,v);

}
