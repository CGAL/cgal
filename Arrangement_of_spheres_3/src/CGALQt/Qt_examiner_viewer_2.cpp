#include "Qt_examiner_viewer_window_2.moc"
#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/IO/Qt_widget_circular_arc_2.h>

#define CGAL_QTEV_LOCK CGAL_precondition(mutex_.locked());


/*std::auto_ptr<QMutexLocker> lockptr;					\
  if (!mutex_.locked()) lockptr= std::auto_ptr<QMutexLocker>(new QMutexLocker(&mutex_));*/

CGAL_BEGIN_NAMESPACE
Qt_examiner_viewer_2 *qt_debug_examiner_viewer_2__=NULL;

/*struct Qt_examiner_viewer_2::QTEV_layer::Conditional_mutex {
  Conditional_mutex(bool skip_lock, QMutex *qm) {
    if (!skip_lock) {
      qm->lock();
      qm_=qm;
    } else {
      qm_=NULL;
    }
  }
  ~Conditional_mutex() {
    if (qm_ != NULL) {
      qm_->unlock();
    }
  }
  QMutex *qm_;
  };*/


/*void Qt_examiner_viewer_2::clear_layer() {
  layers_[cur_layer_]->clear();
  }*/

template <class V, class VC>
CGAL::Color Qt_examiner_viewer_2::QTEV_layer::draw(CGAL::Color cc, const V &v, const VC &vc) {
  for (unsigned int i=0; i< v.size(); ++i){
    if (vc[i] != cc) {
      *P::widget << vc[i];
      cc= vc[i];
    }
    *P::widget << v[i];
  }
  return cc;
}

template <class V, class VC, class VS>
CGAL::Color Qt_examiner_viewer_2::QTEV_layer::draw(CGAL::Color cc, const V &v, const VC &vc, const VS &vs) {
  typename VS::value_type sc;
  if (!vs.empty()) sc= vs[0];
  *P::widget << sc;
  for (unsigned int i=0; i< v.size(); ++i){
    if (vc[i] != cc) {
      *P::widget << vc[i];
      cc= vc[i];
    }
    if (vs[i] != sc){
      *P::widget << vs[i];
      sc= vs[i];
    }
    *P::widget << v[i];
  }
  return cc;
}


/*void Qt_examiner_viewer_2::QTEV_layer::set_is_editing(bool tf) {
  if (tf && !is_editing_) {
    mutex_.lock();
  } else if (!tf && is_editing_) {
    mutex_.unlock();
  }
  is_editing_=tf;
  }*/

Qt_examiner_viewer_2::QTEV_layer::QTEV_layer(NT scale): bbox_(std::numeric_limits<double>::max(),
							      std::numeric_limits<double>::max(),
							      -std::numeric_limits<double>::max(),
							      -std::numeric_limits<double>::max()),
							cur_color_(CGAL::BLACK),
							scale_(scale){
  ubb_=true;
  cur_point_style_= CGAL::DISC;
  
  circles_.reserve(1000);
  circle_colors_.reserve(1000);
  points_.reserve(1000);
  point_colors_.reserve(1000);
  point_styles_.reserve(1000);
  segments_.reserve(1000);
  segment_colors_.reserve(1000);
  lines_.reserve(1000);
  line_colors_.reserve(1000);
  arcs_.reserve(1000);
  arc_colors_.reserve(1000);
  labels_.reserve(1000);
  label_colors_.reserve(1000);
}

void Qt_examiner_viewer_2::QTEV_layer::clear() {
  CGAL_QTEV_LOCK;
  circles_.clear();
  circle_colors_.clear();
  points_.clear();
  point_colors_.clear();
  point_styles_.clear();
  segments_.clear();
  segment_colors_.clear();
  lines_.clear();
  line_colors_.clear();
  arcs_.clear();
  arc_colors_.clear();
  labels_.clear();
  label_colors_.clear();
}



CGAL::Bbox_2 Qt_examiner_viewer_2::QTEV_layer::bounding_box(){
  return bbox_;
}


void Qt_examiner_viewer_2::QTEV_layer::draw(){
  QMutexLocker lock(&mutex_);
  //::cout << "Drawing layer " << std::endl;
  CGAL::Qt_widget *w= P::widget;
  w->lock();
  CGAL::Color cc= cur_color_;
  *w << cur_color_;
  cc=draw(cc, lines_, line_colors_);
  cc=draw(cc, segments_, segment_colors_);
  // draw them last since they have problems. This ordering covers errors. 
  cc=draw(cc, circles_, circle_colors_);
  cc=draw(cc, arcs_, arc_colors_);
  cc=draw(cc, points_, point_colors_, point_styles_);
  {
    CGAL::Color cc;
    if (!label_colors_.empty()) cc= label_colors_[0];
    *w << cc;
    QFont f= w->get_painter().font();
    f.setPixelSize(9);
    w->get_painter().setFont(f);
    for (unsigned int i=0; i< labels_.size(); ++i){
      if (label_colors_[i] != cc) {
	*w << label_colors_[i];
	cc= label_colors_[i];
      }
      int ind= labels_[i].first;
      w->get_painter().drawText(w->x_pixel(CGAL::to_double(points_[ind].x()))+3,
				w->y_pixel(CGAL::to_double(points_[ind].y()))-3,
				QString(labels_[i].second));
    }
      
  }
  w->unlock();
}

void Qt_examiner_viewer_2::QTEV_layer::new_circle(const Circle &c) {
  CGAL_QTEV_LOCK;
  circles_.push_back(Circle(rescale(c.center()),
			    scale_*scale_*NT(1.0)*c.squared_radius()));
  circle_colors_.push_back(cur_color_);
  if (ubb_) bbox_=bbox_+ circles_.back().bbox();
}
void Qt_examiner_viewer_2::QTEV_layer::new_point(const Point &c) {
  CGAL_QTEV_LOCK;
  points_.push_back(rescale(c));
  point_colors_.push_back(cur_color_);
  point_styles_.push_back(cur_point_style_);
  if (ubb_) bbox_=bbox_+ points_.back().bbox();
}
void Qt_examiner_viewer_2::QTEV_layer::new_segment(const Segment &c) {
  CGAL_QTEV_LOCK;
  segments_.push_back(Segment(rescale(c.source()),
					    rescale(c.target())));
  segment_colors_.push_back(cur_color_);
  if (ubb_) bbox_=bbox_+ segments_.back().bbox();
}

void Qt_examiner_viewer_2::QTEV_layer::new_circular_arc(const Circle &c,
							const Point &s,
							const Point &t) {
  CGAL_QTEV_LOCK;
  Point rcc= rescale(c.center());
  Point rs= rescale(s);
  Point rt= rescale(t);
  NT sr= scale_*scale_*c.squared_radius();
  typedef Circular_k::Circular_arc_2 CA;
  typedef Circular_k::Circle_2 CC;
  typedef Circular_k::Point_2 CP;
  typedef Circular_k::Circular_arc_point_2 CAP;
  arcs_.push_back(CA(CC(CP(rcc.x(), rcc.y()), sr),
		     CAP(CP(rs.x(), rs.y())),
		     CAP(CP(rt.x(), rt.y()))));
  arc_colors_.push_back(cur_color_);
  if (ubb_) bbox_=bbox_+ arcs_.back().bbox();
}
    
void Qt_examiner_viewer_2::QTEV_layer::new_line(const Line &c) {
  CGAL_QTEV_LOCK;
  lines_.push_back(Line(rescale(c.point()),
				      c.to_vector()));
  line_colors_.push_back(cur_color_);
  //bbox_=bbox_+ lines_.back().bbox();
}

void Qt_examiner_viewer_2::QTEV_layer::new_label(const std::string str) {
  CGAL_QTEV_LOCK;
  CGAL_precondition(!points_.empty());
  typedef std::pair<int, std::string> SP;
  labels_.push_back(SP(points_.size()-1, str));
  label_colors_.push_back(cur_color_);
}
    
void Qt_examiner_viewer_2::QTEV_layer::set_line_color(CGAL::Color c){
  cur_color_=c;
}
void Qt_examiner_viewer_2::QTEV_layer::set_point_style(CGAL::PointStyle p){
  cur_point_style_=p;
}
  
bool Qt_examiner_viewer_2::QTEV_layer::updating_box() const {
  return ubb_;
}
void Qt_examiner_viewer_2::QTEV_layer::set_updating_box(bool v) {
  ubb_=v;
}

Qt_examiner_viewer_2::Point Qt_examiner_viewer_2::QTEV_layer::rescale(const Point &p) const {
  return Point(scale_*p.x(), scale_*p.y());
}

void  Qt_examiner_viewer_2::set_layer(unsigned int li) {
  std::vector<QTEV_layer *> new_layers;
  {
    QMutexLocker lock(&mutex_);
    while (li >= layers_.size()){
      layers_.push_back(new QTEV_layer(scale_));
      new_layers.push_back(layers_.back());
    }
    cur_layer_=li;
    //layers_[cur_layer_]->clear();
  } 
  // not really safe, things could go away here
  while (!new_layers.empty()) {
    window_->widget()->attach(new_layers.back());
    new_layers.pop_back();
  }
}

void Qt_examiner_viewer_2::set_is_dirty(bool tf) {
  window_->set_is_dirty(tf);
}

void Qt_examiner_viewer_2::clear() {
  {
    QMutexLocker lock(&mutex_);;
    for(unsigned int i=0; i< layers_.size(); ++i){
      {
	//QMutexLocker ilock(layers_[i]->mutex_);
	window_->widget()->detach(layers_[i]);
      }
      delete layers_[i];
    }
    layers_.clear();
  }
  set_layer(0);
}
void Qt_examiner_viewer_2::show_everything() {
  CGAL::Bbox_2 ncb;
  {
    QMutexLocker lock(&mutex_);
    CGAL::Bbox_2 cb=layers_[0]->bounding_box();
    for (unsigned int i=1; i< layers_.size(); ++i){
      cb=cb+ layers_[i]->bounding_box();
    }
    
    double width= cb.xmax()-cb.xmin();
    double height= cb.ymax()-cb.ymin();
    if (cb.xmin() >= cb.xmax() || cb.ymin() >= cb.ymax()) {
      std::cerr << "Nothing to see here folks..." << std::endl;
      ncb= CGAL::Bbox_2(-100,-100,100,100);
    } else {
      double gf=.05;
      ncb= CGAL::Bbox_2(cb.xmin()- gf*width,
			cb.ymin()- gf*height,
			cb.xmax()+ gf*width,
			cb.ymax()+ gf*height);
    }
  }
  //window_->widget()->set_window(ncb.xmin(), ncb.xmax(), ncb.ymin(), ncb.ymax());
  window_->change_view(ncb);
  //window_->widget()->redraw();
  
  // not safe
  //CGAL_assertion(! mutex_.locked());
}

void Qt_examiner_viewer_2::set_viewport(const Bbox_2& bb) {
  window_->change_view(bb);
}

void Qt_examiner_viewer_2::show() {
  qApp->postEvent(window_, new P::Show_event());
  //P::widget()->redraw();
  //window_->set_is_dirty(true);
}

Qt_examiner_viewer_2::~Qt_examiner_viewer_2(){
  for(unsigned int i=0; i< layers_.size(); ++i){
    window_->widget()->detach(layers_[i]);
    delete layers_[i];
  }
}
CGAL_END_NAMESPACE
