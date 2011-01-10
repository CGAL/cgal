#include "Scene_c3t3_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QPainter>

#include <map>
#include <vector>
#include <CGAL/gl.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include <CGAL_demo/Scene_item_with_display_list.h>
#include <CGAL_demo/Scene_interface.h>
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

namespace {
  void CGALglcolor(QColor c, int dv = 0)
  {
    if ( 0 != dv )
    {
// workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
      c = c.darker(dv);
#undef darker
    }
    
    ::glColor4f(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
  }
}

template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value);


enum { DRAW = 0, DRAW_EDGES = 1 };

void draw_triangle(const Kernel::Point_3& pa,
                          const Kernel::Point_3& pb,
                          const Kernel::Point_3& pc) {
  Kernel::Vector_3 n = cross_product(pb - pa, pc -pa);
  n = n / CGAL::sqrt(n*n);

  ::glNormal3d(n.x(),n.y(),n.z());

  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glVertex3d(pc.x(),pc.y(),pc.z());
}

double complex_diag(const Scene_item* item) {
  const Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax-bbox.xmin;
  const double& ydelta = bbox.ymax-bbox.ymin;
  const double& zdelta = bbox.zmax-bbox.zmin;
  const double diag = std::sqrt(xdelta*xdelta +
                                ydelta*ydelta +
                                zdelta*zdelta);
  return diag * 0.7;
}


struct Scene_c3t3_item_priv {
  Scene_c3t3_item_priv() : c3t3() {}
  Scene_c3t3_item_priv(const C3t3& c3t3_) : c3t3(c3t3_) {}

  C3t3 c3t3;
  QVector<QColor> colors;
};


Scene_c3t3_item::
Scene_c3t3_item()
  : d(new Scene_c3t3_item_priv())
  , frame(new ManipulatedFrame())
  , histogram_()
  , data_item_(NULL)
  , indices_()
{
  connect(frame, SIGNAL(modified()), this, SLOT(changed()));
  c3t3_changed();
}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3)
  : d(new Scene_c3t3_item_priv(c3t3)), frame(new ManipulatedFrame())
  , histogram_(), data_item_(NULL), indices_()
{
  connect(frame, SIGNAL(modified()), this, SLOT(changed()));
  c3t3_changed();
}

Scene_c3t3_item::~Scene_c3t3_item()
{
  delete frame;
  delete d;
}

const C3t3& 
Scene_c3t3_item::c3t3() const {
  return d->c3t3;
}

C3t3& 
Scene_c3t3_item::c3t3()
{
  return d->c3t3;
}

Kernel::Plane_3 
Scene_c3t3_item::plane() const {
  const qglviewer::Vec& pos = frame->position();
  const qglviewer::Vec& n =
    frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
}

Scene_item::Bbox 
Scene_c3t3_item::bbox() const {
  if(isEmpty())
    return Bbox();
  else {
    CGAL::Bbox_3 result = c3t3().triangulation().finite_vertices_begin()->point().bbox();
    for(Tr::Finite_vertices_iterator
          vit = ++c3t3().triangulation().finite_vertices_begin(),
          end = c3t3().triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      result = result + vit->point().bbox();
    }
    return Bbox(result.xmin(), result.ymin(), result.zmin(),
                result.xmax(), result.ymax(), result.zmax());
  }
}

QString 
Scene_c3t3_item::toolTip() const {
  int number_of_tets = 0;
  for(Tr::Finite_cells_iterator
        cit = c3t3().triangulation().finite_cells_begin(),
        end = c3t3().triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if( c3t3().is_in_complex(cit) )
      ++number_of_tets;
  }
  return tr("<p>3D complex in a 3D triangulation:<br />"
            "<b>%4</b></p>"
            "<p>Number of vertices: %1<br />"
            "Number of surface facets: %2<br />"
            "Number of volume tetrahedra: %3</p>")
    .arg(c3t3().triangulation().number_of_vertices())
    .arg(c3t3().number_of_facets_in_complex())
    .arg(number_of_tets)
    .arg(this->name());
}

void
Scene_c3t3_item::direct_draw() const {
  direct_draw(DRAW);
}

void
Scene_c3t3_item::direct_draw_edges() const {
  direct_draw(DRAW_EDGES);
}

void
Scene_c3t3_item::direct_draw(int mode) const {
  ::glPushMatrix();
  ::glMultMatrixd(frame->matrix());
  QGLViewer::drawGrid((float)complex_diag(this));
  ::glPopMatrix();

  if(isEmpty())
    return;

  CGALglcolor(QColor(0,0,0));

  //std::cerr << "Direct_draw " << mode << "\n";
  GLboolean lighting = ::glIsEnabled(GL_LIGHTING);
  GLboolean two_side;
  ::glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &two_side);
  ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  if(!lighting)
    ::glDisable(GL_LIGHTING);

  const Kernel::Plane_3& plane = this->plane();

  ::glBegin(GL_TRIANGLES);
  for(C3t3::Facets_in_complex_iterator
        fit = c3t3().facets_in_complex_begin(),
        end = c3t3().facets_in_complex_end();
      fit != end; ++fit)
  {
    const Tr::Cell_handle& cell = fit->first;
    const int& index = fit->second;
    if(cell->subdomain_index() != 0 &&
       cell->neighbor(index)->subdomain_index() != 0)
    {
      continue;
    }

    const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
    const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
    const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
    typedef Kernel::Oriented_side Side;
    using CGAL::ON_ORIENTED_BOUNDARY;
    using CGAL::ON_NEGATIVE_SIDE;
    const Side sa = plane.oriented_side(pa);
    const Side sb = plane.oriented_side(pb);
    const Side sc = plane.oriented_side(pc);
    if(sa == ON_NEGATIVE_SIDE &&
       sb == ON_NEGATIVE_SIDE && 
       sc == ON_NEGATIVE_SIDE)
    {
      if(mode != DRAW_EDGES) {
        if(cell->subdomain_index() == 0) {
          CGALglcolor(d->colors[cell->neighbor(index)->subdomain_index()]);
        }
        else {
          CGALglcolor(d->colors[cell->subdomain_index()]);
        }
      }
      draw_triangle(pa, pb, pc);
    }
  }
  ::glEnd();

  ::glBegin(GL_TRIANGLES);
  for(Tr::Finite_cells_iterator
        cit = c3t3().triangulation().finite_cells_begin(),
        end = c3t3().triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if(! c3t3().is_in_complex(cit) )
      continue;

    const Kernel::Point_3& pa = cit->vertex(0)->point();
    const Kernel::Point_3& pb = cit->vertex(1)->point();
    const Kernel::Point_3& pc = cit->vertex(2)->point();
    const Kernel::Point_3& pd = cit->vertex(3)->point();
    typedef Kernel::Oriented_side Side;
    using CGAL::ON_ORIENTED_BOUNDARY;
    const Side sa = plane.oriented_side(pa);
    const Side sb = plane.oriented_side(pb);
    const Side sc = plane.oriented_side(pc);
    const Side sd = plane.oriented_side(pd);

    if( sa == ON_ORIENTED_BOUNDARY ||
        sb == ON_ORIENTED_BOUNDARY ||
        sc == ON_ORIENTED_BOUNDARY ||
        sd == ON_ORIENTED_BOUNDARY ||
        sb != sa || sc != sa || sd != sa)
    {
      if(mode != DRAW_EDGES) {
        CGALglcolor(d->colors[cit->subdomain_index()],150);
      }
      else
      {
        CGALglcolor(d->colors[cit->subdomain_index()],250);
      }
      draw_triangle(pa, pb, pc);
      draw_triangle(pa, pb, pd);
      draw_triangle(pa, pc, pd);
      draw_triangle(pb, pc, pd);
    }
  }
  ::glEnd();
  
  if(!two_side)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  if(lighting)
    ::glEnable(GL_LIGHTING);
  else
    ::glDisable(GL_LIGHTING);
}


QPixmap
Scene_c3t3_item::graphicalToolTip() const
{
  if ( ! histogram_.isNull() )
  {
    return histogram_;
  }
  else
  {
    const_cast<Scene_c3t3_item&>(*this).build_histogram();
    return histogram_;
  }
}
 
void
Scene_c3t3_item::build_histogram()
{
  // Create an histogram_ and display it
  const int height = 140;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height-top_margin*2;
  const int width = 402;
  const int cell_width = 2;
  const int text_margin = 3;
  const int text_height = 20;
  
  histogram_ = QPixmap(width,height+text_height);
  histogram_.fill(QColor(192,192,192));
  
  QPainter painter(&histogram_);
  painter.setPen(Qt::black);
  painter.setBrush(QColor(128,128,128));
  //painter.setFont(QFont("Arial", 30));
  
  // Build histogram_ data
  double min_value, max_value;
  std::vector<int> histo_data = create_histogram(c3t3(),min_value,max_value);
  
  // Get maximum value (to normalize)
  int max_size = 0;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it )
  {
    max_size = (std::max)(max_size,*it);
  }
  
  // colored histogram
  int j = 0;
  
  // draw
  int i=left_margin;
  for ( std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end() ;
       it != end ; ++it, i+=cell_width )
  {
    int line_height = static_cast<int>( std::ceil(static_cast<double>(drawing_height) *
      static_cast<double>(*it)/static_cast<double>(max_size)) + .5);
    
    painter.fillRect(i,
                     drawing_height+top_margin-line_height,
                     cell_width,
                     line_height,
                     get_histogram_color(j++));
  }
  
  // draw bottom horizontal line
  painter.setPen(Qt::blue);
  
  painter.drawLine(QPoint(left_margin, drawing_height + top_margin),
                   QPoint(left_margin + static_cast<int>(histo_data.size())*cell_width, 
                          drawing_height + top_margin));

  
  // draw min value and max value
  const int min_tr_width = static_cast<int>( 2*(std::floor(min_value)*cell_width + left_margin) );
  const int max_tr_width = static_cast<int>( 
    2*((histo_data.size()-std::floor(max_value))*cell_width + left_margin) );
  const int tr_y = drawing_height + top_margin + text_margin;
  
  painter.setPen(get_histogram_color(min_value));
  QRect min_text_rect (0, tr_y, min_tr_width, text_height);
  painter.drawText(min_text_rect, Qt::AlignCenter, tr("%1").arg(min_value,0,'f',1));
  
  painter.setPen(get_histogram_color(max_value));           
  QRect max_text_rect (width - max_tr_width, tr_y, max_tr_width, text_height);
  painter.drawText(max_text_rect, Qt::AlignCenter, tr("%1").arg(max_value,0,'f',1));
}


template<typename C3t3>
std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value)
{
  typedef typename C3t3::Triangulation::Point Point_3;
  
	std::vector<int> histo(181,0);
  
  min_value = 180.;
  max_value = 0.;
  
	for (typename C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit)
	{
		if( !c3t3.is_in_complex(cit))
			continue;
		
		const Point_3& p0 = cit->vertex(0)->point();
		const Point_3& p1 = cit->vertex(1)->point();
		const Point_3& p2 = cit->vertex(2)->point();
		const Point_3& p3 = cit->vertex(3)->point();
		
		double a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0,p1,p2,p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p2, p1, p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p3, p1, p2)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p2, p0, p3)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p3, p0, p2)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p2, p3, p0, p1)));
		histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);
    
	}
  
	return histo;	
}


QColor
Scene_c3t3_item::get_histogram_color(const double v) const
{
  if ( v < 5 )            { return Qt::red; }
  else if ( v < 10 )      { return QColor(215,108,0); }
  else if ( v < 15 )      { return QColor(138,139,0); }
  else if ( v < 165 )     { return Qt::darkGreen; }
  else if ( v < 170 )     { return QColor(138,139,1); }
  else if ( v < 175 )     { return QColor(215,108,0); }
  else /* 175<v<=180 */   { return Qt::red; }
}


void
Scene_c3t3_item::setColor(QColor c)
{
  color_ = c;
  compute_color_map(c);
}

void
Scene_c3t3_item::c3t3_changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  indices_.clear();
  
  int max = 0;
  for(C3t3::Cells_in_complex_iterator cit = this->c3t3().cells_in_complex_begin(),
      end = this->c3t3().cells_in_complex_end() ; cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    indices_.insert(cit->subdomain_index());
  }

  d->colors.resize(max+1);
  compute_color_map(color_);
  
  // Rebuild histogram
  build_histogram();
}

void
Scene_c3t3_item::compute_color_map(const QColor& c)
{
  typedef Indices::size_type size_type;

  size_type nb_domains = indices_.size();
  size_type i = 0;
  for(Indices::iterator it = indices_.begin(), end = indices_.end();
      it != end; ++it, ++i)
  {
    double hue = c.hueF() + 1./nb_domains * i;
    if ( hue > 1 ) { hue -= 1.; }
    d->colors[*it] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
  }
}

#include "Scene_c3t3_item.moc"
