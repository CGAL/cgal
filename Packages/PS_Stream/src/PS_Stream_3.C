#include <CGAL/IO/PS_Stream_3.h>

CGAL_BEGIN_NAMESPACE

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3,const Direction& d,
			 ostream& os, OutputMode mode):
 PS_Stream(os,mode),_bbox3(bb3),_dir(d)
{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_t=Transformation(yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_t,v[i]);

_bbox=PS_BBox(search_xmin(v),search_ymin(v),search_xmax(v),search_ymax(v));
 set_scale(_bbox);

}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3, const Direction& d, const
			 char* fname, OutputMode mode):
 PS_Stream(fname,mode),_bbox3(bb3),_dir(d)

{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_t=Transformation(yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_t,v[i]);

_bbox=PS_BBox(search_xmin(v),search_ymin(v),search_xmax(v),search_ymax(v));
 set_scale(_bbox);

}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3, const Direction& d, float H, const
	    char* fname, OutputMode mode):
 PS_Stream(H,fname,mode),_bbox3(bb3),_dir(d)
{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_t=Transformation(yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_t,v[i]);

_bbox=PS_BBox(search_xmin(v),search_ymin(v),search_xmax(v),search_ymax(v));
set_window(_bbox,H); 
set_scale(_bbox);
insert_catalogue(); 
}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3,const Direction& d, float H, ostream& os,
	    OutputMode mode ) :
 PS_Stream(H,os,mode),_bbox3(bb3),_dir(d)
{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_t=Transformation(yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_t,v[i]);

_bbox=PS_BBox(search_xmin(v),search_ymin(v),search_xmax(v),search_ymax(v));
set_window(_bbox,H); 
set_scale(_bbox);
insert_catalogue();
}


Point2 PS_Stream_3::transform(Transformation t,const Point3& p)
{
Point2 p2;
Point3 q=t(p);
p2=Point2(q.x(),q.y());
return p2;
}


double PS_Stream_3::search_xmin(Point3 v[])
{
double xmin=v[0].x();
for(int i=1;i<8;i++)
  {if (xmin>v[i].x())
    xmin=v[i].x();
  }
return xmin;
}


double PS_Stream_3::search_ymin(Point3 v[])
{
double ymin=v[0].y();
for(int i=1;i<8;i++)
  {if (ymin>v[i].y())
    ymin=v[i].y();
  }
return ymin;
}

double PS_Stream_3::search_xmax(Point3 v[])
{
double xmax=v[0].x();
for(int i=1;i<8;i++)
  {if (xmax<v[i].x())
    xmax=v[i].x();
  }
return xmax;
}


double PS_Stream_3::search_ymax(Point3 v[])
{
double ymax=v[0].y();
for(int i=1;i<8;i++)
  {if (ymax<v[i].y())
    ymax=v[i].y();
  }
return ymax;
}

Plane PS_Stream_3::compute_plane_equations(const Facet& f)
{
Halfedge_handle h=f.halfedge();
Facet f1;
f1.plane()=f.plane();
return (f1.plane()=Plane(h->opposite()->vertex()->point(),h->vertex()->point(),h->next()->vertex()->point()));

}

CGAL_END_NAMESPACE
