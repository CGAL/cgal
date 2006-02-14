
#ifndef CGAL_POLYAP_TRAITS
#define CGAL_POLYAP_TRAITS

#include <CGAL/basic.h>

#include <CGAL/Point_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>


CGAL_BEGIN_NAMESPACE

struct Max_tag{};
struct Sum_tag{};
struct Line_tag{};
struct Segment_tag{};






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//				P O L Y G O N A L   A P P R O X I M A T I O N   T R A I T S   C L A S S E S 

//////////////////////////////////////////////////////
//
//  Squared Euclidean Distance Error

template<class KT,class DistCumul,class DistM,int>
class Squared_euclidean_error;


//  2D MAX Squared Euclidean error measured to the line support
	
template<class KT>
class Squared_euclidean_error<KT,Max_tag,Line_tag,2>
{		
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_y,dx,dy_n,max;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)		// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
					max=tr_y;
			}

			return max;
		}
		
		FT d,v;

		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=v+cr->x()*dy_n+cr->y()*dx;
			tr_y=tr_y*tr_y/d;
						
			if(tr_y>max)
				max=tr_y;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_y,dx,dy_n,max;
		InputIterator cr,cr_split_pt;
		InputIterator end;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)		// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
				{
					max=tr_y;
					cr_split_pt=cr;
				}
			}
			
			split_pt=cr_split_pt;

			return max;
		}
		
		FT d,v;

		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=v+cr->x()*dy_n+cr->y()*dx;
			tr_y=tr_y*tr_y/d;
						
			if(tr_y>max)
			{
				max=tr_y;
				cr_split_pt=cr;
			}
		}

		split_pt=cr_split_pt;

		return max;
	}
};


//  2D MAX Squared euclidean error measured to the segment

template<class KT>
class Squared_euclidean_error<KT,Max_tag,Segment_tag,2>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy_n,max,tr_y;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
					max=tr_y;
			}

			return max;
		}
		
		FT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
				tr_y=dx1*dx1+dy1*dy1;
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
				{
					tr_y=v+cr->x()*dy_n+cr->y()*dx;
					tr_y=tr_y*tr_y/d;
				}
				
			if(tr_y>max)
				max=tr_y;
		}

		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy_n,max,tr_y;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
				{
					max=tr_y;
					cr_split_pt=cr;
				}
			}

			split_pt=cr_split_pt;

			return max;
		}
		
		FT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
				tr_y=dx1*dx1+dy1*dy1;
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
				{
					tr_y=v+cr->x()*dy_n+cr->y()*dx;
					tr_y=tr_y*tr_y/d;
				}
				
			if(tr_y>max)
			{
				max=tr_y;
				cr_split_pt=cr;
			}
		}
		
		split_pt=cr_split_pt;

		return max;
	}
};


//  2D SUM Euclidean error measured to the line support

template<class KT>
class Squared_euclidean_error<KT,Sum_tag,Line_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;
		InputIterator end;

		cr_b=begin;
		cr_b++;
		
		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					dy_n=(begin->y()-cr_b->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					dy_n=(end->y()-cr_e->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=(end->x()-cr_e->x());
			dy_n=(end->y()-cr_e->y());
			tr_y=dx*dx+dy_n*dy_n;
			
			sum_e+=tr_y;

			return sum_b+sum_e;
		}
		
		FT d,v;
		
		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=v+cr_b->x()*dy_n+cr_b->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
		tr_y=tr_y*tr_y/d;
		sum_e+=tr_y;
	
		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;
		InputIterator end;

		cr_b=begin;
		cr_b++;
		
		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					dy_n=(begin->y()-cr_b->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					dy_n=(end->y()-cr_e->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=(end->x()-cr_e->x());
			dy_n=(end->y()-cr_e->y());
			tr_y=dx*dx+dy_n*dy_n;
			
			sum_e+=tr_y;
				
			split_pt=cr_b;

			return sum_b+sum_e;
		}
		
		FT d,v;
		
		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=v+cr_b->x()*dy_n+cr_b->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
		tr_y=tr_y*tr_y/d;
		sum_e+=tr_y;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  2D SUM Euclidean error measured to the segment

template<class KT>
class Squared_euclidean_error<KT,Sum_tag,Segment_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					dy_n=(begin->y()-cr_b->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_b+=tr_y;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					dy_n=(end->y()-cr_e->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_e+=tr_y;
					cr_e--;
				}
			}
			
			dx=(end->x()-cr_e->x());
			dy_n=(end->y()-cr_e->y());
			tr_y=dx*dx+dy_n*dy_n;
			
			sum_e+=tr_y;

			return sum_b+sum_e;
		}

		FT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
						tr_y=dx1*dx1+dy1*dy1;
					}
					else
					{
						tr_y=v+cr_b->x()*dy_n+cr_b->y()*dx;
						tr_y=tr_y*tr_y/d;
					}
				
				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
						tr_y=dx1*dx1+dy1*dy1;
					}
					else
					{
						tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
						tr_y=tr_y*tr_y/d;
					}
				
				sum_e+=tr_y;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
			tr_y=dx1*dx1+dy1*dy1;
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
				tr_y=dx1*dx1+dy1*dy1;
			}
			else
			{
				tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
				tr_y=tr_y*tr_y/d;
			}
			
		sum_e+=tr_y;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					dy_n=(begin->y()-cr_b->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_b+=tr_y;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					dy_n=(end->y()-cr_e->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_e+=tr_y;
					cr_e--;
				}
			}
			
			dx=(end->x()-cr_e->x());
			dy_n=(end->y()-cr_e->y());
			tr_y=dx*dx+dy_n*dy_n;
			
			sum_e+=tr_y;
			
			split_pt=cr_b;

			return sum_b+sum_e;
		}

		FT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
						tr_y=dx1*dx1+dy1*dy1;
					}
					else
					{
						tr_y=v+cr_b->x()*dy_n+cr_b->y()*dx;
						tr_y=tr_y*tr_y/d;
					}
				
				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
					tr_y=dx1*dx1+dy1*dy1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
						tr_y=dx1*dx1+dy1*dy1;
					}
					else
					{
						tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
						tr_y=tr_y*tr_y/d;
					}
				
				sum_e+=tr_y;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
			tr_y=dx1*dx1+dy1*dy1;
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
				tr_y=dx1*dx1+dy1*dy1;
			}
			else
			{
				tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
				tr_y=tr_y*tr_y/d;
			}
			
		sum_e+=tr_y;

		split_pt=cr_e;

		return sum_b+sum_e;
	}
};


//  3D MAX Euclidean error measured to the line support

template<class KT>
class Squared_euclidean_error<KT,Max_tag,Line_tag,3>
{		
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;

		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				dy=(begin->y()-cr->y());
				dz=(begin->z()-cr->z());

				tr_y=dx*dx+dy*dy+dz*dz;
			
				if(tr_y>max)
					max=tr_y;
			}
			return max;
		}
		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;
		d=d*d;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz;
			tr_y=vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz;
			tr_z=vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz;

			tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
						
			if(tr_y>max)
				max=tr_y;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;

		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				dy=(begin->y()-cr->y());
				dz=(begin->z()-cr->z());

				tr_y=dx*dx+dy*dy+dz*dz;
			
				if(tr_y>max)
				{
					max=tr_y;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;
		d=d*d;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz;
			tr_y=vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz;
			tr_z=vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz;

			tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
						
			if(tr_y>max)
			{
				max=tr_y;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D MAX Euclidean error measured to the segment

template<class KT>
class Squared_euclidean_error<KT,Max_tag,Segment_tag,3>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx1,dy1,dz1,dx2,dy2,dz2,vb,ve;
		InputIterator cr,end;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				dy=begin->y()-cr->y();
				dz=begin->z()-cr->z();

				tr_y=dx*dx+dy*dy+dz*dz;
			
				if(tr_y>max)
					max=tr_y;
			}

			return max;
		}

		max=FT(0);

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;
		d=d*d;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
				dz1=cr->z()-begin->z();
				tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
					dz1=cr->z()-end->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
				{
					tr_x=vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz;
					tr_y=vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz;
					tr_z=vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz;

					tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
				}
				
			if(tr_y>max)
				max=tr_y;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx1,dy1,dz1,dx2,dy2,dz2,vb,ve;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				dy=begin->y()-cr->y();
				dz=begin->z()-cr->z();

				tr_y=dx*dx+dy*dy+dz*dz;
			
				if(tr_y>max)
				{
					max=tr_y;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}

		max=FT(0);
		cr_split_pt=begin;

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;
		d=d*d;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
				dz1=cr->z()-begin->z();
				tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
					dz1=cr->z()-end->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
				{
					tr_x=vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz;
					tr_y=vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz;
					tr_z=vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz;

					tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
				}
				
			if(tr_y>max)
			{
				max=tr_y;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D SUM Euclidean error measured to the line support

template<class KT>
class Squared_euclidean_error<KT,Sum_tag,Line_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					dy=begin->y()-cr_b->y();
					dz=begin->z()-cr_b->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					dy=end->y()-cr_e->y();
					dz=end->z()-cr_e->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			dy=end->y()-cr_e->y();
			dz=end->z()-cr_e->z();

			tr_y=dx*dx+dy*dy+dz*dz;
			
			sum_e+=tr_y;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;
		d=d*d;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz;
				tr_y=vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz;
				tr_z=vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;

				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
				tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
				tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
		tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
		tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

		tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;

		sum_e+=tr_y;
	
		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					dy=begin->y()-cr_b->y();
					dz=begin->z()-cr_b->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					dy=end->y()-cr_e->y();
					dz=end->z()-cr_e->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			dy=end->y()-cr_e->y();
			dz=end->z()-cr_e->z();

			tr_y=dx*dx+dy*dy+dz*dz;
			
			sum_e+=tr_y;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;
		d=d*d;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz;
				tr_y=vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz;
				tr_z=vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;

				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
				tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
				tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
		tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
		tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

		tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;

		sum_e+=tr_y;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  3D SUM Euclidean error measured to the segment

template<class KT>
class Squared_euclidean_error<KT,Sum_tag,Segment_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd,dx1,dy1,dz1;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					dy=begin->y()-cr_b->y();
					dz=begin->z()-cr_b->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					dy=end->y()-cr_e->y();
					dz=end->z()-cr_e->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			dy=end->y()-cr_e->y();
			dz=end->z()-cr_e->z();

			tr_y=dx*dx+dy*dy+dz*dz;
			
			sum_e+=tr_y;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;
		d=d*d;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
					dz1=cr_b->z()-begin->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
						dz1=cr_b->z()-end->z();
						tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
					}
					else
					{ 
						tr_x=vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz;
						tr_y=vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz;
						tr_z=vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz;

						tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
					}

				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
					dz1=cr_e->z()-begin->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
						dz1=cr_e->z()-end->z();
						tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
					}
					else
					{
						tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
						tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
						tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

						tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
					}
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
			dz1=cr_e->z()-begin->z();
			tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
				dz1=cr_e->z()-end->z();
				tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
			}
			else
			{
				tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
				tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
				tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
			}

		sum_e+=tr_y;
	
		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd,dx1,dy1,dz1;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					dy=begin->y()-cr_b->y();
					dz=begin->z()-cr_b->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_b+=tr_y;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					dy=end->y()-cr_e->y();
					dz=end->z()-cr_e->z();

					tr_y=dx*dx+dy*dy+dz*dz;
			
					sum_e+=tr_y;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			dy=end->y()-cr_e->y();
			dz=end->z()-cr_e->z();

			tr_y=dx*dx+dy*dy+dz*dz;
			
			sum_e+=tr_y;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;
		d=d*d;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
					dz1=cr_b->z()-begin->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
						dz1=cr_b->z()-end->z();
						tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
					}
					else
					{ 
						tr_x=vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz;
						tr_y=vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz;
						tr_z=vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz;

						tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
					}

				sum_b+=tr_y;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
					dz1=cr_e->z()-begin->z();
					tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
						dz1=cr_e->z()-end->z();
						tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
					}
					else
					{
						tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
						tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
						tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

						tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
					}
						
				sum_e+=tr_y;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
			dz1=cr_e->z()-begin->z();
			tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
				dz1=cr_e->z()-end->z();
				tr_y=dx1*dx1+dy1*dy1+dz1*dz1;
			}
			else
			{
				tr_x=vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz;
				tr_y=vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz;
				tr_z=vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz;

				tr_y=(tr_x*tr_x+tr_y*tr_y+tr_z*tr_z)/d;
			}

		sum_e+=tr_y;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//////////////////////////////////////////////////////////
//
// L-infinity Distance Error

template<class KT,class DistCumul,class DistM,int>
class L_infinity_error;


//  2D MAX L-infinity error measured to the line support

template<class KT>
class L_infinity_error<KT,Max_tag,Line_tag,2>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT max,dx1,dy1;
		InputIterator end,cr;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
					max=dx1;
			}
			return max;
		}

		FT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
								
			dx1=dy_n*tr_y;
			if(dx1<FT(0))
				dx1=-dx1;
	
			dy1=dx*tr_y;
			if(dy1<FT(0))
				dy1=-dy1;
	
			if(dy1<dx1)
				dy1=dx1;
								
			if(dy1>max)
				max=dy1;	
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT max,dx1,dy1;
		InputIterator end,cr,cr_split_pt;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
				{
					max=dx1;
					cr_split_pt=cr;
				}
			}
				
			split_pt=cr_split_pt;

			return max;
		}

		FT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
								
			dx1=dy_n*tr_y;
			if(dx1<FT(0))
				dx1=-dx1;
	
			dy1=dx*tr_y;
			if(dy1<FT(0))
				dy1=-dy1;
	
			if(dy1<dx1)
				dy1=dx1;
								
			if(dy1>max)
			{
				max=dy1;	
				cr_split_pt=cr;
			}
		}
			
		split_pt=cr_split_pt;

		return max;
	}
};


//  2D MAX L-infinity error measured to the segment

template<class KT>
class L_infinity_error<KT,Max_tag,Segment_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT max,dx1,dy1;
		InputIterator end,cr;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;
 
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
					max=dx1;
			}
			return max;
		}

		FT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
				}
				else
				{
					tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
					dx1=dy_n*tr_y;
					dy1=dx*tr_y;
				}
				
			if(dx1<FT(0))
				dx1=-dx1;
			if(dy1<FT(0))
				dy1=-dy1;

			if(dy1<dx1)
				dy1=dx1;

			if(dy1>max)
				max=dy1;
		}

		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT max,dx1,dy1;
		InputIterator end,cr,cr_split_pt;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;
 
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
				{
					max=dx1;
					cr_split_pt=cr;
				}
			}
				
			split_pt=cr_split_pt;

			return max;
		}

		FT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
				}
				else
				{
					tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
					dx1=dy_n*tr_y;
					dy1=dx*tr_y;
				}
				
			if(dx1<FT(0))
				dx1=-dx1;
			if(dy1<FT(0))
				dy1=-dy1;

			if(dy1<dx1)
				dy1=dx1;

			if(dy1>max)
			{
				max=dy1;
				cr_split_pt=cr;
			}
		}
			
		split_pt=cr_split_pt;

		return max;
	}
};


//  2D SUM L-infinity error measured to the line support

template<class KT>
class L_infinity_error<KT,Sum_tag,Line_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
 	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;
			
			return sum_b+sum_e;
		}

		FT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;
		
				if(dy1<dx1)
					sum_b+=dx1;
				else
					sum_b+=dy1;

				cr_b++;
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;
		
				if(dy1<dx1)
					sum_e+=dx1;
				else
					sum_e+=dy1;

				cr_e--;
			}
		}

		tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
		dx1=dy_n*tr_y;
		if(dx1<FT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<FT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;
			
			split_pt=cr_b;

			return sum_b+sum_e;
		}

		FT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;
		
				if(dy1<dx1)
					sum_b+=dx1;
				else
					sum_b+=dy1;

				cr_b++;
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;
		
				if(dy1<dx1)
					sum_e+=dx1;
				else
					sum_e+=dy1;

				cr_e--;
			}
		}

		tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
		dx1=dy_n*tr_y;
		if(dx1<FT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<FT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  2D SUM L-infinity error measured to the segment

template<class KT>
class L_infinity_error<KT,Sum_tag,Segment_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;

			return sum_b+sum_e;
		}
				
		FT d,tr_y,v,dx1,dy1,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				
				if(dy1<dx1)
					sum_b+=dx1;
				else
					sum_b+=dy1;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				
				if(dy1<dx1)
					sum_e+=dx1;
				else
					sum_e+=dy1;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
				dx1=dy_n*tr_y;
				dy1=dx*tr_y;
			}

		if(dx1<FT(0))
			dx1=-dx1;
		if(dy1<FT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}
				
		FT d,tr_y,v,dx1,dy1,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				
				if(dy1<dx1)
					sum_b+=dx1;
				else
					sum_b+=dy1;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				
				if(dy1<dx1)
					sum_e+=dx1;
				else
					sum_e+=dy1;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
				dx1=dy_n*tr_y;
				dy1=dx*tr_y;
			}

		if(dx1<FT(0))
			dx1=-dx1;
		if(dy1<FT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  3D MAX L-infinity error measured to the line support

template<class KT>
class L_infinity_error<KT,Max_tag,Line_tag,3>
{		
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)	// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				if(dx<0)
					dx=-dx;

				dy=(begin->y()-cr->y());
				if(dy<0)
					dy=-dy;

				dz=(begin->z()-cr->z());
				if(dz<0)
					dz=-dz;
				if(dx<dy)
					dx=dy;
				if(dx<dz)
					dx=dz;

				if(dx>max)
					max=dx;
			}
			return max;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
			if(tr_x<0)
				tr_x=-tr_x;

			tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
			if(tr_y<0)
				tr_y=-tr_y;

			tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
			if(tr_z<0)
				tr_z=-tr_z;

			if(tr_x<tr_y)
				tr_x=tr_y;
			if(tr_x<tr_z)
				tr_x=tr_z;

			if(tr_x>max)
				max=tr_x;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)	// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				if(dx<0)
					dx=-dx;

				dy=(begin->y()-cr->y());
				if(dy<0)
					dy=-dy;

				dz=(begin->z()-cr->z());
				if(dz<0)
					dz=-dz;
				if(dx<dy)
					dx=dy;
				if(dx<dz)
					dx=dz;

				if(dx>max)
				{
					max=dx;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
			if(tr_x<0)
				tr_x=-tr_x;

			tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
			if(tr_y<0)
				tr_y=-tr_y;

			tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
			if(tr_z<0)
				tr_z=-tr_z;

			if(tr_x<tr_y)
				tr_x=tr_y;
			if(tr_x<tr_z)
				tr_x=tr_z;

			if(tr_x>max)
			{
				max=tr_x;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D MAX L-infinity error measured to the segment

template<class KT>
class L_infinity_error<KT,Max_tag,Segment_tag,3>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx2,dy2,dz2,vb,ve;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				if(dx<0)
					dx=-dx;

				dy=begin->y()-cr->y();
				if(dy<0)
					dy=-dy;

				dz=begin->z()-cr->z();
				if(dz<0)
					dz=-dz;

				if(dx<dy)
					dx=dy;
				if(dx<dz)
					dx=dz;

				if(dx>max)
					max=dx;
			}
			return max;
		}

		max=FT(0);

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				tr_x=cr->x()-begin->x();
				tr_y=cr->y()-begin->y();
				tr_z=cr->z()-begin->z();
			}
			else
				if(dd>ve)
				{
					tr_x=cr->x()-end->x();
					tr_y=cr->y()-end->y();
					tr_z=cr->z()-end->z();
				}
				else
				{
					tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
					tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
					tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
				}

			if(tr_x<0)
				tr_x=-tr_x;
			if(tr_y<0)
				tr_y=-tr_y;
			if(tr_z<0)
				tr_z=-tr_z;

			if(tr_x<tr_y)
				tr_x=tr_y;
			if(tr_x<tr_z)
				tr_x=tr_z;

			if(tr_x>max)
				max=tr_x;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx2,dy2,dz2,vb,ve;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				if(dx<0)
					dx=-dx;

				dy=begin->y()-cr->y();
				if(dy<0)
					dy=-dy;

				dz=begin->z()-cr->z();
				if(dz<0)
					dz=-dz;

				if(dx<dy)
					dx=dy;
				if(dx<dz)
					dx=dz;

				if(dx>max)
				{
					max=dx;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}

		max=FT(0);
		cr_split_pt=begin;

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				tr_x=cr->x()-begin->x();
				tr_y=cr->y()-begin->y();
				tr_z=cr->z()-begin->z();
			}
			else
				if(dd>ve)
				{
					tr_x=cr->x()-end->x();
					tr_y=cr->y()-end->y();
					tr_z=cr->z()-end->z();
				}
				else
				{
					tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
					tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
					tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
				}

			if(tr_x<0)
				tr_x=-tr_x;
			if(tr_y<0)
				tr_y=-tr_y;
			if(tr_z<0)
				tr_z=-tr_z;

			if(tr_x<tr_y)
				tr_x=tr_y;
			if(tr_x<tr_z)
				tr_x=tr_z;

			if(tr_x>max)
			{
				max=tr_x;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D SUM L-infinity error measured to the line support

template<class KT>
class L_infinity_error<KT,Sum_tag,Line_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;
	
		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;
					
					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;

					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			if(dx<dy)
				dx=dy;
			if(dx<dz)
				dx=dz;

			sum_e+=dx;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;
						
				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
		if(tr_x<0)
			tr_x=-tr_x;
		tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
		if(tr_y<0)
			tr_y=-tr_y;
		tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
		if(tr_z<0)
			tr_z=-tr_z;

		if(tr_x<tr_y)
			tr_x=tr_y;
		if(tr_x<tr_z)
			tr_x=tr_z;

		sum_e+=tr_x;
	
		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;
	
		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;
					
					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;

					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			if(dx<dy)
				dx=dy;
			if(dx<dz)
				dx=dz;

			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;
						
				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
		if(tr_x<0)
			tr_x=-tr_x;
		tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
		if(tr_y<0)
			tr_y=-tr_y;
		tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
		if(tr_z<0)
			tr_z=-tr_z;

		if(tr_x<tr_y)
			tr_x=tr_y;
		if(tr_x<tr_z)
			tr_x=tr_z;

		sum_e+=tr_x;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  3D SUM L-infinity error measured to the segment

template<class KT>
class L_infinity_error<KT,Sum_tag,Segment_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			if(dx<dy)
				dx=dy;
			if(dx<dz)
				dx=dz;
			
			sum_e+=dx;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_b->x()-begin->x();
					tr_y=cr_b->y()-begin->y();
					tr_z=cr_b->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_b->x()-end->x();
						tr_y=cr_b->y()-end->y();
						tr_z=cr_b->z()-end->z();
					}
					else
					{ 
						tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
						tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
						tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
					}

				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_e->x()-begin->x();
					tr_y=cr_e->y()-begin->y();
					tr_z=cr_e->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_e->x()-end->x();
						tr_y=cr_e->y()-end->y();
						tr_z=cr_e->z()-end->z();
					}
					else
					{
						tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
						tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
						tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
					}
						
				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			tr_x=cr_e->x()-begin->x();
			tr_y=cr_e->y()-begin->y();
			tr_z=cr_e->z()-begin->z();
		}
		else
			if(dd>ve)
			{
				tr_x=cr_e->x()-end->x();
				tr_y=cr_e->y()-end->y();
				tr_z=cr_e->z()-end->z();
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
			}

		if(tr_x<0)
			tr_x=-tr_x;
		if(tr_y<0)
			tr_y=-tr_y;
		if(tr_z<0)
			tr_z=-tr_z;
			
		if(tr_x<tr_y)
			tr_x=tr_y;
		if(tr_x<tr_z)
			tr_x=tr_z;

		sum_e+=tr_x;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)			// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					if(dx<dy)
						dx=dy;
					if(dx<dz)
						dx=dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			if(dx<dy)
				dx=dy;
			if(dx<dz)
				dx=dz;
			
			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_b->x()-begin->x();
					tr_y=cr_b->y()-begin->y();
					tr_z=cr_b->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_b->x()-end->x();
						tr_y=cr_b->y()-end->y();
						tr_z=cr_b->z()-end->z();
					}
					else
					{ 
						tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
						tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
						tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
					}

				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_e->x()-begin->x();
					tr_y=cr_e->y()-begin->y();
					tr_z=cr_e->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_e->x()-end->x();
						tr_y=cr_e->y()-end->y();
						tr_z=cr_e->z()-end->z();
					}
					else
					{
						tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
						tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
						tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
					}
						
				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				if(tr_x<tr_y)
					tr_x=tr_y;
				if(tr_x<tr_z)
					tr_x=tr_z;

				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			tr_x=cr_e->x()-begin->x();
			tr_y=cr_e->y()-begin->y();
			tr_z=cr_e->z()-begin->z();
		}
		else
			if(dd>ve)
			{
				tr_x=cr_e->x()-end->x();
				tr_y=cr_e->y()-end->y();
				tr_z=cr_e->z()-end->z();
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
			}

		if(tr_x<0)
			tr_x=-tr_x;
		if(tr_y<0)
			tr_y=-tr_y;
		if(tr_z<0)
			tr_z=-tr_z;
			
		if(tr_x<tr_y)
			tr_x=tr_y;
		if(tr_x<tr_z)
			tr_x=tr_z;

		sum_e+=tr_x;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


/////////////////////////////////////////////////////////
//
//  Manhattan Distance Error

template<class KT,class DistCumul,class DistM,int>
class Manhattan_error;


//  2D MAX Manhattan error measured to the line support

template<class KT>
class Manhattan_error<KT,Max_tag,Line_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT max,dx1,dy1;
		InputIterator end,cr;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);
		
		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
					max=dx1;
			}
			return max;
		}
		
		FT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
									
			dx1=dy_n*tr_y;
			if(dx1<FT(0))
				dx1=-dx1;
		
			dy1=dx*tr_y;
			if(dy1<FT(0))
				dy1=-dy1;
		
			dy1+=dx1;
									
			if(dy1>max)
				max=dy1;
		}

		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT max,dx1,dy1;
		InputIterator end,cr,cr_split_pt;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}
		
		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
				{
					max=dx1;
					cr_split_pt=cr;
				}
			}
				
			split_pt=cr_split_pt;

			return max;
		}
		
		FT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
									
			dx1=dy_n*tr_y;
			if(dx1<FT(0))
				dx1=-dx1;
		
			dy1=dx*tr_y;
			if(dy1<FT(0))
				dy1=-dy1;
		
			dy1+=dx1;
									
			if(dy1>max)
			{
				max=dy1;
				cr_split_pt=cr;
			}
		}

		split_pt=cr_split_pt;

		return max;
	}
};


//  2D MAX Manhattan error measured to the segment

template<class KT>
class Manhattan_error<KT,Max_tag,Segment_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT max,dx1,dy1;
		InputIterator end,cr;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
					max=dx1;
			}

			return max;
		}

		FT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
				}
				else
				{
					tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
					dx1=dy_n*tr_y;
					dy1=dx*tr_y;
				}

			if(dx1<FT(0))
				dx1=-dx1;
			if(dy1<FT(0))
				dy1=-dy1;
				
			dy1+=dx1;

			if(dy1>max)
				max=dy1;
		}

		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT max,dx1,dy1;
		InputIterator end,cr,cr_split_pt;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<FT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<FT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
				{
					max=dx1;
					cr_split_pt=cr;
				}
			}
				
			split_pt=cr_split_pt;

			return max;
		}

		FT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=FT(0);
		cr_split_pt=begin;

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()-dy_n*cr->y();

			if(dd<vb)
			{
				dx1=cr->x()-begin->x();
				dy1=cr->y()-begin->y();
			}
			else
				if(dd>ve)
				{
					dx1=cr->x()-end->x();
					dy1=cr->y()-end->y();
				}
				else
				{
					tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
					dx1=dy_n*tr_y;
					dy1=dx*tr_y;
				}

			if(dx1<FT(0))
				dx1=-dx1;
			if(dy1<FT(0))
				dy1=-dy1;
				
			dy1+=dx1;

			if(dy1>max)
			{
				max=dy1;
				cr_split_pt=cr;
			}
		}
			
		split_pt=cr_split_pt;

		return max;
	}
};


//  2D SUM Manhattan error measured to the line support

template<class KT>
class Manhattan_error<KT,Sum_tag,Line_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			return sum_b+sum_e;
		}

		FT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;

				sum_b+=dy1+dx1;

				cr_b++;
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;

				sum_e+=dy1+dx1;

				cr_e--;
			}
		}
		tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
					
		dx1=dy_n*tr_y;
		if(dx1<FT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<FT(0))
			dy1=-dy1;

		sum_e+=dy1+dx1;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			split_pt=cr_b;

			return sum_b+sum_e;
		}

		FT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;

				sum_b+=dy1+dx1;

				cr_b++;
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<FT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<FT(0))
					dy1=-dy1;

				sum_e+=dy1+dx1;

				cr_e--;
			}
		}
		tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
					
		dx1=dy_n*tr_y;
		if(dx1<FT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<FT(0))
			dy1=-dy1;

		sum_e+=dy1+dx1;

		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  2D SUM Manhattan error measured to the segment

template<class KT>
class Manhattan_error<KT,Sum_tag,Segment_tag,2>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;

					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;

					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			return sum_b+sum_e;
		}
		
		FT d,tr_y,v,dx1,dy1,dd,vb,ve;
		
		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				sum_b+=dx1+dy1;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				sum_e+=dx1+dy1;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
				dx1=dy_n*tr_y;
				dy1=dx*tr_y;
			}

		if(dx1<FT(0))
			dx1=-dx1;
		if(dy1<FT(0))
			dy1=-dy1;
		sum_e+=dx1+dy1;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy_n,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;

					dy_n=(begin->y()-cr_b->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;

					dy_n=(end->y()-cr_e->y());
					if(dy_n<FT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<FT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			split_pt=cr_b;

			return sum_b+sum_e;
		}
		
		FT d,tr_y,v,dx1,dy1,dd,vb,ve;
		
		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()-dy_n*cr_b->y();

				if(dd<vb)
				{
					dx1=cr_b->x()-begin->x();
					dy1=cr_b->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_b->x()-end->x();
						dy1=cr_b->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				sum_b+=dx1+dy1;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()-dy_n*cr_e->y();

				if(dd<vb)
				{
					dx1=cr_e->x()-begin->x();
					dy1=cr_e->y()-begin->y();
				}
				else
					if(dd>ve)
					{
						dx1=cr_e->x()-end->x();
						dy1=cr_e->y()-end->y();
					}
					else
					{
						tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
						dx1=dy_n*tr_y;
						dy1=dx*tr_y;
					}

				if(dx1<FT(0))
					dx1=-dx1;
				if(dy1<FT(0))
					dy1=-dy1;
				sum_e+=dx1+dy1;

				cr_e--;
			}
		}
		
		dd=dx*cr_e->x()-dy_n*cr_e->y();

		if(dd<vb)
		{
			dx1=cr_e->x()-begin->x();
			dy1=cr_e->y()-begin->y();
		}
		else
			if(dd>ve)
			{
				dx1=cr_e->x()-end->x();
				dy1=cr_e->y()-end->y();
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
				dx1=dy_n*tr_y;
				dy1=dx*tr_y;
			}

		if(dx1<FT(0))
			dx1=-dx1;
		if(dy1<FT(0))
			dy1=-dy1;
		sum_e+=dx1+dy1;

		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  3D MAX Manhattan error measured to the line support

template<class KT>
class Manhattan_error<KT,Max_tag,Line_tag,3>
{		
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)	// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				if(dx<0)
					dx=-dx;

				dy=(begin->y()-cr->y());
				if(dy<0)
					dy=-dy;

				dz=(begin->z()-cr->z());
				if(dz<0)
					dz=-dz;
					
				dx=dx+dy+dz;

				if(dx>max)
					max=dx;
			}
			return max;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
			if(tr_x<0)
				tr_x=-tr_x;

			tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
			if(tr_y<0)
				tr_y=-tr_y;

			tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
			if(tr_z<0)
				tr_z=-tr_z;

			tr_x=tr_x+tr_y+tr_z;

			if(tr_x>max)
				max=tr_x;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,max,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;
		
		if(begin==end || cr==end)	// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)				// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x()); 
				if(dx<0)
					dx=-dx;

				dy=(begin->y()-cr->y());
				if(dy<0)
					dy=-dy;

				dz=(begin->z()-cr->z());
				if(dz<0)
					dz=-dz;
					
				dx=dx+dy+dz;

				if(dx>max)
				{
					max=dx;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
			if(tr_x<0)
				tr_x=-tr_x;

			tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
			if(tr_y<0)
				tr_y=-tr_y;

			tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
			if(tr_z<0)
				tr_z=-tr_z;

			tr_x=tr_x+tr_y+tr_z;

			if(tr_x>max)
			{
				max=tr_x;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D MAX Manhattan error measured to the segment

template<class KT>
class Manhattan_error<KT,Max_tag,Segment_tag,3>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx2,dy2,dz2,vb,ve;
		InputIterator end,cr;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				if(dx<0)
					dx=-dx;

				dy=begin->y()-cr->y();
				if(dy<0)
					dy=-dy;

				dz=begin->z()-cr->z();
				if(dz<0)
					dz=-dz;

				dx=dx+dy+dz;

				if(dx>max)
					max=dx;
			}
			return max;
		}
		
		max=FT(0);

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				tr_x=cr->x()-begin->x();
				tr_y=cr->y()-begin->y();
				tr_z=cr->z()-begin->z();
			}
			else
				if(dd>ve)
				{
					tr_x=cr->x()-end->x();
					tr_y=cr->y()-end->y();
					tr_z=cr->z()-end->z();
				}
				else
				{
					tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
					tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
					tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
				}

			if(tr_x<0)
				tr_x=-tr_x;
			if(tr_y<0)
				tr_y=-tr_y;
			if(tr_z<0)
				tr_z=-tr_z;

			tr_x=tr_x+tr_y+tr_z;

			if(tr_x>max)
				max=tr_x;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,dy,dz,max,tr_x,tr_y,tr_z,p_xy,p_xz,p_yz,s_xy,s_xz,s_yz;
		FT d,vx,vy,vz,dd,dx2,dy2,dz2,vb,ve;
		InputIterator end,cr,cr_split_pt;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=begin->x()-cr->x();
				if(dx<0)
					dx=-dx;

				dy=begin->y()-cr->y();
				if(dy<0)
					dy=-dy;

				dz=begin->z()-cr->z();
				if(dz<0)
					dz=-dz;

				dx=dx+dy+dz;

				if(dx>max)
				{
					max=dx;
					cr_split_pt=cr;
				}
			}
			split_pt=cr_split_pt;

			return max;
		}
		
		max=FT(0);
		cr_split_pt=begin;

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			dd=dx*cr->x()+dy*cr->y()+cr->z()*dz;

			if(dd<vb)
			{
				tr_x=cr->x()-begin->x();
				tr_y=cr->y()-begin->y();
				tr_z=cr->z()-begin->z();
			}
			else
				if(dd>ve)
				{
					tr_x=cr->x()-end->x();
					tr_y=cr->y()-end->y();
					tr_z=cr->z()-end->z();
				}
				else
				{
					tr_x=(vx+cr->x()*s_yz-cr->y()*p_xy-cr->z()*p_xz)/d;
					tr_y=(vy+cr->y()*s_xz-cr->x()*p_xy-cr->z()*p_yz)/d;
					tr_z=(vz+cr->z()*s_xy-cr->x()*p_xz-cr->y()*p_yz)/d;
				}

			if(tr_x<0)
				tr_x=-tr_x;
			if(tr_y<0)
				tr_y=-tr_y;
			if(tr_z<0)
				tr_z=-tr_z;

			tr_x=tr_x+tr_y+tr_z;

			if(tr_x>max)
			{
				max=tr_x;
				cr_split_pt=cr;
			}
		}
		split_pt=cr_split_pt;

		return max;
	}
};


//  3D SUM Manhattan error measured to the line support

template<class KT>
class Manhattan_error<KT,Sum_tag,Line_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;
					
					dx=dx+dy+dz;

					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			dx=dx+dy+dz;

			sum_e+=dx;

			return sum_b+sum_e;
		}
		
		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);
		
		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				tr_x=tr_x+tr_y+tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				tr_x=tr_x+tr_y+tr_z;
						
				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
		if(tr_x<0)
			tr_x=-tr_x;
		tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
		if(tr_y<0)
			tr_y=-tr_y;
		tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
		if(tr_z<0)
			tr_z=-tr_z;

		tr_x=tr_x+tr_y+tr_z;

		sum_e+=tr_x;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;
					
					dx=dx+dy+dz;

					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			dx=dx+dy+dz;

			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}
		
		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();
		
		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx=dx*dx;
		dy=dy*dy;
		dz=dz*dz;

		d=dx+dy+dz;

		s_xz=dx+dz;
		s_xy=dx+dy;
		s_yz=dy+dz;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);
		
		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				tr_x=tr_x+tr_y+tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				if(tr_x<0)
					tr_x=-tr_x;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				if(tr_y<0)
					tr_y=-tr_y;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
				if(tr_z<0)
					tr_z=-tr_z;

				tr_x=tr_x+tr_y+tr_z;
						
				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
		if(tr_x<0)
			tr_x=-tr_x;
		tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
		if(tr_y<0)
			tr_y=-tr_y;
		tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
		if(tr_z<0)
			tr_z=-tr_z;

		tr_x=tr_x+tr_y+tr_z;

		sum_e+=tr_x;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


//  3D SUM Manhattan error measured to the segment

template<class KT>
class Manhattan_error<KT,Sum_tag,Segment_tag,3>
{												
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)		// Error=0, curve segment contains 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			dx=dx+dy+dz;
			
			sum_e+=dx;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_b->x()-begin->x();
					tr_y=cr_b->y()-begin->y();
					tr_z=cr_b->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_b->x()-end->x();
						tr_y=cr_b->y()-end->y();
						tr_z=cr_b->z()-end->z();
					}
					else
					{ 
						tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
						tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
						tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
					}

				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				tr_x=tr_x+tr_y+tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_e->x()-begin->x();
					tr_y=cr_e->y()-begin->y();
					tr_z=cr_e->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_e->x()-end->x();
						tr_y=cr_e->y()-end->y();
						tr_z=cr_e->z()-end->z();
					}
					else
					{
						tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
						tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
						tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
					}
						
				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				tr_x=tr_x+tr_y+tr_z;

				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			tr_x=cr_e->x()-begin->x();
			tr_y=cr_e->y()-begin->y();
			tr_z=cr_e->z()-begin->z();
		}
		else
			if(dd>ve)
			{
				tr_x=cr_e->x()-end->x();
				tr_y=cr_e->y()-end->y();
				tr_z=cr_e->z()-end->z();
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
			}

		if(tr_x<0)
			tr_x=-tr_x;
		if(tr_y<0)
			tr_y=-tr_y;
		if(tr_z<0)
			tr_z=-tr_z;
			
		tr_x=tr_x+tr_y+tr_z;

		sum_e+=tr_x;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT tr_x,tr_y,tr_z,dx,dy,dz,sum_b,sum_e,d;
		FT p_xy,p_xz,p_yz,s_xy,s_xz,s_yz,vx,vy,vz;
		FT vb,ve,dx2,dy2,dz2,dd;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)		// Error=0, curve segment contains 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=begin->x()-cr_b->x();
					if(dx<0)
						dx=-dx;
					dy=begin->y()-cr_b->y();
					if(dy<0)
						dy=-dy;
					dz=begin->z()-cr_b->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_b+=dx;
					
					cr_b++;
				}
				else
				{
					dx=end->x()-cr_e->x();
					if(dx<0)
						dx=-dx;
					dy=end->y()-cr_e->y();
					if(dy<0)
						dy=-dy;
					dz=end->z()-cr_e->z();
					if(dz<0)
						dz=-dz;

					dx=dx+dy+dz;
			
					sum_e+=dx;
					
					cr_e--;
				}
			}
			
			dx=end->x()-cr_e->x();
			if(dx<0)
				dx=-dx;
			dy=end->y()-cr_e->y();
			if(dy<0)
				dy=-dy;
			dz=end->z()-cr_e->z();
			if(dz<0)
				dz=-dz;

			dx=dx+dy+dz;
			
			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();
		dz=end->z()-begin->z();

		vb=begin->x()*dx+begin->y()*dy+begin->z()*dz;
		ve=end->x()*dx+end->y()*dy+end->z()*dz;

		p_xy=dx*dy;
		p_xz=dx*dz;
		p_yz=dy*dz;

		dx2=dx*dx;
		dy2=dy*dy;
		dz2=dz*dz;

		d=dx2+dy2+dz2;

		s_xz=dx2+dz2;
		s_xy=dx2+dy2;
		s_yz=dy2+dz2;

		vx=-begin->x()*s_yz+begin->y()*p_xy+begin->z()*p_xz;
		vy=-begin->y()*s_xz+begin->x()*p_xy+begin->z()*p_yz;
		vz=-begin->z()*s_xy+begin->x()*p_xz+begin->y()*p_yz;

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				dd=dx*cr_b->x()+dy*cr_b->y()+cr_b->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_b->x()-begin->x();
					tr_y=cr_b->y()-begin->y();
					tr_z=cr_b->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_b->x()-end->x();
						tr_y=cr_b->y()-end->y();
						tr_z=cr_b->z()-end->z();
					}
					else
					{ 
						tr_x=(vx+cr_b->x()*s_yz-cr_b->y()*p_xy-cr_b->z()*p_xz)/d;
						tr_y=(vy+cr_b->y()*s_xz-cr_b->x()*p_xy-cr_b->z()*p_yz)/d;
						tr_z=(vz+cr_b->z()*s_xy-cr_b->x()*p_xz-cr_b->y()*p_yz)/d;
					}

				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				tr_x=tr_x+tr_y+tr_z;

				sum_b+=tr_x;

				cr_b++;
			}
			else
			{
				dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

				if(dd<vb)
				{
					tr_x=cr_e->x()-begin->x();
					tr_y=cr_e->y()-begin->y();
					tr_z=cr_e->z()-begin->z();
				}
				else
					if(dd>ve)
					{
						tr_x=cr_e->x()-end->x();
						tr_y=cr_e->y()-end->y();
						tr_z=cr_e->z()-end->z();
					}
					else
					{
						tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
						tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
						tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
					}
						
				if(tr_x<0)
					tr_x=-tr_x;
				if(tr_y<0)
					tr_y=-tr_y;
				if(tr_z<0)
					tr_z=-tr_z;
				
				tr_x=tr_x+tr_y+tr_z;

				sum_e+=tr_x;
	
				cr_e--;
			}
		}
			
		dd=dx*cr_e->x()+dy*cr_e->y()+cr_e->z()*dz;

		if(dd<vb)
		{
			tr_x=cr_e->x()-begin->x();
			tr_y=cr_e->y()-begin->y();
			tr_z=cr_e->z()-begin->z();
		}
		else
			if(dd>ve)
			{
				tr_x=cr_e->x()-end->x();
				tr_y=cr_e->y()-end->y();
				tr_z=cr_e->z()-end->z();
			}
			else
			{
				tr_x=(vx+cr_e->x()*s_yz-cr_e->y()*p_xy-cr_e->z()*p_xz)/d;
				tr_y=(vy+cr_e->y()*s_xz-cr_e->x()*p_xy-cr_e->z()*p_yz)/d;
				tr_z=(vz+cr_e->z()*s_xy-cr_e->x()*p_xz-cr_e->y()*p_yz)/d;
			}

		if(tr_x<0)
			tr_x=-tr_x;
		if(tr_y<0)
			tr_y=-tr_y;
		if(tr_z<0)
			tr_z=-tr_z;
			
		tr_x=tr_x+tr_y+tr_z;

		sum_e+=tr_x;
	
		split_pt=cr_b;

		return sum_b+sum_e;
	}
};


///////////////////////////////
//
//  Vertical Distance Error

template<class KT,class DistCumul>
class Vertical_error;


//  2D MAX Vertical error

template<class KT>
class Vertical_error<KT,Max_tag>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT max,v;
		InputIterator end,cr;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);

			for(cr=begin,cr++;cr!=end;cr++)
			{
				v=(begin->y()-cr->y());
				if(v<FT(0))
					v=-v;

				if(v>max)
					max=v;
			}

			return max;
		}
		
		FT dx,a,b;
		
		dx=end->x()-begin->x();
		if(dx==FT(0))
			return FT(typename KT::RT(INT_MAX));

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		max=FT(0);

		for(cr=begin,cr++;cr!=end;cr++)
		{
			v=cr->y()-a*cr->x()-b;
			if(v<FT(0))
				v=-v;
			
			if(max<v)
				max=v;
		}
		return max;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT max,v;
		InputIterator end,cr,cr_split_pt;
		
		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=FT(0);
			cr_split_pt=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				v=(begin->y()-cr->y());
				if(v<FT(0))
					v=-v;

				if(v>max)
				{
					max=v;
					cr_split_pt=cr;
				}
			}
			
			split_pt=cr_split_pt;

			return max;
		}
		
		FT dx,a,b;
		
		dx=end->x()-begin->x();
		if(dx==FT(0))
		{
			split_pt=cr;
			return FT(typename KT::RT(INT_MAX));
		}

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		max=FT(0);
		cr_split_pt=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			v=cr->y()-a*cr->x()-b;
			if(v<FT(0))
				v=-v;
			
			if(max<v)
			{
				max=v;
				cr_split_pt=cr;
			}
		}

		split_pt=cr_split_pt;
		return max;
	}
};


//  2D SUM Vertical error

template<class KT>
class Vertical_error<KT,Sum_tag>
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT dx,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;

			sum_e+=dx;

			return sum_b+sum_e;
		}

		FT a,b,v;
		
		dx=end->x()-begin->x();
		 
		if(dx==FT(0))
		{
			cr_b=begin;
			cr_b++;

			return FT(typename KT::RT(INT_MAX));
		}

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				v=cr_b->y()-a*cr_b->x()-b;
				if(v<FT(0))
					sum_b-=v;
				else
					sum_b+=v;

				cr_b++;
			}
			else
			{
				v=cr_e->y()-a*cr_e->x()-b;
				if(v<FT(0))
					sum_e-=v;
				else
					sum_e+=v;

				cr_e--;
			}
		}
				
		v=cr_e->y()-a*cr_e->x()-b;

		if(v<FT(0))
			sum_e-=v;
		else
			sum_e+=v;

		return sum_b+sum_e;
	}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		FT dx,sum_b,sum_e;
		InputIterator end,cr_b,cr_e;

		cr_b=begin;
		cr_b++;

		end=beyond;
		end--;

		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=FT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<FT(0))
						dx=-dx;
					
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<FT(0))
						dx=-dx;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<FT(0))
				dx=-dx;

			sum_e+=dx;

			split_pt=cr_b;

			return sum_b+sum_e;
		}

		FT a,b,v;
		
		dx=end->x()-begin->x();
		 
		if(dx==FT(0))
		{
			cr_b=begin;
			cr_b++;

			split_pt=cr_b;

			return FT(typename KT::RT(INT_MAX));
		}

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		sum_b=sum_e=FT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				v=cr_b->y()-a*cr_b->x()-b;
				if(v<FT(0))
					sum_b-=v;
				else
					sum_b+=v;

				cr_b++;
			}
			else
			{
				v=cr_e->y()-a*cr_e->x()-b;
				if(v<FT(0))
					sum_e-=v;
				else
					sum_e+=v;

				cr_e--;
			}
		}
				
		v=cr_e->y()-a*cr_e->x()-b;

		if(v<FT(0))
			sum_e-=v;
		else
			sum_e+=v;

		split_pt=cr_e;

		return sum_b+sum_e;
	}
};


/////////////////////////////////////////
//
// Bounded Volumes Distance Error

template<class KT>
class Bounded_volume_error		
{
  public:
	typedef typename KT::FT FT;		// base number type

  public:
	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		InputIterator end,cr;
	
		FT sem,dx,dy,dx1,dy1,t,a1,b1,a2,b2;
		FT ccx_cr,ccy_cr,rc_cr,ccx_p,ccy_p,rc_p,ccx_n,ccy_n,rc_n;
		FT dst,temp;
		int poz_cr,poz_cc;
		bool cd1;
		bool same_side_p, same_side_n;
		bool c_inf_p,c_inf_n;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
			return FT(0);

		if(*begin==*end)								// Closed curve case
		{
			typedef CGAL::Min_circle_2_traits_2<KT> C_Traits;
			typedef CGAL::Min_circle_2<C_Traits> Min_circle;
			typedef CGAL::Point_2<KT> Point_2;

			int npc;
			Point_2* cp;
			FT d01,d12,d02;
			Point_2 tt;

			Min_circle  mc2(begin,end,true);     

			npc=mc2.number_of_support_points( );

			cp=(Point_2*)mc2.support_points_begin();

			dx=cp[0].x()-cp[1].x();
			dy=cp[0].y()-cp[1].y();
			d01=dx*dx+dy*dy;

			if(npc==3)
			{
				dx=cp[0].x()-cp[2].x();
				dy=cp[0].y()-cp[2].y();
				d02=dx*dx+dy*dy;

				dx=cp[2].x()-cp[1].x();
				dy=cp[2].y()-cp[1].y();
				d12=dx*dx+dy*dy;
			
				if(d02>d01 && d02>d12)
				{
					tt=cp[2];
					cp[2]=cp[1];
					cp[1]=tt;
					d01=d02;
				}
				else
					if(d12>d01 && d12>d02)
					{
						tt=cp[2];
						cp[2]=cp[1];
						cp[1]=cp[0];
						cp[0]=tt;
						d01=d12;
					}
			}

			InputIterator cr1;
			while(cp[0]!=*begin && cp[1]!=*begin)
			{
				*end=*begin;
				for(cr=cr1=begin,cr1++;cr!=end;cr++,cr1++)
					*cr=*cr1;
			}
			
			*end=*begin;

			return d01;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();

		sem=begin->x()*dy-begin->y()*dx;

		cd1=false;

		t=begin->y()-end->y();
		if(t!=FT(0))
		{
			a1=(end->x()-begin->x())/t;
			b1=(begin->y()+end->y())/FT(2)-a1*(begin->x()+end->x())/FT(2);
		}
		else
		{
			ccx_cr=(begin->x()+end->x())/FT(2);
			cd1=true;
		}

		same_side_p=same_side_n=false;
		c_inf_p=c_inf_n=true;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			temp=cr->y()*dx-cr->x()*dy+sem;

			if(temp==FT(0))
				continue;

			if(temp>FT(0))
				poz_cr=1;
			else 
				if(temp<FT(0))
					poz_cr=-1;

			if(poz_cr>0)
			{
				if(!c_inf_p)
				{
					dx1=cr->x()-ccx_p;
					dy1=cr->y()-ccy_p;
					dst=dx1*dx1+dy1*dy1;

					if(dst<rc_p)
						continue;
				}
			}
			else
				if(!c_inf_n)
				{
					dx1=cr->x()-ccx_n;
					dy1=cr->y()-ccy_n;
					dst=dx1*dx1+dy1*dy1;

					if(dst<rc_n)
						continue;
				}

			t=cr->y()-end->y();
			if(t!=FT(0))
			{
				a2=(end->x()-cr->x())/t;
				b2=(cr->y()+end->y())/FT(2)-a2*(cr->x()+end->x())/FT(2);
				if(cd1)
					ccy_cr=a2*ccx_cr+b2;
				else
				{
					ccx_cr=(b2-b1)/(a1-a2);
					ccy_cr=a1*ccx_cr+b1;
				}
			}
			else
			{
				ccx_cr=(cr->x()+end->x())/FT(2);
				ccy_cr=a1*ccx_cr+b1;
			}
				
			a2=begin->x()-ccx_cr;
			b2=begin->y()-ccy_cr;
			rc_cr=a2*a2+b2*b2;

			temp=ccy_cr*dx-ccx_cr*dy+sem;
			if(temp>FT(0))
				poz_cc=1;
			else 
				if(temp<FT(0))
					poz_cc=-1;
				else
					poz_cc=poz_cr;

			if(poz_cr>0)
				if(c_inf_p)
				{
					c_inf_p=false;
					if(poz_cc>0)
						same_side_p=true;

					rc_p=rc_cr;
					ccx_p=ccx_cr;
					ccy_p=ccy_cr;
				}
				else
					if(same_side_p)
					{
						if(poz_cc>0)
							if(rc_p<rc_cr)
							{
								rc_p=rc_cr;
								ccx_p=ccx_cr;
								ccy_p=ccy_cr;
							}
					}
					else
						if(poz_cc<0)
						{
							if(rc_p>rc_cr)
							{
								rc_p=rc_cr;
								ccx_p=ccx_cr;
								ccy_p=ccy_cr;
							}
						}
						else
						{
							same_side_p=true;
							rc_p=rc_cr;
							ccx_p=ccx_cr;
							ccy_p=ccy_cr;
						}
			else
				if(c_inf_n)
				{
					c_inf_n=false;
					if(poz_cc<0)
						same_side_n=true;

					rc_n=rc_cr;
					ccx_n=ccx_cr;
					ccy_n=ccy_cr;
				}
				else
					if(same_side_n)
					{
						if(poz_cc<0)
							if(rc_n<rc_cr)
							{
								rc_n=rc_cr;
								ccx_n=ccx_cr;
								ccy_n=ccy_cr;
							}
					}
					else
						if(poz_cc>0)
						{
							if(rc_n>rc_cr)
							{
								rc_n=rc_cr;
								ccx_n=ccx_cr;
								ccy_n=ccy_cr;
							}
						}
						else
						{
							same_side_n=true;
							rc_n=rc_cr;
							ccx_n=ccx_cr;
							ccy_n=ccy_cr;
						}
		}

		dst=(dx*dx+dy*dy)/FT(4);

		if(c_inf_p)
			if(c_inf_n)
				return 0;
			else
			{
				if(same_side_n)
					return FT(4)*rc_n-FT(2)*dst;
				else
					return dst*dst/(FT(4)*rc_n-FT(3)*dst);
			}
		else
			if(c_inf_n)
			{
				if(same_side_p)
					return FT(4)*rc_p-FT(2)*dst;
				else
					return dst*dst/(FT(4)*rc_p-FT(3)*dst);
			}
			else
				if(same_side_p)
					if(same_side_n)
						if(rc_p>rc_n)
							return FT(4)*rc_p-FT(2)*dst;
						else
							return FT(4)*rc_n-FT(2)*dst;
					else
						return FT(4)*rc_p-FT(2)*dst;
				else
					if(same_side_n)
						return FT(4)*rc_n-FT(2)*dst;
					else
						if(rc_p<rc_n)
							return dst*dst/(FT(4)*rc_p-FT(3)*dst);
						else
							return dst*dst/(FT(4)*rc_n-FT(3)*dst);
}

	template<class InputIterator>
	FT operator()(InputIterator begin,InputIterator beyond,InputIterator &split_pt)
	{
		InputIterator end,cr,cr_split_pt;
		InputIterator p_max_p,p_max_n;
	
		FT sem,dx,dy,dx1,dy1,t,a1,b1,a2,b2;
		FT ccx_cr,ccy_cr,rc_cr,ccx_p,ccy_p,rc_p,ccx_n,ccy_n,rc_n;
		FT dst,temp;
		int poz_cr,poz_cc;
		bool cd1;
		bool same_side_p, same_side_n;
		bool c_inf_p,c_inf_n;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			split_pt=begin;
			return FT(0);
		}

		if(*begin==*end)								// Closed curve case
		{
			typedef CGAL::Min_circle_2_traits_2<KT>  C_Traits;
			typedef CGAL::Min_circle_2<C_Traits>     Min_circle;
			typedef CGAL::Point_2<KT> Point_2;

			int npc;
			Point_2* cp;
			FT d01,d12,d02;
			Point_2 tt;

			Min_circle  mc2(begin,end,true);     

			npc=mc2.number_of_support_points( );

			cp=(Point_2*)mc2.support_points_begin();

			dx=cp[0].x()-cp[1].x();
			dy=cp[0].y()-cp[1].y();
			d01=dx*dx+dy*dy;

			if(npc==3)
			{
				dx=cp[0].x()-cp[2].x();
				dy=cp[0].y()-cp[2].y();
				d02=dx*dx+dy*dy;

				dx=cp[2].x()-cp[1].x();
				dy=cp[2].y()-cp[1].y();
				d12=dx*dx+dy*dy;
			
				if(d02>d01 && d02>d12)
				{
					tt=cp[2];
					cp[2]=cp[1];
					cp[1]=tt;
					d01=d02;
				}
				else
					if(d12>d01 && d12>d02)
					{
						tt=cp[2];
						cp[2]=cp[1];
						cp[1]=cp[0];
						cp[0]=tt;
						d01=d12;
					}
			}

			InputIterator cr1;
			while(cp[0]!=*begin && cp[1]!=*begin)
			{
				*end=*begin;
				for(cr=cr1=begin,cr1++;cr!=end;cr++,cr1++)
					*cr=*cr1;
			}
			
			*end=*begin;

			if(*begin==cp[0])
			{
				for(cr=begin,cr++;*cr!=cp[1];cr++);
				cr_split_pt=cr;
			}
			else
			{
				for(cr=begin,cr++;*cr!=cp[0];cr++);
				cr_split_pt=cr;
			}
				
			split_pt=cr_split_pt;

			return d01;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();

		sem=begin->x()*dy-begin->y()*dx;

		cd1=false;

		t=begin->y()-end->y();
		if(t!=FT(0))
		{
			a1=(end->x()-begin->x())/t;
			b1=(begin->y()+end->y())/FT(2)-a1*(begin->x()+end->x())/FT(2);
		}
		else
		{
			ccx_cr=(begin->x()+end->x())/FT(2);
			cd1=true;
		}

		same_side_p=same_side_n=false;
		c_inf_p=c_inf_n=true;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			temp=cr->y()*dx-cr->x()*dy+sem;

			if(temp==FT(0))
				continue;

			if(temp>FT(0))
				poz_cr=1;
			else 
				if(temp<FT(0))
					poz_cr=-1;

			if(poz_cr>0)
			{
				if(!c_inf_p)
				{
					dx1=cr->x()-ccx_p;
					dy1=cr->y()-ccy_p;
					dst=dx1*dx1+dy1*dy1;

					if(dst<rc_p)
						continue;
				}
			}
			else
				if(!c_inf_n)
				{
					dx1=cr->x()-ccx_n;
					dy1=cr->y()-ccy_n;
					dst=dx1*dx1+dy1*dy1;

					if(dst<rc_n)
						continue;
				}

			t=cr->y()-end->y();
			if(t!=FT(0))
			{
				a2=(end->x()-cr->x())/t;
				b2=(cr->y()+end->y())/FT(2)-a2*(cr->x()+end->x())/FT(2);
				if(cd1)
					ccy_cr=a2*ccx_cr+b2;
				else
				{
					ccx_cr=(b2-b1)/(a1-a2);
					ccy_cr=a1*ccx_cr+b1;
				}
			}
			else
			{
				ccx_cr=(cr->x()+end->x())/FT(2);
				ccy_cr=a1*ccx_cr+b1;
			}
				
			a2=begin->x()-ccx_cr;
			b2=begin->y()-ccy_cr;
			rc_cr=a2*a2+b2*b2;

			temp=ccy_cr*dx-ccx_cr*dy+sem;
			if(temp>FT(0))
				poz_cc=1;
			else 
				if(temp<FT(0))
					poz_cc=-1;
				else
					poz_cc=poz_cr;

			if(poz_cr>0)
				if(c_inf_p)
				{
					c_inf_p=false;
					if(poz_cc>0)
						same_side_p=true;

					rc_p=rc_cr;
					ccx_p=ccx_cr;
					ccy_p=ccy_cr;
					p_max_p=cr;
				}
				else
					if(same_side_p)
					{
						if(poz_cc>0)
							if(rc_p<rc_cr)
							{
								rc_p=rc_cr;
								ccx_p=ccx_cr;
								ccy_p=ccy_cr;
								p_max_p=cr;
							}
					}
					else
						if(poz_cc<0)
						{
							if(rc_p>rc_cr)
							{
								rc_p=rc_cr;
								ccx_p=ccx_cr;
								ccy_p=ccy_cr;
								p_max_p=cr;
							}
						}
						else
						{
							same_side_p=true;
							rc_p=rc_cr;
							ccx_p=ccx_cr;
							ccy_p=ccy_cr;
							p_max_p=cr;
						}
			else
				if(c_inf_n)
				{
					c_inf_n=false;
					if(poz_cc<0)
						same_side_n=true;

					rc_n=rc_cr;
					ccx_n=ccx_cr;
					ccy_n=ccy_cr;
					p_max_n=cr;
				}
				else
					if(same_side_n)
					{
						if(poz_cc<0)
							if(rc_n<rc_cr)
							{
								rc_n=rc_cr;
								ccx_n=ccx_cr;
								ccy_n=ccy_cr;
								p_max_n=cr;
							}
					}
					else
						if(poz_cc>0)
						{
							if(rc_n>rc_cr)
							{
								rc_n=rc_cr;
								ccx_n=ccx_cr;
								ccy_n=ccy_cr;
								p_max_n=cr;
							}
						}
						else
						{
							same_side_n=true;
							rc_n=rc_cr;
							ccx_n=ccx_cr;
							ccy_n=ccy_cr;
							p_max_n=cr;
						}
		}

		dst=(dx*dx+dy*dy)/FT(4);

		if(c_inf_p)
			if(c_inf_n)
			{
				split_pt=begin;
				return 0;
			}
			else
			{
				split_pt=p_max_n;

				if(same_side_n)
					return FT(4)*rc_n-FT(2)*dst;
				else
					return dst*dst/(FT(4)*rc_n-FT(3)*dst);
			}
		else
			if(c_inf_n)
			{
				split_pt=p_max_p;

				if(same_side_p)
					return FT(4)*rc_p-FT(2)*dst;
				else
					return dst*dst/(FT(4)*rc_p-FT(3)*dst);
			}
			else
				if(same_side_p)
					if(same_side_n)
						if(rc_p>rc_n)
						{
							split_pt=p_max_p;
							return FT(4)*rc_p-FT(2)*dst;
						}
						else
						{
							split_pt=p_max_n;
							return FT(4)*rc_n-FT(2)*dst;
						}
					else
					{
						split_pt=p_max_p;
						return FT(4)*rc_p-FT(2)*dst;
					}
				else
					if(same_side_n)
					{
						split_pt=p_max_n;
						return FT(4)*rc_n-FT(2)*dst;
					}
					else
						if(rc_p<rc_n)
						{
							split_pt=p_max_p;
							return dst*dst/(FT(4)*rc_p-FT(3)*dst);
						}
						else
						{
							split_pt=p_max_n;
							return dst*dst/(FT(4)*rc_n-FT(3)*dst);
						}
	}
};


////////////////////////////////////////////////////
//
//  2D SUM Euclidean Distance - measured to the line support; 
//								Perez & Vidal algorithm (incremental alg.); 
//								It works only for Dynamic Prog.and Graph Search Methods
template<class KT>
class Sum_squared_euclidean_error_inc
{
  public:
	typedef typename KT::FT FT;		// base number type

	FT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	FT nr_p;
	void* beginPoint;
	void* endPoint;

  public:
	Sum_squared_euclidean_error_inc()
	{
		x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
		nr_p=0;
		beginPoint=endPoint=NULL;
	}

	template<class InputIterator>				
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT ssvdErr,ssedErr,xk,yk,yValY,dir_coef,dif;
	
		InputIterator cr,end;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)
		{	x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
			nr_p=FT(0);
			beginPoint=&begin;
			endPoint=&end;
			return FT(0);
		}

		cr++;
		if(cr==end)
		{
			cr--;
			x_sum=y_sum=x2_sum=y2_sum=xy_sum=FT(0);
			nr_p=FT(1);
			xk=cr->x();	 
			yk=cr->y();
			beginPoint=&begin;
			endPoint=&end;
		}
		else 
			if (begin==*(InputIterator*)beginPoint)			//First endpoint is the same as in the last call.
			{												//Take the coordinates of the point before the last      				 
				cr=end;
				cr--;
				xk=cr->x();				//endpoint.
				yk=cr->y();
				nr_p=nr_p+FT(1);
			}
			else 
				if (end==*(InputIterator*)endPoint)			//Second endpoint is the same as in the last call.
				{							//Take the coordinates of the point after the first
					cr=begin;
					cr++;
					xk=cr->x();			//endpoint.
					yk=cr->y();
					nr_p=nr_p+FT(1);
				}
				else						//The two endpoint are both new, which shouldn't be possible.
					return FT(typename KT::RT(INT_MAX));
    
		x_sum+=xk;
		y_sum+=yk;
		x2_sum+=(xk*xk);
		y2_sum+=(yk*yk); 
		xy_sum+=(xk*yk);

		//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		dif=(end->x()-begin->x());
		if(dif!=FT(0))
		{
			dir_coef=(end->y()-begin->y())/dif;
			yValY=begin->y()-dir_coef*begin->x();

			ssvdErr=yValY*yValY*nr_p+FT(2)*yValY*dir_coef*x_sum-FT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-FT(2)*dir_coef*xy_sum;
			ssedErr=(FT(1)/(dir_coef*dir_coef+FT(1)))*ssvdErr;
		}
		else
		{
			if((end->y()-begin->y())==FT(0))				// Closed curve case (*begin==*end)
			{	
				xx=begin->x();
				yValY=begin->y();

				ssedErr=xx*xx*nr_p+yValY*yValY*nr_p-FT(2)*yValY*y_sum+y2_sum+x2_sum-FT(2)*xx*x_sum;
			}
			else
			{
				xx=begin->x();
				return x2_sum-FT(2)*xx*x_sum+nr_p*xx*xx;
			}
		}

		if(ssedErr<FT(0))
			return -ssedErr;
		return ssedErr;
	}
};


////////////////////////////////////////////////////
//
//  SUM Verical Distance -	measured to the line support; 
//							Perez & Vidal algorithm (incremental alg.); 
//							It works only for Dynamic Prog. and Graph Search Methods
template<class KT>
class Sum_squared_vertical_error_inc
{
  public:
	typedef typename KT::FT FT;		// base number type

	FT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	FT nr_p;
	void* beginPoint;
	void* endPoint;

  public:
	Sum_squared_vertical_error_inc()
	{
		x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
		nr_p=0;
		beginPoint=endPoint=NULL;
	}

	template<class InputIterator>					
	FT operator()(InputIterator begin,InputIterator beyond)
	{
		FT ssvdErr,xk,yk,yValY,dir_coef,dif;
	
		InputIterator cr,end;

		cr=begin;
		cr++;

		end=beyond;
		end--;

		if(begin==end || cr==end)
		{	x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
			nr_p=FT(0);
			beginPoint=&begin;
			endPoint=&end;
			return FT(0);
		}

		cr++;
		if(cr==end)
		{
			cr--;
			x_sum=y_sum=x2_sum=y2_sum=xy_sum=FT(0);
			nr_p=FT(1);
			xk=cr->x();	 
			yk=cr->y();
			beginPoint=&begin;
			endPoint=&end;
		}
		else 
			if (begin==*(InputIterator*)beginPoint)			//First endpoint is the same as in the last call.
			{								//Take the coordinates of the point before the last      				 
				cr=end;
				cr--;
				xk=cr->x();				//endpoint.
				yk=cr->y();
				nr_p=nr_p+FT(1);
			}
			else 
				if (end==*(InputIterator*)endPoint)			//Second endpoint is the same as in the last call.
				{							//Take the coordinates of the point after the first
					cr=begin;
					cr++;
					xk=cr->x();			//endpoint.
					yk=cr->y();
					nr_p=nr_p+FT(1);
				}
				else	//The two endpoint are both new, which shouldn't be possible.
					return FT(typename KT::RT(INT_MAX));
    
		x_sum+=xk;
		y_sum+=yk;
		x2_sum+=(xk*xk);
		y2_sum+=(yk*yk); 
		xy_sum+=(xk*yk);

		//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		dif=(end->x()-begin->x());
		if(dif!=FT(0))
		{
			dir_coef=(end->y()-begin->y())/dif;
			yValY=begin->y()-dir_coef*begin->x();

			ssvdErr=yValY*yValY*nr_p+FT(2)*yValY*dir_coef*x_sum-FT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-FT(2)*dir_coef*xy_sum;
		}
		else											// Closed curve case (*begin==*end)
		{
			if((end->y()-begin->y())==FT(0))   
				ssvdErr=nr_p*begin->x()*begin->x()+x2_sum-FT(2)*begin->x()*x_sum;
			else
				return FT(typename KT::RT(INT_MAX));
		}

		if(ssvdErr<FT(0))
			return -ssvdErr;
		return ssvdErr;
	}
};



CGAL_END_NAMESPACE

#endif // CGAL_POLYAP_TRAITS

