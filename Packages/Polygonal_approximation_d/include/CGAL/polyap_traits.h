
#ifndef PLOYAP_TRAITS
#define PLOYAP_TRAITS

#include <assert.h>
#include <CGAL/basic.h>

#include <CGAL/Point_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>



CGAL_BEGIN_NAMESPACE



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//				P O L Y G O N A L   A P P R O X I M A T I O N   T R A I T S   C L A S S E S 

///////////////
//
//  Base class

template<class KT>
class DistTypeTraits		
{	
  public:

	typedef typename KT::FT DataT;						// base number type
	typedef typename KT::RT DT;							// base data type

  private:
	
	bool isIntType;

	DT cmmdc(DT n1,DT n2)
	{	DT rest;

		rest=n1-n1/n2*n2;
		
		if(rest==0)
			return n2;
		else
			return cmmdc(n2,rest);
	}

  public:

	DistTypeTraits()
	{
		isIntType=((DT)5./(DT)2.*2.)==4;
	}
	
	virtual DataT cumulate(DataT curveError, DataT deltaError)
	{
		if(deltaError>curveError)
			return deltaError;
		else
			return curveError;
	}

	void rescale(DataT &nr)
	{	
		if(isIntType)
		{
			DT num,den;

			num=KT::FT_numerator(nr);
			den=KT::FT_denominator(nr);

			if(num!=1 && den!=1 && num!=0 && den!=0)
			{
				DT dc=cmmdc(num,den);

				if(dc!=1)
				{
					num=num/dc;
					den=den/dc;

					nr=KT::make_FT(DT(num),DT(den));
				}
			}
		}
	}
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//				MAX Distance Functions


//////////////////////////////////////////////////////
//
//  MAX Euclidean Distance - measured to the line support
template<class KT>
class MaxSquaredEuclideanDistanceError_ls: public DistTypeTraits<KT>
{		
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT tr_y,dx,dy_n,max;
		InputIterator cr,p_max;

		cr=begin;
		cr++;
		
		if(begin==end || cr==end)		// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)				// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
				{
					max=tr_y;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;

			return max;
		}
		
		typename DistTypeTraits<KT>::DataT d,v;

		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=v+cr->x()*dy_n+cr->y()*dx;
			tr_y=tr_y*tr_y/d;
						
			if(tr_y>max)
			{
				max=tr_y;
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;
		return max;
	}
};


//////////////////////////////////////////////////////
//
//  MAX Euclidean Distance - measured to the segment


template<class KT>
class MaxSquaredEuclideanDistanceError_seg: public DistTypeTraits<KT>
{
  public:
	
	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,dy_n,max,tr_y;
		InputIterator cr,p_max;

		cr=begin;
		cr++;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx=(begin->x()-cr->x());
				dy_n=(begin->y()-cr->y());

				tr_y=dx*dx+dy_n*dy_n;
			
				if(tr_y>max)
				{
					max=tr_y;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;

			return max;
		}
		
		typename DistTypeTraits<KT>::DataT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

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
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;
		return max;
	}
};


//////////////////////////////////////////////////////////
//
// MAX L-infinity Distance - measured to the line support

template<class KT>
class MaxLinfinityError_ls: public DistTypeTraits<KT>
{
  public:
	
	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT max,dx1,dy1;
		InputIterator cr,p_max;
		
		cr=begin;
		cr++;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<DataT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<DataT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
				{
					max=dx1;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;
			return max;
		}

		typename DistTypeTraits<KT>::DataT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
								
			dx1=dy_n*tr_y;
			if(dx1<DataT(0))
				dx1=-dx1;
	
			dy1=dx*tr_y;
			if(dy1<DataT(0))
				dy1=-dy1;
	
			if(dy1<dx1)
				dy1=dx1;
								
			if(dy1>max)
			{
				max=dy1;	
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;
		return max;
	}
};


//////////////////////////////////////////////////////////
//
// MAX L-infinity Distance - measured to the segment

template<class KT>
class MaxLinfinityError_seg: public DistTypeTraits<KT>				
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT max,dx1,dy1;
		InputIterator cr,p_max;
		
		cr=begin;
		cr++;
 
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<DataT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<DataT(0))
					dy1=-dy1;

				if(dx1<dy1)
					dx1=dy1;

				if(dx1>max)
				{
					max=dx1;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;

			return max;
		}

		typename DistTypeTraits<KT>::DataT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

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
				
			if(dx1<DataT(0))
				dx1=-dx1;
			if(dy1<DataT(0))
				dy1=-dy1;

			if(dy1<dx1)
				dy1=dx1;

			if(dy1>max)
			{
				max=dy1;
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;

		return max;
	}
};


/////////////////////////////////////////////////////////
//
//  MAX Manhattan Distance - measured to the line support

template<class KT>
class MaxManhattanError_ls: public DistTypeTraits<KT>		
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT max,dx1,dy1;
		InputIterator cr,p_max;
		
		cr=begin;
		cr++;
		
		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}
		
		if(*begin==*end)						// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<DataT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<DataT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
				{
					max=dx1;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;
			return max;
		}
		
		typename DistTypeTraits<KT>::DataT d,dx,dy_n,tr_y,v;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			tr_y=(v+cr->x()*dy_n+cr->y()*dx)/d;
									
			dx1=dy_n*tr_y;
			if(dx1<DataT(0))
				dx1=-dx1;
		
			dy1=dx*tr_y;
			if(dy1<DataT(0))
				dy1=-dy1;
		
			dy1+=dx1;
									
			if(dy1>max)
			{
				max=dy1;
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;

		return max;
	}
};


/////////////////////////////////////////////////////////
//
//  MAX Manhattan Distance - measured to the segment

template<class KT>
class MaxManhattanError_seg: public DistTypeTraits<KT>						
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT max,dx1,dy1;
		InputIterator cr,p_max;
		
		cr=begin;
		cr++;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				dx1=(begin->x()-cr->x());
				if(dx1<DataT(0))
					dx1=-dx1;

				dy1=(begin->y()-cr->y());
				if(dy1<DataT(0))
					dy1=-dy1;

				dx1+=dy1;

				if(dx1>max)
				{
					max=dx1;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;

			return max;
		}

		typename DistTypeTraits<KT>::DataT d,dx,dy_n,tr_y,v,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		max=DataT(0);
		p_max=begin;

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

			if(dx1<DataT(0))
				dx1=-dx1;
			if(dy1<DataT(0))
				dy1=-dy1;
				
			dy1+=dx1;

			if(dy1>max)
			{
				max=dy1;
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;

		return max;
	}
};


///////////////////////////////
//
//  MAX Vertical Distance 

template<class KT>
class MaxVerticalError: public DistTypeTraits<KT>
{
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT max,v;
		InputIterator cr,p_max;
		
		cr=begin;
		cr++;

		if(begin==end || cr==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			max=DataT(0);
			p_max=begin;

			for(cr=begin,cr++;cr!=end;cr++)
			{
				v=(begin->y()-cr->y());
				if(v<DataT(0))
					v=-v;

				if(v>max)
				{
					max=v;
					p_max=cr;
				}
			}
			if(p_MAX)
				*p_MAX=p_max;

			return max;
		}
		
		typename DistTypeTraits<KT>::DataT dx,a,b;
		
		dx=end->x()-begin->x();
		if(dx==DataT(0))
		{
			if(p_MAX)
				*p_MAX=cr;

			return DataT((typename KT::RT)(INT_MAX));
		}

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		max=DataT(0);
		p_max=begin;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			v=cr->y()-a*cr->x()-b;
			if(v<DataT(0))
				v=-v;
			
			if(max<v)
			{
				max=v;
				p_max=cr;
			}
		}
		if(p_MAX)
			*p_MAX=p_max;
		return max;
	}
};


/////////////////////////////////////////
//
// Bounding Volumes Distance 
//

template<class KT>
class MaxBoundVol: public DistTypeTraits<KT>		
{
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		InputIterator cr,p_max;
		InputIterator p_max_p,p_max_n;
	
		typename DistTypeTraits<KT>::DataT sem,dx,dy,dx1,dy1,t,a1,b1,a2,b2;
		typename DistTypeTraits<KT>::DataT ccx_cr,ccy_cr,rc_cr,ccx_p,ccy_p,rc_p,ccx_n,ccy_n,rc_n;
		typename DistTypeTraits<KT>::DataT dst,temp;
		int poz_cr,poz_cc;
		bool cd1;
		bool same_side_p, same_side_n;
		bool c_inf_p,c_inf_n;

		cr=begin;
		cr++;

		if(begin==end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)								// Closed curve case
		{
			typedef CGAL::Min_circle_2_traits_2<KT>  Traits;
			typedef CGAL::Min_circle_2<Traits>      Min_circle;
			typedef CGAL::Point_2<KT> Point2;

			int npc;
			Point2* cp;
			typename DistTypeTraits<KT>::DataT d01,d12,d02;
			Point2 tt;

			Min_circle  mc2(begin,end,true);     

			npc=mc2.number_of_support_points( );

			cp=(Point2*)mc2.support_points_begin();

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
				p_max=cr;
			}
			else
			{
				for(cr=begin,cr++;*cr!=cp[0];cr++);
				p_max=cr;
			}
			if(p_MAX)
				*p_MAX=p_max;
			return d01;
		}

		dx=end->x()-begin->x();
		dy=end->y()-begin->y();

		sem=begin->x()*dy-begin->y()*dx;

		cd1=false;

		t=begin->y()-end->y();
		if(t!=DataT(0))
		{
			a1=(end->x()-begin->x())/t;
			b1=(begin->y()+end->y())/DataT(2)-a1*(begin->x()+end->x())/DataT(2);
		}
		else
		{
			ccx_cr=(begin->x()+end->x())/DataT(2);
			cd1=true;
		}

		same_side_p=same_side_n=false;
		c_inf_p=c_inf_n=true;

		for(cr=begin,cr++;cr!=end;cr++)
		{
			temp=cr->y()*dx-cr->x()*dy+sem;

			if(temp==DataT(0))
				continue;

			if(temp>DataT(0))
				poz_cr=1;
			else 
				if(temp<DataT(0))
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
			if(t!=DataT(0))
			{
				a2=(end->x()-cr->x())/t;
				b2=(cr->y()+end->y())/DataT(2)-a2*(cr->x()+end->x())/DataT(2);
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
				ccx_cr=(cr->x()+end->x())/DataT(2);
				ccy_cr=a1*ccx_cr+b1;
			}
				
			a2=begin->x()-ccx_cr;
			b2=begin->y()-ccy_cr;
			rc_cr=a2*a2+b2*b2;

			temp=ccy_cr*dx-ccx_cr*dy+sem;
			if(temp>DataT(0))
				poz_cc=1;
			else 
				if(temp<DataT(0))
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

		dst=(dx*dx+dy*dy)/DataT(4);

		if(c_inf_p)
			if(c_inf_n)
			{
				if(p_MAX)
					*p_MAX=begin;
				return 0;
			}
			else
			{
				if(p_MAX)
					*p_MAX=p_max_n;

				if(same_side_n)
					return DataT(4)*rc_n-DataT(2)*dst;
				else
					return dst*dst/(DataT(4)*rc_n-DataT(3)*dst);
			}
		else
			if(c_inf_n)
			{
				if(p_MAX)
					*p_MAX=p_max_p;

				if(same_side_p)
					return DataT(4)*rc_p-DataT(2)*dst;
				else
					return dst*dst/(DataT(4)*rc_p-DataT(3)*dst);
			}
			else
				if(same_side_p)
					if(same_side_n)
						if(rc_p>rc_n)
						{
							if(p_MAX)
								*p_MAX=p_max_p;
							return DataT(4)*rc_p-DataT(2)*dst;
						}
						else
						{
							if(p_MAX)
								*p_MAX=p_max_n;
							return DataT(4)*rc_n-DataT(2)*dst;
						}
					else
					{
						if(p_MAX)
							*p_MAX=p_max_p;
						return DataT(4)*rc_p-DataT(2)*dst;
					}
				else
					if(same_side_n)
					{
						if(p_MAX)
							*p_MAX=p_max_n;
						return DataT(4)*rc_n-DataT(2)*dst;
					}
					else
						if(rc_p<rc_n)
						{
							if(p_MAX)
								*p_MAX=p_max_p;
							return dst*dst/(DataT(4)*rc_p-DataT(3)*dst);
						}
						else
						{
							if(p_MAX)
								*p_MAX=p_max_n;
							return dst*dst/(DataT(4)*rc_n-DataT(3)*dst);
						}
	}
};




/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//
//   SUM Distance Functions
//

////////////////////////////////////////////////////
//
//  SUM Euclidean Distance - measured to the line support; 
//						     Perez & Vidal algorithm (incremental alg.); It works only for Dynamic Prog. Methods

template<class KT>
class SumSquaredEuclideanDistanceError_inc: public DistTypeTraits<KT>
{
	typename DistTypeTraits<KT>::DataT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	typename DistTypeTraits<KT>::DataT nr_p;
	void* beginPoint;
	void* endPoint;


  public:

	SumSquaredEuclideanDistanceError_inc()
	{
		x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
		nr_p=0;
		beginPoint=endPoint=NULL;
	}

	template<class InputIterator>				
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT ssvdErr,ssedErr,xk,yk,yValY,dir_coef,dif;
	
		InputIterator curr=begin;
		curr++;
		if(begin==end || curr==end)
		{	x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
			nr_p=DataT(0);
			beginPoint=&begin;
			endPoint=&end;
			return DataT(0);
		}

		curr++;
		if(curr==end)
		{
			curr--;
			x_sum=y_sum=x2_sum=y2_sum=xy_sum=DataT(0);
			nr_p=DataT(1);
			xk=curr->x();	 
			yk=curr->y();
			beginPoint=&begin;
			endPoint=&end;
		}
		else 
			if (begin==*(InputIterator*)beginPoint)			//First endpoint is the same as in the last call.
			{								//Take the coordinates of the point before the last      				 
				curr=end;
				curr--;
				xk=curr->x();				//endpoint.
				yk=curr->y();
				nr_p=nr_p+DataT(1);
			}
			else 
				if (end==*(InputIterator*)endPoint)			//Second endpoint is the same as in the last call.
				{							//Take the coordinates of the point after the first
					curr=begin;
					curr++;
					xk=curr->x();			//endpoint.
					yk=curr->y();
					nr_p=nr_p+DataT(1);
				}
				else						//The two endpoint are both new, which shouldn't be possible.
				{
					return DataT((typename KT::RT)(INT_MAX));
				}
    
		x_sum += xk;
		y_sum += yk;
		x2_sum += (xk * xk);
		y2_sum += (yk * yk); 
		xy_sum += (xk * yk);


		//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		dif=(end->x()-begin->x());
		if(dif!=DataT(0))
		{
			dir_coef=(end->y()-begin->y())/dif;
			yValY=begin->y()-dir_coef*begin->x();

			ssvdErr=yValY*yValY*nr_p+DataT(2)*yValY*dir_coef*x_sum-DataT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-DataT(2)*dir_coef*xy_sum;
			ssedErr=(DataT(1)/(dir_coef*dir_coef+DataT(1)))*ssvdErr;
		}
		else
		{
			if((end->y()-begin->y())==DataT(0))				// Closed curve case (*begin==*end)
			{	
				xx=begin->x();
				yValY=begin->y();

				ssedErr=xx*xx*nr_p+yValY*yValY*nr_p-DataT(2)*yValY*y_sum+y2_sum+x2_sum-DataT(2)*xx*x_sum;
			}
			else
			{
				xx=begin->x();
				return x2_sum-DataT(2)*xx*x_sum+nr_p*xx*xx;
			}
		}

		if(ssedErr<DataT(0))
			return -ssedErr;
		return ssedErr;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return curveError+deltaError;
	}	
};


////////////////////////////////////////////////////
//
//  SUM Verical Distance - measured to the line support; 
//						   Perez & Vidal algorithm (incremental alg.); It works only for Dynamic Prog. Methods

template<class KT>
class SumSquaredVerticalDistanceError_inc: public DistTypeTraits<KT>
{
	typename DistTypeTraits<KT>::DataT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	typename DistTypeTraits<KT>::DataT nr_p;
	void* beginPoint;
	void* endPoint;

  public:

	SumSquaredVerticalDistanceError_inc()
	{
		x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
		nr_p=0;
		beginPoint=endPoint=NULL;
	}

	template<class InputIterator>					
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT ssvdErr,xk,yk,yValY,dir_coef,dif;
	
		InputIterator curr=begin;
		curr++;
		if(begin==end || curr==end)
		{	x_sum=y_sum=x2_sum=y2_sum=xy_sum=0;
			nr_p=DataT(0);
			beginPoint=&begin;
			endPoint=&end;
			return DataT(0);
		}

		curr++;
		if(curr==end)
		{
			curr--;
			x_sum=y_sum=x2_sum=y2_sum=xy_sum=DataT(0);
			nr_p=DataT(1);
			xk=curr->x();	 
			yk=curr->y();
			beginPoint=&begin;
			endPoint=&end;
		}
		else 
			if (begin==*(InputIterator*)beginPoint)			//First endpoint is the same as in the last call.
			{								//Take the coordinates of the point before the last      				 
				curr=end;
				curr--;
				xk=curr->x();				//endpoint.
				yk=curr->y();
				nr_p=nr_p+DataT(1);
			}
			else 
				if (end==*(InputIterator*)endPoint)			//Second endpoint is the same as in the last call.
				{							//Take the coordinates of the point after the first
					curr=begin;
					curr++;
					xk=curr->x();			//endpoint.
					yk=curr->y();
					nr_p=nr_p+DataT(1);
				}
				else	//The two endpoint are both new, which shouldn't be possible.
				{
					return DataT((typename KT::RT)(INT_MAX));
				}
    
		x_sum += xk;
		y_sum += yk;
		x2_sum += (xk * xk);
		y2_sum += (yk * yk); 
		xy_sum += (xk * yk);

		//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		dif=(end->x()-begin->x());
		if(dif!=DataT(0))
		{
			dir_coef=(end->y()-begin->y())/dif;
			yValY=begin->y()-dir_coef*begin->x();

			ssvdErr=yValY*yValY*nr_p+DataT(2)*yValY*dir_coef*x_sum-DataT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-DataT(2)*dir_coef*xy_sum;
		}
		else											// Closed curve case (*begin==*end)
		{
			if((end->y()-begin->y())==DataT(0))   
				ssvdErr=nr_p*begin->x()*begin->x()+x2_sum-DataT(2)*begin->x()*x_sum;
			else
				return DataT((typename KT::RT)(INT_MAX));
		}

		if(ssvdErr<DataT(0))
			return -ssvdErr;
		return ssvdErr;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return curveError+deltaError;
	}
};


////////////////////////////////////////////////////
//
//  SUM Euclidean Distance - measured to the line support; 

template<class KT>
class SumSquaredEuclideanDistanceError_ls:public DistTypeTraits<KT>
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;//p_max;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=DataT(0);

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
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}
		
		typename DistTypeTraits<KT>::DataT d,v;
		
		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=DataT(0);

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
	
		if(p_MAX)
			*p_MAX=cr_b;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


//////////////////////////////////////////////////////////////////////
//
//  SUM Euclidean Distance - measured to the line support; uses number rescale

template<class KT>
class SumSquaredEuclideanDistanceError_ls_rescale: public DistTypeTraits<KT>			
{												
  public:

	template<class InputIterator>						
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					dy_n=(begin->y()-cr_b->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_b+=tr_y;
					rescale(sum_b);
					
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					dy_n=(end->y()-cr_e->y());

					tr_y=dx*dx+dy_n*dy_n;
			
					sum_e+=tr_y;
					rescale(sum_e);
					cr_e--;
				}
			}
			
			dx=(end->x()-cr_e->x());
			dy_n=(end->y()-cr_e->y());
			tr_y=dx*dx+dy_n*dy_n;
			
			sum_e+=tr_y;
			rescale(sum_e);

			sum_e+=sum_b;
			rescale(sum_e);
			
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_e;
		}
		
		typename DistTypeTraits<KT>::DataT d,v;
		
		dx=end->x()-begin->x();
		dy_n=begin->y()-end->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=DataT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=v+cr_b->x()*dy_n+cr_b->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_b+=tr_y;
				rescale(sum_b);

				cr_b++;
			}
			else
			{
				tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
				tr_y=tr_y*tr_y/d;
						
				sum_e+=tr_y;
				rescale(sum_e);
	
				cr_e--;
			}
		}
			
		tr_y=v+cr_e->x()*dy_n+cr_e->y()*dx;
		tr_y=tr_y*tr_y/d;
		sum_e+=tr_y;
		rescale(sum_e);

		if(p_MAX)
			*p_MAX=cr_b;

		sum_e+=sum_b;
		rescale(sum_e);
		return sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		deltaError+=curveError;
		rescale(deltaError);

		return deltaError;
	}
};


////////////////////////////////////////////////////
//
//  SUM Euclidean Distance - measured to the segment

template<class KT>
class SumSquaredEuclideanDistanceError_seg:public DistTypeTraits<KT>
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT tr_y,dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=DataT(0);

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
			
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}

		typename DistTypeTraits<KT>::DataT d,v,dd,dx1,dy1,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=DataT(0);

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
		if(p_MAX)
			*p_MAX=cr_e;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


////////////////////////////////////////////////////////
//
// SUM L-infinity Distance - measured to the line support

template<class KT>
class SumLinfinityError_ls:public DistTypeTraits<KT>				
{												
  public:
 
	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<DataT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<DataT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;
			
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}

		typename DistTypeTraits<KT>::DataT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=DataT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<DataT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<DataT(0))
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
				if(dx1<DataT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<DataT(0))
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
		if(dx1<DataT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<DataT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		if(p_MAX)
			*p_MAX=cr_b;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


////////////////////////////////////////////////////////
//
// SUM L-infinity Distance - measured to the segment

template<class KT>
class SumLinfinityError_seg:public DistTypeTraits<KT>
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					if(dx<dy_n)
						dx=dy_n;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<DataT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<DataT(0))
				dy_n=-dy_n;

			if(dx<dy_n)
				dx=dy_n;

			sum_e+=dx;

			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}
				
		typename DistTypeTraits<KT>::DataT d,tr_y,v,dx1,dy1,dd,vb,ve;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=DataT(0);

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

				if(dx1<DataT(0))
					dx1=-dx1;
				if(dy1<DataT(0))
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

				if(dx1<DataT(0))
					dx1=-dx1;
				if(dy1<DataT(0))
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

		if(dx1<DataT(0))
			dx1=-dx1;
		if(dy1<DataT(0))
			dy1=-dy1;
		
		if(dy1<dx1)
			sum_e+=dx1;
		else
			sum_e+=dy1;

		if(p_MAX)
			*p_MAX=cr_b;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


////////////////////////////////////////////////////////
//
// SUM Manhattan Distance - measured to the line support

template<class KT>
class SumManhattanError_ls:public DistTypeTraits<KT>
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)				// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)						// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(begin->y()-cr_b->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<DataT(0))
						dx=-dx;
					dy_n=(end->y()-cr_e->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<DataT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<DataT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}

		typename DistTypeTraits<KT>::DataT d,tr_y,v,dx1,dy1;

		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		sum_b=sum_e=DataT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				tr_y=(v+cr_b->x()*dy_n+cr_b->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<DataT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<DataT(0))
					dy1=-dy1;

				sum_b+=dy1+dx1;

				cr_b++;
			}
			else
			{
				tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
									
				dx1=dy_n*tr_y;
				if(dx1<DataT(0))
					dx1=-dx1;
		
				dy1=dx*tr_y;
				if(dy1<DataT(0))
					dy1=-dy1;

				sum_e+=dy1+dx1;

				cr_e--;
			}
		}
		tr_y=(v+cr_e->x()*dy_n+cr_e->y()*dx)/d;
					
		dx1=dy_n*tr_y;
		if(dx1<DataT(0))
			dx1=-dx1;
		
		dy1=dx*tr_y;
		if(dy1<DataT(0))
			dy1=-dy1;

		sum_e+=dy1+dx1;

		if(p_MAX)
			*p_MAX=cr_b;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


////////////////////////////////////////////////////////
//
// SUM Manhattan Distance - measured to the segment

template<class KT>
class SumManhattanError_seg:public DistTypeTraits<KT>						
{												
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,dy_n,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<DataT(0))
						dx=-dx;

					dy_n=(begin->y()-cr_b->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;
					
					sum_b+=dx+dy_n;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<DataT(0))
						dx=-dx;

					dy_n=(end->y()-cr_e->y());
					if(dy_n<DataT(0))
						dy_n=-dy_n;

					sum_e+=dx+dy_n;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<DataT(0))
				dx=-dx;
			dy_n=(end->y()-cr_e->y());
			if(dy_n<DataT(0))
				dy_n=-dy_n;

			sum_e+=dx+dy_n;
			
			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}
		
		typename DistTypeTraits<KT>::DataT d,tr_y,v,dx1,dy1,dd,vb,ve;
		
		dx=end->x()-begin->x();
		dy_n=-end->y()+begin->y();
		d=dx*dx+dy_n*dy_n;
		v=begin->x()*end->y()-begin->y()*end->x();

		vb=begin->x()*dx-begin->y()*dy_n;
		ve=end->x()*dx-end->y()*dy_n;

		sum_b=sum_e=DataT(0);

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

				if(dx1<DataT(0))
					dx1=-dx1;
				if(dy1<DataT(0))
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

				if(dx1<DataT(0))
					dx1=-dx1;
				if(dy1<DataT(0))
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

		if(dx1<DataT(0))
			dx1=-dx1;
		if(dy1<DataT(0))
			dy1=-dy1;
		sum_e+=dx1+dy1;

		if(p_MAX)
			*p_MAX=cr_b;

		return sum_b+sum_e;
	}

	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};


//////////////////////////
//
// SUM Vertical Distance

template<class KT>
class SumVerticalError:public DistTypeTraits<KT>
{
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator* p_MAX=NULL)
	{
		typename DistTypeTraits<KT>::DataT dx,sum_b,sum_e;
		InputIterator cr_b,cr_e;

		cr_b=begin;
		cr_b++;
		if(begin==end || cr_b==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			return DataT(0);
		}

		if(*begin==*end)							// Closed curve case
		{
			sum_b=sum_e=DataT(0);

			for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
			{
				if(sum_b<=sum_e)
				{
					dx=(begin->x()-cr_b->x());
					if(dx<DataT(0))
						dx=-dx;
					
					sum_b+=dx;
					cr_b++;
				}
				else
				{
					dx=(end->x()-cr_e->x());
					if(dx<DataT(0))
						dx=-dx;

					sum_e+=dx;
					cr_e--;
				}
			}

			dx=(end->x()-cr_e->x());
			if(dx<DataT(0))
				dx=-dx;

			sum_e+=dx;

			if(p_MAX)
				*p_MAX=cr_b;

			return sum_b+sum_e;
		}

		typename DistTypeTraits<KT>::DataT a,b,v;
		
		dx=end->x()-begin->x();
		 
		if(dx==DataT(0))
		{
			cr_b=begin;
			cr_b++;
			if(p_MAX)
				*p_MAX=cr_b;
			return DataT((typename KT::RT)(INT_MAX));
		}

		a=(end->y()-begin->y())/dx;
		b=begin->y()-a*begin->x();

		sum_b=sum_e=DataT(0);

		for(cr_b=begin,cr_b++, cr_e=end,cr_e--;cr_b!=cr_e;)
		{
			if(sum_b<=sum_e)
			{
				v=cr_b->y()-a*cr_b->x()-b;
				if(v<DataT(0))
					sum_b-=v;
				else
					sum_b+=v;

				cr_b++;
			}
			else
			{
				v=cr_e->y()-a*cr_e->x()-b;
				if(v<DataT(0))
					sum_e-=v;
				else
					sum_e+=v;

				cr_e--;
			}
		}
				
		v=cr_e->y()-a*cr_e->x()-b;

		if(v<DataT(0))
			sum_e-=v;
		else
			sum_e+=v;

		if(p_MAX)
			*p_MAX=cr_e;

		return sum_b+sum_e;
	}
	
	typename DistTypeTraits<KT>::DataT cumulate(typename DistTypeTraits<KT>::DataT curveError, typename DistTypeTraits<KT>::DataT deltaError)
	{
		return deltaError+curveError;
	}
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//       Douglas Peucker Hull Algorithm  -  Only for Recursive Split (recursive implementation), MIN E Method !!!
//

#define PUSH_OP 0			 //Operation names saved in history stack 
#define TOP_OP 1
#define BOT_OP 2


///////////////////////////
//
//		PATH HULL 

template<class KT,class InputIterator>
class PATH_HULL
{
  private:

	typedef typename KT::FT DataT;

	
	long 			top;			//Top end of elt.
	long			bot;			//Bottom end of elt.
	long			hp;				//The stack pointer.
	int*			op;				//History stack for operations.				dim=THRICE_HULL_MAX
	InputIterator*	elt;			//Double ended queue storing a convex hull. dim=TWICE_HULL_MAX
	InputIterator*	helt;			//History stack for points.					dim=THRICE_HULL_MAX

	long HULL_MAX;
	long TWICE_HULL_MAX;
	long THRICE_HULL_MAX;


  public:

	PATH_HULL()
	{
		op=NULL;
		elt=NULL;
		helt=NULL;
	}

	PATH_HULL(long dim_h1)
	{
		HULL_MAX=dim_h1;
		TWICE_HULL_MAX=2*HULL_MAX;
		THRICE_HULL_MAX=3*HULL_MAX;

		assert(op=new int[THRICE_HULL_MAX]);
		assert(elt=new InputIterator[TWICE_HULL_MAX]);
		assert(helt=new InputIterator[THRICE_HULL_MAX]);

		top=0;
		bot=0;
		hp=0;
		for (long i=0;i<THRICE_HULL_MAX;i++) 
			op[i]=-1;
	}

	~PATH_HULL()
	{
		delete[] op;
		delete[] elt;
		delete[] helt;
	}


	void push(InputIterator e)
	{
		top++;
		bot--;
		hp++;
		elt[top]=e;
		elt[bot]=e;
		helt[hp]=e;
		op[hp]=PUSH_OP;
	}


	void pop_Top()
	{
		hp++;
		helt[hp]=elt[top];
		top--;
		op[hp]=TOP_OP;
	}


	void pop_Bot()
	{
		hp++;
		helt[hp]=elt[bot];
		bot++;
		op[hp]=BOT_OP;
	}


	void hull_Init(InputIterator e1,InputIterator e2)
	{
		elt[HULL_MAX]=e1;
		top=HULL_MAX+1;
		elt[top]=e2;
		bot=HULL_MAX-1;
		elt[bot]=e2;
		hp=0;
		helt[hp]=e2; 
		op[0]=PUSH_OP;
	}


	bool SLOPE_SIGN(int p,int q,InputIterator begin,InputIterator end)  
	{ 
		DataT sHomogX,sHomogY,deltaXpq,deltaYpq,slopeVal;
 
		sHomogX=begin->y()-end->y();

		sHomogY=end->x()-begin->x();

		if(sHomogX==DataT(0) && sHomogY==DataT(0))
			sHomogX=DataT(1);

		deltaXpq=elt[q]->x()-elt[p]->x();
		deltaYpq=elt[q]->y()-elt[p]->y();
		slopeVal=sHomogX*deltaXpq+sHomogY*deltaYpq;
 
		return (slopeVal>=DataT(0) );
	}


typename DistTypeTraits<KT>::DataT squared_distance(InputIterator begin,InputIterator end,InputIterator cr)		
{
	typename DistTypeTraits<KT>::DataT d,dx,dy,tr_y;
	typename DistTypeTraits<KT>::DataT dbx,dby,dex,dey;
	
	if(end->x()==begin->x() && end->y()==begin->y())
		return (cr->x()-begin->x())*(cr->x()-begin->x())+(cr->y()-begin->y())*(cr->y()-begin->y());


	if((dx=end->x()-begin->x())*(dbx=cr->x()-begin->x())+(dy=end->y()-begin->y())*(dby=cr->y()-begin->y())<DataT(0))
	{
		return dbx*dbx+dby*dby;
	}
	else
		if(dy*(dey=cr->y()-end->y())+dx*(dex=cr->x()-end->x())>DataT(0))
		{
			return dex*dex+dey*dey;
		}
		else
		{
			d=dx*dx+dy*dy;
			tr_y=dby*dx-dbx*dy;
			return tr_y*tr_y/d;
		}
}


	InputIterator Find_Extreme(InputIterator begin,InputIterator end, typename DistTypeTraits<KT>::DataT *dist)
	{			
		long mid,lo,m1,brk,m2,hi;
		bool sbase, sbrk;
		typename DistTypeTraits<KT>::DataT d1, d2;
		InputIterator e;

		if((top-bot)>6)						//If there are > 6  points on the hull 
		{									//(otherwise we will just look at them all.
			lo=bot; 
			hi=top-1;	
			sbase = SLOPE_SIGN(hi,lo,begin, end); //The sign of the base edge.    
			
			do
			{ 
				brk = (lo + hi) / 2;		//Binary search for an edge with opposite sign.
				sbrk = SLOPE_SIGN(brk, brk+1, begin, end);
				
				if(sbase == sbrk )
					if(sbase == (SLOPE_SIGN(lo, brk+1, begin, end)) ) 
						lo = brk + 1;
					else 
						hi = brk;
			}
			while(sbase==sbrk);
  
			m1=brk;							//Now, the sign changes between the base edge and
			while(lo<m1)					// brk+1. Binary search for the extreme point.
			{
				mid=(lo+ m1)/2;
				
				if(sbase==(SLOPE_SIGN(mid,mid+1,begin,end)) ) 
					lo=mid+1;
				else 
					m1=mid;
			}  

			m2=brk;							//The sign also chagnes between brk and the base.
			while(m2<hi)					//Binary search again.
			{
				mid=(m2+hi)/2;
				
				if(sbase==SLOPE_SIGN(mid,mid+1,begin,end) ) 
					hi = mid;
				else 
					m2=mid+1;
			}					

			d1=squared_distance(begin,end,elt[lo]);
			d2=squared_distance(begin,end,elt[m2]);


			if(d1>d2)
			{ 
				e=elt[lo]; 
				*dist=d1; 
			}
			else
			{ 
				e=elt[m2]; 
				*dist=d2; 
			}
		}
		else 								//Few points in hull--search by brute force.
		{
			*dist=DataT(0);
			e=elt[bot];
			for(mid=bot;mid<top;mid++)
			{
				d1=squared_distance(begin,end,elt[mid]); 
				if(d1>*dist)
				{
					*dist=d1; 
					e=elt[mid];
				}
			}
		}
		
		return e;
	} 

	void Hull_Add(InputIterator p)
	{
		register bool topflag,botflag;
  
		topflag=((elt[top]->x() - p->x()) * (elt[top-1]->y() - p->y())) >= ((elt[top-1]->x() - p->x()) * (elt[top]->y() - p->y()));		// left_of
		botflag=((elt[bot+1]->x() - p->x()) * (elt[bot]->y() - p->y())) >= ((elt[bot]->x() - p->x()) * (elt[bot+1]->y() - p->y()));

		if(topflag || botflag)							//If the new point is outside the hull.
		{
			while(topflag && (top>HULL_MAX))			//Pop points until convexity is ensured.
			{
				pop_Top();
				topflag=((elt[top]->x() - p->x()) * (elt[top-1]->y() - p->y())) >= ((elt[top-1]->x() - p->x()) * (elt[top]->y() - p->y()));
			}
			
			while(botflag && (bot<HULL_MAX))
			{
				pop_Bot();
				botflag=((elt[bot+1]->x() - p->x()) * (elt[bot]->y() - p->y())) >= ((elt[bot]->x() - p->x()) * (elt[bot+1]->y() - p->y()));
			}
			
			push(p);									//Then push the new point on the top and bottom 
		}												//of the queue.  
	}  


	void Split(InputIterator e)
	{ 
		register InputIterator tmpe;
		register long tmpo;
														//Loop until we reach the PUSH_OP for e.
		while( (hp>= 0) && ((tmpo = op[hp]),((tmpe =helt[hp])!= e) || (tmpo != PUSH_OP)) )
		{
			hp--;
			
			switch (tmpo)								//Undo the operation.
			{
				case PUSH_OP:	
								top--; 
								bot++; 
								break;
				case TOP_OP:	
								top++; 
								elt[top]=tmpe; 
								break;
				case BOT_OP:
								bot--;
								elt[bot]=tmpe; 
								break;
			} 
		}
	}
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Douglas Peucker Hull Traits Class  -  Only for Recursive Split (recursive implementation), MIN E Method !!!
//
//	Split the HULL version

template<class KT>
class DouglasPeuckerHullTraits_Split: public DistTypeTraits<KT>
{
  protected:

	static long PHtag;
	static void *left;
	static void *right;
	static long dim_h;
	static bool first_st;
	
	long PHtag_r;
	void *left_r;
	void *right_r;
	DouglasPeuckerHullTraits_Split *ant;

	bool fSplitL;
	bool fFirstIt;

	template<class InputIterator>
	void Build(InputIterator begin,InputIterator end)
	{
		register InputIterator cr;
		InputIterator PHtag_it;
		
		long nn,N=0;

		for(cr=begin;cr!=end;cr++,N++);
	
		if(N<2)
			return;

		N=N/2;
		for(cr=begin,nn=0;nn<N;cr++,nn++);
/*
		if(left)
		{
			delete (PATH_HULL<KT,InputIterator>*)left;
			delete (PATH_HULL<KT,InputIterator>*)right;
		}
		left=new PATH_HULL<KT,InputIterator>(N+2); 
		assert(left);
		right=new PATH_HULL<KT,InputIterator>(N+2);
		assert(right);
*/

		if(!left)
		{
			left=new PATH_HULL<KT,InputIterator>(N+2); 
			assert(left);
			right=new PATH_HULL<KT,InputIterator>(N+2);
			assert(right);
		}

		PHtag=N;

		PHtag_it=cr;
		cr--;

		((PATH_HULL<KT,InputIterator>*)left)->hull_Init(PHtag_it,cr);	 			//Build left hull.
	
		if(cr!=begin)
		{
			for(cr--; cr!=begin; cr--)
				((PATH_HULL<KT,InputIterator>*)left)->Hull_Add(cr);  
			((PATH_HULL<KT,InputIterator>*)left)->Hull_Add(cr);  
		}

		cr=PHtag_it;
		cr++;
		((PATH_HULL<KT,InputIterator>*)right)->hull_Init(PHtag_it,cr); 			//Build right hull.

		if(cr!=end)
		{
			for(cr++;cr!=end;cr++)
				((PATH_HULL<KT,InputIterator>*)right)->Hull_Add(cr);  
			((PATH_HULL<KT,InputIterator>*)right)->Hull_Add(cr);  
		}
	}

  public:

	bool preproc;
	int fl_proc2;


	DouglasPeuckerHullTraits_Split()
	{
		preproc=true;
		fFirstIt=false;
		first_st=true;
		fl_proc2=0;
		left_r=right_r=NULL;
	}

	DouglasPeuckerHullTraits_Split(DouglasPeuckerHullTraits_Split& t)
	{
		preproc=true;
		fFirstIt=false;
		left_r=right_r=NULL;

		t.fl_proc2++;
		fl_proc2=0;

		ant=&t;
	}

	virtual ~DouglasPeuckerHullTraits_Split()
	{
		if(fFirstIt)
		{
//			delete left;
//			delete right;
			left=right=NULL;
			first_st=true;
		}

	}

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator *p_MAX=NULL)
	{
		InputIterator lextr,rextr;
		typename DistTypeTraits<KT>::DataT ldist, rdist;
		InputIterator cr;


		if(first_st)
		{
			first_st=false;
			InputIterator cr;
			ant->fl_proc2=0;
			ant->fSplitL=0;

			left=right=NULL;
			Build(begin,end);
		}

/////////////////////////////////////////////////////// Preproc2
		if(ant->fl_proc2==1)
			if(ant->fSplitL)
			{
				left=right=NULL;
				Build(begin,end);
			}

		if(ant->fl_proc2==2)
		{
			if(ant->fSplitL)
			{
				if(left)
				{
					delete (PATH_HULL<KT,InputIterator>*)left;
					delete (PATH_HULL<KT,InputIterator>*)right;
					left=right=NULL;
				}
				left=ant->left_r;
				right=ant->right_r;
				ant->left_r=ant->right_r=NULL;
				PHtag=ant->PHtag_r;

				((PATH_HULL<KT,InputIterator>*)left)->Split(begin);
			}
			else
				Build(begin,end);
		}
													
////////////////////////////////////////////////////// Current error computation
		cr=begin;
		cr++;
		if (begin== end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			rdist=DataT(0);
		}
		else
		{
			lextr=((PATH_HULL<KT,InputIterator>*)left)->Find_Extreme(begin,end,&ldist);
			rextr=((PATH_HULL<KT,InputIterator>*)right)->Find_Extreme(begin,end,&rdist);

			if (ldist>rdist)
			{
				if(p_MAX)
					*p_MAX=lextr;
				rdist=ldist;
			}
			else
			{
   				if(p_MAX)
					*p_MAX=rextr;
			}
		}
	
//////////////////////////////////////////////////////////// Preproc 1
		long n;

		fSplitL=false;

		for(cr=begin,n=0;n!=PHtag;cr++,n++)
			if(cr==*p_MAX)
			{
				fSplitL=true;
				break;
			}

		if(fSplitL)
		{
			left_r=left;
			right_r=right;
			PHtag_r=PHtag-n;
		}
		else
			if(cr==*p_MAX)
				Build(begin,*p_MAX);
			else
				((PATH_HULL<KT,InputIterator>*)right)->Split(*p_MAX); 
	
		return rdist;	
	}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Douglas Peucker Hull Traits Class  -  Only for Recursive Split, MIN E Method !!!
//
//	Always REBUILD the HULL version

template<class KT>
class DouglasPeuckerHullTraits_Build: public DouglasPeuckerHullTraits_Split<KT>
{
  public:

	template<class InputIterator>
	typename DistTypeTraits<KT>::DataT operator()(InputIterator begin,InputIterator end,InputIterator *p_MAX=NULL)
	{
		InputIterator lextr,rextr;
		typename DistTypeTraits<KT>::DataT ldist, rdist;
		InputIterator cr;

		if(first_st)
		{
			first_st=false;
			InputIterator cr;
			ant->fl_proc2=0;
			fFirstIt=true;

			left=right=NULL;
			Build(begin,end);
		}

///////////////////////////////////////////////////////// Preproc2
		if(ant->fl_proc2==2)
			Build(begin,end);
													
/////////////////////////////////////////////////////// Current error computation
		cr=begin;
		cr++;
		if (begin== end || cr==end)					// Error=0, curve segment contain 1 or 2 points
		{
			if(p_MAX)
				*p_MAX=begin;
			rdist=DataT(0);
		}
		else
		{
			lextr=((PATH_HULL<KT,InputIterator>*)left)->Find_Extreme(begin,end,&ldist);
			rextr=((PATH_HULL<KT,InputIterator>*)right)->Find_Extreme(begin,end,&rdist);

			if (ldist>rdist)
			{
				if(p_MAX)
					*p_MAX=lextr;
				rdist=ldist;
			}
			else
			{
   				if(p_MAX)
					*p_MAX=rextr;

			}
		}

//////////////////////////////////////////////////////////// Preproc 1
			
		Build(begin,*p_MAX);
		
		return rdist;	
	
	}

};


template<class KT>
long DouglasPeuckerHullTraits_Split<KT>::PHtag;

template<class KT>
void* DouglasPeuckerHullTraits_Split<KT>::left;

template<class KT>
void* DouglasPeuckerHullTraits_Split<KT>::right;

template<class KT>
long DouglasPeuckerHullTraits_Split<KT>::dim_h;

template<class KT>
bool DouglasPeuckerHullTraits_Split<KT>::first_st;



CGAL_END_NAMESPACE

#endif

