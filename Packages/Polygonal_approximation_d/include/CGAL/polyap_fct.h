
#ifndef CGAL_POLYAP_FCT
#define CGAL_POLYAP_FCT


#include <assert.h>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE


//////////////////////////////////////////////////////////////////////////////////////////
//
//			Dynamic Programming

//  Bounded-# (min-E) local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_DP_bnp_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i;
	
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;
	
	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_i=begin;c_i!=beyond;c_i++)
			*result++=*c_i;

		return result;
	}

	typename std::size_t N=n+1;
	
	FT **ApproxError;
	typename std::size_t **SplitPoints;
	FT **Err;

	ApproxError=new FT*[N];
	assert(ApproxError);
	for(m=0;m<N;m++)
	{
		ApproxError[m]=new FT[n_pt_bound];
		assert(ApproxError[m]);
	}

	SplitPoints=new typename std::size_t*[N];
	assert(SplitPoints);
	for(m=0;m<N;m++)
	{
		SplitPoints[m]=new typename std::size_t[n_pt_bound];
		assert(SplitPoints[m]);
	}

	Err=new FT*[N];
	assert(Err);
	for(m=0;m<N;m++)
	{
		Err[m]=new FT[N];
		assert(Err[m]);
	}

	Err[0][0]=FT(0);
	Err[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=FT(0);
	
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Err[m][n]=Err[n][m]=error(c_n,c_m);

	if(n_pt_bound==2)
	{
		approx_error=Err[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=1;n<N;n++)
		ApproxError[n][1]=Err[n][0];

	FT err_crt;
	
	for(m=2;m<n_pt_bound;m++)
	{
		ApproxError[m][m]=FT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				if(ApproxError[i][m-1]>Err[i][n])
					err_crt=ApproxError[i][m-1]; 
				else
					err_crt=Err[i][n];

				if(err_crt<ApproxError[n][m])
				{
					ApproxError[n][m]=err_crt;
					SplitPoints[n][m]=i;
				}
			}
			err_crt=Err[m-1][n];
			if(err_crt<ApproxError[n][m])
			{
				ApproxError[n][m]=err_crt;
				SplitPoints[n][m]=m-1;
			}
			
		}
	}
	
	approx_error=ApproxError[N-1][n_pt_bound-1];

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;

	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=n_pt_bound-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;
		nn=n;

		local_container.push_front(c_i);
	}

	local_container.push_front(begin);
	
	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	for(m=0;m<N;m++)
		delete[] ApproxError[m];
	delete[] ApproxError;

	for(m=0;m<N;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Err[m];
	delete[] Err;

	return result;
};


//  Bounded-E (min-#) local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_DP_be_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i;
	FT err_crt;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;
	
	FT **ApproxError;
	typename std::size_t **SplitPoints;
	FT **Err;

	ApproxError=new FT*[N];
	assert(ApproxError);
	for(m=0;m<N;m++)
	{
		ApproxError[m]=new FT[N];
		assert(ApproxError[m]);
	}

	Err=new FT*[N];
	assert(Err);
	for(m=0;m<N;m++)
	{
		Err[m]=new FT[N];
		assert(Err[m]);
	}

	Err[0][0]=Err[0][1]=Err[1][0]=FT(0);

	if(N>2)
		for(n=2,c_n=begin,c_n++,c_n++,c_n++;n<N;n++,c_n++)
			ApproxError[n][1]=Err[n][0]=Err[0][n]=error(begin,c_n);

	if(ApproxError[N-1][1]<=error_bound)
	{
		approx_n_pt=2;
	
		c_i=beyond;
		c_i--;

		*result++=*begin;
		*result=*c_i;
	
		for(m=0;m<N;m++)
			delete[] ApproxError[m];
		delete[] ApproxError;

		for(m=0;m<N;m++)
			delete[] Err[m];
		delete[] Err;

		return result;
	}

	SplitPoints=new typename std::size_t*[N];
	assert(SplitPoints);
	for(m=0;m<N;m++)
	{
		SplitPoints[m]=new typename std::size_t[N];
		assert(SplitPoints[m]);
	}

	Err[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=FT(0);

	for(n=1,c_n=begin,c_n++;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Err[m][n]=Err[n][m]=error(c_n,c_m);

	for(m=2;m<N;m++)
	{
		ApproxError[m][m]=FT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				if(ApproxError[i][m-1]>Err[i][n])
					err_crt=ApproxError[i][m-1]; 
				else
					err_crt=Err[i][n];

				if(err_crt<ApproxError[n][m])
				{
					ApproxError[n][m]=err_crt;
					SplitPoints[n][m]=i;
				}
			}
			err_crt=Err[m-1][n];
			if(err_crt<ApproxError[n][m])
			{
				ApproxError[n][m]=err_crt;
				SplitPoints[n][m]=m-1;
			}
		}

		if(ApproxError[N-1][m]<=error_bound)
		{
			approx_n_pt=m+1;
			break;
		}
	}
	
	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=approx_n_pt-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;
	
	for(m=0;m<N;m++)
		delete[] ApproxError[m];
	delete[] ApproxError;

	for(m=0;m<N;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Err[m];
	delete[] Err;
 
	return result;
};


//  Bounded-# (min-E) global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_DP_bnp_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i;
	
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_i=begin;c_i!=beyond;c_i++)
			*result++=*c_i;

		return result;
	}

	typename std::size_t N=n+1;
	
	FT **ApproxError;
	typename std::size_t **SplitPoints;
	FT **Err;

	ApproxError=new FT*[N];
	assert(ApproxError);
	for(m=0;m<N;m++)
	{
		ApproxError[m]=new FT[n_pt_bound];
		assert(ApproxError[m]);
	}

	SplitPoints=new typename std::size_t*[N];
	assert(SplitPoints);
	for(m=0;m<N;m++)
	{
		SplitPoints[m]=new typename std::size_t[n_pt_bound];
		assert(SplitPoints[m]);
	}

	Err=new FT*[N];
	assert(Err);
	for(m=0;m<N;m++)
	{
		Err[m]=new FT[N];
		assert(Err[m]);
	}

	Err[0][0]=FT(0);
	Err[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=FT(0);
	
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Err[m][n]=Err[n][m]=error(c_n,c_m);;

	if(n_pt_bound==2)
	{
		approx_error=Err[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=1;n<N;n++)
		ApproxError[n][1]=Err[n][0];

	FT err_crt;
	
	for(m=2;m<n_pt_bound;m++)
	{
		ApproxError[m][m]=FT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				err_crt=ApproxError[i][m-1]+Err[i][n]; 
				if(err_crt<ApproxError[n][m])
				{
					ApproxError[n][m]=err_crt;
					SplitPoints[n][m]=i;
				}
			}
			err_crt=Err[m-1][n];
			if(err_crt<ApproxError[n][m])
			{
				ApproxError[n][m]=err_crt;
				SplitPoints[n][m]=m-1;
			}
		}
	}
	
	approx_error=ApproxError[N-1][n_pt_bound-1];

	c_i=beyond;
	c_i--;


	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;

	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=n_pt_bound-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;
		nn=n;

		local_container.push_front(c_i);
	}

	local_container.push_front(begin);
	
	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	for(m=0;m<N;m++)
		delete[] ApproxError[m];
	delete[] ApproxError;

	for(m=0;m<N;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Err[m];
	delete[] Err;

	return result;
};



//  Bounded-E (min-#) global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_DP_be_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i;
	FT err_crt;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;
	
	FT **ApproxError;
	typename std::size_t **SplitPoints;
	FT **Err;

	ApproxError=new FT*[N];
	assert(ApproxError);
	for(m=0;m<N;m++)
	{
		ApproxError[m]=new FT[N];
		assert(ApproxError[m]);
	}

	Err=new FT*[N];
	assert(Err);
	for(m=0;m<N;m++)
	{
		Err[m]=new FT[N];
		assert(Err[m]);
	}

	Err[0][0]=Err[0][1]=Err[1][0]=FT(0);

	if(N>2)
		for(n=2,c_n=begin,c_n++,c_n++,c_n++;n<N;n++,c_n++)
			ApproxError[n][1]=Err[n][0]=Err[0][n]=error(begin,c_n);

	if(ApproxError[N-1][1]<=error_bound)
	{
		approx_n_pt=2;
	
		c_i=beyond;
		c_i--;

		*result++=*begin;
		*result=*c_i;
	
		for(m=0;m<N;m++)
			delete[] ApproxError[m];
		delete[] ApproxError;

		for(m=0;m<N;m++)
			delete[] Err[m];
		delete[] Err;

		return result;
	}

	SplitPoints=new typename std::size_t*[N];
	assert(SplitPoints);
	for(m=0;m<N;m++)
	{
		SplitPoints[m]=new typename std::size_t[N];
		assert(SplitPoints[m]);
	}

	Err[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=FT(0);

	for(n=1,c_n=begin,c_n++;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Err[m][n]=Err[n][m]=error(c_n,c_m);

	for(m=2;m<N;m++)
	{
		ApproxError[m][m]=FT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				err_crt=ApproxError[i][m-1]+Err[i][n]; 
				if(err_crt<ApproxError[n][m])
				{
					ApproxError[n][m]=err_crt;
					SplitPoints[n][m]=i;
				}
			}
			err_crt=Err[m-1][n];
			if(err_crt<ApproxError[n][m])
			{
				ApproxError[n][m]=err_crt;
				SplitPoints[n][m]=m-1;
			}
		}

		if(ApproxError[N-1][m]<=error_bound)
		{
			approx_n_pt=m+1;
			break;
		}
	}
	
	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=approx_n_pt-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;
	
	for(m=0;m<N;m++)
		delete[] ApproxError[m];
	delete[] ApproxError;

	for(m=0;m<N;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Err[m];
	delete[] Err;
 
	return result;
};


//////////////////////////////////////////////////////////////////////////////////////////
//
//				Recursive Split

//  Bounded-E local error assessment recursive implementation

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_RSr_be_lea(	InputIterator begin, 
													InputIterator beyond,
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	static bool fl;
	FT crt_err;
	InputIterator split_pt,temp_it;

	if(fl==0)
	{
		if(begin==beyond)
		{
			approx_n_pt=0;
			return result;
		}

		temp_it=begin;
		temp_it++;
		if(temp_it==beyond)
		{
			approx_n_pt=1;
			*result++=*begin;
			return result;
		}
		
		temp_it++;
		if(temp_it==beyond)
		{
			approx_n_pt=2;
			*result++=*begin;
			temp_it--;
			*result++=*temp_it;
			return result;
		}

		fl=1;
		approx_n_pt=1;
		*result++=*begin;		// save begin as the first point of the result
	}

	crt_err=error(begin,beyond,split_pt);
	
	if(crt_err>error_bound)
	{
		temp_it=split_pt;
		temp_it++;
		polygonal_approximation_RSr_be_lea(begin,temp_it,approx_n_pt,error_bound,result,error);
		polygonal_approximation_RSr_be_lea(split_pt,beyond,approx_n_pt,error_bound,result,error);
	}	
	else
	{
		temp_it=beyond;
		temp_it--;
		*result++=*temp_it;
		approx_n_pt++;
	}

	return result;
}


template<class InputIterator>
struct Sp_List
{	InputIterator begin;
	InputIterator beyond;
	InputIterator split_pt;
	Sp_List *next;
};	

template<class InputIterator,class FTp>		
struct WaitList
{	FTp err;
	Sp_List<InputIterator> *adr_el;
	WaitList *next;
};


//  Bounded-E local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_RS_be_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	InputIterator split_pt,c_n;
	FT crt_err;

	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	l_sp=new Sp_List<InputIterator>;
	assert(l_sp);

	l_sp->next=NULL;
	l_sp->begin=begin;
	l_sp->beyond=beyond;

	crt_err=error(begin,beyond,split_pt);

	l_sp->split_pt=split_pt;

	sp_cr=l_sp;

	approx_n_pt=2;

	while(crt_err>error_bound)
	{
		sp_new=new Sp_List<InputIterator>;
		assert(sp_new);

		approx_n_pt++;

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->beyond=sp_cr->beyond;
		sp_new->begin=sp_cr->split_pt;
		sp_cr->beyond=sp_cr->split_pt;
		sp_cr->beyond++;

		crt_err=error(sp_cr->begin,sp_cr->beyond,split_pt);

		while(sp_cr->next && crt_err<=error_bound)
		{
			sp_cr=sp_cr->next;
			crt_err=error(sp_cr->begin,sp_cr->beyond,split_pt);
		}
		
		sp_cr->split_pt=split_pt;
	}

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		l_sp->beyond--;
		*result++=*l_sp->beyond;
		l_sp=l_sp->next;
		delete sp_cr;
	}

	return result;
};


//  Bounded-# local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_RS_bnp_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	InputIterator c_n;
	typename std::size_t n=0;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	if(n_pt_bound==2)
	{
		approx_error=error(begin,beyond,c_n);
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}


	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	WaitList<InputIterator,FT> *l_wt,*wt_cr,*wt_new;
	InputIterator split_pt;

	l_sp=new Sp_List<InputIterator>;
	assert(l_sp);
	l_wt=new WaitList<InputIterator,FT>;
	assert(l_wt);

	l_sp->next=NULL;
	l_sp->begin=begin;
	l_sp->beyond=beyond;

	l_wt->next=NULL;
	l_wt->adr_el=l_sp;
	l_wt->err=error(begin,beyond,split_pt);

	l_sp->split_pt=split_pt;

	typename std::size_t m=2;

	while(m<n_pt_bound)
	{
		sp_cr=l_wt->adr_el;
		
		sp_new=new Sp_List<InputIterator>;
		assert(sp_new);

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->beyond=sp_cr->beyond;
		sp_new->begin=sp_cr->split_pt;
		sp_cr->beyond=sp_cr->split_pt;
		sp_cr->beyond++;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_cr->begin,sp_cr->beyond,split_pt);
		sp_cr->split_pt=split_pt;
		wt_new->adr_el=sp_cr;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_new->begin,sp_new->beyond,split_pt);
		wt_new->adr_el=sp_new;
		sp_new->split_pt=split_pt;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;

		m++;
	}

	approx_error=l_wt->err;

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		l_sp->beyond--;
		*result++=*l_sp->beyond;
		l_sp=l_sp->next;
		delete sp_cr;
	}

	while(l_wt!=NULL)
	{
		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;
	}

	return result;
};


//  Bounded-E global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_RS_be_gea(	InputIterator begin, 
													InputIterator beyond,
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	WaitList<InputIterator,FT> *l_wt,*wt_cr,*wt_new;
	InputIterator split_pt,c_n;
	FT crt_err;

	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	l_sp=new Sp_List<InputIterator>;
	assert(l_sp);
	l_wt=new WaitList<InputIterator,FT>;
	assert(l_wt);

	l_sp->next=NULL;
	l_sp->begin=begin;
	l_sp->beyond=beyond;

	l_wt->next=NULL;
	l_wt->adr_el=l_sp;
	l_wt->err=error(begin,beyond,split_pt);

	l_sp->split_pt=split_pt;

	typename std::size_t m=2;

	crt_err=l_wt->err;

	while(crt_err>error_bound)
	{
		sp_cr=l_wt->adr_el;

		crt_err-=l_wt->err; 
		
		sp_new=new Sp_List<InputIterator>;
		assert(sp_new);

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->beyond=sp_cr->beyond;
		sp_new->begin=sp_cr->split_pt;
		sp_cr->beyond=sp_cr->split_pt;
		sp_cr->beyond++;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_cr->begin,sp_cr->beyond,split_pt);
		crt_err+=wt_new->err;
		sp_cr->split_pt=split_pt;
		wt_new->adr_el=sp_cr;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_new->begin,sp_new->beyond,split_pt);
		crt_err+=wt_new->err;
		wt_new->adr_el=sp_new;
		sp_new->split_pt=split_pt;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;

		m++;
	}

	approx_n_pt=m;

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		l_sp->beyond--;
		*result++=*l_sp->beyond;
		l_sp=l_sp->next;
		delete sp_cr;
	}

	while(l_wt!=NULL)
	{
		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;
	}

	return result;
};


//  Bounded-# global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_RS_bnp_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	InputIterator c_n;
	typename std::size_t n=0;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	if(n_pt_bound==2)
	{
		approx_error=error(begin,beyond,c_n);
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}


	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	WaitList<InputIterator,FT> *l_wt,*wt_cr,*wt_new;
	InputIterator split_pt;
	FT crt_err; 

	l_sp=new Sp_List<InputIterator>;
	assert(l_sp);
	l_wt=new WaitList<InputIterator,FT>;
	assert(l_wt);

	l_sp->next=NULL;
	l_sp->begin=begin;
	l_sp->beyond=beyond;

	l_wt->next=NULL;
	l_wt->adr_el=l_sp;
	l_wt->err=error(begin,beyond,split_pt);

	l_sp->split_pt=split_pt;

	typename std::size_t m=2;

	crt_err=l_wt->err; 

	while(m<n_pt_bound)
	{
		sp_cr=l_wt->adr_el;

		crt_err-=l_wt->err;
		
		sp_new=new Sp_List<InputIterator>;
		assert(sp_new);

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->beyond=sp_cr->beyond;
		sp_new->begin=sp_cr->split_pt;
		sp_cr->beyond=sp_cr->split_pt;
		sp_cr->beyond++;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_cr->begin,sp_cr->beyond,split_pt);
		crt_err+=wt_new->err;
		sp_cr->split_pt=split_pt;
		wt_new->adr_el=sp_cr;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_new=new WaitList<InputIterator,FT>;
		assert(wt_new);

		wt_new->err=error(sp_new->begin,sp_new->beyond,split_pt);
		crt_err+=wt_new->err;
		wt_new->adr_el=sp_new;
		sp_new->split_pt=split_pt;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;

		m++;
	}

	approx_error=crt_err;

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		l_sp->beyond--;
		*result++=*l_sp->beyond;
		l_sp=l_sp->next;
		delete sp_cr;
	}

	while(l_wt!=NULL)
	{
		wt_cr=l_wt;
		l_wt=l_wt->next;
		delete wt_cr;
	}

	return result;
};


///////////////////////////////////////////////////////////////////////////////////
//
//		Graph Search

//  Bounded-E local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_GS_be_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t path;
	typename std::size_t n,m,m_pt,i,n_vis,m_crt;
	typename DistTraits::FT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	
	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;

	typedef struct cel
	{	typename std::size_t pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

	m_crt=N+1;
	err_crt=error_bound+1;
	
	FT **Graph;
	CEL **visited;

	visited=new CEL*[N];
	assert(visited);
	for(m=0;m<N;m++)
		visited[m]=NULL;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);
	Graph[N-1][N-1]=FT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Graph[n][m]=Graph[m][n]=error(c_n,c_m);

	if(Graph[0][N-1]<=error_bound)
	{
		approx_n_pt=2;

		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(Graph[0][N-1]<=error_bound)
		Graph[0][N-1]=Graph[N-1][0]=error_bound+1;
	
	queue_b=queue_e=new CEL;
	assert(queue_b);
	queue_b->next=NULL;
	queue_b->pt_crt=0;

	breadth_crt=new CEL;
	assert(breadth_crt);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++)
		{
				if(Graph[queue_b->pt_crt][i]<=error_bound)
				{
					if(i==N-1)
					{
						max=Graph[queue_b->pt_crt][N-1];
						breadth_crt=visited[queue_b->pt_crt];
						m_pt=2;
						while(breadth_crt->next)
						{
							m_pt++;
							vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
							if(vm>max)
								max=vm;
							breadth_crt=breadth_crt->next;
						}
						if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
						{
							m_crt=m_pt;
							err_max=max;
							path=queue_b->pt_crt;
						}
					}

					if(visited[i]==0)
					{
						breadth_crt=new CEL;
						assert(breadth_crt);

						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;

						queue_e->next=new CEL;
						assert(queue_e->next);

						queue_e=queue_e->next;
						queue_e->next=NULL;
						queue_e->pt_crt=i;
					}
				}
		}

		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	while(queue_b)
	{
		if(queue_b->pt_crt!=N-1)
			if(Graph[queue_b->pt_crt][N-1]<=error_bound)
			{
				max=Graph[queue_b->pt_crt][N-1];
				breadth_crt=visited[queue_b->pt_crt];
				m_pt=2;
				while(breadth_crt->next)
				{
					m_pt++;
					vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
					if(vm>max)
						max=vm;
					breadth_crt=breadth_crt->next;
				}
				if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
				{
					m_crt=m_pt;
					err_max=max;
					path=queue_b->pt_crt;
				}
			}
		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	approx_n_pt=m_crt;

	delete visited[N-1];
	visited[N-1]=NULL;

	breadth_crt=visited[path]->next;

	delete visited[path];
	visited[path]=NULL;

	while(breadth_crt)
	{
		path=breadth_crt->pt_crt;
		breadth_crt=breadth_crt->next;

		delete visited[path];
		visited[path]=NULL;
	}

	for(c_n=begin,i=0;c_n!=beyond;c_n++,i++)
		if(visited[i]==NULL)
			*result++=*c_n;

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	delete[] visited;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
};



template<class FT>
void QuickSort(FT* p,typename std::size_t n)  //quick sort
{
	FT tmp,pivot;

	if(n<2)
		return;

	if(n==2)
	{
		if(p[0]>p[1])
		{
			tmp=p[0];
			p[0]=p[1];
			p[1]=tmp;
		}
		return;
	}

    pivot=p[n/2];
	std::size_t i=0,j=n-1;

	while(i<j)
	{
		while(p[i]<pivot)
			i++;
		while(p[j]>pivot)
			j--;
		if(i<j)
		{
			tmp=p[i];
			p[i]=p[j];
			p[j]=tmp;
			if(p[i]==p[j])
			{
				i++;
				j--;
			}
		}
	}

	if(i==j)
	{	if(p[i]>pivot)
			i=i-1;
	}
	else if(j<i)
		i=j;
	QuickSort<FT>(p,i+1);
	j=n-i-1;
	QuickSort<FT>(&p[i+1],j);
}


//  Bounded-# local error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_GS_bnp_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{

	typedef typename DistTraits::FT FT;

	typename std::size_t path;
	typename std::size_t n,m,m_pt,i,n_vis,m_crt;
	typename DistTraits::FT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}

	typename std::size_t N=n+1;

	typedef struct cel
	{	typename std::size_t pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

	FT **Graph;
	CEL **visited;

	visited=new CEL*[N];
	assert(visited);
	for(m=0;m<N;m++)
		visited[m]=NULL;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);
	Graph[N-1][N-1]=FT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);

	typename std::size_t N_err=(N-1)*(N-2)/2;
	FT *List_err;
	List_err=new FT[N_err];
	assert(List_err);

	typename std::size_t idx=0;

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Graph[n][m]=Graph[m][n]=List_err[idx++]=error(c_n,c_m);

	if(n_pt_bound==2)
	{
		approx_error=Graph[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	FT cr_err_bound;
	bool fl_circ=false;

	typename std::size_t N_crt,N_mf;
	double delta_N,nc;
	
	if(Graph[0][N-1]<=0)
		fl_circ=true;

	QuickSort<FT>(List_err,N_err);

	delta_N=N_err/2.;
	nc=0;
	m=n_pt_bound+1;

	while(delta_N>0.4)
	{
		if(m>n_pt_bound)
			nc+=delta_N;
		else
			nc-=delta_N;

		N_crt=(int)(nc+0.5);
		if(N_crt<0)
			N_crt=0;
		if(N_crt>=N_err)
			N_crt=N_err;

		cr_err_bound=List_err[N_crt];
		m_crt=N+1;
		err_crt=cr_err_bound+1;

		if(fl_circ)									
			Graph[0][N-1]=Graph[N-1][0]=cr_err_bound+1;

		for(m=0;m<N;m++)
			visited[m]=NULL;
		
		queue_b=queue_e=new CEL;
		assert(queue_b);

		queue_b->next=NULL;
		queue_b->pt_crt=0;

		breadth_crt=new CEL;
		assert(breadth_crt);

		breadth_crt->pt_crt=0;
		breadth_crt->next=NULL;

		visited[0]=breadth_crt;
		n_vis=1;

		while(n_vis<N)
		{
			for(i=(queue_b->pt_crt+1);i<N;i++)
			{
					if(Graph[queue_b->pt_crt][i]<=cr_err_bound)
					{
						if(i==N-1)
						{
							max=Graph[queue_b->pt_crt][N-1];
							breadth_crt=visited[queue_b->pt_crt];
							m_pt=2;
							while(breadth_crt->next)
							{
								m_pt++;
								vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
								if(vm>max)
									max=vm;
								breadth_crt=breadth_crt->next;
							}
							if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
							{
								m_crt=m_pt;
								err_max=max;
								path=queue_b->pt_crt;
							}
						}

						if(visited[i]==0)
						{
							breadth_crt=new CEL;
							assert(breadth_crt);
							breadth_crt->pt_crt=i;
							breadth_crt->next=visited[queue_b->pt_crt];
							visited[i]=breadth_crt;
							n_vis++;
							queue_e->next=new CEL;
							assert(queue_e->next);
							queue_e=queue_e->next;
							queue_e->next=NULL;
							queue_e->pt_crt=i;
						}
					}
			}

			queue_crt=queue_b;
			queue_b=queue_b->next;
			delete queue_crt;
		}

		while(queue_b)
		{
			if(queue_b->pt_crt!=N-1)
				if(Graph[queue_b->pt_crt][N-1]<=cr_err_bound)
				{
					max=Graph[queue_b->pt_crt][N-1];
					breadth_crt=visited[queue_b->pt_crt];
					m_pt=2;
					while(breadth_crt->next)
					{
						m_pt++;
						vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
						if(vm>max)
							max=vm;
						breadth_crt=breadth_crt->next;
					}
					if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
					{
						m_crt=m_pt;
						err_max=max;
						path=queue_b->pt_crt;
					}
				}
			queue_crt=queue_b;
			queue_b=queue_b->next;
			delete queue_crt;
		}

		m=m_crt;

		delete visited[N-1];
		visited[N-1]=NULL;

		breadth_crt=visited[path]->next;

		delete visited[path];
		visited[path]=NULL;

		while(breadth_crt)
		{
			path=breadth_crt->pt_crt;
			breadth_crt=breadth_crt->next;

			delete visited[path];
			visited[path]=NULL;
		}

		for(i=0;i<N;i++)
			if(visited[i])
				delete visited[i];

		approx_error=cr_err_bound;
	
		delta_N/=2.;

		if(m<=n_pt_bound)
			N_mf=N_crt;
	}

	if(m!=n_pt_bound)
		N_crt=N_mf;
	
	cr_err_bound=List_err[N_crt];
	m_crt=N+1;
	err_crt=cr_err_bound+1;

	if(fl_circ)									
		Graph[0][N-1]=Graph[N-1][0]=cr_err_bound+1;

	for(m=0;m<N;m++)
		visited[m]=NULL;
		
	queue_b=queue_e=new CEL;
	assert(queue_b);

	queue_b->next=NULL;
	queue_b->pt_crt=0;

	breadth_crt=new CEL;
	assert(breadth_crt);

	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++) 
		{
				if(Graph[queue_b->pt_crt][i]<=cr_err_bound)
				{
					if(i==N-1)
					{
						max=Graph[queue_b->pt_crt][N-1];
						breadth_crt=visited[queue_b->pt_crt];
						m_pt=2;
						while(breadth_crt->next)
						{
							m_pt++;
							vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
							if(vm>max)
								max=vm;
							breadth_crt=breadth_crt->next;
						}
						if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
						{
							m_crt=m_pt;
							err_max=max;
							path=queue_b->pt_crt;
						}
					}

					if(visited[i]==0)
					{
						breadth_crt=new CEL;
						assert(breadth_crt);

						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;

						queue_e->next=new CEL;
						assert(queue_e->next);

						queue_e=queue_e->next;
						queue_e->next=NULL;
						queue_e->pt_crt=i;
					}
				}
		}

		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	while(queue_b)
	{
		if(queue_b->pt_crt!=N-1)
			if(Graph[queue_b->pt_crt][N-1]<=cr_err_bound)
			{
				max=Graph[queue_b->pt_crt][N-1];
				breadth_crt=visited[queue_b->pt_crt];
				m_pt=2;
				while(breadth_crt->next)
				{
					m_pt++;
					vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
					
					if(vm>max)
						max=vm;
					breadth_crt=breadth_crt->next;
				}
				if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
				{
					m_crt=m_pt;
					err_max=max;
					path=queue_b->pt_crt;
				}
			}
		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	m=m_crt;
	n_pt_bound=m;

	delete visited[N-1];
	visited[N-1]=NULL;

	breadth_crt=visited[path]->next;

	delete visited[path];
	visited[path]=NULL;

	while(breadth_crt)
	{
		path=breadth_crt->pt_crt;
		breadth_crt=breadth_crt->next;

		delete visited[path];
		visited[path]=NULL;
	}

	for(c_n=begin,i=0;c_n!=beyond;c_n++,i++)
		if(visited[i]==NULL)
			*result++=*c_n;

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	approx_error=cr_err_bound;

	delete[] visited;

	delete[] List_err;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;

	return result;
}


//  Bounded-E global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_GS_be_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													typename DistTraits::FT error_bound, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i,j;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;

	FT **Graph;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);
	Graph[N-1][N-1]=FT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Graph[n][m]=Graph[m][n]=error(c_n,c_m);

	if(Graph[0][N-1]<=error_bound)
	{
		approx_n_pt=2;

		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(Graph[0][N-1]<=error_bound)
		Graph[0][N-1]=Graph[N-1][0]=error_bound+1;
	
	bool fl_find=false;
	FT err_crt;
	FT* ApproxError;
	typename std::size_t **SplitPoints;

	SplitPoints=new typename std::size_t*[N+1];
	assert(SplitPoints);

	ApproxError=new FT[N+1];
	assert(ApproxError);

	m=2;

	for(i=0;i<=N-1;i++)
		ApproxError[i]=Graph[i][0];

	while(!fl_find)
	{
		m++;
		SplitPoints[m]=new typename std::size_t[N];
		assert(SplitPoints[m]);

		ApproxError[N-1]=ApproxError[N-2];
		SplitPoints[m][N-1]=N-2;
		if(ApproxError[N-1]<=error_bound)
			fl_find=true;
		else
			for(j=N-3;j>=m-2;j--)
			{
				if(ApproxError[j]<=error_bound && Graph[j][N-1]<=error_bound)
				{
					err_crt=ApproxError[j]+Graph[j][N-1];

					if(err_crt<ApproxError[N-1])
					{
						ApproxError[N-1]=err_crt;
						SplitPoints[m][N-1]=j;
						if(ApproxError[N-1]<=error_bound)
						{
							fl_find=true;
							break;
						}
					}
				}
			}

		if(!fl_find)
			for(i=N-2;i>=m-1;i--)
			{
				ApproxError[i]=ApproxError[i-1];
				SplitPoints[m][i]=i-1;
				for(j=i-2;j>=m-2;j--)
				{
					if(ApproxError[j]<=error_bound && Graph[j][i]<=error_bound)
					{
						err_crt=ApproxError[j]+Graph[j][i];
						if(err_crt<ApproxError[i])
						{
							ApproxError[i]=err_crt;
							SplitPoints[m][i]=j;
						}
					}
				}
			}
	}

	approx_n_pt=m;

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=approx_n_pt;m>2;m--)
	{
		n=SplitPoints[m][n];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	delete[] ApproxError;

	for(m=3;m<=approx_n_pt;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
};


//  Bounded-# global error assessment 

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator polygonal_approximation_GS_bnp_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t n_pt_bound, 
													typename DistTraits::FT &approx_error, 
													OutputIterator result,
													DistTraits error )
{
	typedef typename DistTraits::FT FT;

	typename std::size_t n,m,i,j;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}

	typename std::size_t N=n+1;

	FT **Graph;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);
	Graph[N-1][N-1]=FT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++,c_m++;m<N;m++,c_m++)
			Graph[n][m]=Graph[m][n]=error(c_n,c_m);

	if(n_pt_bound==2)
	{
		approx_error=Graph[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	typename std::size_t mm;

	FT err_crt;
	FT* ApproxError;
	typename std::size_t **SplitPoints;

	SplitPoints=new typename std::size_t*[N+1];
	assert(SplitPoints);

	ApproxError=new FT[N+1];
	assert(ApproxError);

	for(i=0;i<=N-1;i++)
		ApproxError[i]=Graph[i][0];

	FT epsi_crt;

	mm=2;

	while(mm<n_pt_bound)
	{
		mm++;
		SplitPoints[mm]=new typename std::size_t[N];
		assert(SplitPoints[mm]);

		epsi_crt=ApproxError[N-1]=ApproxError[N-2];
		SplitPoints[mm][N-1]=N-2;
		
		for(j=N-3;j>=mm-2;j--)
		{
			if(ApproxError[j]<=epsi_crt && Graph[j][N-1]<=epsi_crt)
			{
				err_crt=ApproxError[j]+Graph[j][N-1];

				if(err_crt<ApproxError[N-1])
				{
					epsi_crt=ApproxError[N-1]=err_crt;
					SplitPoints[mm][N-1]=j;
				}
			}
		}
	
		if(mm<n_pt_bound)
			for(i=N-2;i>=mm-1;i--)
			{
				epsi_crt=ApproxError[i]=ApproxError[i-1];
				SplitPoints[mm][i]=i-1;
				for(j=i-2;j>=mm-2;j--)
				{
					if(ApproxError[j]<=epsi_crt && Graph[j][i]<=epsi_crt)
					{
						err_crt=ApproxError[j]+Graph[j][i];
						if(err_crt<ApproxError[i])
						{
							epsi_crt=ApproxError[i]=err_crt;
							SplitPoints[mm][i]=j;
						}
					}
				}
			}
	}

	approx_error=ApproxError[N-1];

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=n_pt_bound;m>2;m--)
	{
		n=SplitPoints[m][n];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	delete[] ApproxError;

	for(m=3;m<=mm;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
}


///////////////////////////////////////////////////////////////////////////////////
//
//		Optimal algorithms for Squared Euclidean distance

template<class InputIterator, class FT>
bool slope_sign(InputIterator p,InputIterator q,InputIterator begin,InputIterator end)  
{ 
	FT sHomogX,sHomogY,deltaXpq,deltaYpq,slopeVal;
 
	sHomogX=begin->y()-end->y();
	sHomogY=end->x()-begin->x();

	if(sHomogX==FT(0) && sHomogY==FT(0))
		sHomogX=FT(1);

	deltaXpq=q->x()-p->x();
	deltaYpq=q->y()-p->y();
	slopeVal=sHomogX*deltaXpq+sHomogY*deltaYpq;
 
	return (slopeVal>=FT(0) );
}


//  Bounded-E local error assessment; Based on the online convex hull (Toussaint)

template<class InputIterator,class OutputIterator, class FT>
OutputIterator polygonal_approximation_CHGS_be_lea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													FT error_bound, 
													OutputIterator result )
{
	typename std::size_t path;

	typename std::size_t n,m,m_pt,i,n_vis,m_crt;
	FT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	
	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;


	// build the double stack

	typename std::size_t top;			//Top end of elt.
	typename std::size_t bot;			//Bottom end of elt.
	InputIterator*	elt;			//Double ended queue storing a convex hull. dim=TWICE_HULL_MAX

	typename std::size_t HULL_MAX;
	typename std::size_t TWICE_HULL_MAX;
	typename std::size_t THRICE_HULL_MAX;

	bool topflag,botflag;

	HULL_MAX=N;
	TWICE_HULL_MAX=2*HULL_MAX;
	THRICE_HULL_MAX=3*HULL_MAX;

	elt=new InputIterator[TWICE_HULL_MAX];
	assert(elt);

	for (typename std::size_t j=0;j<TWICE_HULL_MAX;j++)
		elt[j]=NULL;

	typedef struct cel
	{	typename std::size_t pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

	m_crt=N+1;
	err_crt=error_bound+1;
	
	FT **Graph;
	CEL **visited;

	visited=new CEL*[N];
	assert(visited);

	for(m=0;m<N;m++)
		visited[m]=NULL;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);

	Graph[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);
	
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		top=0;
		bot=0;
		c_i=c_n;
		c_i++;

		elt[HULL_MAX]=c_n;
		top=HULL_MAX+1;
		elt[top]=c_i;
		bot=HULL_MAX-1;
		elt[bot]=c_i;

		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{
			topflag=((elt[top]->x()-c_m->x())*(elt[top-1]->y()-c_m->y())) >= ((elt[top-1]->x()-c_m->x())*(elt[top]->y()-c_m->y()));		// left_of
			botflag=((elt[bot+1]->x()-c_m->x())*(elt[bot]->y()-c_m->y())) >= ((elt[bot]->x()-c_m->x())*(elt[bot+1]->y()-c_m->y()));

			if(topflag || botflag)							//If the new point is outside the hull.
			{
				while(topflag && (top>HULL_MAX))			//Pop points until convexity is ensured.
				{
					top--;
					topflag=((elt[top]->x()-c_m->x())*(elt[top-1]->y()-c_m->y())) >= ((elt[top-1]->x()-c_m->x())*(elt[top]->y()-c_m->y()));
				}
				
				while(botflag && (bot<HULL_MAX))
				{
					bot++;
					botflag=((elt[bot+1]->x()-c_m->x())*(elt[bot]->y()-c_m->y())) >= ((elt[bot]->x()-c_m->x())*(elt[bot+1]->y()-c_m->y()));
				}
			
				top++;					//Then push the new point on the top and bottom 
				bot--;					//of the queue. 
				elt[top]=c_m;
				elt[bot]=c_m;
			}												 

			// find the furthest point

			typename std::size_t mid,lo,m1,brk,m2,hi;
			bool sbase, sbrk;
			FT d1,d2,dx,dy,d,tr_y;

			if((top-bot)>6)						//If there are > 6  points on the hull 
			{									//(otherwise we will just look at them all.
				lo=bot; 
				hi=top-1;	
				sbase = slope_sign<InputIterator,FT>(elt[hi],elt[lo],c_n,c_m); //The sign of the base edge.    

				do
				{ 
					brk = (lo + hi) / 2;		//Binary search for an edge with opposite sign.
					sbrk = slope_sign<InputIterator,FT>(elt[brk], elt[brk+1], c_n,c_m);
				
					if(sbase == sbrk )
						if(sbase == (slope_sign<InputIterator,FT>(elt[lo], elt[brk+1], c_n,c_m)) ) 
							lo = brk + 1;
						else 
							hi = brk;
				}
				while(sbase==sbrk);
  
				m1=brk;							//Now, the sign changes between the base edge and
				while(lo<m1)					// brk+1. Binary search for the extreme point.
				{
					mid=(lo+ m1)/2;
				
					if(sbase==(slope_sign<InputIterator,FT>(elt[mid],elt[mid+1],c_n,c_m)) ) 
						lo=mid+1;
					else 
						m1=mid;
				}  

				m2=brk;							//The sign also chagnes between brk and the base.
				while(m2<hi)					//Binary search again.
				{
					mid=(m2+hi)/2;
				
					if(sbase==slope_sign<InputIterator,FT>(elt[mid],elt[mid+1],c_n,c_m) ) 
						hi = mid;
					else 
						m2=mid+1;
				}					

				dx=c_m->x()-c_n->x();
				dy=c_m->y()-c_n->y();
				d=dx*dx+dy*dy;
				if(d==0)
				{			
					d1=(elt[lo]->x()-c_n->x())*(elt[lo]->x()-c_n->x())+(elt[lo]->y()-c_n->y())*(elt[lo]->y()-c_n->y());
					d2=(elt[hi]->x()-c_n->x())*(elt[hi]->x()-c_n->x())+(elt[hi]->y()-c_n->y())*(elt[hi]->y()-c_n->y());
				}
				else
				{					
					//d1=squared_distance(begin,end,elt[lo]);
					tr_y=(elt[lo]->y()-c_n->y())*dx-(elt[lo]->x()-c_n->x())*dy;
					d1=tr_y*tr_y/d;

					//d2=squared_distance(begin,end,elt[hi]);
					tr_y=(elt[hi]->y()-c_n->y())*dx-(elt[hi]->x()-c_n->x())*dy;
					d2=tr_y*tr_y/d;
				}

				if(d1>d2)
					d2=d1; 
			}
			else 								//Few points in hull--search by brute force.
			{
				d2=FT(0);
				for(mid=bot;mid<top;mid++)
				{
					//d1=squared_distance(begin,end,elt[mid]); 
					dx=c_m->x()-c_n->x();
					dy=c_m->y()-c_n->y();
					d=dx*dx+dy*dy;
					if(d==0)
					{
						d1=(elt[mid]->x()-c_n->x())*(elt[mid]->x()-c_n->x())+(elt[mid]->y()-c_n->y())*(elt[mid]->y()-c_n->y());
					}
					else
					{
						tr_y=(elt[mid]->y()-c_n->y())*dx-(elt[mid]->x()-c_n->x())*dy;
						d1=tr_y*tr_y/d;
					}

					if(d1>d2)
						d2=d1; 
				}
			}
			Graph[n][m]=Graph[m][n]=d2;
		}
	}

	delete[] elt;
	elt=NULL;

	if(Graph[0][N-1]<=error_bound)
	{
		approx_n_pt=2;

		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(Graph[0][N-1]<=error_bound)
		Graph[0][N-1]=Graph[N-1][0]=error_bound+1;
	
	queue_b=queue_e=new CEL;
	assert(queue_b);

	queue_b->next=NULL;
	queue_b->pt_crt=0;

	breadth_crt=new CEL;
	assert(breadth_crt);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++)
		{
				if(Graph[queue_b->pt_crt][i]<=error_bound)
				{
					if(i==N-1)
					{
						max=Graph[queue_b->pt_crt][N-1];
						breadth_crt=visited[queue_b->pt_crt];
						m_pt=2;
						while(breadth_crt->next)
						{
							m_pt++;
							vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
							if(vm>max)
								max=vm;
							breadth_crt=breadth_crt->next;
						}
						if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
						{
							m_crt=m_pt;
							err_max=max;
							path=queue_b->pt_crt;
						}
					}

					if(visited[i]==0)
					{
						breadth_crt=new CEL;
						assert(breadth_crt);

						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;

						queue_e->next=new CEL;
						assert(queue_e->next);

						queue_e=queue_e->next;
						queue_e->next=NULL;
						queue_e->pt_crt=i;
					}
				}
		}

		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	while(queue_b)
	{
		if(queue_b->pt_crt!=N-1)
			if(Graph[queue_b->pt_crt][N-1]<=error_bound)
			{
				max=Graph[queue_b->pt_crt][N-1];
				breadth_crt=visited[queue_b->pt_crt];
				m_pt=2;
				while(breadth_crt->next)
				{
					m_pt++;
					vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
					if(vm>max)
						max=vm;
					breadth_crt=breadth_crt->next;
				}
				if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
				{
					m_crt=m_pt;
					err_max=max;
					path=queue_b->pt_crt;
				}
			}
		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	approx_n_pt=m_crt;

	delete visited[N-1];
	visited[N-1]=NULL;

	breadth_crt=visited[path]->next;

	delete visited[path];
	visited[path]=NULL;

	while(breadth_crt)
	{
		path=breadth_crt->pt_crt;
		breadth_crt=breadth_crt->next;

		delete visited[path];
		visited[path]=NULL;
	}

	for(c_n=begin,i=0;c_n!=beyond;c_n++,i++)
		if(visited[i]==NULL)
			*result++=*c_n;

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	delete[] visited;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
};


//  Bounded-# local error assessment; Based on the online convex hull (Toussaint)

template<class InputIterator,class OutputIterator, class FT>
OutputIterator polygonal_approximation_CHGS_bnp_lea(	InputIterator begin, 
														InputIterator beyond, 
														typename std::size_t n_pt_bound, 
														FT &approx_error, 
														OutputIterator result )
{
	typename std::size_t path;
	typename std::size_t n,m,i;
	typename std::size_t m_pt,n_vis,m_crt;
	FT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	FT cr_err_bound;
	bool fl_circ=false;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}


	typename std::size_t N=n+1;
	typename std::size_t N_err=(N-1)*(N-2)/2;
	typename std::size_t N_crt,N_mf;
	double delta_N,nc;


	// build the double stack

	typename std::size_t top;			//Top end of elt.
	typename std::size_t bot;			//Bottom end of elt.
	InputIterator*	elt;			//Double ended queue storing a convex hull. dim=TWICE_HULL_MAX

	typename std::size_t HULL_MAX;
	typename std::size_t TWICE_HULL_MAX;
	typename std::size_t THRICE_HULL_MAX;

	bool topflag,botflag;

	HULL_MAX=N;
	TWICE_HULL_MAX=2*HULL_MAX;
	THRICE_HULL_MAX=3*HULL_MAX;

	elt=new InputIterator[TWICE_HULL_MAX];
	assert(elt);

	for (typename std::size_t j=0;j<TWICE_HULL_MAX;j++)
		elt[j]=NULL;
	
	typedef struct cel
	{	typename std::size_t pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;
	
	CEL **visited;

	visited=new CEL*[N];
	assert(visited);

	for(m=0;m<N;m++)
		visited[m]=NULL;

	FT *List_err;
	FT **Graph;

	List_err=new FT[N_err];
	assert(List_err);

	Graph=new FT*[N];
	assert(Graph);
    for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);
	Graph[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);

	typename std::size_t index1=0;

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		top=0;
		bot=0;
		c_i=c_n;
		c_i++;

		elt[HULL_MAX]=c_n;
		top=HULL_MAX+1;
		elt[top]=c_i;
		bot=HULL_MAX-1;
		elt[bot]=c_i;

		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{
			topflag=((elt[top]->x()-c_m->x())*(elt[top-1]->y()-c_m->y())) >= ((elt[top-1]->x()-c_m->x())*(elt[top]->y()-c_m->y()));		// left_of
			botflag=((elt[bot+1]->x()-c_m->x())*(elt[bot]->y()-c_m->y())) >= ((elt[bot]->x()-c_m->x())*(elt[bot+1]->y()-c_m->y()));

			if(topflag || botflag)							//If the new point is outside the hull.
			{
				while(topflag && (top>HULL_MAX))			//Pop points until convexity is ensured.
				{
					top--;
					topflag=((elt[top]->x()-c_m->x())*(elt[top-1]->y()-c_m->y())) >= ((elt[top-1]->x()-c_m->x())*(elt[top]->y()-c_m->y()));
				}
				
				while(botflag && (bot<HULL_MAX))
				{
					bot++;
					botflag=((elt[bot+1]->x()-c_m->x())*(elt[bot]->y()-c_m->y())) >= ((elt[bot]->x()-c_m->x())*(elt[bot+1]->y()-c_m->y()));
				}
			
				top++;					//Then push the new point on the top and bottom 
				bot--;					//of the queue. 
				elt[top]=c_m;
				elt[bot]=c_m;
			}												 

			// find the furthest point

			typename std::size_t mid,lo,m1,brk,m2,hi;
			bool sbase, sbrk;
			FT d1,d2,dx,dy,d,tr_y;

			if((top-bot)>6)						//If there are > 6  points on the hull 
			{									//(otherwise we will just look at them all.
				lo=bot; 
				hi=top-1;	
				sbase = slope_sign<InputIterator,FT>(elt[hi],elt[lo],c_n,c_m); //The sign of the base edge.    
			
				do
				{ 
					brk = (lo + hi) / 2;		//Binary search for an edge with opposite sign.
					sbrk = slope_sign<InputIterator,FT>(elt[brk], elt[brk+1], c_n,c_m);
				
					if(sbase == sbrk )
						if(sbase == (slope_sign<InputIterator,FT>(elt[lo], elt[brk+1], c_n,c_m)) ) 
							lo = brk + 1;
						else 
							hi = brk;
				}
				while(sbase==sbrk);
  
				m1=brk;							//Now, the sign changes between the base edge and
				while(lo<m1)					// brk+1. Binary search for the extreme point.
				{
					mid=(lo+ m1)/2;
				
					if(sbase==(slope_sign<InputIterator,FT>(elt[mid],elt[mid+1],c_n,c_m)) ) 
						lo=mid+1;
					else 
						m1=mid;
				}  

				m2=brk;							//The sign also chagnes between brk and the base.
				while(m2<hi)					//Binary search again.
				{
					mid=(m2+hi)/2;
				
					if(sbase==slope_sign<InputIterator,FT>(elt[mid],elt[mid+1],c_n,c_m) ) 
						hi = mid;
					else 
						m2=mid+1;
				}					

				//d1=squared_distance(begin,end,elt[lo]);
				dx=c_m->x()-c_n->x();
				dy=c_m->y()-c_n->y();
				d=dx*dx+dy*dy;
				if(d==0)
				{
					d1=(elt[lo]->x()-c_n->x())*(elt[lo]->x()-c_n->x())+(elt[lo]->y()-c_n->y())*(elt[lo]->y()-c_n->y());
					d2=(elt[hi]->x()-c_n->x())*(elt[hi]->x()-c_n->x())+(elt[hi]->y()-c_n->y())*(elt[hi]->y()-c_n->y());
				}
				else
				{
					tr_y=(elt[lo]->y()-c_n->y())*dx-(elt[lo]->x()-c_n->x())*dy;
					d1=tr_y*tr_y/d;

					//d2=squared_distance(begin,end,elt[hi]);
					tr_y=(elt[hi]->y()-c_n->y())*dx-(elt[hi]->x()-c_n->x())*dy;
					d2=tr_y*tr_y/d;
				}

				if(d1>d2)
					d2=d1; 
			}
			else 								//Few points in hull--search by brute force.
			{
				d2=FT(0);
				for(mid=bot;mid<top;mid++)
				{
					//d1=squared_distance(begin,end,elt[mid]); 
					dx=c_m->x()-c_n->x();
					dy=c_m->y()-c_n->y();
					d=dx*dx+dy*dy;
					if(d==0)
						d1=(elt[mid]->x()-c_n->x())*(elt[mid]->x()-c_n->x())+(elt[mid]->y()-c_n->y())*(elt[mid]->y()-c_n->y());
					else
						tr_y=(elt[mid]->y()-c_n->y())*dx-(elt[mid]->x()-c_n->x())*dy;
					d1=tr_y*tr_y/d;

					if(d1>d2)
						d2=d1; 
				}
			}
			
			List_err[index1]=d2;
			index1++;

			Graph[n][m]=Graph[m][n]=d2;
		}
	}

	delete[] elt;
	elt=NULL;

	if(n_pt_bound==2)
	{
		approx_error=Graph[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(Graph[0][N-1]<=0)									
		fl_circ=true;

	QuickSort<FT>(List_err,N_err);

	delta_N=N_err/2.;
	nc=0;
	m=n_pt_bound+1;

	while(delta_N>0.4)
	{
		if(m>n_pt_bound)
			nc+=delta_N;
		else
			nc-=delta_N;

		N_crt=(int)(nc+0.5);
		if(N_crt<0)
			N_crt=0;
		if(N_crt>=N_err)
			N_crt=N_err;

		cr_err_bound=List_err[N_crt];
		m_crt=N+1;
		err_crt=cr_err_bound+1;

		if(fl_circ)									
			Graph[0][N-1]=Graph[N-1][0]=cr_err_bound+1;

		for(m=0;m<N;m++)
			visited[m]=NULL;

		queue_b=queue_e=new CEL;
		assert(queue_b);
		queue_b->next=NULL;
		queue_b->pt_crt=0;

		breadth_crt=new CEL;
		assert(breadth_crt);

		breadth_crt->pt_crt=0;
		breadth_crt->next=NULL;

		visited[0]=breadth_crt;
		n_vis=1;

		while(n_vis<N)
		{
			for(i=(queue_b->pt_crt+1);i<N;i++)
			{
					if(Graph[queue_b->pt_crt][i]<=cr_err_bound)
					{
						if(i==N-1)
						{
							max=Graph[queue_b->pt_crt][N-1];
							breadth_crt=visited[queue_b->pt_crt];
							m_pt=2;
							while(breadth_crt->next)
							{
								m_pt++;
								vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
								if(vm>max)
									max=vm;
								breadth_crt=breadth_crt->next;
							}
							if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
							{
								m_crt=m_pt;
								err_max=max;
								path=queue_b->pt_crt;
							}
						}

						if(visited[i]==0)
						{
							breadth_crt=new CEL;
							assert(breadth_crt);

							breadth_crt->pt_crt=i;
							breadth_crt->next=visited[queue_b->pt_crt];
							visited[i]=breadth_crt;
							n_vis++;
							
							queue_e->next=new CEL;
							assert(queue_e->next);

							queue_e=queue_e->next;
							queue_e->next=NULL;
							queue_e->pt_crt=i;
						}
					}
			}

			queue_crt=queue_b;
			queue_b=queue_b->next;
			delete queue_crt;
		}

		while(queue_b)
		{
			if(queue_b->pt_crt!=N-1)
				if(Graph[queue_b->pt_crt][N-1]<=cr_err_bound)
				{
					max=Graph[queue_b->pt_crt][N-1];
					breadth_crt=visited[queue_b->pt_crt];
					m_pt=2;
					while(breadth_crt->next)
					{
						m_pt++;
						vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
						if(vm>max)
							max=vm;
						breadth_crt=breadth_crt->next;
					}
					if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
					{
						m_crt=m_pt;
						err_max=max;
						path=queue_b->pt_crt;
					}
				}
			queue_crt=queue_b;
			queue_b=queue_b->next;
			delete queue_crt;
		}

		m=m_crt;

		delete visited[N-1];
		visited[N-1]=NULL;

		breadth_crt=visited[path]->next;

		delete visited[path];
		visited[path]=NULL;

		while(breadth_crt)
		{
			path=breadth_crt->pt_crt;
			breadth_crt=breadth_crt->next;

			delete visited[path];
			visited[path]=NULL;
		}

		for(i=0;i<N;i++)
			if(visited[i])
				delete visited[i];

		approx_error=cr_err_bound;
	
		delta_N/=2.;

		if(m<=n_pt_bound)
			N_mf=N_crt;
	}

	if(m!=n_pt_bound)
		N_crt=N_mf;
	
	cr_err_bound=List_err[N_crt];
	m_crt=N+1;
	err_crt=cr_err_bound+1;

	if(fl_circ)									
		Graph[0][N-1]=Graph[N-1][0]=cr_err_bound+1;

	for(m=0;m<N;m++)
		visited[m]=NULL;
	
	queue_b=queue_e=new CEL;
	assert(queue_b);

	queue_b->next=NULL;
	queue_b->pt_crt=0;

	breadth_crt=new CEL;
	assert(breadth_crt);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++) 
		{
				if(Graph[queue_b->pt_crt][i]<=cr_err_bound)
				{
					if(i==N-1)
					{
						max=Graph[queue_b->pt_crt][N-1];
						breadth_crt=visited[queue_b->pt_crt];
						m_pt=2;
						while(breadth_crt->next)
						{
							m_pt++;
							vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
							if(vm>max)
								max=vm;
							breadth_crt=breadth_crt->next;
						}
						if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
						{
							m_crt=m_pt;
							err_max=max;
							path=queue_b->pt_crt;
						}
					}

					if(visited[i]==0)
					{
						breadth_crt=new CEL;
						assert(breadth_crt);

						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;
					
						queue_e->next=new CEL;
						assert(queue_e->next);

						queue_e=queue_e->next;
						queue_e->next=NULL;
						queue_e->pt_crt=i;
					}
				}
		}

		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	while(queue_b)
	{
		if(queue_b->pt_crt!=N-1)
			if(Graph[queue_b->pt_crt][N-1]<=cr_err_bound)
			{
				max=Graph[queue_b->pt_crt][N-1];
				breadth_crt=visited[queue_b->pt_crt];
				m_pt=2;
				while(breadth_crt->next)
				{
					m_pt++;
					vm=Graph[breadth_crt->pt_crt][breadth_crt->next->pt_crt];
					
					if(vm>max)
						max=vm;
					breadth_crt=breadth_crt->next;
				}
				if(m_pt<m_crt || (m_pt==m_crt && err_max>max))
				{
					m_crt=m_pt;
					err_max=max;
					path=queue_b->pt_crt;
				}
			}
		queue_crt=queue_b;
		queue_b=queue_b->next;
		delete queue_crt;
	}

	m=m_crt;
	n_pt_bound=m;

	delete visited[N-1];
	visited[N-1]=NULL;

	breadth_crt=visited[path]->next;

	delete visited[path];
	visited[path]=NULL;

	while(breadth_crt)
	{
		path=breadth_crt->pt_crt;
		breadth_crt=breadth_crt->next;

		delete visited[path];
		visited[path]=NULL;
	}

	for(c_n=begin,i=0;c_n!=beyond;c_n++,i++)
		if(visited[i]==NULL)
			*result++=*c_n;

	approx_error=cr_err_bound;

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	delete[] visited;

	delete[] List_err;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;

	return result;
}


//  Bounded-E global error assessment; Based on the incremental technique (Perez and Vidal)

template<class InputIterator,class OutputIterator, class FT>
OutputIterator polygonal_approximation_ITGS_be_gea(	InputIterator begin, 
													InputIterator beyond, 
													typename std::size_t &approx_n_pt, 
													FT error_bound, 
													OutputIterator result )
{
	typename std::size_t n,m,i,j;
	InputIterator c_n,c_i,c_m;
	
	if(begin==beyond)
	{
		approx_n_pt=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=1;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_n_pt=2;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	for(n=0,c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	typename std::size_t N=n+1;

	FT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	FT nr_p;
	FT ssvdErr,ssedErr,xk,yk,yValY,dir_coef,dif;

	FT **Graph;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);

	Graph[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);
	
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		x_sum=0;
		y_sum=0;
		x2_sum=0;
		y2_sum=0;
		xy_sum=0;
		nr_p=FT(0);

		for(m=n+2,c_m=c_n,c_m++;m<N;m++) 
		{
			nr_p+=FT(1);
			xk=c_m->x();
			yk=c_m->y();

			c_m++;

			x_sum += xk;
			y_sum += yk;
			x2_sum += (xk * xk);
			y2_sum += (yk * yk); 
			xy_sum += (xk * yk);

			//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		
			dif=(c_m->x()-c_n->x());

			if(dif!=FT(0))
			{
				dir_coef=(c_m->y()-c_n->y())/dif;
				yValY=c_n->y()-dir_coef*c_n->x();

				ssvdErr=yValY*yValY*nr_p+FT(2)*yValY*dir_coef*x_sum-FT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-FT(2)*dir_coef*xy_sum;
				ssedErr=(FT(1)/(dir_coef*dir_coef+FT(1)))*ssvdErr;
			}
			else
			{
				if((c_m->y()-c_n->y())==FT(0))				// Closed curve case (*begin==*end)
				{	
					xx=c_n->x();
					yValY=c_n->y();

					ssedErr=xx*xx*nr_p+yValY*yValY*nr_p-FT(2)*yValY*y_sum+y2_sum+x2_sum-FT(2)*xx*x_sum;
				}
				else
				{
					xx=c_n->x();
					ssedErr=x2_sum-FT(2)*xx*x_sum+nr_p*xx*xx;
				}
			}

			if(ssedErr<FT(0))
				ssedErr=-ssedErr;

			Graph[n][m]=Graph[m][n]=ssedErr;
		}
	}
	
	if(Graph[0][N-1]<=error_bound)
	{
		approx_n_pt=2;

		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}

	bool fl_find=false;
	FT err_crt;
	FT* ApproxError;
	typename std::size_t **SplitPoints;

	SplitPoints=new typename std::size_t*[N+1];
	assert(SplitPoints);

	ApproxError=new FT[N+1];
	assert(ApproxError);

	m=2;

	for(i=0;i<=N-1;i++)
		ApproxError[i]=Graph[i][0];

	while(!fl_find)
	{
		m++;
		SplitPoints[m]=new typename std::size_t[N];
		assert(SplitPoints[m]);

		ApproxError[N-1]=ApproxError[N-2];
		SplitPoints[m][N-1]=N-2;
		if(ApproxError[N-1]<=error_bound)
			fl_find=true;
		else
			for(j=N-3;j>=m-2;j--)
			{
				if(ApproxError[j]<=error_bound && Graph[j][N-1]<=error_bound)
				{
					err_crt=ApproxError[j]+Graph[j][N-1];

					if(err_crt<ApproxError[N-1])
					{
						ApproxError[N-1]=err_crt;
						SplitPoints[m][N-1]=j;
						if(ApproxError[N-1]<=error_bound)
						{
							fl_find=true;
							break;
						}
					}
				}
			}

		if(!fl_find)
			for(i=N-2;i>=m-1;i--)
			{
				ApproxError[i]=ApproxError[i-1];
				SplitPoints[m][i]=i-1;
				for(j=i-2;j>=m-2;j--)
				{
					if(ApproxError[j]<=error_bound && Graph[j][i]<=error_bound)
					{
						err_crt=ApproxError[j]+Graph[j][i];
						if(err_crt<ApproxError[i])
						{
							ApproxError[i]=err_crt;
							SplitPoints[m][i]=j;
						}
					}
				}
			}
	}

	approx_n_pt=m;

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=approx_n_pt;m>2;m--)
	{
		n=SplitPoints[m][n];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	delete[] ApproxError;

	for(m=3;m<=approx_n_pt;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;

 
	return result;
};


//  Bounded-# global error assessment; Based on the incremental technique (Perez and Vidal)

template<class InputIterator,class OutputIterator, class FT>
OutputIterator polygonal_approximation_ITGS_bnp_gea(	InputIterator begin, 
														InputIterator beyond, 
														typename std::size_t n_pt_bound, 
														FT &approx_error, 
														OutputIterator result )
{
	typename std::size_t n,m,i,j;
	InputIterator c_n,c_i,c_m;

	if(begin==beyond)
	{
		approx_error=0;
		return result;
	}

	c_n=begin;
	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		return result;
	}

	c_n++;
	if(c_n==beyond)
	{
		approx_error=0;
		*result++=*begin;
		c_n--;
		*result++=*c_n;
		return result;
	}

	if(n_pt_bound<2)
		return result;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	if(n_pt_bound>n)
	{
		approx_error=0;

		for(c_n=begin;c_n!=beyond;c_n++)
			*result++=*c_n;

		return result;
	}

	typename std::size_t N=n+1;

	FT x_sum,y_sum,x2_sum,y2_sum,xy_sum,xx;
	FT nr_p;
	FT ssvdErr,ssedErr,xk,yk,yValY,dir_coef,dif;

	FT **Graph;

	Graph=new FT*[N];
	assert(Graph);
	for(m=0;m<N;m++)
	{
		Graph[m]=new FT[N];
		assert(Graph[m]);
	}

	Graph[0][0]=Graph[0][1]=Graph[1][0]=FT(0);

	Graph[N-1][N-1]=FT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=FT(0);
	
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		x_sum=0;
		y_sum=0;
		x2_sum=0;
		y2_sum=0;
		xy_sum=0;
		nr_p=FT(0);

		for(m=n+2,c_m=c_n,c_m++;m<N;m++) 
		{
			nr_p+=FT(1);
			xk=c_m->x();
			yk=c_m->y();

			c_m++;

			x_sum += xk;
			y_sum += yk;
			x2_sum += (xk * xk);
			y2_sum += (yk * yk); 
			xy_sum += (xk * yk);

			//Compute the direction coefficient and the y-value at the y-axis of approx. curve.
		
			dif=(c_m->x()-c_n->x());

			if(dif!=FT(0))
			{
				dir_coef=(c_m->y()-c_n->y())/dif;
				yValY=c_n->y()-dir_coef*c_n->x();

				ssvdErr=yValY*yValY*nr_p+FT(2)*yValY*dir_coef*x_sum-FT(2)*yValY*y_sum+dir_coef*dir_coef*x2_sum+y2_sum-FT(2)*dir_coef*xy_sum;
				ssedErr=(FT(1)/(dir_coef*dir_coef+FT(1)))*ssvdErr;
			}
			else
			{
				if((c_m->y()-c_n->y())==FT(0))				// Closed curve case (*begin==*end)
				{	
					xx=c_n->x();
					yValY=c_n->y();

					ssedErr=xx*xx*nr_p+yValY*yValY*nr_p-FT(2)*yValY*y_sum+y2_sum+x2_sum-FT(2)*xx*x_sum;
				}
				else
				{
					xx=c_n->x();
					ssedErr=x2_sum-FT(2)*xx*x_sum+nr_p*xx*xx;
				}
			}

			if(ssedErr<FT(0))
				ssedErr=-ssedErr;

			Graph[n][m]=Graph[m][n]=ssedErr;
		}
	}
	
	if(n_pt_bound==2)
	{
		approx_error=Graph[0][N-1];
		*result++=*begin;
		c_n=beyond;
		c_n--;
		*result++=*c_n;
		return result;
	}


	typename std::size_t mm;

	FT err_crt;
	FT* ApproxError;
	typename std::size_t **SplitPoints;

	SplitPoints=new typename std::size_t*[N+1];
	assert(SplitPoints);

	ApproxError=new FT[N+1];
	assert(ApproxError);

	for(i=0;i<=N-1;i++)
		ApproxError[i]=Graph[i][0];

	FT epsi_crt;

	mm=2;

	while(mm<n_pt_bound)
	{
		mm++;
		SplitPoints[mm]=new typename std::size_t[N];
		assert(SplitPoints[mm]);

		epsi_crt=ApproxError[N-1]=ApproxError[N-2];
		SplitPoints[mm][N-1]=N-2;
		
		for(j=N-3;j>=mm-2;j--)
		{
			if(ApproxError[j]<=epsi_crt && Graph[j][N-1]<=epsi_crt)
			{
				err_crt=ApproxError[j]+Graph[j][N-1];

				if(err_crt<ApproxError[N-1])
				{
					epsi_crt=ApproxError[N-1]=err_crt;
					SplitPoints[mm][N-1]=j;
				}
			}
		}
	
		if(mm<n_pt_bound)
			for(i=N-2;i>=mm-1;i--)
			{
				epsi_crt=ApproxError[i]=ApproxError[i-1];
				SplitPoints[mm][i]=i-1;
				for(j=i-2;j>=mm-2;j--)
				{
					if(ApproxError[j]<=epsi_crt && Graph[j][i]<=epsi_crt)
					{
						err_crt=ApproxError[j]+Graph[j][i];
						if(err_crt<ApproxError[i])
						{
							epsi_crt=ApproxError[i]=err_crt;
							SplitPoints[mm][i]=j;
						}
					}
				}
			}
	}

	approx_error=ApproxError[N-1];

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> local_container;
	typename std::list<InputIterator>::iterator cri;
	
	local_container.push_front(c_i);

	n=N-1;
	
	typename std::size_t nn=N-1;

	for(m=n_pt_bound;m>2;m--)
	{
		n=SplitPoints[m][n];

		for(typename std::size_t ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		local_container.push_front(c_i);
	}

	local_container.push_front(begin);

	for(cri=local_container.begin();cri!=local_container.end();cri++)
		*result++=**cri;

	delete[] ApproxError;

	for(m=3;m<=mm;m++)
		delete[] SplitPoints[m];
	delete[] SplitPoints;


	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
}


CGAL_END_NAMESPACE

#endif // CGAL_POLYAP_FCT
