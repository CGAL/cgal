// Header



#ifndef CGAL_POLYAP_FCT
#define CGAL_POLYAP_FCT

#include <assert.h>
#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//
//			Dynamic Programming
//
/////////////
//
//  MIN E
//

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator dynProgApprox(	InputIterator begin, 
								InputIterator beyond, 
								unsigned long nr_pct, 
								typename DistTraits::DataT *SmallestError, 
								OutputIterator result,
								DistTraits error)
{
	typedef typename DistTraits::DataT DataT;
	
	unsigned long n,m,i;
	
	InputIterator c_n,c_i,c_m;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;
	if(nr_pct>n)
		nr_pct=n;

	unsigned long N=n+1;
	
	DataT **ApproxError;
	unsigned long **SplitPoints;
	DataT **Err;

	assert(ApproxError=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(ApproxError[m]=new DataT[nr_pct]);

	assert(SplitPoints=new unsigned long*[N]);
	for(m=0;m<N;m++)
		assert(SplitPoints[m]=new unsigned long[nr_pct]);

	assert(Err=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Err[m]=new DataT[N]);

	Err[0][0]=DataT(0);
	Err[N-1][N-1]=DataT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=DataT(0);
	for(n=0,c_n=begin;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{	
			Err[n][m]=error(c_n,c_m);
			Err[m][n]=Err[n][m];
		}

	for(n=1;n<N;n++)
		ApproxError[n][1]=Err[n][0];

	DataT epsilon;
	
	for(m=2;m<nr_pct;m++)
	{
		ApproxError[m][m]=DataT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				epsilon=error.cumulate(ApproxError[i][m-1],Err[i][n]); 
				if(epsilon<ApproxError[n][m])
				{
					ApproxError[n][m]=epsilon;
					SplitPoints[n][m]=i;
				}
			}
			epsilon=Err[m-1][n];
			if(epsilon<ApproxError[n][m])
			{
				ApproxError[n][m]=epsilon;
				SplitPoints[n][m]=m-1;
			}
			
		}
	}
	
	*SmallestError=ApproxError[N-1][nr_pct-1];

	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> result1;
	typename std::list<InputIterator>::iterator cri;

	result1.push_front(c_i);

	n=N-1;
	
	unsigned long nn=N-1;

	for(m=nr_pct-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(unsigned long ii=nn;ii>n;ii--)
			c_i--;

		nn=n;

		result1.push_front(c_i);
	}

	result1.push_front(begin);

	for(cri=result1.begin();cri!=result1.end();cri++)
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



//////////////
//
//  MIN #
//
template<class InputIterator,class OutputIterator, class DistTraits>

OutputIterator dynProgApprox(	InputIterator begin, 
								InputIterator beyond, 
								unsigned long *nr_pct, 
								typename DistTraits::DataT epsi, 
								OutputIterator result,
								DistTraits error)
{
	typedef typename DistTraits::DataT DataT;

	unsigned long n,m,i;
	DataT epsilon;
	InputIterator c_n,c_i,c_m;

	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	unsigned long N=n+1;
	
	DataT **ApproxError;
	unsigned long **SplitPoints;
	DataT **Err;

	assert(ApproxError=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(ApproxError[m]=new DataT[N]);

	assert(Err=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Err[m]=new DataT[N]);

	Err[0][0]=Err[0][1]=Err[1][0]=DataT(0);

	if(N>2)
		for(n=2,c_n=begin,c_n++,c_n++;n<N;n++,c_n++)
			ApproxError[n][1]=Err[n][0]=Err[0][n]=error(begin,c_n);

	if(ApproxError[N-1][1]<=epsi)
	{
		*nr_pct=2;
	
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

	assert(SplitPoints=new unsigned long*[N]);
	for(m=0;m<N;m++)
		assert(SplitPoints[m]=new unsigned long[N]);

	Err[N-1][N-1]=DataT(0);
	for(n=1;n<N-1;n++)
		Err[n][n]=Err[n][n-1]=Err[n][n+1]=Err[n-1][n]=Err[n+1][n]=DataT(0);

	for(n=1,c_n=begin,c_n++;n<N-2;n++,c_n++)
		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{	Err[n][m]=error(c_n,c_m);
			Err[m][n]=Err[n][m];
		}

	for(m=2;m<N;m++)
	{
		ApproxError[m][m]=DataT(0);
		SplitPoints[m][m]=m-1;

		for(n=m+1;n<N;n++)
		{
			ApproxError[n][m]=ApproxError[n-1][m-1];
			SplitPoints[n][m]=n-1;
			for(i=n-2;i>=m;i--)
			{
				epsilon=error.cumulate(ApproxError[i][m-1],Err[i][n]); 
				if(epsilon<ApproxError[n][m])
				{
					ApproxError[n][m]=epsilon;
					SplitPoints[n][m]=i;
				}
			}
			epsilon=Err[m-1][n];
			if(epsilon<ApproxError[n][m])
			{
				ApproxError[n][m]=epsilon;
				SplitPoints[n][m]=m-1;
			}
		}

		if(ApproxError[N-1][m]<=epsi)
		{
			*nr_pct=m+1;
			break;
		}
	}
	
	c_i=beyond;
	c_i--;

	typename std::list<InputIterator> result1;
	typename std::list<InputIterator>::iterator cri;
	
	result1.push_front(c_i);

	n=N-1;
	
	unsigned long nn=N-1;

	for(m=*nr_pct-1;m>=2;m--)
	{
		n=SplitPoints[n][m];

		for(unsigned long ii=nn;ii>n;ii--)
			c_i--;

		nn=n;
		result1.push_front(c_i);
	}

	result1.push_front(begin);

	for(cri=result1.begin();cri!=result1.end();cri++)
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
//////////////////////////////////////////////////////////////////////////////////////////
//
//				Recursive Split
//
/////////////////
//
//  Min # recursive implementation
//

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator recSplitApprox(	InputIterator begin, 
								InputIterator end, 
								unsigned long *nr_pct, 
								typename DistTraits::DataT epsi, 
								OutputIterator result,
								DistTraits error)
{
	typedef typename DistTraits::DataT DataT;

	static unsigned long fl;
	DataT d;
	InputIterator b,e,p_max;

	fl++;

	if(fl==1)
	{
		end--;					// put the correct value of end point
		*nr_pct=1;

		d=error(begin,end,&p_max);
		*result++=*begin;		// save begin as the first point of the result
	}
	else
		d=error(begin,end,&p_max);
	
	if(d>epsi)
	{
		recSplitApprox(begin,p_max,nr_pct,epsi,result,error);
		recSplitApprox(p_max,end,nr_pct,epsi,result,error);
	}	
	else
	{
		*result++=*end;
		(*nr_pct)++;
	}

	fl--;

	return result;
}




/////////////////
//
//  Min E
//
template<class InputIterator>
struct Sp_List
{	InputIterator bg;
	InputIterator ed;
	InputIterator pm;
	Sp_List *next;
};	


template<class InputIterator,class DataTp>		
struct WaitList
{	DataTp err;
	Sp_List<InputIterator> *adr_el;
	WaitList *next;
};


template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator recSplitApprox(	InputIterator begin, 
								InputIterator end, 
								unsigned long nr_pct, 
								typename DistTraits::DataT* SmallestError, 
								OutputIterator result,
								DistTraits error)
{
	typedef typename DistTraits::DataT DataT;

	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	WaitList<InputIterator,DataT> *l_wt,*wt_cr,*wt_new;
	InputIterator cr,p_max;

	l_sp=new Sp_List<InputIterator>;
	assert(l_sp);
	l_wt=new WaitList<InputIterator,DataT>;
	assert(l_wt);

	cr=end;
	cr--;
	l_sp->next=NULL;
	l_sp->bg=begin;
	l_sp->ed=cr;

	l_wt->next=NULL;
	l_wt->adr_el=l_sp;
	l_wt->err=error(begin,cr,&p_max);

	l_sp->pm=p_max;

	unsigned long m=2;

	while(m<nr_pct)
	{
		sp_cr=l_wt->adr_el;
		
		assert(sp_new=new Sp_List<InputIterator>);

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->ed=sp_cr->ed;
		sp_new->bg=sp_cr->pm;
		sp_cr->ed=sp_cr->pm;

		wt_new=new WaitList<InputIterator,DataT>;
		assert(wt_new);

		wt_new->err=error(sp_cr->bg,sp_cr->ed,&p_max);
		sp_cr->pm=p_max;
		wt_new->adr_el=sp_cr;

		wt_cr=l_wt;
		while(wt_cr->next!=NULL && wt_cr->next->err>wt_new->err)
			wt_cr=wt_cr->next;
		wt_new->next=wt_cr->next;
		wt_cr->next=wt_new;

		wt_new=new WaitList<InputIterator,DataT>;
		assert(wt_new);

		wt_new->err=error(sp_new->bg,sp_new->ed,&p_max);
		wt_new->adr_el=sp_new;
		sp_new->pm=p_max;

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

	*SmallestError=l_wt->err;

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		*result++=*l_sp->ed;
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



////////////////////////
//
//	 Min # iterativ implementation

template<class InputIterator,class OutputIterator, class DistTraits>
OutputIterator recSplitApprox_itr(	InputIterator begin, 
									InputIterator beyond, 
									unsigned long* nr_pct, 
									typename DistTraits::DataT epsi, 
									OutputIterator result,
									DistTraits error)
{
	typedef typename DistTraits::DataT DataT;

	Sp_List<InputIterator> *l_sp,*sp_cr,*sp_new;
	InputIterator cr,p_max;
	DataT err;

	assert(l_sp=new Sp_List<InputIterator>);

	cr=beyond;
	cr--;
	l_sp->next=NULL;
	l_sp->bg=begin;
	l_sp->ed=cr;

	err=error(begin,cr,&p_max);

	l_sp->pm=p_max;

	sp_cr=l_sp;

	*nr_pct=2;

	while(err>epsi)
	{
		sp_new=new Sp_List<InputIterator>;
		(*nr_pct)++;

		sp_new->next=sp_cr->next;
		sp_cr->next=sp_new;

		sp_new->ed=sp_cr->ed;
		sp_new->bg=sp_cr->pm;
		sp_cr->ed=sp_cr->pm;

		err=error(sp_cr->bg,sp_cr->ed,&p_max);

		while(sp_cr->next && err<=epsi)
		{
			sp_cr=sp_cr->next;
			err=error(sp_cr->bg,sp_cr->ed,&p_max);
		}
		
		sp_cr->pm=p_max;
	}

	*result++=*begin;

	while(l_sp!=NULL)
	{
		sp_cr=l_sp;
		*result++=*l_sp->ed;
		l_sp=l_sp->next;
		delete sp_cr;
	}

	return result;
};





////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
//		Graph Search
//

//////////////
//
//  MIN # 
//
template<class DistTraits,class InputIterator,class OutputIterator>
OutputIterator GraphSearchApprox(	InputIterator begin, 
									InputIterator beyond, 
									unsigned long *nr_pct, 
									typename DistTraits::DataT epsi, 
									OutputIterator result,
									DistTraits error
								)
{
	typedef typename DistTraits::DataT DataT;
	unsigned long path;

	unsigned long n,m,m_pt,i,n_vis,m_crt;
	typename DistTraits::DataT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	
	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	unsigned long N=n+1;

	typedef struct cel
	{	unsigned long pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

	m_crt=N+1;
	err_crt=epsi+1;
	
	DataT **Graph;
	CEL **visited;

	assert(visited=new CEL*[N]);
	for(m=0;m<N;m++)
		visited[m]=NULL;

	assert(Graph=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Graph[m]=new DataT[N]);

	Graph[0][0]=Graph[0][1]=Graph[1][0]=DataT(0);
	Graph[N-1][N-1]=DataT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=DataT(0);

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{
			Graph[n][m]=Graph[m][n]=error(c_n,c_m);
		}
	}

	if(Graph[0][N-1]<=epsi)
		Graph[0][N-1]=Graph[N-1][0]=epsi+1;
	
	assert(queue_b=queue_e=new CEL);
	queue_b->next=NULL;
	queue_b->pt_crt=0;

	assert(breadth_crt=new CEL);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++)
		{
				if(Graph[queue_b->pt_crt][i]<=epsi)
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
						assert(breadth_crt=new CEL);
						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;
						assert(queue_e->next=new CEL);
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
			if(Graph[queue_b->pt_crt][N-1]<=epsi)
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

	*nr_pct=m_crt;

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
	{
		if(visited[i]==NULL)
		{
			*result++=*c_n;

		}
	}

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	delete[] visited;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
};



template<class DistTraits>

void ordonare_qs(typename DistTraits::DataT* p,unsigned long n)  //quick sort
{
	typename DistTraits::DataT tmp,pivot;

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
	int i=0,j=n-1;

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
//	if(i>0)
		ordonare_qs<DistTraits>(p,i+1);
	j=n-i-1;
//	if(j>1)
		ordonare_qs<DistTraits>(&p[i+1],j);
}




//////////////
//
//  MIN E - based on Binary search and MIN # 
//
template<class DistTraits,class InputIterator,class OutputIterator>
OutputIterator GraphSearchApprox(	InputIterator begin, 
									InputIterator beyond, 
									unsigned long nr_pct, 
									typename DistTraits::DataT *SmallestError, 
									OutputIterator result,
									DistTraits error)
{

	typedef typename DistTraits::DataT DataT;
	unsigned long path;

	unsigned long n,m,m_pt,i,n_vis,m_crt;
	typename DistTraits::DataT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	
	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	unsigned long N=n+1;

	typedef struct cel
	{	unsigned long pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

//	m_crt=N+1;
//	err_crt=epsi+1;
	
	DataT **Graph;
	CEL **visited;

	assert(visited=new CEL*[N]);
	for(m=0;m<N;m++)
		visited[m]=NULL;

	assert(Graph=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Graph[m]=new DataT[N]);

	Graph[0][0]=Graph[0][1]=Graph[1][0]=DataT(0);
	Graph[N-1][N-1]=DataT(0);

	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=DataT(0);

	long index1=0;
	unsigned long N_err=(N-1)*(N-2)/2;
	DataT *List_err;
	assert(List_err=new DataT[N_err]);

	for(n=0,c_n=begin;n<N-2;n++,c_n++)
	{
		for(m=n+2,c_m=c_n,c_m++,c_m++;c_m!=beyond;m++,c_m++)
		{
			Graph[n][m]=Graph[m][n]=error(c_n,c_m);
			List_err[index1]=Graph[m][n];
			index1++;
		}
	}
	
	
	DataT epsi;
	bool fl_circ=false;

	unsigned long N_crt,N_mf;
	double delta_N,nc;
	

	if(Graph[0][N-1]<=0)
		fl_circ=true;

	ordonare_qs<DistTraits>(List_err,N_err);

	delta_N=N_err/2.;
	nc=0;
	m=nr_pct+1;

	while(delta_N>0.4)
	{
		if(m>nr_pct)
			nc+=delta_N;
		else
			nc-=delta_N;

		N_crt=(int)(nc+0.5);
		if(N_crt<0)
			N_crt=0;
		if(N_crt>=N_err)
			N_crt=N_err;

		epsi=List_err[N_crt];
		m_crt=N+1;
		err_crt=epsi+1;

		if(fl_circ)									
			Graph[0][N-1]=Graph[N-1][0]=epsi+1;

		for(m=0;m<N;m++)
			visited[m]=NULL;
		
		assert(queue_b=queue_e=new CEL);
		queue_b->next=NULL;
		queue_b->pt_crt=0;

		assert(breadth_crt=new CEL);
		breadth_crt->pt_crt=0;
		breadth_crt->next=NULL;

		visited[0]=breadth_crt;
		n_vis=1;

		while(n_vis<N)
		{
			for(i=(queue_b->pt_crt+1);i<N;i++)
			{
					if(Graph[queue_b->pt_crt][i]<=epsi)
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
							assert(breadth_crt=new CEL);
							breadth_crt->pt_crt=i;
							breadth_crt->next=visited[queue_b->pt_crt];
							visited[i]=breadth_crt;
							n_vis++;
							assert(queue_e->next=new CEL);
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
				if(Graph[queue_b->pt_crt][N-1]<=epsi)
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

		*SmallestError=epsi;
	
		delta_N/=2.;

		if(m<=nr_pct)
			N_mf=N_crt;
	}

	if(m!=nr_pct)
		N_crt=N_mf;
	
	epsi=List_err[N_crt];
	m_crt=N+1;
	err_crt=epsi+1;

	if(fl_circ)									
		Graph[0][N-1]=Graph[N-1][0]=epsi+1;

	for(m=0;m<N;m++)
		visited[m]=NULL;
		
	assert(queue_b=queue_e=new CEL);
	queue_b->next=NULL;
	queue_b->pt_crt=0;

	assert(breadth_crt=new CEL);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++) 
		{
				if(Graph[queue_b->pt_crt][i]<=epsi)
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
						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;
						assert(queue_e->next=new CEL);
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
			if(Graph[queue_b->pt_crt][N-1]<=epsi)
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
	nr_pct=m;

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
	{
		if(visited[i]==NULL)
		{
			*result++=*c_n;
		}
	}

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	*SmallestError=epsi;

	delete[] visited;

	delete[] List_err;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;

	return result;
}


























////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
//		Optimal algorithms for Max Euclidean distance
//

template<class InputIterator, class DistTraits>
bool SLOPE_SIGN1(InputIterator p,InputIterator q,InputIterator begin,InputIterator end)  
{ 
	typedef typename DistTraits::DataT DataT;
	
	DataT sHomogX,sHomogY,deltaXpq,deltaYpq,slopeVal;
 
	sHomogX=begin->y()-end->y();
	sHomogY=end->x()-begin->x();

	if(sHomogX==DataT(0) && sHomogY==DataT(0))
		sHomogX=DataT(1);

	deltaXpq=q->x()-p->x();
	deltaYpq=q->y()-p->y();
	slopeVal=sHomogX*deltaXpq+sHomogY*deltaYpq;
 
	return (slopeVal>=DataT(0) );
}


//////////////
//
//  MIN # - Based on the online convex hull (Toussaint)
//
template<class DistTraits,class InputIterator,class OutputIterator>
OutputIterator GraphToussaintApprox(	InputIterator begin, 
										InputIterator beyond, 
										unsigned long *nr_pct, 
										typename DistTraits::DataT epsi, 
										OutputIterator result
										//,DistTraits error
										)
{
	typedef typename DistTraits::DataT DataT;

	unsigned long path;

	unsigned long n,m,m_pt,i,n_vis,m_crt;
	typename DistTraits::DataT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	
	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	unsigned long N=n+1;


	// build the double stack


	long 			top;			//Top end of elt.
	long			bot;			//Bottom end of elt.
	InputIterator*	elt;			//Double ended queue storing a convex hull. dim=TWICE_HULL_MAX

	long HULL_MAX;
	long TWICE_HULL_MAX;
	long THRICE_HULL_MAX;

	bool topflag,botflag;

	HULL_MAX=N;
	TWICE_HULL_MAX=2*HULL_MAX;
	THRICE_HULL_MAX=3*HULL_MAX;

	assert(elt=new InputIterator[TWICE_HULL_MAX]);

	for (long j=0;j<TWICE_HULL_MAX;j++)
		elt[j]=NULL;

	typedef struct cel
	{	unsigned long pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;

	m_crt=N+1;
	err_crt=epsi+1;
	
	DataT **Graph;
	CEL **visited;

	visited=new CEL*[N];
	for(m=0;m<N;m++)
		visited[m]=NULL;

	assert(Graph=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Graph[m]=new DataT[N]);

	Graph[0][0]=Graph[0][1]=Graph[1][0]=DataT(0);

	Graph[N-1][N-1]=DataT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=DataT(0);

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

			long mid,lo,m1,brk,m2,hi;
			bool sbase, sbrk;
			DataT d1,d2,dx,dy,d,tr_y;

			if((top-bot)>6)						//If there are > 6  points on the hull 
			{									//(otherwise we will just look at them all.
				lo=bot; 
				hi=top-1;	
				sbase = SLOPE_SIGN1<InputIterator,DistTraits>(elt[hi],elt[lo],c_n,c_m); //The sign of the base edge.    
			
				do
				{ 
					brk = (lo + hi) / 2;		//Binary search for an edge with opposite sign.
					sbrk = SLOPE_SIGN1<InputIterator,DistTraits>(elt[brk], elt[brk+1], c_n,c_m);
				
					if(sbase == sbrk )
						if(sbase == (SLOPE_SIGN1<InputIterator,DistTraits>(elt[lo], elt[brk+1], c_n,c_m)) ) 
							lo = brk + 1;
						else 
							hi = brk;
				}
				while(sbase==sbrk);
  
				m1=brk;							//Now, the sign changes between the base edge and
				while(lo<m1)					// brk+1. Binary search for the extreme point.
				{
					mid=(lo+ m1)/2;
				
					if(sbase==(SLOPE_SIGN1<InputIterator,DistTraits>(elt[mid],elt[mid+1],c_n,c_m)) ) 
						lo=mid+1;
					else 
						m1=mid;
				}  

				m2=brk;							//The sign also chagnes between brk and the base.
				while(m2<hi)					//Binary search again.
				{
					mid=(m2+hi)/2;
				
					if(sbase==SLOPE_SIGN1<InputIterator,DistTraits>(elt[mid],elt[mid+1],c_n,c_m) ) 
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
				d2=DataT(0);
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

	if(Graph[0][N-1]<=epsi)
		Graph[0][N-1]=Graph[N-1][0]=epsi+1;
	
	assert(queue_b=queue_e=new CEL);
	queue_b->next=NULL;
	queue_b->pt_crt=0;

	assert(breadth_crt=new CEL);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++)
		{
				if(Graph[queue_b->pt_crt][i]<=epsi)
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
						assert(breadth_crt=new CEL);
						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;
						assert(queue_e->next=new CEL);
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
			if(Graph[queue_b->pt_crt][N-1]<=epsi)
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

	*nr_pct=m_crt;

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
	{
		if(visited[i]==NULL)
		{
			*result++=*c_n;

		}
	}

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	delete[] visited;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;
 
	return result;
};

/*
template<class DistTraits>
void ordonare_qs(typename DistTraits::DataT* p,unsigned long n)  //quick sort
{
	typename DistTraits::DataT tmp,pivot;

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
	int i=0,j=n-1;

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
//	if(i>0)
		ordonare_qs<DistTraits>(p,i+1);
	j=n-i-1;
//	if(j>1)
		ordonare_qs<DistTraits>(&p[i+1],j);
}

*/


//////////////
//
//  MIN E - based on Binary search and MIN # (Toussaint)
//
template<class DistTraits,class InputIterator,class OutputIterator>
OutputIterator BinarySearchApprox(	InputIterator begin, 
									InputIterator beyond, 
									unsigned long nr_pct, 
									typename DistTraits::DataT *SmallestError, 
									OutputIterator result)
//									DistTraits error)
{
	typedef typename DistTraits::DataT DataT;
	unsigned long path;

	unsigned long n,m;
	unsigned long i;
	unsigned long m_pt,n_vis,m_crt;
	DataT err_crt,max,vm,err_max;
	InputIterator c_n,c_i,c_m;
	DataT epsi;
	bool fl_circ=false;
	
	n=0;
	for(c_n=begin,c_n++;c_n!=beyond;c_n++)
		n++;

	unsigned long N=n+1;
	unsigned long N_err=(N-1)*(N-2)/2;
	unsigned long N_crt,N_mf;
	double delta_N,nc;


	// build the double stack

	long 			top;			//Top end of elt.
	long			bot;			//Bottom end of elt.
	InputIterator*	elt;			//Double ended queue storing a convex hull. dim=TWICE_HULL_MAX

	long HULL_MAX;
	long TWICE_HULL_MAX;
	long THRICE_HULL_MAX;

	bool topflag,botflag;

	HULL_MAX=N;
	TWICE_HULL_MAX=2*HULL_MAX;
	THRICE_HULL_MAX=3*HULL_MAX;

	assert(elt=new InputIterator[TWICE_HULL_MAX]);

	for (long j=0;j<TWICE_HULL_MAX;j++)
		elt[j]=NULL;

	
	typedef struct cel
	{	unsigned long pt_crt;
		struct cel *next;
	} CEL;

	CEL *queue_b,*queue_e,*queue_crt;
	queue_e=queue_b=queue_crt=NULL;

	CEL *breadth_crt;
	breadth_crt=NULL;
	
	CEL **visited;

	assert(visited=new CEL*[N]);
	for(m=0;m<N;m++)
		visited[m]=NULL;

	DataT *List_err;
	DataT **Graph;

	assert(List_err=new DataT[N_err]);

	assert(Graph=new DataT*[N]);
	for(m=0;m<N;m++)
		assert(Graph[m]=new DataT[N]);

	Graph[0][0]=Graph[0][1]=Graph[1][0]=DataT(0);
	Graph[N-1][N-1]=DataT(0);
	for(n=1;n<N-1;n++)
		Graph[n][n]=Graph[n][n-1]=Graph[n][n+1]=Graph[n-1][n]=Graph[n+1][n]=DataT(0);

	long index1=0;

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

			long mid,lo,m1,brk,m2,hi;
			bool sbase, sbrk;
			DataT d1,d2,dx,dy,d,tr_y;

			if((top-bot)>6)						//If there are > 6  points on the hull 
			{									//(otherwise we will just look at them all.
				lo=bot; 
				hi=top-1;	
				sbase = SLOPE_SIGN1<InputIterator,DistTraits>(elt[hi],elt[lo],c_n,c_m); //The sign of the base edge.    
			
				do
				{ 
					brk = (lo + hi) / 2;		//Binary search for an edge with opposite sign.
					sbrk = SLOPE_SIGN1<InputIterator,DistTraits>(elt[brk], elt[brk+1], c_n,c_m);
				
					if(sbase == sbrk )
						if(sbase == (SLOPE_SIGN1<InputIterator,DistTraits>(elt[lo], elt[brk+1], c_n,c_m)) ) 
							lo = brk + 1;
						else 
							hi = brk;
				}
				while(sbase==sbrk);
  
				m1=brk;							//Now, the sign changes between the base edge and
				while(lo<m1)					// brk+1. Binary search for the extreme point.
				{
					mid=(lo+ m1)/2;
				
					if(sbase==(SLOPE_SIGN1<InputIterator,DistTraits>(elt[mid],elt[mid+1],c_n,c_m)) ) 
						lo=mid+1;
					else 
						m1=mid;
				}  

				m2=brk;							//The sign also chagnes between brk and the base.
				while(m2<hi)					//Binary search again.
				{
					mid=(m2+hi)/2;
				
					if(sbase==SLOPE_SIGN1<InputIterator,DistTraits>(elt[mid],elt[mid+1],c_n,c_m) ) 
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
				d2=DataT(0);
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

	if(Graph[0][N-1]<=0)
		fl_circ=true;

	ordonare_qs<DistTraits>(List_err,N_err);

	delta_N=N_err/2.;
	nc=0;
	m=nr_pct+1;

	while(delta_N>0.4)
	{
		if(m>nr_pct)
			nc+=delta_N;
		else
			nc-=delta_N;

		N_crt=(int)(nc+0.5);
		if(N_crt<0)
			N_crt=0;
		if(N_crt>=N_err)
			N_crt=N_err;

		epsi=List_err[N_crt];
		m_crt=N+1;
		err_crt=epsi+1;

		if(fl_circ)									
			Graph[0][N-1]=Graph[N-1][0]=epsi+1;

		for(m=0;m<N;m++)
			visited[m]=NULL;
		
		assert(queue_b=queue_e=new CEL);
		queue_b->next=NULL;
		queue_b->pt_crt=0;

		assert(breadth_crt=new CEL);
		breadth_crt->pt_crt=0;
		breadth_crt->next=NULL;

		visited[0]=breadth_crt;
		n_vis=1;

		while(n_vis<N)
		{
			for(i=(queue_b->pt_crt+1);i<N;i++)
			{
					if(Graph[queue_b->pt_crt][i]<=epsi)
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
							assert(breadth_crt=new CEL);
							breadth_crt->pt_crt=i;
							breadth_crt->next=visited[queue_b->pt_crt];
							visited[i]=breadth_crt;
							n_vis++;
							assert(queue_e->next=new CEL);
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
				if(Graph[queue_b->pt_crt][N-1]<=epsi)
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

		*SmallestError=epsi;
	
		delta_N/=2.;

		if(m<=nr_pct)
			N_mf=N_crt;
	}

	if(m!=nr_pct)
		N_crt=N_mf;
	
	epsi=List_err[N_crt];
	m_crt=N+1;
	err_crt=epsi+1;

	if(fl_circ)									
		Graph[0][N-1]=Graph[N-1][0]=epsi+1;

	for(m=0;m<N;m++)
		visited[m]=NULL;
		
	assert(queue_b=queue_e=new CEL);
	queue_b->next=NULL;
	queue_b->pt_crt=0;

	assert(breadth_crt=new CEL);
	breadth_crt->pt_crt=0;
	breadth_crt->next=NULL;

	visited[0]=breadth_crt;
	n_vis=1;

	while(n_vis<N)
	{
		for(i=queue_b->pt_crt+1;i<N;i++) 
		{
				if(Graph[queue_b->pt_crt][i]<=epsi)
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
						breadth_crt->pt_crt=i;
						breadth_crt->next=visited[queue_b->pt_crt];
						visited[i]=breadth_crt;
						n_vis++;
						assert(queue_e->next=new CEL);
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
			if(Graph[queue_b->pt_crt][N-1]<=epsi)
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
	nr_pct=m;

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
	{
		if(visited[i]==NULL)
		{
			*result++=*c_n;
		}
	}

	for(i=0;i<N;i++)
		if(visited[i])
			delete visited[i];

	*SmallestError=epsi;

	delete[] visited;

	delete[] List_err;

	for(m=0;m<N;m++)
		delete[] Graph[m];
	delete[] Graph;

	return result;
}


CGAL_END_NAMESPACE



#endif
