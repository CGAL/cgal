#ifndef _INTER_H
#define _INTER_H


#define EPSILON   (0.000001)

template<class number_type>
class Lines
{
public:
    class pmpoint;
	class pmcurve;

	Lines(int fp = 1)
	{
		use_fp = fp;
	}

	int use_fp;

	void Intersect(pmpoint* pnts, int num_points, 
			       pmcurve* cvs, int num_curves, 
				   pmpoint* &intr_pnts, int &intr_np, 
				   pmcurve* &intr_cvs, int &intr_nc)
	{
		int *cvs2 = new int[num_curves];
		int i;
		
		for (i = 0; i < num_curves; i++) // if cvs[i] is false then it was 
			cvs2[i] = 1;    // merged with another curve - hence, should
		                    // be ignored.
		
		// upper boundary for the number of 
		// intersection points + endpoints.
		int PNUM = num_points * (num_points + 1); 
		pmpoint * pnt_buff = new pmpoint[PNUM];
		
		// initializing with the endpoints.
		for (i = 0; i < num_points; i++)
		{
			pnt_buff[i] = pnts[i];
			pnt_buff[i].n  = i;
		} 
		PNUM = num_points;     
		
		// upper boundary on the number of segments.
		int CNUM = num_curves * num_curves;
		pmcurve * crv_buff = new pmcurve[CNUM];
		CNUM = 0;    
		
		int k, n;
		// an auxilary points buffer for keeping the intersection
		// points of one curve.
		pmpoint * pnt_aux_buff = new pmpoint[num_curves];
		
		for (i = 0; i < num_curves; i++)
		{
			// n = number of intersections of the line
			// pnt_aux_buff will contain the n intersection points.
			// Intersect is the only place where cvs2[i] can turn to false.
			n = Intersect(i, num_curves, cvs, cvs2, pnt_aux_buff);
			
			// sort the n points from source to target.
			// if cvs[i] intersects two curves in one point-
			// the sort process will eliminate one.
			n = Sort(n, pnt_aux_buff, cvs[i]);
			
			
			if (cvs2[i])
			{
				// adding the new segments
				crv_buff[CNUM].s = cvs[i].s;
				
				for (k = 0; k < n; k++)
				{
					// add the intersection points to the buffer 
					pnt_aux_buff[k].n = PNUM;
					pnt_buff[PNUM] = pnt_aux_buff[k];
					PNUM++;
					
					crv_buff[CNUM].t = pnt_aux_buff[k];   
					CNUM++;
					
					crv_buff[CNUM].s = pnt_aux_buff[k];
				}
				crv_buff[CNUM].t = cvs[i].t;
				CNUM++;
			}
		}
		
		intr_np = PNUM;
		intr_nc = CNUM;
		intr_pnts = pnt_buff;
		intr_cvs = crv_buff;
	}
	
	int Intersect(int i, int num_curves, pmcurve * cvs, 
		int* cvs2, pmpoint *pnt_aux_buff)
	{
		int j, n = 0;
		pmcurve cv;
		pmpoint p;
		
		for (j = 0; j < num_curves; j++)  // checking intersection with every segment
		{
			if ((j != i) && (cvs2[j]))	  
			{
				if ( (!is_vertical(cvs[i])) &&
					(!is_vertical(cvs[j])) &&
					is_same(a(cvs[i]),a(cvs[j]))	)
				{
					// case 1:the segments have the same derivitive 
					// (however they are not vertical)
					if ( is_same(b(cvs[i]), b(cvs[j])))
						// the segments are on the same line
						if (is_in_x_range(cvs[i],cvs[j].s) ||
							is_in_x_range(cvs[i],cvs[j].t) ||
							is_in_x_range(cvs[j],cvs[i].s) )
						{
							// the segments overlap
							if (j > i)  // always true (?)
							{
								cv.s = leftmost( leftmost(cvs[i].s, cvs[j].s), 
									leftmost(cvs[i].t, cvs[j].t) );
								cv.t = rightmost(rightmost(cvs[i].s, cvs[j].s), 
									rightmost(cvs[i].t, cvs[j].t));
								cvs[j] = cv;
								cvs2[i] = 0;
								return 0;
							}
						}
						// else - the segments don't intersect
						// else - the segments don't intersect (parralel)
				} // end of case 1.
				else
				{
					if (is_vertical(cvs[i]) && is_vertical(cvs[j]))
					{ 	
						// case 2: both of the segments are vertical
						if ( is_same(cvs[i].s.x, cvs[j].s.x) &&
							(is_in_y_range(cvs[i], cvs[j].s) ||
							is_in_y_range(cvs[i], cvs[j].t) ||
							is_in_y_range(cvs[j],	cvs[i].s)) )
						{	
							// the segments overlap
							if (j > i) // always true (?)
							{
								cv.s = lowest( lowest(cvs[i].s, cvs[j].s), 
									lowest(cvs[i].t, cvs[j].t) );
								cv.t = highest(highest(cvs[i].s, cvs[j].s), 
									highest(cvs[i].t, cvs[j].t));
								cvs[j] = cv;
								cvs2[i] = 0;
								return 0;
							}
						} 
						// if the segments are parralel - do nothing
					}
					else
					{
						if (is_vertical(cvs[i]))
						{
							// case 3: cvs[i] is vertical but cvs[j] isn't
							p.x = cvs[i].s.x;
							p.y = a(cvs[j]) * p.x + b(cvs[j]);
						}
						else
							if (is_vertical(cvs[j]))
							{    
								// case 3: cvs[j] is vertical but cvs[i] isn't
								p.x = cvs[j].s.x;
								p.y = a(cvs[i]) * p.x + b(cvs[i]);
							}
							else
							{
								p.x = (b(cvs[i]) - b(cvs[j]))/(a(cvs[j]) - a(cvs[i]));
								p.y = a(cvs[i]) * p.x + b(cvs[i]);
							}
							
							if (is_in_x_range(cvs[i],p) && is_in_x_range(cvs[j],p) &&
								is_in_y_range(cvs[i],p) && is_in_y_range(cvs[j],p) &&
								(!is_same(p,cvs[i].s)) && (!is_same(p,cvs[i].t)))
							{
								pnt_aux_buff[n] = p;
								n++;
							}
					}
				}
			}
		}
		return n;
	}
	
	
	int Sort(int n, pmpoint* buff, pmcurve &cv)
	{
		if (n == 0) return 0;
		
		int flag1, flag2;
		flag1 = flag2 = 0;
		if (is_same_x(cv.s, cv.t))
		{
			flag1 = 1;    // indicates that we have to compare on the y-axes
			if (is_higher(cv.s, cv.t))
				flag2 = 1; // we have to sort from high to low
		}
		else
			if (is_right(cv.s, cv.t))
				flag2 = 1;
			
			//cases: flag1 = 0:   flag2 = 0  from left to right
			//					 flag2 = 1  from right to left
			//       flag1 = 1:   flag2 = 0  from low to up
			//					 flag2 = 1  from up	to low
			
			pmpoint tmp;
			int i, j;
			
			for (i = 0; i < n - 1; i++) // buff[i] contains the min in the i+1'th iteration
				for (j = i+1; j < n; j++)
				{
					if ( ((!flag1) && (!flag2) && (is_left(buff[j], buff[i])))  ||
						((!flag1) && (flag2) && (is_right(buff[j], buff[i])))	 ||
						((flag1) && (!flag2) && (is_lower(buff[j], buff[i])))	 ||
						((flag1) && (flag2) && (is_higher(buff[j], buff[i])))	  )
					{
						tmp = buff[i];
						buff[i] = buff[j];
						buff[j] = tmp;
					}
				}
				
				int n2 = n;
				for (i = n-1; i > 0; i--)
				{
					if (is_same(buff[i], buff[i-1]))
					{
						for (j = i; j < n2 - 1; j++)
							buff[j] = buff[j+1];
						n2--;
					}
				}
				
				return n2;
	}       
	

	
	class pmpoint
	{
	public:
		number_type x;
		number_type y;
		int n;
	};
	
	
	int is_left(const pmpoint &p1, const pmpoint &p2)
	{ 
		return (p1.x < p2.x); 
	}
	
	int is_right(const pmpoint &p1, const pmpoint &p2) 
	{ 
		return (p1.x > p2.x); 
	}
	
	int is_same_x(const pmpoint &p1, const pmpoint &p2)
	{ 
		return (p1.x == p2.x); 
	}
	
	int is_lower(const pmpoint &p1, const pmpoint &p2) 
	{ 
		return (p1.y < p2.y); 
	}
	
	int is_higher(const pmpoint &p1, const pmpoint &p2) 
	{ 
		return (p1.y > p2.y); 
	}
	
	int is_same_y(const pmpoint &p1, const pmpoint &p2) 
	{ 
		return (p1.y == p2.y); 
	}
	
	int is_same(number_type d1, number_type d2)
	{
		return (d1 == d2);
	}
	
	int is_same(const pmpoint &p1, const pmpoint &p2) 
	{
		if (!use_fp)
		{
			return is_same_x(p1,p2) && is_same_y(p1,p2);
		}
		else
		{
			number_type d = (p1.x - p2.x)*(p1.x - p2.x) + 
				(p1.y - p2.y)*(p1.y - p2.y);
			if (d < EPSILON * EPSILON)
				return 1;
			else
				return 0;
			// return (check_left(p1,p2) == 0) && (check_lower(p1,p2) == 0); 
		}
	}
	
	pmpoint leftmost(const pmpoint &p1, const pmpoint &p2)
	{ return (is_left(p1, p2) ? p1 : p2); }
	
	pmpoint rightmost(const pmpoint &p1, const pmpoint &p2)
	{ return (is_right(p1, p2) ? p1 : p2); }
	
	pmpoint lowest(const pmpoint &p1, const pmpoint &p2)
	{ return (is_lower(p1, p2) ? p1 : p2); }
	
	pmpoint highest(const pmpoint &p1, const pmpoint &p2)
	{ return (is_higher(p1, p2) ? p1 : p2); }
	
	
	number_type a(pmcurve &c) 
	{
		if (is_same_x(c.s ,c.t))
			return 0;
		return (c.t.y - c.s.y)/(c.t.x - c.s.x);
	}
	
	number_type b(pmcurve &c)
	{
		if (is_same_x(c.s ,c.t))
			return 0;
		return (c.t.y * c.s.x - c.s.y * c.t.x )/(c.s.x - c.t.x);
	}
	
	int is_vertical(pmcurve &c)
	{
		return (c.s.x == c.t.x);
	}
	
	
	class pmcurve
	{
	public:
		pmpoint s;
		pmpoint t;
	};
	
	
	int is_in_x_range(pmcurve &cv, pmpoint &q) 
	{ 
		int r = !( is_right(q, rightmost(cv.t, cv.s)) ||
			is_left(q, leftmost(cv.t, cv.s))	 );
		return r;
	}
	
	int is_in_y_range(pmcurve &cv, pmpoint &q) 
	{ 
		int r = !( is_higher(q, highest(cv.t, cv.s)) ||
			is_lower(q, lowest(cv.t, cv.s))	 );
		return r;
	}
	
};
	

	
#endif
