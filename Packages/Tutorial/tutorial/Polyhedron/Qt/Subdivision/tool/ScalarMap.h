//**************************************
// ScalarMap
//**************************************
// Pierre Alliez
// Mathieu Desbrun
// Created:  08 02 01
// Modified: 09 19 01
////////////////////////////////////////

#ifndef _MAP_
#define _MAP_

#include<assert.h>

template<class T>
class CScalarSeg
{

private :

	T *m_pData;
	int m_w;

public :

	int w() { return m_w; }

  // life cycle
  CScalarSeg()
	{
		m_pData = NULL;
		m_w = 0;
	}
  CScalarSeg(int w)
	{
		m_pData = NULL;
		m_w = 0;
		alloc(w);
	}

  ~CScalarSeg()
	{
		free();
	}

	// alloc
	int alloc(int w,
		        int fill = 0,
						T v = 0)
	{
		free();
		m_pData = new T[w];
		assert(m_pData);
		if(!m_pData)
			return 0;
		if(fill)
			Fill(v);
		m_w = w;
		return 1;
	}

	void SetAt(int i,T v)
	{
		if(i<m_w)
			m_pData[i] = v;
	}

  T operator[](int i)
	{
		return m_pData[i];
	}

  T GetAt(int i)
	{
		return m_pData[i];
	}

	// free
	void free()
	{
		delete [] m_pData;
		m_pData = NULL;


		m_w = 0;
	}

	// fill
	void fill(T v)
	{
		for(int j=0;j<m_w;j++)
			m_pData[j] = v;
	}

	int nb(T v)
	{
		int nb = 0;
		for(int j=0;j<m_w;j++)
			nb += (m_pData[j] == v);
		return nb;
	}

	void error_diffusion(T threshold,
		                   T low,
											 T high)
	{
		for(int i=0;i<m_w;i++)
		{
			// set value
			T value = m_pData[i];
			if(value > threshold)
				m_pData[i] = high;
			else
				m_pData[i] = low;
			
			// error
			T error = value-m_pData[i];
			if(i<(m_w-1))
				m_pData[i+1] += error;
		}
	}
};

template<class T>
class CScalarMap
{

private :

	T **m_ppData;
	int m_w;
	int m_h;

public :

  // life cycle
  CScalarMap()
	{
		m_ppData = NULL;
		m_w = 0;
		m_h = 0;
	}
  CScalarMap(int w,int h)
	{
		m_ppData = NULL;
		m_w = 0;
		m_h = 0;
		alloc(w,h);
	}
  ~CScalarMap()
	{
		free();
	}

	// geometry
	int w() { return m_w; }
	int h() { return m_h; }
	int area() { return m_w*m_h; }

	// data
	T **get_fata() { return m_ppData; }

	T value(int x,int y)
	{
		assert(x >= 0 && x < m_w);
		assert(y >= 0 && y < m_h);
		return m_ppData[y][x];
	}
	/*
  T operator[][](int x,int y)
  {
		assert(x >= 0 && x < m_w);
		assert(y >= 0 && y < m_h);
		return m_ppData[y][x];
	} */
	void set(int x,int y,T v)
	{
		if(x >= 0  &&
			 x < m_w &&
			 y >= 0  &&
			 y < m_h)
			m_ppData[y][x] = v;
	}

	// free
	void free()
	{
		for(int i=0;i<m_h;i++)
			delete [] m_ppData[i];
		delete [] m_ppData;
		m_ppData = NULL;
		m_w = 0;
		m_h = 0;
	}
	
	// alloc
	int alloc(int w,int h,int fill = 0,T v = 0)
	{
		free();
		m_ppData = new T*[h];
		assert(m_ppData);
		if(!m_ppData)
			return 0;
		for(int i=0;i<h;i++)
		{
			m_ppData[i] = new T[w];
			assert(m_ppData[i]);
			if(!m_ppData[i])
				return 0;
		}
		m_w = w;
		m_h = h;
		if(fill)
			Fill(v);
		return 1;
	}


	// min

	T min()
	{
		if(area() == 0)
		{
			//TRACE("** CScalarMap::min -> empty image\n");
			return 0;
		}
		T min = m_ppData[0][0];
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				min = min(min,m_ppData[i][j]);
		return min;
	}

	T max()
	{
		if(area() == 0)
		{
			//TRACE("** CScalarMap::max -> empty image\n");
			return 0;
		}
		T max = m_ppData[0][0];
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				max = max(max,m_ppData[i][j]);
		return max;
	}

	T average()
	{
	  if(area() == 0)
	    return 0;
		return sum()/(T)area();
	}
	
	T sum()
	{
		T sum = 0;
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				sum += m_ppData[i][j];
		return sum;
	}

	T sum(int left,
		    int right,
				int top,
				int bottom)
	{
		T sum = 0;
		for(int i=top;i<=bottom;i++)
			for(int j=left;j<=right;j++)
				sum += m_ppData[i][j];
		return sum;
	}

	T sum_sqrt(int left,
		         int right,
						 int top,
						 int bottom)
	{
		T sum = 0;
		for(int i=top;i<=bottom;i++)
			for(int j=left;j<=right;j++)
				sum += sqrt(m_ppData[i][j]);
		return sum;
	}

	// stretch
	void stretch(T a,T b)
	{
		T min = Min();
		T max = Max();
		T range = max-min;
		if(range == (T)0 || a == b)
		{

			Fill(min);
			return;
		}
		if(a < b)
		{
			T nrange = b-a;
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
					m_ppData[i][j] = a+(m_ppData[i][j]-min)/range*nrange;
		}
	}

	// stretch
	void rescale(T m)
	{
		T max = max();
		if(max == (T)0)
		{
			Fill(m);
			return;
		}
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				m_ppData[i][j] = m_ppData[i][j]/max*m;
	}

	void fill(T v)
	{
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				m_ppData[i][j] = v;
	}

	void mult(CScalarMap<T> *pMap)
	{
		if(m_w == pMap->w() && m_h == pMap->h())
		{
			T ** ppData = pMap->GetData();
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
					m_ppData[i][j] *= ppData[i][j];
		}
	}
	void mult(T v)
	{
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				m_ppData[i][j] *= v;
	}

	void mask(CScalarMap<T> *pMap,
		        T src,
						T target)
	{
		if(m_w == pMap->w() && m_h == pMap->h())
		{
			T ** ppData = pMap->GetData();
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
					if(ppData[i][j] == src)
						m_ppData[i][j] = target;
		}
	}

	void pow(T power)
	{
		if(power != 1)
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
					m_ppData[i][j] = pow(m_ppData[i][j],power);
	}
	void log()
	{
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				m_ppData[i][j] = log(m_ppData[i][j]);
	}
	void add(T v)
	{
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
				m_ppData[i][j] += v;
	}

	void add(CScalarMap<T> *pMap)
	{
		if(m_w == pMap->w() && m_h == pMap->h())
		{
			T ** ppData = pMap->GetData();
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
					m_ppData[i][j] += ppData[i][j];
		}
	}

	void integrate(CScalarMap<T> *pMapIndex,
		             T *pArray)
	{
		if(m_w == pMapIndex->w() && m_h == pMapIndex->h())
		{
			T ** ppDataIndex = pMapIndex->GetData();
			for(int i=0;i<m_h;i++)
				for(int j=0;j<m_w;j++)
				{
					int index = (int)ppDataIndex[i][j];
					pArray[index] += m_ppData[i][j];
				}
		}
	}

	void filter_gaussian(int iter = 1)
	{
		CScalarMap<T> copy;
		copy.alloc(m_w,m_h);
		T ** ppData = copy.GetData();

		for(int k=0;k<iter;k++)
		{

			// filter...
			T d = 1.0/16.0;
			int i,j;
			for(i=1;i<(m_h-1);i++)
				for(j=1;j<(m_w-1);j++)
					ppData[i][j] = d*(m_ppData[i-1][j-1]+m_ppData[i-1][j+1]+
														m_ppData[i+1][j-1]+m_ppData[i+1][j+1]+
														2*(m_ppData[i][j-1]+m_ppData[i-1][j]+
															 m_ppData[i][j+1]+m_ppData[i+1][j])+
														4*m_ppData[i][j]);
			// v borders
			d = 1.0/12.0;
			for(i=1,j=0;i<(m_h-1);i++)
				ppData[i][j] = d*(m_ppData[i-1][j+1]+m_ppData[i+1][j+1]+
													2*(m_ppData[i-1][j]+m_ppData[i][j+1]+
														 m_ppData[i+1][j])+4*m_ppData[i][j]);
			for(i=1,j=m_w-1;i<(m_h-1);i++)
				ppData[i][j] = d*(m_ppData[i-1][j-1]+m_ppData[i+1][j-1]+
													2*(m_ppData[i-1][j]+m_ppData[i][j-1]+
														 m_ppData[i+1][j])+4*m_ppData[i][j]);
			// v borders
			for(j=1,i=0;j<(m_w-1);j++)
				ppData[i][j] = d*(m_ppData[i+1][j-1]+m_ppData[i+1][j+1]+
													2*(m_ppData[i][j-1]+m_ppData[i][j+1]+
														 m_ppData[i+1][j])+4*m_ppData[i][j]);
			for(j=1,i=m_h-1;j<(m_w-1);j++)
				ppData[i][j] = d*(m_ppData[i-1][j-1]+m_ppData[i-1][j+1]+
													2*(m_ppData[i][j-1]+m_ppData[i][j+1]+
														 m_ppData[i-1][j])+4*m_ppData[i][j]);
			// corners
			d = 1.0/9.0;
			ppData[0][0] = d*(m_ppData[1][1]+2*m_ppData[1][0]+
												2*m_ppData[0][1]+4*m_ppData[0][0]);
			ppData[0][m_w-1] = d*(m_ppData[1][m_w-2]+2*m_ppData[0][m_w-2]+
																2*m_ppData[1][m_w-1]+4*m_ppData[0][m_w-1]);
			ppData[m_h-1][0] = d*(m_ppData[m_h-2][1]+2*m_ppData[m_h-2][0]+
																2*m_ppData[m_h-1][1]+4*m_ppData[m_h-1][0]);
			ppData[m_h-1][m_w-1] = d*(m_ppData[m_h-2][m_w-2]+2*m_ppData[m_h-2][m_w-1]+
																2*m_ppData[m_h-1][m_w-2]+4*m_ppData[m_h-1][m_w-1]);


			copy(&copy);

		}
	}
	void copy(CScalarMap<T> *pMap)
	{
		// conditional alloc
		if(m_w != pMap->w() ||
			 m_h != pMap->h())
			alloc(pMap->w(),pMap->h());
		T ** ppData = pMap->GetData();
		int i;
		for(i=0;i<m_h;i++)
			memcpy(m_ppData[i],ppData[i],m_w*sizeof(T));
	}

	void copy_mirror(CScalarMap<T> *pMap)
	{
		int w = pMap->w();
		int h = pMap->h();

		// conditional alloc
		if(m_w != w ||
			 m_h != 2*h) // twice the height
			alloc(w,2*h);
		T ** ppData = pMap->GetData();

		// copy upper half
		int i,j;
		for(i=0;i<h;i++)
			for(j=0;j<m_w;j++)
				m_ppData[i+h][j] = ppData[i][j];
		for(i=0;i<h;i++)
			for(j=0;j<m_w;j++)
				m_ppData[i][j] = ppData[h-1-i][j];
	}

	void copy_upper_half(CScalarMap<T> *pMap)
	{
		int w = pMap->w();
		int h = pMap->h();

		// conditional alloc
		if(m_h == h/2) // half the height
		{
			T ** ppData = pMap->GetData();
			// copy upper half
			int i,j;
			for(i=0;i<m_h;i++)
				for(j=0;j<m_w;j++)
					m_ppData[i][j] = ppData[i+m_h][j];
		}
	}

	// add ramp
	void add_ramp(int size,
		            T min,
							  T max)
	{
		if(area() == 0)
			return;
		T range = max-min;
		int i,j;
		for(i=0;i<m_h;i++)
		{
			T v = min+(T)i/(T)(m_h-1)*range;
			for(j=(m_w-size);j<m_w;j++)
				m_ppData[i][j] = v;
		}
	}

	void saturate_high(T v)
	{
	  int i,j;
		for(i=0;i<m_h;i++)
			for(j=0;j<m_w;j++)
				if(m_ppData[i][j]>v)
					m_ppData[i][j] = v;
	}

	void inflate_mean(T v)
	{
	  int i,j,k,l;
		for(i=0;i<m_h;i++)
			for(j=0;j<m_w;j++)
				if(m_ppData[i][j] == v)
				{
					T sum = 0;
					int nb = 0;
					for(k=(i-1);k<=(i+1);k++)
						for(l=(j-1);l<=(j+1);l++)
							if(k>=0 && k<m_h && l>=0 && l<m_w)
								if(m_ppData[k][l] == v)
								{
									sum += m_ppData[k][l];
									nb++;
								}
					if(nb > 0)
					{
						sum /= (T)nb;
						m_ppData[i][j] = sum;
					}
				}
	}

	void error_diffusion(T threshold,
		                   T low,
											 T high)
	{
		T ratio = 1.0/16.0;
		int i,j;
		for(i=0;i<m_h;i++)
			for(j=0;j<m_w;j++)
			{
				// set value
				T value = m_ppData[i][j];
				if(value > threshold)
					m_ppData[i][j] = high;
				else
					m_ppData[i][j] = low;

				// error
				T error = value-m_ppData[i][j];

				// diffuse error
				T tmp = ratio*error;
				if(j>=3 &&
					 j<=(m_w-2) &&
					 i<=(m_h-2))
				{
					m_ppData[i][j+1]   +=  8*tmp;
					m_ppData[i+1][j-3] +=    tmp;
					m_ppData[i+1][j-2] +=    tmp;
					m_ppData[i+1][j-1] +=  2*tmp;
					m_ppData[i+1][j]   +=  4*tmp;
				}
			}
	}

	T threshold(int i,
		          int j,
							T threshold,
							T low,
							T high,
							T &value)
	{
		value = m_ppData[i][j];
		if(value > threshold)
			m_ppData[i][j] = high;
		else
			m_ppData[i][j] = low;
		return (value-m_ppData[i][j]); // return error
	}

	void error_diffusion_victor(T threshold,
		                          T low,
														  T high)
	{
		static float vctab[256*3] = {
		0.722222,	0,	0.277778,	
		0.722222,	0,	0.277778,	
		0.677419,	0,	0.322581,	
		0.636364,	0,	0.363636,	
		0.615385,	0,	0.384615,	
		0.602564,	0.0384615,	0.358974,	
		0.589744,	0.0769231,	0.333333,	
		0.576923,	0.115385,	0.307692,	
		0.564103,	0.153846,	0.282051,	
		0.551282,	0.192308,	0.25641,	
		0.538462,	0.230769,	0.230769,	
		0.535256,	0.239316,	0.225427,	
		0.532051,	0.247863,	0.220085,	
		0.528846,	0.25641,	0.214744,	
		0.525641,	0.264957,	0.209402,	
		0.522436,	0.273504,	0.20406,	
		0.519231,	0.282051,	0.198718,	
		0.516026,	0.290598,	0.193376,	
		0.512821,	0.299145,	0.188034,	
		0.509615,	0.307692,	0.182692,	
		0.50641,	0.316239,	0.17735,	
		0.503205,	0.324786,	0.172009,	
		0.5,	0.333333,	0.166667,	
		0.496753,	0.329004,	0.174242,	
		0.493506,	0.324675,	0.181818,	
		0.49026,	0.320346,	0.189394,	
		0.487013,	0.316017,	0.19697,	
		0.483766,	0.311688,	0.204545,	
		0.480519,	0.307359,	0.212121,	
		0.477273,	0.30303,	0.219697,	
		0.474026,	0.298701,	0.227273,	
		0.470779,	0.294372,	0.234848,	
		0.467532,	0.290043,	0.242424,	
		0.464286,	0.285714,	0.25,	
		0.461039,	0.281385,	0.257576,	
		0.457792,	0.277056,	0.265152,	
		0.454545,	0.272727,	0.272727,	
		0.456169,	0.280844,	0.262987,	
		0.457792,	0.288961,	0.253247,	
		0.459416,	0.297078,	0.243506,	
		0.461039,	0.305195,	0.233766,	
		0.462662,	0.313312,	0.224026,	
		0.464286,	0.321429,	0.214286,	
		0.465909,	0.329545,	0.204545,	
		0.467532,	0.337662,	0.194805,	
		0.469156,	0.345779,	0.185065,	
		0.470779,	0.353896,	0.175325,	
		0.472403,	0.362013,	0.165584,	
		0.474026,	0.37013,	0.155844,	
		0.475649,	0.378247,	0.146104,	
		0.477273,	0.386364,	0.136364,	
		0.478896,	0.394481,	0.126623,	
		0.480519,	0.402597,	0.116883,	
		0.482143,	0.410714,	0.107143,	
		0.483766,	0.418831,	0.0974026,	
		0.48539,	0.426948,	0.0876623,	
		0.487013,	0.435065,	0.0779221,	
		0.488636,	0.443182,	0.0681818,	
		0.49026,	0.451299,	0.0584416,	
		0.491883,	0.459416,	0.0487013,	
		0.493506,	0.467532,	0.038961,	
		0.49513,	0.475649,	0.0292208,	
		0.496753,	0.483766,	0.0194805,	
		0.498377,	0.491883,	0.00974026,	
		0.5,	0.5,	0,	
		0.485577,	0.504808,	0.00961538,	
		0.471154,	0.509615,	0.0192308,	
		0.456731,	0.514423,	0.0288462,	
		0.442308,	0.519231,	0.0384615,	
		0.427885,	0.524038,	0.0480769,	
		0.413462,	0.528846,	0.0576923,	
		0.399038,	0.533654,	0.0673077,	
		0.384615,	0.538462,	0.0769231,	
		0.441026,	0.464103,	0.0948718,	
		0.497436,	0.389744,	0.112821,	
		0.553846,	0.315385,	0.130769,	
		0.610256,	0.241026,	0.148718,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.65,	0.18,	0.17,	
		0.633333,	0.193333,	0.173333,	
		0.616667,	0.206667,	0.176667,	
		0.6,	0.22,	0.18,	
		0.583333,	0.233333,	0.183333,	
		0.566667,	0.246667,	0.186667,	
		0.55,	0.26,	0.19,	
		0.533333,	0.273333,	0.193333,	
		0.516667,	0.286667,	0.196667,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.508333,	0.293333,	0.198333,	
		0.516667,	0.286667,	0.196667,	
		0.525,	0.28,	0.195,	
		0.533333,	0.273333,	0.193333,	
		0.541667,	0.266667,	0.191667,	
		0.55,	0.26,	0.19,	
		0.558333,	0.253333,	0.188333,	
		0.566667,	0.246667,	0.186667,	
		0.575,	0.24,	0.185,	
		0.583333,	0.233333,	0.183333,	
		0.591667,	0.226667,	0.181667,	
		0.6,	0.22,	0.18,	
		0.608333,	0.213333,	0.178333,	
		0.616667,	0.206667,	0.176667,	
		0.625,	0.2,	0.175,	
		0.633333,	0.193333,	0.173333,	
		0.641667,	0.186667,	0.171667,	
		0.65,	0.18,	0.17,	
		0.658333,	0.173333,	0.168333,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.658333,	0.173333,	0.168333,	
		0.65,	0.18,	0.17,	
		0.641667,	0.186667,	0.171667,	
		0.633333,	0.193333,	0.173333,	
		0.625,	0.2,	0.175,	
		0.616667,	0.206667,	0.176667,	
		0.608333,	0.213333,	0.178333,	
		0.6,	0.22,	0.18,	
		0.591667,	0.226667,	0.181667,	
		0.583333,	0.233333,	0.183333,	
		0.575,	0.24,	0.185,	

		0.566667,	0.246667,	0.186667,	
		0.558333,	0.253333,	0.188333,	
		0.55,	0.26,	0.19,	
		0.541667,	0.266667,	0.191667,	
		0.533333,	0.273333,	0.193333,	
		0.525,	0.28,	0.195,	
		0.516667,	0.286667,	0.196667,	
		0.508333,	0.293333,	0.198333,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.5,	0.3,	0.2,	
		0.516667,	0.286667,	0.196667,	
		0.533333,	0.273333,	0.193333,	
		0.55,	0.26,	0.19,	
		0.566667,	0.246667,	0.186667,	
		0.583333,	0.233333,	0.183333,	
		0.6,	0.22,	0.18,	
		0.616667,	0.206667,	0.176667,	
		0.633333,	0.193333,	0.173333,	
		0.65,	0.18,	0.17,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.666667,	0.166667,	0.166667,	
		0.610256,	0.241026,	0.148718,	
		0.553846,	0.315385,	0.130769,	
		0.497436,	0.389744,	0.112821,	
		0.441026,	0.464103,	0.0948718,	
		0.384615,	0.538462,	0.0769231,	
		0.399038,	0.533654,	0.0673077,	
		0.413462,	0.528846,	0.0576923,	
		0.427885,	0.524038,	0.0480769,	
		0.442308,	0.519231,	0.0384615,	
		0.456731,	0.514423,	0.0288462,	
		0.471154,	0.509615,	0.0192308,	
		0.485577,	0.504808,	0.00961538,	
		0.5,	0.5,	0,	
		0.498377,	0.491883,	0.00974026,	
		0.496753,	0.483766,	0.0194805,	
		0.49513,	0.475649,	0.0292208,	
		0.493506,	0.467532,	0.038961,	
		0.491883,	0.459416,	0.0487013,	
		0.49026,	0.451299,	0.0584416,	
		0.488636,	0.443182,	0.0681818,	
		0.487013,	0.435065,	0.0779221,	
		0.48539,	0.426948,	0.0876623,	
		0.483766,	0.418831,	0.0974026,	
		0.482143,	0.410714,	0.107143,	
		0.480519,	0.402597,	0.116883,	
		0.478896,	0.394481,	0.126623,	
		0.477273,	0.386364,	0.136364,	
		0.475649,	0.378247,	0.146104,	
		0.474026,	0.37013,	0.155844,	
		0.472403,	0.362013,	0.165584,	
		0.470779,	0.353896,	0.175325,	
		0.469156,	0.345779,	0.185065,	
		0.467532,	0.337662,	0.194805,	
		0.465909,	0.329545,	0.204545,	
		0.464286,	0.321429,	0.214286,	
		0.462662,	0.313312,	0.224026,	
		0.461039,	0.305195,	0.233766,	

		0.459416,	0.297078,	0.243506,	
		0.457792,	0.288961,	0.253247,	
		0.456169,	0.280844,	0.262987,	
		0.454545,	0.272727,	0.272727,	
		0.457792,	0.277056,	0.265152,	
		0.461039,	0.281385,	0.257576,	
		0.464286,	0.285714,	0.25,	
		0.467532,	0.290043,	0.242424,	
		0.470779,	0.294372,	0.234848,	
		0.474026,	0.298701,	0.227273,	
		0.477273,	0.30303,	0.219697,	
		0.480519,	0.307359,	0.212121,	
		0.483766,	0.311688,	0.204545,	
		0.487013,	0.316017,	0.19697,	
		0.49026,	0.320346,	0.189394,	
		0.493506,	0.324675,	0.181818,	
		0.496753,	0.329004,	0.174242,	
		0.5,	0.333333,	0.166667,	
		0.503205,	0.324786,	0.172009,	
		0.50641,	0.316239,	0.17735,	
		0.509615,	0.307692,	0.182692,	
		0.512821,	0.299145,	0.188034,	
		0.516026,	0.290598,	0.193376,	
		0.519231,	0.282051,	0.198718,	
		0.522436,	0.273504,	0.20406,	
		0.525641,	0.264957,	0.209402,	
		0.528846,	0.25641,	0.214744,	
		0.532051,	0.247863,	0.220085,	
		0.535256,	0.239316,	0.225427,	
		0.538462,	0.230769,	0.230769,	
		0.551282,	0.192308,	0.25641,	
		0.564103,	0.153846,	0.282051,	
		0.576923,	0.115385,	0.307692,	
		0.589744,	0.0769231,	0.333333,	
		0.602564,	0.0384615,	0.358974,	
		0.615385,	0,	0.384615,	
		0.636364,	0,	0.363636,	
		0.677419,	0,	0.322581,	
		0.722222,	0,	0.277778,	
		0.722222,	0,	0.277778};

		int i,j;
		i = j = 0;
		T value = 0;
		T error = 0;
		for(i=0;i<(m_h-1);i++)
		{
			if(i%2 == 0) // even lines
			{
				// first pixel
				j = 0;
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i][j+1] += error;

				for(j=1;j<(m_w-1);j++)
				{
					error = threshold(i,j,threshold,low,high,value);
					int index = (int)(value*255.0f);
					index = (index<0) ? 0 : (index>255) ? 255 : index;
					m_ppData[i][j+1]   += error*vctab[3*index];
					m_ppData[i+1][j-1] += error*vctab[3*index+1];
					m_ppData[i+1][j]   += error*vctab[3*index+2];
				}

				// last pixel
				j = m_w-1;
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i+1][j] += error;
			}
			else // odd lines
			{
				// last pixel
				j = m_w-1;
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i][j-1] += error;

				for(j=(m_w-2);j>=1;j--)
				{
					error = threshold(i,j,threshold,low,high,value);
					// diffuse error (backward)
					int index = (int)(value*255.0f);
					index = (index<0) ? 0 : (index>255) ? 255 : index;
					m_ppData[i][j-1]   += error*vctab[3*index];
					m_ppData[i+1][j+1] += error*vctab[3*index+1];
					m_ppData[i+1][j]   += error*vctab[3*index+2];

				}

				// first pixel
				j = 0;
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i+1][j] += error;
			}
		}

		// last line
		i = m_h-1;
		if(i%2 == 0) // even line
		{
			for(j=0;j<(m_w-1);j++)
			{
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i][j+1]   += error;
			}

			// last pixel (no error diffusion)
			j = m_w-1;
			error = threshold(i,j,threshold,low,high,value);
		}
		else // odd line
		{
			for(j=(m_w-1);j>=1;j--)
			{

				// set value
				error = threshold(i,j,threshold,low,high,value);
				m_ppData[i][j-1] += error;
			}

			// first pixel (no error diffusion)
			j = 0;
			error = threshold(i,j,threshold,low,high,value);
		}
	}

	int nb_points(T v)
	{
		int nb = 0;
		int i,j;
		for(i=0;i<m_h;i++)
			for(j=0;j<m_w;j++)
				nb += (m_ppData[i][j] == v);
		return nb;
	}

	void decimate()
	{
		CScalarMap<T> map;
		map.Copy(this);
		T **ppData = map.GetData();
		int w = map.w();
		int h = map.h();
		alloc(w/2,h/2);
		int i,j;
		for(i=0;i<m_h;i++)
			for(j=0;j<m_w;j++)
				m_ppData[i][j] = ppData[2*i][2*j];
	}

	void fill(int l,
		        int r,
						int t,
						int b,
						T v)

	{
		if(l < 0) return;
		if(r < 0) return;
		if(t < 0) return;
		if(b < 0) return;

		if(l >= m_w) return;
		if(r >= m_w) return;
		if(t >= m_h) return;
		if(b >= m_h) return;

		int i,j;
		for(i=l;i<=r;i++)
			for(j=t;j<=b;j++)
				m_ppData[j][i] = v;
	}

	// stroke
	void stroke(int Ax,
		          int Ay,
							int Bx,

							int By,
							T v)
	{
		if(Ax < 0) return;
		if(Bx < 0) return;
		if(Ay < 0) return;
		if(By < 0) return;

		if(Ax >= m_w) return;
		if(Bx >= m_w) return;
		if(Ay >= m_h) return;
		if(By >= m_h) return;

		// store the change in X and Y 
		// of the line endpoints
		int dX = abs(Bx-Ax);	
		int dY = abs(By-Ay);
			
		int Xincr, Yincr;
		if (Ax > Bx) { Xincr=-1; } else { Xincr=1; }	// which direction in X?
		if (Ay > By) { Yincr=-1; } else { Yincr=1; }	// which direction in Y?
			
		if(dX >= dY)	// if X is the independent variable
		{           
			int dPr 	= dY<<1;         // amount to increment decision if right is chosen (always)
			int dPru 	= dPr - (dX<<1); // amount to increment decision if up is chosen
			int P 		= dPr - dX;      // decision variable start value

			// process each point in the line 
			// one at a time (just use dX)
			for (; dX>=0; dX--)            
			{
				m_ppData[Ay][Ax] = v; // plot the pixel
				if(P > 0)               // is the pixel going right AND up?
				{ 
					Ax+=Xincr;	       // increment independent variable

					Ay+=Yincr;         // increment dependent variable
					P+=dPru;           // increment decision (for up)
				}
				else                 // is the pixel just going right?
				{
					Ax+=Xincr;         // increment independent variable
					P+=dPr;            // increment decision (for right)
				}
			}		
		}
		else              // if Y is the independent variable
		{
			int dPr 	= dX<<1;           // amount to increment decision if right is chosen (always)
			int dPru 	= dPr - (dY<<1);   // amount to increment decision if up is chosen
			int P 		= dPr - dY;  // decision variable start value

			for(; dY>=0; dY--)            // process each point in the line one at a time (just use dY)
			{
				m_ppData[Ay][Ax] = v; // plot the pixel
				if (P > 0)               // is the pixel going up AND right?
				{ 
					Ax+=Xincr;         // increment dependent variable
					Ay+=Yincr;         // increment independent variable
					P+=dPru;           // increment decision (for up)
				}
				else                     // is the pixel just going up?
				{
					Ay+=Yincr;         // increment independent variable
					P+=dPr;            // increment decision (for right)
				}
			}		
		}		
	}


	/*
	void Transfert(CTransfertFunction *pTf)
	{
		float max = pTf->Size()-2;
		for(int i=0;i<m_h;i++)
			for(int j=0;j<m_w;j++)
			{
				float v = m_ppData[i][j]*max;
				int start = (int)floor(v);
				int end = start+1;
				float r = v-(float)start;
				m_ppData[i][j] = r*pTf->f(start)+(1-r)*pTf->f(end);
			}
	}*/
};

#endif // _MAP_

