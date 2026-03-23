#include<iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<stdint.h>
#include<vector>
#include<map>
#include<sstream>

using namespace std;

#define SIZE 4

typedef struct
{
	vector<string> name;
	vector<string> strand;
} sequence;

typedef map<string,float> map_string;

/****************************************************************************
 *
 * Construct the random number generator based on Mersenne Twister algorithm
 *
 ****************************************************************************/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

void write_seqs(char* file,sequence seqs)
{
	char* csv=".csv";
	char* filename=strcat(file,csv);
	ofstream out(filename);//,ios::app);
	string temp;
	cout<<seqs.name.size()<<endl;
	{
		for(int j=0;j<seqs.name.size();j++)
		{
			out.write(seqs.name[j].c_str(),seqs.name[j].length());
			out<<endl;
			out.write(seqs.strand[j].c_str(),seqs.strand[j].length());
			out<<endl;
		}
	}
	out.close();
	return;
}


void load_data(char* file,sequence &seqs)
{
	ifstream in(file);
	istringstream istr;
	string s,s_q;
	int i=1;    
    while(getline(in,s))
    {
        
        i=i+1;
	if (i%2==0)
	{
        	seqs.name.push_back(s);
		cout<<s<<endl;
		istr.clear();
	}
	else
	{
        	seqs.strand.push_back(s);
        	istr.clear();
	}
    }
    return;
}













/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
		(1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
		+ init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
		- i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
	
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t genrand_int32(void) {
	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */
	
	if (mti >= N) { /* generate N words at one time */
		int kk;
		
		if (mti == N+1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */
		
		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
		
		mti = 0;
	}
	
	y = mt[mti++];
	
	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);
	
	return y;
}

uint64_t genrand_int64(void) {
	uint64_t x, y;
	
	x = genrand_int32(); y = genrand_int32();
	return (x<<32)|y;
}

void init_mersenne(void) {
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	init_by_array(init, length);
}

float genrand_real2(void)
{
	return genrand_int64()*(1.0/(4294967296.0*4294967296.0)); 
}

float ranf(void)
{    
	init_genrand((unsigned)time(NULL)+rand());
	return (genrand_real2());
}

float genunf(float low,float high)
{
	static float genunf;
    genunf = low+(high-low)*ranf(); 
	
    return genunf;
}

/****************************************************************************
 *
 * Function: partition
 * 1. Used in quickSort algorithm
 *
 ****************************************************************************/

template <class T>
int partition(T a[], int start, int stop, int id[])
{
        int temp_id, up=start, down=stop-1;
        T temp_value, part=a[stop];
        if(stop<=start) return start;

        while(true)
        {    
                while(a[up]<part) up++;
                while(part<a[down] && (up<down)) down--;

                if(up>=down) break;

                temp_value=a[up];  a[up]=a[down]; a[down]=temp_value;
                temp_id=id[up]; id[up]=id[down]; id[down]=temp_id;

                up++; down--;
        }    

        temp_value=a[up]; a[up]=a[stop]; a[stop]=temp_value;
        temp_id=id[up]; id[up]=id[stop]; id[stop]=temp_id;
        return up;
}

/****************************************************************************
 *
 * Function: quickSort
 * 1. Sort the integers saved in an array in increasing order
 * 2. Call the function partition
 *
 ****************************************************************************/

template <class T>
void quickSort(T a[], int start, int stop, int id[])
{
        int i;
        if(stop<=start) return;

        i=partition(a,start,stop,id);
        quickSort(a,start,i-1,id);
        quickSort(a,i+1,stop,id);
}

/****************************************************************************
 *
 * Function: genmulone
 * 1. Generate a random variable following the multinomial distribution
 *
 ****************************************************************************/

int genmulone(float *p_input, long ncat)
{
	int i;
	float *p = NULL;
	p = new float [ncat];
	for(i=0;i<ncat;i++) p[i]=p_input[i];
	
	int*q = NULL;
	q = new int [ncat];
	for(i=0;i<ncat;i++) q[i]=i;
	
	quickSort(p, 0, ncat - 1, q);
	
	float prob=genunf(0, 1);
	int outcome=ncat - 1;
	while (outcome>=0  && prob > p[outcome])
	{
		prob -= p[outcome];
		outcome --;
	}
	
	return q[outcome];
}

/******************************************************************************************************
 *
 * Function: markov function
 * Used to generate background sequences
 *
 ******************************************************************************************************/

void markov(sequence& f_in, sequence& f_bg, int flag)
{

/******************************************************************************************************
 *
 * Construct the hash table for various k-mers (k=1, 2, 3, 4)
 *
 ******************************************************************************************************/

	char* alphabet;	// the alphabet for DNA nucleotides A, C, G, T
	alphabet=new char [SIZE];

	alphabet[0]='A';
	alphabet[1]='C';
	alphabet[2]='G';
	alphabet[3]='T';

	map_string hash_onebase, hash_twobase, hash_threebase, hash_fourbase;

	if(flag == 0 || flag == 1 || flag == 2 || flag == 3)
	{
		string onebase;	// the string to contain A, C, G, T

		for(int i = 0; i < SIZE; i ++)
		{
			onebase.clear();
			onebase.push_back(alphabet[i]);
			hash_onebase[onebase] = 0;
		}	// initialization for the hash table
	}
	else
	{
		cerr<<"Error: The input flag parameter is wrong! Please check it again!"<<endl;
		exit(0);
	}
	
	if(flag == 1 || flag == 2 || flag == 3)
	{
		string twobase;

		for(int i = 0; i < SIZE; i ++)
		{
			for(int j = 0; j < SIZE; j ++)
			{
				twobase.clear();
				twobase.push_back(alphabet[i]);
				twobase.push_back(alphabet[j]);
				hash_twobase[twobase] = 0;
			}
		}	// initialization for the hash table
	}
	
	if(flag == 2 || flag == 3)
	{
		string threebase;

		for(int i = 0; i < SIZE; i ++)
		{
			for(int j = 0; j < SIZE; j ++)
			{
				for(int k = 0; k < SIZE; k ++)
				{
					threebase.clear();
					threebase.push_back(alphabet[i]);
					threebase.push_back(alphabet[j]);
					threebase.push_back(alphabet[k]);
					hash_threebase[threebase]=0;
				}
			}
		}	// initialization for the hash table
	}
	
	if(flag == 3)
	{
		string fourbase;

		for(int i = 0; i < SIZE; i ++)
		{
			for(int j = 0; j < SIZE; j ++)
			{
				for(int k = 0; k < SIZE; k ++)
				{
					for(int l = 0; l < SIZE; l ++)
					{
						fourbase.clear();
						fourbase.push_back(alphabet[i]);
						fourbase.push_back(alphabet[j]);
						fourbase.push_back(alphabet[k]);
						fourbase.push_back(alphabet[l]);
						hash_fourbase[fourbase]=0;
					}
				}
			}
		}
	}

/******************************************************************************************************
 *
 * Count the occurrence numbers of various k-mers (k=1, 2, 3, 4)
 *
 ******************************************************************************************************/

	if(flag == 0)
	{
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			for(int j = 0; j < f_in.strand[i].size(); j ++)
			{
				if(f_in.strand[i][j] != 'N')
				{
					hash_onebase[f_in.strand[i].substr(j,1)]++;
				}
			}
		}
	}
	else if(flag == 1)
	{
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			for(int j = 0; j < f_in.strand[i].size(); j ++)
			{
				if(f_in.strand[i][j] != 'N')
				{
					hash_onebase[f_in.strand[i].substr(j,1)] ++;
				}

				if(j < f_in.strand[i].size() - 1 && f_in.strand[i].substr(j,2).find('N') 
						== f_in.strand[i].substr(j,2).npos)
				{
					hash_twobase[f_in.strand[i].substr(j,2)] ++;
				}
			}
		}
	}
	else if(flag == 2)
	{
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			for(int j = 0; j < f_in.strand[i].size(); j ++)
			{
				if(f_in.strand[i][j] != 'N')
				{
					hash_onebase[f_in.strand[i].substr(j,1)] ++;
				}

				if(j < f_in.strand[i].size() - 1 && f_in.strand[i].substr(j,2).find('N') 
						== f_in.strand[i].substr(j,2).npos)
				{
					hash_twobase[f_in.strand[i].substr(j,2)] ++;
				}

				if(j < f_in.strand[i].size() - 2 && f_in.strand[i].substr(j,3).find('N') 
						== f_in.strand[i].substr(j,3).npos)
				{
					hash_threebase[f_in.strand[i].substr(j,3)] ++;
				}
			}
		}
	}
	else if(flag == 3)
	{
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			for(int j = 0; j < f_in.strand[i].size(); j ++)
			{
				if(f_in.strand[i][j] != 'N')
				{
					hash_onebase[f_in.strand[i].substr(j,1)] ++;
				}
				if(j < f_in.strand[i].size() - 1 && f_in.strand[i].substr(j,2).find('N') 
						== f_in.strand[i].substr(j,2).npos)
				{
					hash_twobase[f_in.strand[i].substr(j,2)] ++;
				}

				if(j < f_in.strand[i].size() - 2 && f_in.strand[i].substr(j,3).find('N') 
						== f_in.strand[i].substr(j,3).npos)
				{
					hash_threebase[f_in.strand[i].substr(j,3)] ++;
				}

				if(j < f_in.strand[i].size() - 3 && f_in.strand[i].substr(j,4).find('N') 
						== f_in.strand[i].substr(j,4).npos)
				{
					hash_fourbase[f_in.strand[i].substr(j,4)]++;
				}
			}
		}
	}

	cout<<"k-mer counting has been completed.\n";

/******************************************************************************************************
 *
 * Transform k-mer occurrences into k-mer frequencies (k = 0, 1, 2, 3)
 *
 ******************************************************************************************************/

	map<string, float>::iterator iter;
	float sum = 0;	// the sum of all nucleotides in sequences

	if(flag == 3)
	{
		for(iter = hash_fourbase.begin(); iter != hash_fourbase.end(); iter ++)
		{
			iter->second /= hash_threebase[iter->first.substr(0, 3)];
		}

		for(iter = hash_threebase.begin(); iter != hash_threebase.end(); iter ++)
		{
			iter->second /= hash_twobase[iter->first.substr(0, 2)];
		}

		for(iter = hash_twobase.begin(); iter != hash_twobase.end(); iter ++)
		{
			iter->second /= hash_onebase[iter->first.substr(0, 1)];
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			sum += iter->second;
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			iter->second /= sum;
		}
	}
	else if(flag == 2)
	{
		for(iter = hash_threebase.begin(); iter != hash_threebase.end(); iter ++)
		{
			iter->second /= hash_twobase[iter->first.substr(0, 2)];
		}

		for(iter = hash_twobase.begin(); iter != hash_twobase.end(); iter ++)
		{
			iter->second /= hash_onebase[iter->first.substr(0, 1)];
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			sum += iter->second;
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			iter->second /= sum;
		}
	}
	else if(flag == 1)
	{
		for(iter = hash_twobase.begin(); iter != hash_twobase.end(); iter ++)
		{
			iter->second /= hash_onebase[iter->first.substr(0, 1)];
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			sum += iter->second;
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			iter->second /= sum;
		}
	}
	else if(flag == 0)
	{
		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			sum += iter->second;
		}

		for(iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
		{
			iter->second /= sum;
		}
	}

	cout << "Finished claculating the k-mer frequencies.\n";

/******************************************************************************************************
 *
 * Generate the simulated sequences with the same number and lengths
 *
 ******************************************************************************************************/

	float* values0, *values1, *values2, *values3;
	values0 = new float [SIZE];
	values1 = new float [SIZE];
	values2 = new float [SIZE];
	values3 = new float [SIZE];
	string tempstr;

	for(int i = 0; i < f_in.strand.size(); i ++)
	{
		ostringstream segment_fd;
		segment_fd << (i + 1);
		//f_bg.name.push_back(segment_fd.str());
		f_bg.name.push_back(f_in.name[i]);
		f_bg.strand.push_back(tempstr);
	}

	if(flag == 0)
	{
		for(int i = 0; i < SIZE; i++)
		{
			tempstr.clear();
			tempstr.push_back(alphabet[i]);
			values0[i] = hash_onebase[tempstr];
		}
		for(int i = 0;i < f_in.strand.size(); i ++)
		{
			for(int j = 0; j < f_in.strand[i].size(); j ++)
			{
				f_bg.strand[i].push_back(alphabet[genmulone(values0, SIZE)]);
			}
		}
	}
	else if(flag == 1)
	{
		for(int i = 0; i < SIZE; i ++)
		{
			tempstr.clear();
			tempstr.push_back(alphabet[i]);
			values0[i] = hash_onebase[tempstr];
		}
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			f_bg.strand[i].push_back(alphabet[genmulone(values0, SIZE)]);

			for(int j = 1; j < f_in.strand[i].size(); j ++)
			{
				for(int k = 0; k < SIZE; k ++)
				{
					tempstr.clear();
					tempstr.push_back(f_bg.strand[i][j - 1]);
					tempstr.push_back(alphabet[k]);
					values1[k] = hash_twobase[tempstr];
				}

				f_bg.strand[i].push_back(alphabet[genmulone(values1, SIZE)]);
			}
		}
	}
	else if(flag == 2)
	{
		for(int i = 0; i < SIZE; i ++)
		{
			tempstr.clear();
			tempstr.push_back(alphabet[i]);
			values0[i] = hash_onebase[tempstr];
		}
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			f_bg.strand[i].push_back(alphabet[genmulone(values0, SIZE)]);

			for(int k = 0; k < SIZE; k ++)
			{
				tempstr.clear();
				tempstr.push_back(f_bg.strand[i][0]);
				tempstr.push_back(alphabet[k]);
				values1[k] = hash_twobase[tempstr];
			}

			f_bg.strand[i].push_back(alphabet[genmulone(values1, SIZE)]);

			for(int j = 2; j < f_in.strand[i].size(); j ++)
			{
				tempstr.clear();
				tempstr.push_back(f_bg.strand[i][j - 2]);
				tempstr.push_back(f_bg.strand[i][j - 1]);
				for(int k = 0; k < SIZE; k ++)
				{
					tempstr.push_back(alphabet[k]);
					values2[k] = hash_threebase[tempstr];
				}
				f_bg.strand[i].push_back(alphabet[genmulone(values2, SIZE)]);
			}
		}
	}
	else if(flag == 3)
	{
		for(int i = 0; i < SIZE; i ++)
		{
			tempstr.clear();
			tempstr.push_back(alphabet[i]);
			values0[i] = hash_onebase[tempstr];
		}
		for(int i = 0; i < f_in.strand.size(); i ++)
		{
			f_bg.strand[i].push_back(alphabet[genmulone(values0,SIZE)]);
			for(int k = 0; k < SIZE; k ++)
			{
				tempstr.clear();
				tempstr.push_back(f_bg.strand[i][0]);
				tempstr.push_back(alphabet[k]);
				values1[k] = hash_twobase[tempstr];
			}
			f_bg.strand[i].push_back(alphabet[genmulone(values1, SIZE)]);

			for(int k = 0; k < SIZE; k ++)
			{
				tempstr.clear();
				tempstr.push_back(f_bg.strand[i][0]);
				tempstr.push_back(f_bg.strand[i][1]);
				tempstr.push_back(alphabet[k]);
				values2[k] = hash_twobase[tempstr];
			}
			f_bg.strand[i].push_back(alphabet[genmulone(values2, SIZE)]);

			for(int j = 3; j < f_in.strand[i].size(); j ++)
			{
				for(int k = 0; k < SIZE; k ++)
				{
					tempstr.clear();
					tempstr.push_back(f_bg.strand[i][j - 3]);
					tempstr.push_back(f_bg.strand[i][j - 2]);
					tempstr.push_back(f_bg.strand[i][j - 1]);
					tempstr.push_back(alphabet[k]);
					values3[k] = hash_fourbase[tempstr];
				}
				f_bg.strand[i].push_back(alphabet[genmulone(values3, SIZE)]);
			}
		}
	}

	cout << "Finished generating background sequences.\n";
}
