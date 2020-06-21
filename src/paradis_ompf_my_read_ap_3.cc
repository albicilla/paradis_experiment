#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <sstream>
#include <thread>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
using namespace std;

#define FOR(i,a,b) for(int i=a;i<b;i++)
#define rep(i,b) FOR(i,0,b)
#define INF 1e9
#define dump(x) cerr<<#x<<"="<<x<<endl
#define ALL(a) (a).begin(),(a).end()
#define EACH(e,v) for(auto& e:v)
#define SORT(v) sort(ALL(v))
#define PERM(v) SORT(v);for(bool c##p=1;c##p;c##p=next_permutation(ALL(v)))
#define printArr(x) for(auto itr:x)cerr<<"="<<itr<<endl

long long DATASIZE=1000LL*1000LL*1000LL*10LL;
long long SIZE=DATASIZE/400LL;
string inputfile="";

#include <omp.h>


class CppClock{
    public:
     std::chrono::system_clock::time_point start;
     CppClock(){
         start = std::chrono::system_clock::now();
     }
     double get(){
         auto end = std::chrono::system_clock::now();
         double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());
         return elapsed;
     }
 };


template <class T> void chmin(T & a, T const & b) { if (b < a) a = b; }
template <class T> void chmax(T & a, T const & b) { if (b > a) a = b; }
template <typename T>
string to_string(const T &n){ostringstream stm;stm << n;return stm.str();}
inline int toInt(string s) { int v; istringstream sin(s); sin >> v; return v; }
using ll =long long;
const ll mod = 1000'000'007;

// n次元配列の初期化。第２引数の型のサイズごとに初期化していく。
template<typename A, size_t N, typename T>
void Fill(A (&array)[N], const T &val){
    std::fill( (T*)array, (T*)(array+N), val );
}

template< class _Type>
//gcc can't inline this funciton for some reson 
inline void _swap(_Type &a, _Type&b)  {
  _Type temp = b;
  b = a;
  a = temp;
}




const int kisuu=256;
const int MaxThreadNum=224;
const long long MaxDataSize=10000000000;
const long long MaxDataNum=4294967295;
const int MaxKisuu=256;
int threadNum;
int NumRange=0;
ll unsortCnt=0;

//std::vector<ll> Dataset;
//concurrent_vector<int> v;

ll Datasize;



//Sort 対象データセット
void makeDataset(int* Dataset){    

  int fd =open(inputfile.c_str(),O_RDONLY);
    if(fd==-1){printf("cannot open file.\n");exit(1);}
    else{printf("successfully open file.\n");}
   
    for(long long i=0;i<400LL;i++){
      int ret=read(fd,&Dataset[i*SIZE],sizeof(int)*SIZE);
      if(ret==-1){
	cout<<"read failure."<<endl;
	exit(1);
      }else{
	//	cout<<"ret:"<<ret<<endl;
      }
    }
    cout<<"read done"<<endl;

    for (long long i = 0; i < DATASIZE; i++) {
     if (Dataset[i] == 0) {
       printf("Error: %lld\n", i);
       exit(1);
     }
    }
  
    close(fd);
    cout << "close done" << endl;
}


static const ll kRadixBits = 8;
static const size_t kInsertSortThreshold = 0;
static const ll kRadixMask = (1LL << kRadixBits) - 1LL;
static const ll kRadixBin = 1LL << kRadixBits;



template<class D>
inline int determineDigitBucket(int stage,D num){
  return ((num>>(8*stage))&kRadixMask);
}


void report_num_threads(int level)
{
    #pragma omp single
    {
        printf("Level %d: number of threads in the team - %d\n",
               level, omp_get_num_threads());
    }
}

template<class T>
bool compare(const T &x,const T &y){
  return x < y;
}
template <class RandomIt>
inline void insert_sort_core_(RandomIt s, RandomIt e)
{
    for (RandomIt i = s + 1; i < e; ++i) {
        if (compare(*i, *(i - 1))) {
            RandomIt j;
            auto tmp = *i;
            *i = *(i - 1);
            for (j = i - 1; j > s && compare(tmp, *(j - 1)); --j) {
                *j = *(j - 1);
            }
            *j = tmp;
        }
    }
}



double createHistTime=0.0;
double paradisParmutateTime=0.0;
double paradisRepairTime=0.0;
double createGhGt=0.0;
double quickSortTime = 0.0;
double MiTime = 0.0;
int localhistCall=0;
int recum=0;
int CallProcesses=0;
int swapNum=0;
int needRepairNum=0;
int prefer_insert=0;

double repair_ck=0.0;
ofstream writing_file;


template<class D,int kth_byte>
inline void RadixSort(int* arr,ll elenum,ll start,int processes=1){
    ll cnt[MaxKisuu];
     
    rep(i,kisuu){cnt[i]=0LL;}    

    //step1
    ll part=elenum/(ll)processes;
    ll res=elenum%(ll)(processes);

    //threadNum個のヒストグラムを作る
    ll localHists[MaxThreadNum][MaxKisuu];
    ll gh[MaxKisuu],gt[MaxKisuu],starts[MaxKisuu],ends[MaxKisuu];
    ll ph[MaxThreadNum][MaxKisuu];
    ll pt[MaxThreadNum][MaxKisuu];
    
    //indicesをpartitionしてカウント
    ll nvalue = elenum;

    // sum i(Ci) > 0 初めは全要素
    int SumCi=elenum,nth=2;
    int roop=0;
    //for paradis repair
    //    ll pfp[processes+1];
    int var_p=processes;
    double start_ck,end_ck;
#pragma omp parallel num_threads(processes)   
      {

          int th = omp_get_thread_num();
	  #pragma omp for
          rep(i,kisuu){
	    rep(th,processes)localHists[th][i]=0LL;
	  }
	  #pragma omp barrier
	  #pragma omp for
	  for(ll i=start;i<start+elenum;i++){
	    //	    assert(arr[i]>=(D)0);
	    int digit=determineDigitBucket(kth_byte,arr[i]);
	    localHists[th][digit]++;
	  }
	  #pragma omp barrier
          #pragma omp for      
          for(int i=0;i<kisuu;i++){
            for(int j=0;j<processes;j++){
             cnt[i]+=localHists[j][i];
            }
          }
    #pragma omp barrier

    #pragma omp single
    {
    gh[0]=start;
    gt[0]=start+cnt[0];
    starts[0]=gh[0];
    }
    //step2
    #pragma omp single
    for(int i=1;i<kRadixBin;i++){
        //calc ghi
        gh[i]=gh[i-1]+cnt[i-1];
        //calc gti;
        gt[i]=gh[i]+cnt[i];

        starts[i]=gh[i];
    }

    #pragma omp barrier
    //step3
    while(SumCi!=0){
      

#pragma omp for
      for(int ii=0;ii<var_p;ii++)
        {
	    int pID=omp_get_thread_num();

	
	    for(int i=0;i<kisuu;i++){
	      ll part=(ll)(gt[i]-gh[i])/(ll)var_p;
	      ll res=(ll)(gt[i]-gh[i])%(ll)(var_p);
	      
	      if(pID<var_p-1){
		ph[pID][i]=part*pID+gh[i];
		pt[pID][i]=part*(pID+1LL)+gh[i];
	      }else{
		ph[pID][i]=part*pID+gh[i];
		pt[pID][i]=part*(pID+1LL)+gh[i]+res;
	      }
	    }

	    for(int i=0;i<kisuu;i++){
	      ll head=ph[pID][i];
	      while(head<pt[pID][i]){
		int v=arr[head];
		int k=determineDigitBucket(kth_byte,v);
		
		while(k!=i&&ph[pID][k]<pt[pID][k]){
		  //		  int arph=(int)arr[ph[pID][k]++];
		  _swap(v,arr[ph[pID][k]++]);
		  k=determineDigitBucket(kth_byte,v);
		  //swapNum++;
		}
		if(k==i){
		  arr[head++]=arr[ph[pID][i]];
		  arr[ph[pID][i]++]=v;
		}else{
		  arr[head++]=v;
		}
	      }
	    }
	      
	    
	}//end of omp pfp

#pragma omp single
      {
      nth=1;
      SumCi=0;      
      }
#pragma omp barrier
#pragma omp single
      start_ck=omp_get_wtime();
 #pragma omp for
 for(int k=0;k<kRadixBin;k++){
   ll nthAftKey=0LL;
   for(int pID=0;pID<var_p;pID++){
     if(pt[pID][k]-ph[pID][k]>=1e7 && ph[pID][k]-nthAftKey>=pt[pID][k]-ph[pID][k]){
       //             cout<<"more than 1e7!"<<endl;
    #pragma omp parallel for num_threads(thread::hardware_concurrency()/processes)
    for(ll ii=ph[pID][k];ii<pt[pID][k];ii++){
    // if(gh[k]==-1)gh[k]=ph[pID][k];
        _swap(arr[gh[k]+(ii-ph[pID][k])+nthAftKey],arr[ii]);
    }
    }else{
    for(ll ii=ph[pID][k];ii<pt[pID][k];ii++){
    // if(gh[k]==-1)gh[k]=ph[pID][k];
    _swap(arr[gh[k]+(ii-ph[pID][k])+nthAftKey],arr[ii]);
    }
   }
     nthAftKey+=(pt[pID][k]-ph[pID][k]);
   }
   gt[k]=gh[k]+nthAftKey;
   if(gt[k]-gh[k]>0LL)SumCi=1;
 }

 #pragma omp barrier
#pragma omp single
 {
   end_ck=omp_get_wtime();
   repair_ck=(end_ck-start_ck)*1000;
      cout<<"clock \t"<<repair_ck<<"[ms]"<<endl;
      //repairは複数回行われるため
      writing_file<<"repair time "<<repair_ck<<"[ms]"<<endl;

 }
      
    }//end of while
    }//end of omp2

     
    if(kth_byte>0){
      {
#pragma omp parallel  num_threads(processes) 
#pragma omp single
	{

	for(int i=0;i<kisuu;i++){
	  int nextStageThreads=1;
	  nextStageThreads=processes*(cnt[i]*(log(cnt[i])/log(kRadixBin))/(elenum*(log(elenum)/log(kRadixBin))));

     	  if(cnt[i]>64LL){
	  #pragma omp task	   
	    RadixSort<D,(kth_byte > 0 ? (kth_byte - 1) : 0)>(arr,cnt[i],starts[i],max((int)nextStageThreads,(int)1));
	  }
	  else if(cnt[i]>1LL){
	    //if elements less than 64 call insertion sort
	    insert_sort_core_(arr+starts[i],arr+starts[i]+cnt[i]);
	  }
	}
	#pragma omp taskwait
	}
      }
    }
}

signed main(int argc, char** argv){

    if (argc < 2) {
        printf("ERROR! set thread num\n");
        return 1;
    }
   
    threadNum = atoi(argv[1]);

    if(threadNum>MaxThreadNum){
        printf("ERROR! too many thread num\n");
        return 1;
    }
    if (argc < 4) {  
        printf("ERROR! too few argments. (note) ./paradis [threadNum] [file] [datasize] \n (etc)./paradis_ompf_my_read 64 ./1_100_10^9 1000000000\n");
        return 1;
    }
        
    inputfile = to_string(argv[2]);

    DATASIZE = atoll(argv[3]);
    
    cout<<"radix="<<kisuu<<endl;
    cout<<"maximum threads="<<thread::hardware_concurrency()<<endl;

   
    cout<<"creating dataset..."<<flush;

    int *Dataset = (int*)calloc(sizeof(int),DATASIZE);
    makeDataset(Dataset);

    bool flag_one=false;
    for (long long i = 0; i < DATASIZE; i++) {
      //      cout<<"Dataset["<<i<<"]="<<Dataset[i]<<endl;
    if (Dataset[i] == 0) {
      printf("Error: %lld\n", i);
      exit(1);
    }
    if(Dataset[i]==1)flag_one=true;
    }

    //in general case this is not need
    /*
    if(!flag_one){
      cout<<"Error no one"<<endl;
      exit(1);
    }
    */

    
    cout<<" finish!"<<endl;cout<<endl;


    string filename = "log_ap_3.txt";

    
    writing_file.open(filename,ios::app);

     ll Ds = -1;

    NumRange=0;while(Ds){NumRange++;Ds/=kisuu;}
   cout<<"NumRange="<<NumRange<<" "<<endl;
   writing_file<<argv[0]<<endl;
    writing_file<<"threadNum="<<threadNum<<" Datasize="<<DATASIZE<<" NumRange="<<NumRange<<endl;
   
    
       
    cout<<"PARADIS is running..."<<flush;
    auto start = std::chrono::system_clock::now();
    //sortしたい目的の配列,levelの数,次のlevelに渡すindexの配列,levelの深さ
    omp_set_max_active_levels(4);
    RadixSort<ll,3>(Dataset,DATASIZE,0,threadNum);
    auto end = std::chrono::system_clock::now();


    for(ll i=1;i<DATASIZE;i++){
      if(Dataset[i-1]>Dataset[i]){
	cout<<"not sorted at "<<i<<endl;
	cout<<"Dataset["<<i-1<<"]="<<Dataset[i-1]<<" Dataset["<<i<<"]="<<Dataset[i]<<endl;
	exit(1);
      }
    }
    if(!std::is_sorted(Dataset,Dataset+DATASIZE)){
        std::cerr<<"Not sorted"<<std::endl;
    }else{
      //std::cout<<"sorted!! good job!"<<std::endl;
      cout<<" finish!"<<endl;
    }
    

    auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000.0);

    printf("paradis time %lf[ms]\n",elapsed);
    writing_file<<"paradis time "<<elapsed<<"\n"<<endl;
}

