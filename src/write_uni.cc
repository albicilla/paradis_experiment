#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <random>
#include <climits>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#define SIZE (1000LL * 100LL * 25LL * 1LL)
#define FILE "/work/albicilla/uni_10^9"

using namespace std;

using ll = long long;

signed main()
{
  cout << "-1" << endl;
  int *data = (int*)malloc(sizeof(int)*SIZE*400LL);
  if(!data) {
    cout << "malloc fail" << endl;
    exit(1);
  }
  struct timeval begin,end;
  gettimeofday(&begin,NULL);
  cout << "0" << endl;
  
  long long idx=0LL,n=SIZE;
  uniform_int_distribution<int> dist(0.0,INT_MAX);
  random_device seed_gen;
  default_random_engine engine(seed_gen());
  
  for(ll i=0;i<n;i++){
      for(int j=0;j<100;j++){
          data[idx]=dist(engine);
          idx++;
      }
  }
  cout << "1:" <<idx<< endl;
  for(ll i=0;i<n;i++){
      for(int j=0;j<100;j++){
          data[idx]=dist(engine);
          idx++;
      }
  }
  cout << "2:" << idx<<endl;
  for(ll i=0;i<n;i++){
      for(int j=0;j<100;j++){
          data[idx]=dist(engine);
          idx++;
      }
  }
  
  cout << "3" << endl;
  for(ll i=0;i<n;i++){
      for(int j=0;j<100;j++){
          data[idx]=dist(engine);
          idx++;
      }
  }
  cout<<"idx="<<idx<<endl;

  for(int i=0;i<100;i++)
    {
      cout<<"data[i]="<<data[i]<<endl;
    }
  int fd = open(FILE, O_WRONLY|O_TRUNC|O_CREAT|O_APPEND,0644);
  for (long long i = 0; i < 400LL; i++) {
    int ret = write(fd, &data[i*SIZE], sizeof(int)*SIZE);
    if (ret == -1) {
      cout << "write failure" << endl;
      exit(1);
    } //else printf("ret: %d\n", ret);
  }
  close(fd);
  gettimeofday(&end,NULL);
  //printf("%ld usec\n",(end.tv_sec-begin.tv_sec)*1000*1000+(end.tv_usec-begin.tv_usec));
  printf("%ld sec\n",(end.tv_sec-begin.tv_sec));
}



