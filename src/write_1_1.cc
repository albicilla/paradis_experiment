#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#define SIZE (1000LL * 100LL * 25LL)
#define FILE "/data/1_1_10^9"

using namespace std;

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
  for(long long i=0;i<n*100LL;i++){
      data[idx]=1;
      idx++;
  }

  cout << "1:" <<idx<< endl;
  
  for(long long i=0;i<n*100LL;i++){
      data[idx]=2;
      idx++;
  }

  cout << "2:" << idx<<endl;
  for(long long i=0;i<n*100LL;i++){
    data[idx]=1;
    idx++;
  }

  for(long long i=0;i<n*100LL;i++){
    data[idx]=2;
    idx++;
  }
  printf("idx: %lld\n", idx);
  
  cout << "3" << endl;
  for (long long i = 0; i < SIZE * 400LL; i++) {
    if (data[i] == 0) {
      printf("Error: %lld (%lld)\n", i, SIZE*400LL);
      printf("idx: %lld\n", idx);
      exit(1);
    }
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



