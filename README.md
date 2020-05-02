# PARADIS
implemention of PARADIS - fast parallel radix sort algorithm.
http://www.vldb.org/pvldb/vol8/p1518-cho.pdf
## Requirements

* CMake `>= 3.50`
* C++ Compiler `>= C++17`

## How to build

```sh
mkdir build
cd build
cmake ..
make
```
## Run
### Data generation
write.cc generates binary style input file. The default filename is "1_100_10^9".
```sh
./write
```

### Run PARADIS
usage
./paradis_ompf_my_read [thread] [inputfile] [datasize].

```sh
./paradis_ompf_my_read 64 ./1_100_10^9 1000000000
```


