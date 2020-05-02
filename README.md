# paradis
implemention of PARADIS - fast parallel radix sort algorithm.

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
## run
### data generation
write.cc generates binary style input file. The default filename is "1_100_10^9".
```sh
./write
```

### run paradis
usage
./paradis_ompf_my_read [thread] [inputfile] [datasize].

```sh
./paradis_ompf_my_read 64 ./1_100_10^9 1000000000
```


