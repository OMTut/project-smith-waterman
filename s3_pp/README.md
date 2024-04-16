compile
```
g++ -fopenmp -o app s3_pp.cpp
```

## Parallel Version
Currently, only the inner loop takes advantage of parallelization. I haven't found a way to safely parallelize the outer loop. My attempts have led to a race condition that skews the end result quite significantly.

## GPU Version
compile
``` nvcc -o app filname.cu