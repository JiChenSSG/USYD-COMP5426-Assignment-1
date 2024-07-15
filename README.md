# USYD 2024 COMP5426 Assignment1

a parallel design for Gaussian Elimination in c with openmpi

## Usage

```sh
make
chmod +x test.sh
chmod +x test\_omp.sh

# For running and check the serial result, use
# ./gepp N
./gepp 1000

# For running and check the parallel result, use
# ./gepp_omp N T
./gepp_omp 1000 2

# For serial result time testing, use
# sh ./test.sh N
./test.sh 1000

# For serial result time testing, use:
# ./test_omp.sh N T
./test_omp.sh 1000 2
```
