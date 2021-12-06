# Genetic Algorithm
Genetic Algorithm implementation to estimate cardiomyocyte mathematical model parameters.

## Build

To compile the code type the following commands in the terminal:
```C
git clone https://github.com/humanphysiologylab/Genetic-Algorithm.git
cd Genetic-Algorithm
mkdir build-release && cd $_
module load mpi/openmpi-x86_64
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

On MacOS use
```C
git clone https://github.com/humanphysiologylab/Genetic-Algorithm.git
cd Genetic-Algorithm
mkdir build-release && cd $_
cmake -D CMAKE_C_COMPILER=/opt/homebrew/bin/gcc-11 -D CMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-11 ..
make
```

To generate the documentation, please run
```C
doxygen Doxyfile
```
The documentation will be available for your browser at `doxygen/html/index.html` .

Important note: this implementation of genetic algorithm requires MPI.

## Run
The main configuration parameters (number of generations, number of organisms etc) should be specified in `config.json` file.
```C
mpirun -np [num_mpi_processes] ./ga config.json
```  
The folder also contains an example of MPI job script file (`task.sh`) to run the program on high performance computing clusters.

## Citing Genetic Algorithm
<a href="https://link.springer.com/article/10.1007/s10863-018-9775-7#citeas" alt=""><img src="https://img.shields.io/badge/DOI%3A-doi.org%2F10.1007%2Fs10863--018--9775--7-brightgreen.svg"></a>

    Smirnov DN, Belyakova GA, Syunyaev RA, Gusev OA, Deviatiiarov RM, Efimov I “Personalized ionic models of cardiomyocytes on the basis of molecular data and genetic algorithms” in Biomembranes 2018, Dolgoprudny, Russia, Oct 2018.

```
@misc{genetic_algorithm_2018,
  author       = {Smirnov DN, Belyakova GA, Syunyaev RA, Gusev OA, Deviatiiarov RM, Efimov I},
  title        = {Personalized ionic models of cardiomyocytes on the basis of molecular data and genetic algorithms},
  journal      = {J Bioenerg Biomembr}
  year         = 2018,
  doi          = {doi.org/10.1007/s10863-018-9775-7},
  url          = {https://doi.org/10.1007/s10863-018-9775-7}
}
```
