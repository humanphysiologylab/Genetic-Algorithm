# Genetic-Algorithm
Genetic Algorithm implementation to determine the set of cardiomyocyte model parameters.

## Build and compile

To compile the code type the following commands in the terminal:
```C
cmake CMakeLists.txt
make
```
IMPORTANT NOTE: this implementation of genetic algorithm requires MPI!

## Run
The main configuration parameters (number of generations, number of organisms etc) are in the `input_ga.txt` file.
```C
mpirun -np [num_process] ./ga input_ga.txt  
```  
The folder also contains an example of MPI job script file (`snode.sh`) that can be used for high-perfomance calculations.

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
