# Landau-Siegel zero tester

This repository is associated with the paper "Numerical computations concerning Landau-Siegel zeros" by Rick F. Lu, [Asif Zaman](www.math.toronto.edu/zaman/), and Haonan Zhao. This project was part of a 2025 undergraduate research project at the University of Toronto Math Department. We gratefully acknowledge funding by NSERC and UTEA Undergraduate Student Research Awards. This research was also supported by Compute Canada ([alliancecan.ca](https://alliancecan.ca/en)), where most of our computation was performed.

Given an absolute constant $c > 0$ and primitive quadratic Dirichlet character $\chi$ modulo $q$, we provide an algorithm which shows that 

$$L(\sigma,\chi) \neq 0 \quad \text{ for real } \sigma > 1-\frac{c}{\log(q)}.$$  

In other words, we computationally eliminate the existence of exceptional zeros for a fixed modulus. Our code also handles all moduli $q \leq Q$ with some additional efficiencies. 

## Table of Contents
1. [Background](#background)
2. [Requirements](#requirements)
3. [Usage](#usage)
4. [Directory Structure](#directory-structure)

## Background

Let $\chi$ be a Dirichlet character mod $q$,

$$L(s,\chi) := \sum_{n=1}^\infty \frac{\chi(n)}{n^s}, \quad \Re(s) > 1,$$

be the associated Dirichlet $L$-function, and $s = \sigma + it$ be a complex number.

Landau proved that there is a constant $c > 0$ such that for all $q$,
$\prod_\chi L(s,\chi)$ has at most one zero in the region

$$
\sigma > 1 - \frac{c}{\log(q(|t|+4))}
$$

where $\chi$ ranges over all Dirichlet characters of modulus $q$.  
If such a zero exists, then it is necessarily real and the associated character $\chi$ is quadratic.

Such a zero, if it exists, is called a **Landau–Siegel zero** or an **exceptional zero**.  

Running this program on Compute Canada's Nibi clusters, we computationally verified that Landau–Siegel zeros do not 
exist for any primitive quadratic character of modulus $q \le 10^{10}$. In particular, we proved the 
following theorem.

> **Theorem 1.** If $\chi$ is a primitive quadratic Dirichlet character of modulus $q \leq 10^{10}$, then
>
> $$L(\sigma,\chi) \neq 0 \quad \text{ for real } \sigma > 1-\frac{1}{5\log(q)}.$$  


## Requirements

To build and run the computational verification code, you will need the following system dependencies:

- [**FLINT**](http://flintlib.org/) (Fast Library for Number Theory) ≥ 2.9
- **CMake** ≥ 3.20 (for building the project)
- **MPI** (Message Passing Interface, e.g. OpenMPI or MPICH) for parallel execution

Make sure these are installed and available in your system path.

## Usage

After installing the requirements and building the project, you can run the computations as follows.

### 1. Precomputation

Run `misc/precompute_kronecker.c` to generate `chi.txt`. After the file is generated, move it to `/input` dirctory.

We have provided a list of the first 5761455 primes in `input/primes.txt`. If you wish to use more primes than this, 
i.e. if you wish to set $N_0>5761455$, then you must generate your own list of primes up to $N_0$, call the file `primes.txt` 
and move it to the `/input` dirctory. The file `primes.txt` has to contain all primes upto $N_0$, one prime a line, 
in increasing order. One may generate such list using [Primesieve](https://github.com/kimwalisch/primesieve).

### 2. Choose the version

Edit `src/CMakeLists.txt` to choose the version of the main program to run (either main.c, main_by_file.c, or main_test_ind.c).

### 3. Edit parameters

Set the parameters of the program to your use case by changing the values defined by `#define`. In particular, modify `N00`, 
`checkDistance`, `qMax`, `c0`, and `preset`.

### 4. Build the project
From the repository root:  
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

### 5. Run with MPI
From the repository root use the following mpi command to run (change 10 to desired number of MPI processes).
```bash
mpirun -np 10 ./build/src/mpiTest
```

### 6. Using `main_by_file.c`
If one wants to run the verification on a list of $q$ in a file, name the file to `input.txt` and store it in the `/input` 
dirctory. The file `input.txt` has to contain one $q$ in each line and the `lineMax` constant in `main_by_file.c` has to 
match the line number in `input.txt`.

## Directory Structure

All source code is stored in the directory `/src`. All input files are stored in the `/input` directory. 
The `/test` directory contains all unit tests for this project. `/misc` contains helpful scripts that are not used 
during the main process.


By Rick Lu, Asif Zaman, and Haonan Zhao.
