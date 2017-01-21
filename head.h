#include <iostream>
#include <fstream>
#include <random>
#include<iomanip>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
// #include <malloc.h>
#include <stdlib.h>
#include <ctime>
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <bitset>
#include <complex>
#include <map>

#define PI 3.14159265358979

//using namespace std;

int main();
void PrintHam(double JJ);

// precondition.hpp
int SpinBasisConstruct(std::vector<int> &vec);
int HoleBasisConstruct(std::vector<int> &vec);

// hamiltonian.hpp
int Flip(int s, int i, int j);
