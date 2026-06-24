#pragma once
#include"root.h"
#include<random>
#include<iostream>
#include<fstream>
#include<chrono>

pp::vector make_r(double rmax, double dr);
pp::matrix make_H(double rmax, double dr);
void hamiltonian();