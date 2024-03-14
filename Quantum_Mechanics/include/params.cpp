#pragma once;
#include <string>

constexpr double pi = 3.141592653589793;
constexpr long double hbar = 1.054571817e-34;
const std::string RAW_PATH = "../data/raw";

constexpr double L = 1.00; // size of well
constexpr unsigned int SPACESTEPS = 2000; // number of space discretization steps

constexpr unsigned int TIMESTEPS  = 2000;// number of time discretization steps 
constexpr double END_TIME = 1; // maximum time  

constexpr double dx = static_cast<double>(L) / (SPACESTEPS+1);
constexpr double dt = static_cast<double>(END_TIME) / (TIMESTEPS+1);