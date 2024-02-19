#pragma once;
const double pi = 3.141592653589793;



const double L = 1.00; // size of well
const unsigned long int SPACESTEPS = 1000; // number of space discretization steps
const unsigned long int TIMESTEPS  = 1000;// number of time discretization steps 
const double END_TIME = 5; // maximum time  

const double dx = static_cast<double>(L) / (SPACESTEPS+1);
const double dt = static_cast<double>(END_TIME) / (TIMESTEPS+1);