/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include<fstream>

#include"molecular_dynamics.h"



int main(int argc, char *argv[]){
    // Initial conditions
    double m0[N] = {333030.0, 1.0}; // earth mass
    double r0[N] = {109.0, 1.0}; // earth radius
    double x0[N] = {0.0, 1.0}, y0[N] = {0.0, 0.0}; // AU
    double vx0[N] = {0.0, 0.0}, vy0[N] = {0.0, 2*M_PI/365.25}; // AU/day

    double t_min = 0, t_max = 1e3; // days
    int steps = 1000;

    // Implementation
    Body BodySystem[N];
    Collider Newton;

    std::ofstream orbit("earth.txt"), energy("energy.txt"), angular("angular.txt");

    for(int i=0; i<N; i++)
        BodySystem[i].initialize(x0[i], y0[i], vx0[i], vy0[i], m0[i], r0[i]);

    double dt = (t_max - t_min)/(double) steps;
    double E0 = Newton.energy(BodySystem);

    for(int t=0; t<steps; t++){
        orbit << BodySystem[1].get_x() << '\t' << BodySystem[1].get_y() << '\n';
        energy << t << '\t' << Newton.energy(BodySystem)/E0 << '\n';
        Newton.move_with_pefrl(BodySystem, dt);
    }

    orbit.close(); energy.close(); angular.close();

    return 0;
}