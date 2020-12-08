#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

#include"molecular_dynamics.h"


#define str_(text) std::string text
#define print(text) std::cout << text <<std::endl

void calculate_orbits(double **initial, int N, double t_max, str_(name), double t_min=0.0, int steps=1000);

int selected = 1; // Angular momentum calculated by default is earth's (0 is sun, 2 is moon)

int main(int argc, char *argv[]){
    int N, c; double **initial_conditions = NULL;
    double TMAX = 1e3; //days

    fh::load_file("earth_sun.csv", initial_conditions, N, c);
    calculate_orbits(initial_conditions, N, TMAX, "earth_sun");
    fh::clear(initial_conditions, N);

    fh::load_file("earth_moon_sun.csv", initial_conditions, N, c);
    calculate_orbits(initial_conditions, N, TMAX, "earth_moon_sun");
    fh::clear(initial_conditions, N);

    fh::load_file("kozai.csv", initial_conditions, N, c); // This starts with the moon at 90 degrees
    double moon_earth = initial_conditions[2][1];

    for(int i=90; i>=0; i--){
        double theta = i*M_PI/180.0; // Angle from x-axis
        initial_conditions[2][0] = 1.0 + moon_earth*std::cos(theta);
        initial_conditions[2][1] = moon_earth*std::sin(theta);

        std::stringstream i_s; i_s << i;
        str_(name) = "kozai_" + i_s.str();
        calculate_orbits(initial_conditions, N, TMAX, name);
    }

    fh::clear(initial_conditions,  N, false);

    return 0;
}

void calculate_orbits(double **initial, int N, double t_max, str_(name), double t_min, int steps){
    Body *Bodies = new Body[N];
    Collider Newton(N);

    for(int i=0; i<N; i++) Bodies[i] = Body(N);

    std::ofstream orbit("orbits_"+name+".txt"), 
        energy("energy_"+name+".txt"), 
        angular("angular_"+name+".txt");

    orbit << "#x_earth,y_earth,z_earth,x_sun,y_sun,z_sun(,x_moon,y_moon,z_moon) (if the moon is there)\n";
    energy << "#t,E/E0\n";
    angular << "#t,L/L0\n";

    for(int i=0; i<N; i++)
        Bodies[i].initialize(
            initial[i][0], initial[i][1], initial[i][2],
            initial[i][3], initial[i][4], initial[i][5], 
            initial[i][6], initial[i][7]
        );

    double dt = (t_max - t_min)/(double) steps;
    double E0 = Newton.energy(Bodies), L0 = Newton.angular_momentum(Bodies[selected]);

    for(int t=0; t<steps; t++){

        for(int i=0; i<N; i++)
            orbit << Bodies[i].get_x() << ',' << Bodies[i].get_y() << ',' << Bodies[i].get_z() << ',';
        orbit << '\n';

        energy << t << ',' << Newton.energy(Bodies)/E0 << '\n';
        angular << t << ',' << Newton.angular_momentum(Bodies[selected])/L0 << '\n';

        Newton.move_with_pefrl(Bodies, dt);
    }

    orbit.close(); energy.close(); angular.close();

    delete[] Bodies;

    print("All done for " + name);

}