#include<fstream>
#include<string>

#include"molecular_dynamics.h"

#define str(text) std::string text
#define print(text) std::cout << text <<std::endl

void calculate_orbits(double **initial, int N, double t_max, str(name), double t_min=0.0, int steps=1000);

int main(int argc, char *argv[]){
    int N, c; double **initial_conditions = NULL;
    double TMAX = 1e3; //days

    load_file("earth_sun.csv", initial_conditions, N, c);
    calculate_orbits(initial_conditions, N, TMAX, "earth_sun");

    load_file("earth_moon_sun.csv", initial_conditions, N, c);
    calculate_orbits(initial_conditions, N, TMAX, "earth_moon_sun");

    load_file("kozai.csv", initial_conditions, N, c);
    calculate_orbits(initial_conditions, N, TMAX, "kozai");

    for(int i=0; i<N; i++) delete[] initial_conditions[i];
    delete[] initial_conditions;

    return 0;
}

void calculate_orbits(double **initial, int N, double t_max, str(name), double t_min, int steps){
    Body *Bodies = new Body[N];
    Collider Newton(N);

    for(int i=0; i<N; i++) Bodies[i] = Body(N);

    std::ofstream orbit("orbits_"+name+".txt"), 
        energy("energy_"+name+".txt"), 
        angular("angular"+name+".txt");

    orbit << "#x,y,z (for each  body)\n";
    energy << "#t,E/E0\n";
    angular << "#t,L/L0\n";

    for(int i=0; i<N; i++)
        Bodies[i].initialize(
            initial[i][0], initial[i][1], initial[i][2],
            initial[i][3], initial[i][4], initial[i][5], 
            initial[i][6], initial[i][7]
        );

    double dt = (t_max - t_min)/(double) steps;
    double E0 = Newton.energy(Bodies), L0 = Newton.angular_momentum(Bodies[1]);

    for(int t=0; t<steps; t++){

        for(int i=0; i<N; i++)
            orbit << Bodies[i].get_x() << ',' << Bodies[i].get_y() << ',' << Bodies[i].get_z() << ',';
        orbit << '\n';

        energy << t << ',' << Newton.energy(Bodies)/E0 << '\n';
        angular << t << ',' << Newton.angular_momentum(Bodies[1])/L0 << '\n';

        Newton.move_with_pefrl(Bodies, dt);
    }

    orbit.close(); energy.close(); angular.close();

    print("All done for " + name);

}