#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <chrono>

class Particles {


public:

    typedef struct {
        double fx;
        double fy;
    } f_var;

    typedef struct {
        double x_min;
        double x;
        double y;
    } min_var;

    typedef struct {
        double ek;
        double ep;
        double em;
        double ef;
    } e_var;

    Particles(){
        create_particles();
    }

    int create_particles();

    std::vector<std::vector<double>> count_acc(std::vector<std::vector<double>> r);

    min_var min_r(std::vector<double> p1, std::vector<double> p2);

    f_var force(std::vector<double> p1, std::vector<double> s1, std::vector<double> p2, std::vector<double> s2);

    e_var energy();

    int moving();

    //methods of spin

    int create_spins();
    std::vector<double> normalized(std::vector<double> massive);
    int change_spins(std::vector<std::vector<double>> r);

private:
    int n = 25;
    double e = 0.03;
    double dt = 0.01;
    int s = 80;
    int l = 250;
    std::vector<int> tscale;
    std::vector<double> ekscale;
    std::vector<double> epscale;

    std::vector<std::vector<double>> r;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> a;

    //components of spin

    double h = 0.1;
    double spin_l = 2.0;
    double rc = (double)l/4;
    double j_0 = 50;
    std::vector<std::vector<double>> spins;

};