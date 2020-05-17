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
//    M=80002 #число итераций
//    m=20000 #через сколько итераций вывод
//    mcr=50 #через сколько итераций учёт температуры и расчёт корреляций
//    N=16 # число частиц
//    Lx=4.5 #размер области
//    Ly=np.sqrt(3)*Lx/2
//    dt=0.01 #шаг по времени
//            E=0.05 #эпсилон
//            J0=0.2
//    rc=2
//    s=1 #сигма
//            T=0.01
//    v=0.01 #макс скорость
//    DT=1600 #время измерения корреляций
//            t2=0 #t0 для автокорр. ф-ции скорости
    int n = 16;
    double e = 0.05;
    double dt = 0.001;
    double s = 1;
    double lx = 4.5;
    double ly = pow(3,0.5)/2 * lx;
    double vmax=0.01;
    std::vector<int> tscale;
    std::vector<double> ekscale;
    std::vector<double> epscale;

    std::vector<std::vector<double>> r;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> a;

    //components of spin

    double h = 1;
    double spin_l =1;
    double rc = 1;
    double j_0 = 0.0;
    std::vector<std::vector<double>> spins;

};