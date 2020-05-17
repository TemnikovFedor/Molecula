#include "particles.h"

int Particles::create_particles(){
    int k = n;
    while(pow(k,0.5)-(int)pow(k,0.5)!=0){
        k += 1;
    }
    k = pow(k,0.5);

    for(int i=0; i<(int)k; i++){
        for(int j=0; j<(int)k; j++){
            std::vector<double> rr;
            rr.push_back((0.5+i)*l/k+0.25*(rand() % 10000001 / 10000000.0)*l/k);
            rr.push_back((0.5+j)*l/k+0.25*(rand() % 10000001 / 10000000.0)*l/k);
            r.push_back(rr);
            if (r.size()==n) break;
        }
        if (r.size()==n) break;
    }

    std::vector<double> v_sum(2, 0);
    for(int i=0; i<n; i++){
        std::vector<double> vv;
        vv.push_back(((rand() % 20000000-9999999) / 10000000.0)*l/16);
        vv.push_back(((rand() % 20000000-9999999) / 10000000.0)*l/16);
        v.push_back(vv);
        v_sum[0] += v[i][0];
        v_sum[1] += v[i][1];
    }

    for(int i=0; i<n; i++){
        v[i][0] -= v_sum[0]/n;
        v[i][1] -= v_sum[1]/n;
    }

    create_spins();
    count_acc(r);
    return 0;
}

std::vector<std::vector<double>> Particles::count_acc(std::vector<std::vector<double>> r){
    a.clear();
    for (int i = 0; i < n; i++) {
        double ax=0;
        double ay=0;
        for (int j = 0; j < n; j++) {
            if (i != j){
                Particles::f_var Force = force(r[i],spins[i],r[j],spins[j]);
                ax += Force.fx;
                ay += Force.fy;
            }
        }
        std::vector<double> aa;
        aa.push_back(ax);
        aa.push_back(ay);
        a.push_back(aa);
    }
    return a;
}

Particles::min_var Particles::min_r(std::vector<double> p1, std::vector<double> p2){
    Particles::min_var min_distance;
    min_distance.x_min = pow(l,2)*10;
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            double distance = pow(p2[0]+i*l-p1[0],2)+pow(p2[1]+j*l-p1[1],2);
            if (distance < min_distance.x_min){
                min_distance.x_min = distance;
                min_distance.x = p2[0]+i*l;
                min_distance.y = p2[1]+j*l;
            }
        }
    }
    min_distance.x_min = pow(min_distance.x_min,0.5);
    return min_distance;
}

Particles::f_var Particles::force(std::vector<double> p1, std::vector<double> s1, std::vector<double> p2, std::vector<double> s2){
    Particles::min_var min_distance = min_r(p1,p2);
    Particles::f_var F;

    F.fx = 4*e*(12*pow(s,12)*(p1[0]-min_distance.x)/pow(min_distance.x_min,14)-6*pow(s,6)*(p1[0]-min_distance.x)/pow(min_distance.x_min,8));
    F.fy = 4*e*(12*pow(s,12)*(p1[1]-min_distance.y)/pow(min_distance.x_min,14)-6*pow(s,6)*(p1[1]-min_distance.y)/pow(min_distance.x_min,8));
    if (min_distance.x_min < rc) {
        double scalar = 0;
        for (int i = 0; i < 3; i++) {
            scalar += s1[i]*s2[i];
        }
        F.fx += j_0*3*pow(1-min_distance.x_min/rc,2)*(-1/rc)*(1/(2*min_distance.x_min))*(p1[0]-min_distance.x);
        F.fy += j_0*3*pow(1-min_distance.x_min/rc,2)*(-1/rc)*(1/(2*min_distance.x_min))*(p1[1]-min_distance.y);
    }
    return F;
}

Particles::e_var Particles::energy(){
    Particles::e_var E;

    E.ek=0;
    for (int i = 0; i < n; i++) {
        E.ek += 0.5*(pow(v[i][0],2)+pow(v[i][1],2));
    }

    E.ep=0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                Particles::min_var min_distance = min_r(r[i],r[j]);
                E.ep += 4*e*(pow((s/min_distance.x_min),12)-pow((s/min_distance.x_min),6));
            }
        }
    }

    E.em=0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                Particles::min_var min_distance = min_r(r[i],r[j]);
                if(min_distance.x_min < rc) {
                    double scalar = 0;
                    for (int ii = 0; ii < 3; ii++) {
                        scalar += spins[i][ii]*spins[j][ii];
                    }
                    E.em -= j_0*pow(1-min_distance.x_min/rc,3)*scalar;
                }
            }
        }
    }

    E.ep/=2;
    E.ef = E.ek + E.ep + E.em;
    return E;
}

int Particles::moving(){
    auto start = std::chrono::system_clock::now();
    for (int iter = 0; iter < 500; iter++) {
        tscale.push_back(iter);
        std::vector<std::vector<double>> r_d;
        std::vector<std::vector<double>> v_d;
        std::vector<std::vector<double>> a_d = a;

        for (int i = 0; i < n; i++) {
            std::vector<double> rr;
            rr.push_back(r[i][0]+v[i][0]*dt+a[i][0]*pow(dt,2)/2);
            rr.push_back(r[i][1]+v[i][1]*dt+a[i][1]*pow(dt,2)/2);
            r_d.push_back(rr);

            if (r_d[i][0] > l or r_d[i][0] < 0){
                r_d[i][0] -= l*(int)(r_d[i][0]/l);
            }
            if (r_d[i][1] > l or r_d[i][1] < 0){
                r_d[i][1] -= l*(int)(r_d[i][1]/l);
            }
        }

        change_spins(r_d);
        count_acc(r_d);

        for (int i = 0; i < n; i++) {
            std::vector<double> vv;
            vv.push_back(v[i][0]+0.5*(a[i][0]+a_d[i][0])*dt);
            vv.push_back(v[i][1]+0.5*(a[i][1]+a_d[i][1])*dt);
            v_d.push_back(vv);
        }

        Particles::e_var Energy = energy();
        std::cout << "Total energy is: " << Energy.ef << std::endl;
        std::cout << "Magnit energy is: " << Energy.em << std::endl;
        std::cout << std::endl;
        ekscale.push_back(Energy.ek);
        epscale.push_back(Energy.ep);
        r = r_d;
        v = v_d;
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;
    return 0;
}

//methods of spin

int Particles::create_spins(){
    for (int i = 0; i < n; i++) {
        std::vector<double> din_s;
        for (int j = 0; j < 3; j++) {
            din_s.push_back(rand() % 10000000 / 10000000.0);
        }
        din_s = normalized(din_s);
        spins.push_back(din_s);
    }
    return 0;
}

std::vector<double> Particles::normalized(std::vector<double> massive){
    double z = 0;
    for (int i = 0; i < massive.size(); i++) {
        z += pow(massive[i],2);
    }
    for (int i = 0; i < massive.size(); i++) {
        massive[i] *= pow(spin_l/z,0.5);
    }
    return massive;
}

int Particles::change_spins(std::vector<std::vector<double>> r){

    std::vector<std::vector<double>> H_ex;
    for (int i = 0; i < n; i++) {
        std::vector<double> din_vect;
        H_ex.push_back(din_vect);
        for (int j = 0; j < 3; j++) {
            H_ex[i].push_back(0);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                Particles::min_var min_dist = min_r(r[i],r[j]);
                if (min_dist.x_min < rc) {
                    for (int ii = 0; ii < 3; ii++) {
                        H_ex[i][ii] += -j_0*pow(1-min_dist.x_min/rc,3)*spins[j][ii];
                    }
                }
            }
        }
    }

    std::vector<std::vector<double>> ds;
    for (int i = 0; i < n; i++) {
        std::vector<double> din_vect;
        ds.push_back(din_vect);
        for (int j = 0; j < 3; j++) {
            ds[i].push_back(0);
        }
    }

    for (int i = 0; i < n; i++) {
        ds[i][0] = (1/(h*spin_l))*(H_ex[i][1]*spins[i][2]-H_ex[i][2]*spins[i][1]);
        ds[i][1] = -(1/(h*spin_l))*(H_ex[i][0]*spins[i][2]-H_ex[i][2]*spins[i][0]);
        ds[i][2] = (1/(h*spin_l))*(H_ex[i][0]*spins[i][1]-H_ex[i][1]*spins[i][0]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            spins[i][j] += ds[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        std::vector<double> din_s = normalized(spins[i]);
        spins[i] = din_s;
    }
}