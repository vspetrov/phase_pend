#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <math.h>
#include <stdio.h>
#include <numeric>
using namespace std;

class ChainPend {
public:
    ChainPend(int _N) {
        N = _N;
        Phi.resize(N);
        tmpPhi.resize(N);
        tmpX.resize(N);
        X.resize(N);
        std::fill(Phi.begin(),Phi.end(),0.1);
        std::fill(X.begin(), X.end(), 0.1);
        for (int j=0; j<4; j++) {
            for (int i=0; i<2; i++) {
                rk[j][i].resize(N);
                std::fill(rk[j][i].begin(), rk[j][i].end(), 0.0);
            }
        }
        lambda = 0.3;
        gamma.resize(N);
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(1.0,1.5);
        double g0 = 1.01;
        double g_max = g0+0.02;
        for (int i=0; i<N; i++) {
            // gamma[i] = distribution(generator);
            gamma[i] = g0 + (g_max-g0)*i/(N-1);
        }

        d = 0;
        Dt = 0.1;
        freqCalcSkip = 0;
        omega.resize(N);
    }
    double reset_gamma(double g0) {
        double g_max = g0+0.2;
        for (int i=0; i<N; i++) {
            // gamma[i] = distribution(generator);
            gamma[i] = g0 + (g_max-g0)*i/(N-1);
        }

    }
    double rhs_phi(int i) {
        return X[i];
    }

    double rhs_x(int i) {
        double rst = -lambda*X[i]-sin(Phi[i])+gamma[i];
        if (i == 0) {
            rst += d*sin(Phi[1]-Phi[0]);
        } else if (i == N-1) {
            rst += d*sin(Phi[N-2]-Phi[N-1]);
        } else {
            rst += d*sin(Phi[i-1]-Phi[i]) + d*sin(Phi[i+1]-Phi[i]);
        }
        return rst;
    }

    void reset() {
        std::fill(Phi.begin(),Phi.end(),0.1);
        std::fill(X.begin(), X.end(), 0.1);
    }
    void make_step(const double dt) {
        std::copy(Phi.begin(),Phi.end(),tmpPhi.begin());
        std::copy(X.begin(),X.end(),tmpX.begin());
        for (int i=0; i<N; i++) {
            rk[0][0][i] = dt*rhs_phi(i);
            rk[0][1][i] = dt*rhs_x(i);
        }

        for (int i=0; i<N; i++) {
            Phi[i] = tmpPhi[i]+rk[0][0][i]/2.0;
            X[i] = tmpX[i]+rk[0][1][i]/2.0;
        }
        for (int i=0; i<N; i++) {
            rk[1][0][i] = dt*rhs_phi(i);
            rk[1][1][i] = dt*rhs_x(i);
        }

        for (int i=0; i<N; i++) {
            Phi[i] = tmpPhi[i]+rk[1][0][i]/2.0;
            X[i] = tmpX[i]+rk[1][1][i]/2.0;
        }
        for (int i=0; i<N; i++) {
            rk[2][0][i] = dt*rhs_phi(i);
            rk[2][1][i] = dt*rhs_x(i);
        }

        for (int i=0; i<N; i++) {
            Phi[i] = tmpPhi[i]+rk[2][0][i];
            X[i] = tmpX[i]+rk[2][1][i];
        }
        for (int i=0; i<N; i++) {
            rk[3][0][i] = dt*rhs_phi(i);
            rk[3][1][i] = dt*rhs_x(i);
        }

        for (int i=0; i<N; i++) {
            Phi[i] = tmpPhi[i] + (rk[0][0][i] + 2*rk[1][0][i] + 2*rk[2][0][i] + rk[3][0][i])/6.0;
            X[i] = tmpX[i] + (rk[0][1][i] + 2*rk[1][1][i] + 2*rk[2][1][i] + rk[3][1][i])/6.0;
        }
    }

    void solve(double MaxTime) {
        double time = 0;
        bool freqStarted = false;
        vector<double> p1;
        p1.resize(N);

        while (time < MaxTime) {
            if (!freqStarted && time > freqCalcSkip) {
                freqStarted = true;
                std::copy(Phi.begin(),Phi.end(),p1.begin());
            }

#if 0
            rst.push_back(Phi);
            time_v.push_back(time);
#endif
            make_step(Dt);

            time += Dt;
        }
        for (int i=0; i<N; i++) {
            omega[i] = (Phi[i]-p1[i])/(MaxTime - freqCalcSkip);
        }

    }

    vector<vector<double> > rst;
    vector<double> time_v;
    double freqCalcSkip;
    vector<double> omega;
    double d;
    double lambda;
    int N;
private:

    double Dt;
    vector<double> Phi;
    vector<double> tmpPhi;
    vector<double> tmpX;
    vector<double> X;
    vector<double> rk[4][2];
    vector<double> gamma;
};

int main(int argc, char **argv)
{
    ChainPend CP(100);
    CP.freqCalcSkip = 1500;

    CP.d = 5;

    FILE *ofs=fopen("f_lambda2.dat","w");
    for (double lambda=0; lambda < 1; lambda += 0.01) {
        fprintf(stderr,"lambda = %g\n",lambda);
        CP.reset();
        CP.lambda = lambda;
        CP.solve(3000);
        double av_f = std::accumulate(CP.omega.begin(), CP.omega.end(), 0.0) / CP.N;
        bool is_sync = true;

        for (int i=0; i<CP.N; i++) {
            double om = CP.omega[i];
            if (fabs(om-av_f)/av_f > 0.01) {
                is_sync = false;
                break;
            }
        }
        if (is_sync) {
            fprintf(ofs,"%g %g\n",lambda,av_f);
        }
    }

    fclose(ofs);

#if 0
    FILE *ofs=fopen("omegas_l_0.3.dat","w");
    for (double d=0; d<5+1e-5; d+=0.01) {
        printf("calculating d = %g ...",d);fflush(stdout);
        CP.d = d;
        CP.solve(3000);
        fprintf(ofs,"%g ", d);
        for (int i=0; i<CP.N; i++) {
            fprintf(ofs,"%g ",CP.omega[i]);
        }
        fprintf(ofs,"\n");
        printf(" done.\n");fflush(stdout);
        CP.reset();
    }
    fclose(ofs);
#endif

#if 0
    FILE *ofs = fopen("rst.dat","w");
    for (int i=0; i<CP.rst.size(); i++) {
        fprintf(ofs,"%g ",CP.time_v[i]);
        for (int j=0; j<CP.rst[i].size(); j++)
            fprintf(ofs,"%g ", CP.rst[i][j]);
        fprintf(ofs,"\n");
    }
    fclose(ofs);
#endif

    return 0;
}
