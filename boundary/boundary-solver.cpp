#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    const int Nx = 160000;
    const int Ny = 80;
    const double delta0 = 0.05;
    const double L = 0.02;
    const double rho = 1;
    const double mu = 1.95e-5;
    const double nu = mu / rho;
    const double dp_dx = 0.0;
    const double kappa = 0.41;

    const int yp_iter = 100;

    int count = 0;

    vector<double> x(Nx);
    vector<double> y(Ny);
    for (int i = 0; i < Nx; ++i) {
        x[i] = L * i / (Nx-1);
    }

    for (int j = 0; j < Ny; ++j) {
        y[j] = pow(1.0 * j / (Ny-1), 1.3) * delta0;
    }

    double dx = x[1] - x[0];
    vector<double> dy(Ny);
    dy[0] = y[1] - y[0];
    dy[Ny-1] = y[Ny-1] - y[Ny-2];

    for (int j=1; j < Ny-1; ++j) {
        dy[j] = (y[j+1] - y[j-1]) / 2;
    }

    vector<vector<double>> u(Nx, vector<double>(Ny, 0));
    vector<vector<double>> v(Nx, vector<double>(Ny, 0));
    vector<vector<double>> nut(Nx, vector<double>(Ny, 0));
    vector<vector<double>> mut(Nx, vector<double>(Ny, 0));
    vector<vector<double>> tau_v(Nx, vector<double>(Ny, 0));
    vector<vector<double>> tau_t(Nx, vector<double>(Ny, 0));
    vector<vector<double>> dudy(Nx, vector<double>(Ny, 0));

    for (int j = 0; j < Ny; ++j) {
        double eta = y[j] / delta0;
        u[0][j] = (eta < 1.0) ? (1.5 * eta - 0.5 * pow(eta, 3)) : 1.0;
    }

    dudy[0][0] = (u[0][1] - u[0][0]) / dy[0];
    dudy[0][Ny-1] = (u[0][Ny-1] - u[0][Ny-2]) / dy[Ny-1];

    for (int j = 1; j < Ny-1; ++j) {
        dudy[0][j] = (u[0][j+1] - u[0][j-1]) / (2*dy[j]);
    }

    for (int i = 1; i < Nx; ++i) {

        for (int j = 0; j < Ny; ++j) {
            nut[i-1][j] = pow((kappa * y[j]), 2) * fabs(dudy[i-1][j]);
            mut[i-1][j] = rho * nut[i-1][j];
        }

        for (int j = 1; j < Ny-1; ++j) {
            double du_dy = (u[i-1][j+1] - u[i-1][j-1]) / (2*dy[j]);
            dudy[i][j] = du_dy;

            double d2u_dy2 = (u[i-1][j+1] - 2 * u[i-1][j] + u[i-1][j-1]) / pow(dy[j], 2);

            double visc_term = (mu + mut[i-1][j]) * d2u_dy2 + (
                (mut[i-1][j+1] - mut[i-1][j-1]) / (2*dy[j])
            ) * du_dy;

            tau_v[i-1][j] = mu * du_dy;
            tau_t[i-1][j] = mut[i-1][j]*du_dy;

            double du_dx = (fabs(u[i-1][j]) > 1e-8) ? (
                ((visc_term - dp_dx) / rho - v[i-1][j] * du_dy) / u[i-1][j]
            ) : 0;

            u[i][j] = u[i-1][j] + dx * du_dx;
        }

        u[i][0] = 0;        // No slip condition
        u[i][Ny-1] = 1;     // Outer stream flow condition U(x) = 1

        v[i][Ny-1] = 0;

        for (int j = Ny-1; j > 0; --j) {
            double du_dx = (u[i][j] - u[i-1][j]) / dx;
            v[i][j-1] = v[i][j] - dy[j] * du_dx;
        }
    }

    double dudy_w = (u[yp_iter][1] - u[yp_iter][0]) / dy[0];
    double tau_w = mu * dudy_w;
    double u_tau = sqrt(tau_w / rho);
    
    vector<double> y_plus(Ny);
    vector<vector<double>> u_plus(Nx, vector<double>(Ny, 0.0));

//    for (int j = 0; j < Ny; ++j) {
//        y_plus[j] = y[j] * u_tau / nu;
//    }

//    for (int i = 0; i < Nx; ++i) {          // Add to main loop for greater efficiency
//        for (int j = 0; j < Ny; ++j) {
//            u_plus[i][j] = u[i][j] / u_tau;
//        }
//    }

    ofstream outfile("./boundary-layer.csv");
    outfile << "y" << "," << "u[0]" << "," << "u[1000]" << "," << "u[10'000]" << ",";
    outfile << "u[40'000]" << "," << "u[80'000]" << "," << "u[160'000]" << ",";
    outfile << "y+" << "," << "u+" << endl;
    for (int j = 0; j<Ny; ++j) {
        y_plus[j] = y[j] * u_tau / nu;
        outfile << y[j] << "," << u[0][j] << "," << u[1000][j] << "," << u[10000][j] << ",";
        outfile << u[40000][j] << "," << u[80000-1][j] << "," << u[160000-1][j] << ",";
        outfile << y_plus[j] << "," << 0 << "," << endl;
    }
    outfile.close();

    cout << "Simulation complete. Data saved to boundary-layer.csv" << endl;

    system("pause");
}