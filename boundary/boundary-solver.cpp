#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    const int Nx = 160000*4;
    const int Ny = 120;
    const double delta0 = 0.05;
    const double L = 0.02*2;
    const double rho = 1;
    const double mu = 1.95e-5;
    const double nu = mu / rho;
    const double dp_dx = 0.0;
    const double U_x = 1;
    const double kappa = 0.41;

    const int yp_iter = 160000;

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

    cout << "dx = " << dx << endl;
    cout << "dy = " << pow(dy[0], 2) << endl;

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

        dudy[i][0] = (u[i][1] - u[i][0]) / dy[0];       // Wall gradient

        u[i][0] = 0;        // No slip condition
        u[i][Ny-1] = U_x;     // Outer stream flow condition U(x) = 1

        v[i][Ny-1] = 0;

        for (int j = Ny-1; j > 0; --j) {
            double du_dx = (u[i][j] - u[i-1][j]) / dx;
            v[i][j-1] = v[i][j] - dy[j] * du_dx;
        }
    }

    //double dudy_w = (u[yp_iter][1] - u[yp_iter][0]) / dy[0];
    double tau_w = mu * dudy[yp_iter][0];
    double cf = (2 * tau_w) / (rho * pow(U_x, 2));
    cout << "dudy = " << dudy[yp_iter][0] << endl;
    cout << "cf = " << cf << endl;
    double u_tau = sqrt(tau_w / rho);
    double delta_v = nu / u_tau;
    cout << "delta_v = " << delta_v << endl;
    
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
    outfile << "y" << "," << "u[0]" << "," << "u[10000]" << "," << "u[80000]" << ",";
    outfile << "u[160000]" << "," << "u[320000]" << "," << "u[640000]" << ",";
    outfile << "y+" << "," << "u+" << "," << "tau_v[40000]" << "," << "tau_t[40000]" << ",";
    outfile << "tau_v[640000]" << "," << "tau_t[640000]" << endl;
    for (int j = 0; j<Ny; ++j) {
        y_plus[j] = y[j] * u_tau / nu;
        outfile << y[j] << "," << u[0][j] << "," << u[10000][j] << "," << u[80000][j] << ",";
        outfile << u[160000][j] << "," << u[320000-1][j] << "," << u[640000-1][j] << ",";
        outfile << y_plus[j] << "," << 0 << "," << tau_v[40000][j] << "," << tau_t[40000][j] << ",";
        outfile << tau_v[640000-1][j] << "," << tau_t[640000-1][j] << endl;
    }
    outfile.close();

    cout << "Simulation complete. Data saved to boundary-layer.csv" << endl;

    //system("pause");
}