#include <iostream>
#include <fstream> 
#include <vector>
#include <cmath>
#include <random>
using namespace std;

constexpr int dimension = 3; // Dimension of the simulation box (the primary cell)

// Time Variables in femtoseconds
constexpr int steps = 100000; // Total number of timesteps
constexpr double dt = 1E-4; // Time step
constexpr double total_time = steps * dt; // Total time 

// Particle Variables
constexpr int N = 30; // Number of particles
constexpr double mass = 6.6335; // Mass of particles in units of 10**(-26) kg
constexpr double lattice_constant = 5.4691; // Lattice constant for Argon in Armstrong 

// System variables 
constexpr double Temp0 = 273; // Temperature in Kelvin 

// Variables for algorithm
constexpr double boltz_const = 1.380649; // Boltzmann constant in J/K
constexpr double sigma = 3.405; // σ in Armstrong 
constexpr double epsilon = 1.6544; // ε in units of 10^-21 J
constexpr double rcutoff = 2.5 * sigma; // Cutoff radius 
constexpr double box_size = 10 * sigma; // Box size 
const double u_at_cutoff = 4 * epsilon * (pow((sigma / rcutoff), 12) - pow((sigma / rcutoff), 6)); // Value of the potential at the cutoff
const double du_at_cutoff = 24 * epsilon * ((-2 * (pow(sigma, 12) / pow(rcutoff, 13))) + (pow(sigma, 6) / pow(rcutoff, 7))); // Value of the derivative of the potential at the cutoff
const double volume = pow(box_size, 3); // Box volume 
const double density = N / volume; // Particle density
constexpr int input_integral = 10; // File handling
constexpr double v_max = 0.5; // Maximum velocity

// Function 1: Calculates the center-of-mass reference frame
void centerOfMassTranslate(vector<vector<double>>& position) {
    int N = position[0].size();
    int dimensions = position.size();

    vector<double> com(dimensions, 0.0);
    for (int i = 0; i < dimensions; ++i) {
        for (int j = 0; j < N; ++j) {
            com[i] += position[i][j];
        }
        com[i] /= N;
    }
    for (int i = 0; i < dimensions; ++i) {
        for (int j = 0; j < N; ++j) {
            position[i][j] -= com[i];
        }
    }
}

// Function 2: Returns the sign of a real number
double sign(double a, double b) {
    if (b > 0) {
        return a;
    }
    else if (b < 0) {
        return -a;
    }
    else {
        return a;
    }
}

// Function 3: Time evolution loop
void timeEvolution(vector<vector<double>>& position, vector<vector<double>>& velocity) {
    ofstream energy_file("energy.txt");
    ofstream position_file("positions.txt");
    ofstream velocity_file("velocity.txt");

    for (int step = 1; step <= steps; ++step) {
        // Periodic boundary conditions
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < N; ++j) {
                if (position[i][j] > 0.5)
                    position[i][j] -= 1;
                else if (position[i][j] < -0.5)
                    position[i][j] += 1;
            }
        }

        vector<vector<double>> acc(dimension, vector<double>(N, 0.0));
        vector<double> pot(N, 0.0);
        double vrl = 0.0;

        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (j != i) {
                    double r2 = 0.0;
                    vector<double> R(dimension);
                    for (int k = 0; k < dimension; ++k) {
                        R[k] = position[k][i] - position[k][j];
                        if (abs(R[k]) > 0.5)
                            R[k] -= sign(1.0, R[k]);
                        R[k] *= box_size;
                        r2 += R[k] * R[k];
                    }

                    if (r2 < rcutoff * rcutoff) {
                        double r1 = sqrt(r2);
                        double ri2 = 1.0 / r2;
                        double ri6 = ri2 * ri2 * ri2;
                        double ri12 = ri6 * ri6;
                        double sig6 = sigma * sigma * sigma * sigma * sigma * sigma;
                        double sig12 = sig6 * sig6;
                        double u = 4.0 * epsilon * (sig12 * ri12 - sig6 * ri6) - u_at_cutoff - r1 * du_at_cutoff;
                        double du = 24.0 * epsilon * ri2 * (2.0 * sig12 * ri12 - sig6 * ri6) + du_at_cutoff * sqrt(ri2);
                        pot[j] += u;
                        vrl -= du * r2;

                        for (int k = 0; k < dimension; ++k) {
                            acc[k][i] += du * R[k];
                            acc[k][j] -= du * R[k];
                        }
                    }
                }
            }
        }
        // Write particle positions and velocity to file every 100 steps for every tenth particle
        if (step % 100 == 0) {
            position_file << dt * step << " ";
            velocity_file << dt * step << " ";
            for (int j = 0; j < N; j += 10) {
                for (int i = 0; i < dimension; ++i) {
                    position_file << position[i][j] << "  ";
                    velocity_file << velocity[i][j] << "  ";
                }
            }
            position_file << endl;
            velocity_file << endl;

            // Compute kinetic energy
            double kinetic_energy = 0.0;
            double potential_energy = 0.0;
            for (int i = 0; i < N; ++i) {
                kinetic_energy += 0.5 * mass * (velocity[0][i] * velocity[0][i] + velocity[1][i] * velocity[1][i] + velocity[2][i] * velocity[2][i]);
                potential_energy += pot[i];
            }
            // Write energy to file every 100 steps
            energy_file << step * dt << " " << kinetic_energy << " " << potential_energy << " " << kinetic_energy + potential_energy << endl;
        }

        // Update positions
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < N; ++j) {
                position[i][j] += dt * velocity[i][j] + 0.5 * acc[i][j] * dt * dt;
            }
        }

        // Compute temperature and rescale velocities
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < N; ++j) {
                velocity[i][j] += 0.5 * dt * acc[i][j];
            }
        }
    }

    energy_file.close();
    position_file.close();
    velocity_file.close();
}

int main() {
    //--------------------- Initiate positions, velocities, & accelerations -------------------------//
    vector<vector<double>> position(dimension, vector<double>(N));
    vector<vector<double>> velocity(dimension, vector<double>(N));
    vector<vector<double>> acceleration(dimension, vector<double>(N));
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis_pos(0.0, box_size);
    uniform_real_distribution<double> dis_vel(0.0, v_max);

    // Initializing positions and velocities
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < dimension; ++i) {
            position[i][j] = dis_pos(gen);
            velocity[i][j] = dis_vel(gen);
        }
    }
    // Scale positions to the box size
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < N; ++j) {
            position[i][j] /= box_size;
        }
    }
    // Perform center-of-mass translation
    centerOfMassTranslate(position);

    // Time evolution loop
    timeEvolution(position, velocity);

    return 0;
}














