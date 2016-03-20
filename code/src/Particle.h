#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "CanonicalEnsemble.h"

struct Particle2D {
    Particle2D(double q, double x, double y) : q(q), x(x), y(y) {}
    double q = 1.0;
    double x = 0.0;
    double y = 0.0;
};

struct Particle3D {
    Particle3D(double q, double x, double y, double z)
        : q(q), x(x), y(y), z(z) {}
    double q = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

double coulomb_core(const Particle2D &a, const Particle2D &b);
double coulomb_core(const Particle3D &a, const Particle3D &b);

double lennard_jones(const Particle2D &a, const Particle2D &b);
double lennard_jones(const Particle3D &a, const Particle3D &b);

double avg_pair_dist(const CanonicalEnsemble<Particle2D>::State &state);
double avg_pair_dist(const CanonicalEnsemble<Particle3D>::State &state);

#endif // PARTICLE_H_
