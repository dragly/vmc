#include "hamiltonianideal.h"
#include "../matrix.h"
#include "../utils.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"
#include <math.h>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

/*!
  * Constructs an Hamiltonian of the type used in the report
  */
HamiltonianIdeal::HamiltonianIdeal(Config *config) :
    Hamiltonian(config),
    omega(config->omega())
{
}
int ntries = 0;

/*!
  * Calculates the kinetic energy as $-\frac{1}{2} \frac{1}{\Psi} \nabla^2 \Psi$
  */
double HamiltonianIdeal::kineticEnergy(WaveFunction *wave, vec2 r[])
{
    return -0.5*wave->laplace(r);
}

/*!
  * Returns the sum of external potential and interaction potential energy.
  */
double HamiltonianIdeal::potentialEnergy(WaveFunction *wave, vec2 r[])
{
    double ePotential = 0;
    ePotential += externalPotentialEnergy(wave, r);
    if(m_interactionEnabled) {
        ePotential += interactionPotentialEnergy(wave, r);
    }
    return ePotential;
}

/*!
  * Returns the energy from the harmonic oscillator potential
  */
double HamiltonianIdeal::externalPotentialEnergy(WaveFunction *wave, vec2 r[])
{
    (void)wave;
    double externalEnergy = 0;
    double rSingleParticle = 0;
    // contribution from harmonic oscillator potential
    for (int i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (int j = 0; j < m_nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        externalEnergy += rSingleParticle;
    }
    externalEnergy = 0.5 * omega * omega * externalEnergy;
    m_totalExternalPotentialEnergy += externalEnergy;
    return externalEnergy;
}

/*!
  * Returns the energy from the particle interaction.
  */
double HamiltonianIdeal::interactionPotentialEnergy(WaveFunction *wave, vec2 r[])
{
    (void)wave;
    double interactionEnergy = 0;
    //     contribution from electron-electron potential
    // TODO Optimization - store an array of potential energies between particles and only update necessary parts (reduces use of sqrt)
    for (int i = 0; i < m_nParticles-1; i++) {
        for (int j = i+1; j < m_nParticles; j++) {
            double distanceSquared = 0;
            for (int k = 0; k < m_nDimensions; k++) {
                double dist = (r[i][k]-r[j][k]);
                distanceSquared += dist*dist;
            }
            interactionEnergy += 1/sqrt(distanceSquared);
        }
    }
    m_totalInteractionPotentialEnergy += interactionEnergy;
    return interactionEnergy;
}

void HamiltonianIdeal::resetTotalEnergies() {
    Hamiltonian::resetTotalEnergies();
    m_totalExternalPotentialEnergy = 0;
    m_totalInteractionPotentialEnergy = 0;
}

string HamiltonianIdeal::totalsString()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    stringstream myString;
    myString << totalKineticEnergy() / (totalEnergy) << " " << totalInteractionPotentialEnergy() / (totalEnergy) << " " << totalExternalPotentialEnergy() / (totalEnergy);
    return myString.str();
}

void HamiltonianIdeal::outputTotals()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    std::cout << "Kinetic energy: " << setprecision(3) << totalKineticEnergy() / (totalEnergy) << " interactive potential energy: " << setprecision(3) << totalInteractionPotentialEnergy() / (totalEnergy) << " external potential energy: " << setprecision(3) << totalExternalPotentialEnergy() / (totalEnergy) << std::endl;
}
