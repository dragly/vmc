#include <stdlib.h>

#include "densityplotter.h"

#include "wavefunction.h"
#include "inih/cpp/INIReader.h"
#include "config.h"
#include "montecarlostandard.h"

DensityPlotter::DensityPlotter(Config *config) :
    m_config(config)
{

}

void DensityPlotter::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    m_charge = atof(settings->Get("DensityPlotter", "charge", "1.0").c_str());
    m_stepLength = atof(settings->Get("DensityPlotter", "stepLength", "1.0").c_str());
    m_wave = WaveFunction::fromName(settings->Get("Wave", "class", "WaveSimple"), m_config);
    m_hamiltonian = Hamiltonian::fromName(settings->Get("Hamiltonian", "class", "HamiltonianSimple"), m_config, 1.0);
    m_nCycles = settings->GetInteger("DensityPlotter", "nCycles", 1000);
}

void DensityPlotter::makePlot()
{
    double *energies = new double[2];
    double *allEnergies = new double[m_nCycles];
    double wfnew = 0;
    double wfold = 0;
    double energy = 0;
    double energy2 = 0;
    double delta_e = 0;
    // initialisations of variational parameters and energies
    energy = energy2 = 0; delta_e=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < number_particles; i++) {
        for (int j=0; j < dimension; j++) {
            r_old[i][j] = step_length*(ran2(&idum)-0.5);
        }
    }
    wfold = m_wave->wave(r_old);
    // loop over positions
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new position
        for (int i = 0; i < number_particles; i++) {
            for (int j=0; j < dimension; j++) {
                r_new[i][j] = r_old[i][j]+step_length*(ran2(&idum)-0.5);
            }
            // TODO: Optimize this by removing the if-test. Profile first!
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < number_particles; k++) {
                if ( k != i) {
                    for (int l=0; l < dimension; l++) {
                        r_new[k][l] = r_old[k][l];
                    }
                }
            }
            wfnew = m_wave->wave(r_new);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(&idum) <= wfnew*wfnew/wfold/wfold ) {
                for (int l=0; l < dimension; l++) {
                    r_old[i][l]=r_new[i][l];
                }
                wfold = wfnew;
            }
        }  //  end of loop over particles
        // compute local energy
        delta_e = m_hamiltonian->energy(m_wave, r_old);
        // save all energies on last variate
        //        if(variate==max_variations){
        allEnergies[cycle] = delta_e;
        //        }
        // update energies
        energy += delta_e;
        energy2 += delta_e*delta_e;
    }   // end of loop over MC trials
    // return the energy and the energy squared
    energies[0] = energy;
    energies[1] = energy2;
    monteCarlo->sample(m_nCycles, energies, allEnergies);
}
