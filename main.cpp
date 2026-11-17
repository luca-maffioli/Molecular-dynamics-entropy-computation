/*
 * ============================================================================
 *  Molecular dynamics of a WCA fluid: two- and three-particle distribution
 *  functions (g2 and g3)
 * ============================================================================
 *
 *  Equilibrium MD simulation of a three-dimensional fluid of particles that
 *  interact through the Weeks-Chandler-Andersen potential (Lennard-Jones
 *  truncated at r_c = 2^(1/6) sigma). The particles are in a periodic cubic
 *  box.
 *
 *    - Integrator : 5th-order Gear predictor-corrector
 *    - Thermostat : Gaussian isokinetic thermostat with a proportional
 *                   feedback term (or isoenergetic dynamics)
 *    - Forces     : cell lists with ghost cells for periodic boundaries
 *    - Sampling   : histograms of the two-particle (g2) and three-particle
 *                   (g3) distribution functions, also built with cell lists
 *
 *  The program runs for a fixed amount of CPU time (simtime). Every savetime
 *  seconds it saves the histograms, the current microstate and a short
 *  report.
 *
 *  Compile:  g++ -O3 -march=native main.cpp -o g3md
 *
 *  See DOCUMENTATION.md for a description of every function and output file.
 * ============================================================================
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


/* ============================================================================
 *  Physical parameters
 * ========================================================================== */

const int    nside   = 15;                         // lattice cells per side of the box (initial configuration)
const int    npart   = 2 * nside * nside * nside;  // number of particles (2 per lattice cell, bcc)
const double dens    = 1;                          // number density
const double sigma   = 1;                          // particle diameter
const double epsilon = 1;                          // LJ energy parameter
const double mass    = 1;                          // particle mass
const double temp    = 1;                          // temperature

const int    freed  = 3;  // degrees of freedom removed by the Gaussian thermostat (= dimension of the system)
const double fback  = 5;  // coefficient of the Gaussian feedback term

const int constrain = 1;  // 1 = isokinetic (thermostatted), 0 = isoenergetic


/* ============================================================================
 *  Run-time parameters
 * ========================================================================== */

const double simtime  = 86400 * 4.5;  // total simulation time (CPU seconds)
const double savetime = 3600 * 1;     // interval between saves (CPU seconds)
const double dt       = 0.001;        // integration time step


/* ============================================================================
 *  Box and force cell list
 * ========================================================================== */

const double length  = pow(npart / dens, 1. / 3);   // box side length
const double r_cut   = pow(2, 1. / 6) * sigma;       // interaction cutoff (WCA)
const double rc2     = r_cut * r_cut;
const double prefact = 24 * epsilon / (sigma * sigma);  // prefactor of the LJ force
const int    sigma2  = sigma * sigma;

const int    nfcell = length / r_cut;          // force cells per side
const double fdcell = length / nfcell;         // force cell size
const int    flayer = ceil(r_cut / fdcell);    // layers of neighbour cells needed to cover r_cut
const int    tfcell = nfcell + 2 * flayer;     // force cells per side, including ghost layers
const int    upf    = tfcell - flayer - 1;     // index of the last real force cell along a side

const int nfc2 = nfcell * nfcell;
const int nfc3 = nfcell * nfcell * nfcell;
const int tfc2 = tfcell * tfcell;
const int tfc3 = tfcell * tfcell * tfcell;
const int ofc  = tfcell * tfcell * tfcell - nfcell * nfcell * nfcell;  // number of ghost force cells


/* ============================================================================
 *  Distribution-function cell list
 * ========================================================================== */

const int    ndcell   = pow(npart / 8, 1. / 3);  // distribution cells per side
const double ddcell   = length / ndcell;         // distribution cell size
const int    maxlayer = (ndcell - 1) / 2;        // largest stencil that never reaches the same cell twice

const double dthreshold_temp = length / 2;                         // upper limit from the minimum-image convention
const double dthreshold_max  = maxlayer * ddcell;                  // upper limit from the cell stencil
const double dthreshold      = min(dthreshold_max, dthreshold_temp);  // largest distance sampled in g2 and g3
const double dt2             = dthreshold * dthreshold;

const int dlayer_temp = dthreshold_temp / ddcell + 1;
const int dlayer      = min(dlayer_temp, maxlayer);  // layers of neighbour cells for the distribution functions
const int tdcell      = ndcell + 2 * dlayer;         // distribution cells per side, including ghost layers
const int upd         = tdcell - dlayer - 1;         // index of the last real distribution cell along a side

const int ndc2 = ndcell * ndcell;
const int ndc3 = ndcell * ndcell * tdcell;
const int tdc2 = tdcell * tdcell;
const int tdc3 = tdc2 * tdcell;


/* ============================================================================
 *  Histograms
 * ========================================================================== */

const int pbin  = 600;                               // number of bins of g2
const int tbin  = 600;                               // number of bins per distance of g3
const int ttbin = tbin * (tbin + 1) * (tbin + 2) / 6;  // total number of stored g3 bins

// g3 is stored as a 1D vector. Only the region i >= j >= k is kept, because
// g3 is symmetric under any permutation of its three arguments. The full
// 3D histogram is recovered as
//
//     g3(i,j,k) = triad[ i*(i+1)*(i+2)/6 + (j+1)*j/2 + k ],   i >= j >= k
//
// and g3(i,j,k) = g3(j,i,k) = g3(k,i,j) = ... for all 6 permutations.

const double psbin = dthreshold / pbin;  // bin width of g2
const double tsbin = dthreshold / tbin;  // bin width of g3


/* ============================================================================
 *  Gear predictor-corrector coefficients
 * ========================================================================== */

const double PI  = 3.141592653;
const double fo2 = 3.0 / 16.0;
const double f12 = 251.0 / 360.0;
const double f32 = 11.0 / 18.0;
const double f42 = 1.0 / 6.0;
const double f52 = 1.0 / 60.0;
const double fo1 = 251.0 / 720.0;
const double f21 = 11.0 / 12.0;
const double f31 = 1.0 / 3.0;
const double f41 = 1.0 / 24.0;


/* ============================================================================
 *  Data structures
 * ========================================================================== */

// Phase-space state of one particle: position, momentum and their scaled time
// derivatives up to 4th order (Gear algorithm), plus the cell-list links.
typedef struct
{
    double x0, x1, x2, x3, x4;
    double y0, y1, y2, y3, y4;
    double z0, z1, z2, z3, z4;

    double px0, px1, px2, px3, px4;
    double py0, py1, py2, py3, py4;
    double pz0, pz1, pz2, pz3, pz4;

    int findex;  // next particle in the same force cell (-1 = end of list)
    int dindex;  // next particle in the same distribution cell (-1 = end of list)
} Phase;

// Total force acting on one particle
typedef struct
{
    double fx;
    double fy;
    double fz;
} Interact;

// Force cell: head of its particle list and the periodic shift of the cell
// (0 for real cells, +/- length for ghost cells)
typedef struct
{
    int    flist;
    double x;
    double y;
    double z;
} Fcellslist;

// Distribution cell: same layout as Fcellslist
typedef struct
{
    int    dlist;
    double x;
    double y;
    double z;
} Dcellslist;

// Indices of the real force cells lying on the two faces of the box normal
// to one axis
typedef struct
{
    int n[2 * nfc2];
} Borderfcells;

// Pair {ghost cell index, real cell it is a copy of}
typedef struct
{
    int n[2];
} Outerfcells;


/* ============================================================================
 *  Periodic boundary conditions
 * ========================================================================== */

// Wrap a coordinate back into [-length/2, length/2)
void pbcs(double *x)
{
    if (*x <= -length / 2)
    {
        *x += length;
    }
    else if (*x >= length / 2)
    {
        *x -= length;
    }
}

// Map a ghost force-cell coordinate onto the corresponding real cell
void fdpbcs(int *x)
{
    if (*x < flayer)
    {
        *x += nfcell;
    }
    else if (*x > upf)
    {
        *x -= nfcell;
    }
}

// Map a ghost distribution-cell coordinate onto the corresponding real cell
void ddpbcs(int *x)
{
    if (*x < dlayer)
    {
        *x += ndcell;
    }
    else if (*x > upd)
    {
        *x -= ndcell;
    }
}


/* ============================================================================
 *  Initialisation
 * ========================================================================== */

// Seed the random number generator with the current time
void rand0(void)
{
    srand((int)time(NULL));
}

// Place the particles on a bcc lattice and set all position derivatives to 0
void setLattice(Phase *phase)
{
    int i, j, k, c;
    double cell = length / (nside);

    c = 0;
    for (i = 1; i <= nside; i++)
    {
        for (j = 1; j <= nside; j++)
        {
            for (k = 1; k <= nside; k++)
            {
                phase[c].x0 = (i - 0.75) * cell;
                phase[c].y0 = (j - 0.75) * cell;
                phase[c].z0 = (k - 0.75) * cell;
                c++;
                phase[c].x0 = (i - 0.25) * cell;
                phase[c].y0 = (j - 0.25) * cell;
                phase[c].z0 = (k - 0.25) * cell;
                c++;
            }
        }
    }

    for (i = 0; i < npart; i++)
    {
        pbcs(&phase[i].x0);
        pbcs(&phase[i].y0);
        pbcs(&phase[i].z0);
    }

    for (i = 0; i < npart; i++)
    {
        phase[i].x1 = 0;
        phase[i].y1 = 0;
        phase[i].z1 = 0;

        phase[i].x2 = 0;
        phase[i].y2 = 0;
        phase[i].z2 = 0;

        phase[i].x3 = 0;
        phase[i].y3 = 0;
        phase[i].z3 = 0;

        phase[i].x4 = 0;
        phase[i].y4 = 0;
        phase[i].z4 = 0;
    }
}

// Draw random momenta, remove the centre-of-mass momentum, rescale to the
// target temperature and set all momentum derivatives to 0
void randomMomenta(Phase *phase)
{
    double pxsum = 0, pysum = 0, pzsum = 0, kinetictemp = 0, scale;
    int i;

    rand0();

    for (i = 0; i < npart; i++)
    {
        phase[i].px0 = ((double)rand() / (RAND_MAX));
        phase[i].py0 = ((double)rand() / (RAND_MAX));
        phase[i].pz0 = ((double)rand() / (RAND_MAX));
    }

    // Remove the total momentum
    for (i = 0; i < npart; i++)
    {
        pxsum += phase[i].px0;
        pysum += phase[i].py0;
        pzsum += phase[i].pz0;
    }

    pxsum /= npart;
    pysum /= npart;
    pzsum /= npart;

    for (i = 0; i < npart; i++)
    {
        phase[i].px0 -= pxsum;
        phase[i].py0 -= pysum;
        phase[i].pz0 -= pzsum;
    }

    // Rescale to the target temperature
    for (i = 0; i < npart; i++)
    {
        kinetictemp += (phase[i].px0 * phase[i].px0) + (phase[i].py0 * phase[i].py0) + (phase[i].pz0 * phase[i].pz0);
    }

    scale = sqrt(3 * npart * temp / kinetictemp);

    for (i = 0; i < npart; i++)
    {
        phase[i].px0 *= scale;
        phase[i].py0 *= scale;
        phase[i].pz0 *= scale;
    }

    for (i = 0; i < npart; i++)
    {
        phase[i].px1 = 0;
        phase[i].py1 = 0;
        phase[i].pz1 = 0;

        phase[i].px2 = 0;
        phase[i].py2 = 0;
        phase[i].pz2 = 0;

        phase[i].px3 = 0;
        phase[i].py3 = 0;
        phase[i].pz3 = 0;

        phase[i].px4 = 0;
        phase[i].py4 = 0;
        phase[i].pz4 = 0;
    }
}

// Build the static cell-list tables:
//   - periodic shift of every force and distribution cell (0 for real cells)
//   - list of ghost force cells and the real cells they copy (outfcells)
//   - list of real force cells on the box faces (bfc)
//   - spherical truncation of the distribution-function stencil (y, z)
void Bordercells(Dcellslist *dcellslist, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, int *y, int *z)
{
    int cx, cy, cz, ix, iy, iz;
    int i, j;

    // Force cells: periodic shifts and ghost -> real mapping
    int c = 0;
    for (cx = 0; cx < tfcell; cx++)
    {
        for (cy = 0; cy < tfcell; cy++)
        {
            for (cz = 0; cz < tfcell; cz++)
            {
                i = cx * tfc2 + cy * tfcell + cz;

                fcellslist[i].x = 0;
                fcellslist[i].y = 0;
                fcellslist[i].z = 0;

                if (cx < flayer) fcellslist[i].x = -length;
                if (cx > upf)    fcellslist[i].x =  length;
                if (cy < flayer) fcellslist[i].y = -length;
                if (cy > upf)    fcellslist[i].y =  length;
                if (cz < flayer) fcellslist[i].z = -length;
                if (cz > upf)    fcellslist[i].z =  length;

                bool real_cell = (cx >= flayer && cy >= flayer && cz >= flayer && cx <= upf && cy <= upf && cz <= upf);
                if (!real_cell)
                {
                    ix = cx;
                    iy = cy;
                    iz = cz;
                    fdpbcs(&ix);
                    fdpbcs(&iy);
                    fdpbcs(&iz);

                    j = ix * tfc2 + iy * tfcell + iz;
                    outfcells[c].n[0] = i;
                    outfcells[c].n[1] = j;
                    c++;
                }
            }
        }
    }

    // Distribution cells: periodic shifts
    for (cx = 0; cx < tdcell; cx++)
    {
        for (cy = 0; cy < tdcell; cy++)
        {
            for (cz = 0; cz < tdcell; cz++)
            {
                i = cx * tdc2 + cy * tdcell + cz;

                dcellslist[i].x = 0;
                dcellslist[i].y = 0;
                dcellslist[i].z = 0;

                if (cx < dlayer) dcellslist[i].x = -length;
                if (cx > upd)    dcellslist[i].x =  length;
                if (cy < dlayer) dcellslist[i].y = -length;
                if (cy > upd)    dcellslist[i].y =  length;
                if (cz < dlayer) dcellslist[i].z = -length;
                if (cz > upd)    dcellslist[i].z =  length;
            }
        }
    }

    // Real force cells on the faces normal to x (bfc[0]), y (bfc[1]) and z (bfc[2])
    int nx = 0, ny = 0, nz = 0;

    for (cx = flayer; cx <= upf; cx++)
    {
        for (cy = flayer; cy <= upf; cy++)
        {
            for (cz = flayer; cz <= upf; cz++)
            {
                i = cx * tfc2 + cy * tfcell + cz;
                if (cx == flayer || cx == upf)
                {
                    bfc[0].n[nx] = i;
                    nx++;
                }
                if (cy == flayer || cy == upf)
                {
                    bfc[1].n[ny] = i;
                    ny++;
                }
                if (cz == flayer || cz == upf)
                {
                    bfc[2].n[nz] = i;
                    nz++;
                }
            }
        }
    }

    // Spherical stencil: largest |ny| for each nx (y), and largest |nz| for each (nx, |ny|) (z)
    int ymax, zmax, k;

    for (i = 0; i < dlayer + 1; i++)
    {
        j = 0;
        ymax = j;

        while (min((i) * (i), (i - 1) * (i - 1)) + min((j) * (j), (j - 1) * (j - 1)) < dlayer * dlayer)
        {
            ymax = j;
            j++;
        }
        y[i] = ymax;
    }

    for (i = 0; i < dlayer + 1; i++)
    {
        for (j = 0; j < y[i] + 1; j++)
        {
            k = 0;
            zmax = k;
            while (min((i) * (i), (i - 1) * (i - 1)) + min((j) * (j), (j - 1) * (j - 1)) + min((k) * (k), (k - 1) * (k - 1)) < dlayer * dlayer)
            {
                zmax = k;
                k++;
            }
            z[i * (dlayer + 1) + j] = zmax;
        }
    }
}


/* ============================================================================
 *  Forces
 * ========================================================================== */

// Compute the WCA forces on all particles with the cell-list method.
// Each pair is visited once (half-shell stencil) and Newton's third law is
// applied.
void force(Interact *interact, Phase *phase, Fcellslist *fcellslist, Outerfcells *outfcells)
{
    int i, j;
    int cx, cy, cz, index;

    // Reset the linked lists and the forces
    for (i = 0; i < npart; i++)
    {
        phase[i].findex = -1;
    }

    for (i = 0; i < tfc3; i++)
    {
        fcellslist[i].flist = -1;
    }

    // Insert every particle at the head of the list of its cell
    for (i = 0; i < npart; i++)
    {
        interact[i].fx = 0;
        interact[i].fy = 0;
        interact[i].fz = 0;

        cx = ((phase[i].x0 + (length / 2)) / fdcell) + flayer;
        cy = ((phase[i].y0 + (length / 2)) / fdcell) + flayer;
        cz = ((phase[i].z0 + (length / 2)) / fdcell) + flayer;

        index = cx * tfc2 + cy * tfcell + cz;

        if (fcellslist[index].flist == -1)
        {
            fcellslist[index].flist = i;
        }
        else
        {
            phase[i].findex = fcellslist[index].flist;
            fcellslist[index].flist = i;
        }
    }

    // Ghost cells point to the list of the real cell they copy
    for (i = 0; i < ofc; i++)
    {
        fcellslist[outfcells[i].n[0]].flist = fcellslist[outfcells[i].n[1]].flist;
    }

    // Loop over real cells and their half-shell of neighbours
    int diff;
    int nx, ny, nz;
    int p0, p1;
    double rx, ry, rz, d, F;

    for (cx = flayer; cx <= upf; cx++)
    {
        for (cy = flayer; cy <= upf; cy++)
        {
            for (cz = flayer; cz <= upf; cz++)
            {
                i = cx * tfc2 + cy * tfcell + cz;

                for (nx = 0; nx <= flayer; nx++)
                {
                    for (ny = -flayer; ny <= flayer; ny++)
                    {
                        for (nz = -flayer; nz <= flayer; nz++)
                        {
                            j = (cx + nx) * tfc2 + (cy + ny) * tfcell + cz + nz;

                            if (j >= i)
                            {
                                if (j == i)
                                    diff = 0;
                                else
                                    diff = 1;

                                switch (diff)
                                {
                                case 0:  // same cell: each pair once (p1 > p0)
                                    p0 = fcellslist[i].flist;
                                    while (p0 != -1)
                                    {
                                        p1 = fcellslist[j].flist;
                                        while (p1 > p0 && p1 != -1)
                                        {
                                            rx = (phase[p1].x0) - (phase[p0].x0);
                                            ry = (phase[p1].y0) - (phase[p0].y0);
                                            rz = (phase[p1].z0) - (phase[p0].z0);

                                            d = (rx * rx) + (ry * ry) + (rz * rz);

                                            if (d < rc2)
                                            {
                                                d = sigma2 / d;

                                                F = prefact * (2 * d * d * d * d * d * d * d - d * d * d * d);
                                                interact[p1].fx += F * rx;
                                                interact[p1].fy += F * ry;
                                                interact[p1].fz += F * rz;

                                                interact[p0].fx -= F * rx;
                                                interact[p0].fy -= F * ry;
                                                interact[p0].fz -= F * rz;
                                            }

                                            p1 = phase[p1].findex;
                                        }
                                        p0 = phase[p0].findex;
                                    }
                                    break;

                                case 1:  // neighbour cell: all pairs, with the periodic shift of cell j
                                    p0 = fcellslist[i].flist;
                                    while (p0 != -1)
                                    {
                                        p1 = fcellslist[j].flist;
                                        while (p1 != -1)
                                        {
                                            rx = (phase[p1].x0 + fcellslist[j].x) - (phase[p0].x0);
                                            ry = (phase[p1].y0 + fcellslist[j].y) - (phase[p0].y0);
                                            rz = (phase[p1].z0 + fcellslist[j].z) - (phase[p0].z0);

                                            d = (rx * rx) + (ry * ry) + (rz * rz);

                                            if (d < rc2)
                                            {
                                                d = sigma2 / d;

                                                F = prefact * (2 * d * d * d * d * d * d * d - d * d * d * d);
                                                interact[p1].fx += F * rx;
                                                interact[p1].fy += F * ry;
                                                interact[p1].fz += F * rz;

                                                interact[p0].fx -= F * rx;
                                                interact[p0].fy -= F * ry;
                                                interact[p0].fz -= F * rz;
                                            }

                                            p1 = phase[p1].findex;
                                        }
                                        p0 = phase[p0].findex;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


/* ============================================================================
 *  Integrator
 * ========================================================================== */

// Advance the system by one time step with the 5th-order Gear
// predictor-corrector and the Gaussian thermostat.
// Periodic boundary conditions are applied only to particles in the border
// cells, because only those can leave the box in one step.
void gear(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells)
{
    int i;
    double xcor, ycor, zcor, pxcor, pycor, pzcor;
    double zeta, zden, znum;

    // Predictor
    for (i = 0; i < npart; i++)
    {
        phase[i].x0 += phase[i].x1 + phase[i].x2 + phase[i].x3 + phase[i].x4;
        phase[i].y0 += phase[i].y1 + phase[i].y2 + phase[i].y3 + phase[i].y4;
        phase[i].z0 += phase[i].z1 + phase[i].z2 + phase[i].z3 + phase[i].z4;

        phase[i].x1 += 2.0 * phase[i].x2 + 3.0 * phase[i].x3 + 4.0 * phase[i].x4;
        phase[i].y1 += 2.0 * phase[i].y2 + 3.0 * phase[i].y3 + 4.0 * phase[i].y4;
        phase[i].z1 += 2.0 * phase[i].z2 + 3.0 * phase[i].z3 + 4.0 * phase[i].z4;
        phase[i].x2 += 3.0 * phase[i].x3 + 6.0 * phase[i].x4;
        phase[i].y2 += 3.0 * phase[i].y3 + 6.0 * phase[i].y4;
        phase[i].z2 += 3.0 * phase[i].z3 + 6.0 * phase[i].z4;
        phase[i].x3 += 4.0 * phase[i].x4;
        phase[i].y3 += 4.0 * phase[i].y4;
        phase[i].z3 += 4.0 * phase[i].z4;

        phase[i].px0 += phase[i].px1 + phase[i].px2 + phase[i].px3 + phase[i].px4;
        phase[i].py0 += phase[i].py1 + phase[i].py2 + phase[i].py3 + phase[i].py4;
        phase[i].pz0 += phase[i].pz1 + phase[i].pz2 + phase[i].pz3 + phase[i].pz4;

        phase[i].px1 += 2.0 * phase[i].px2 + 3.0 * phase[i].px3 + 4.0 * phase[i].px4;
        phase[i].py1 += 2.0 * phase[i].py2 + 3.0 * phase[i].py3 + 4.0 * phase[i].py4;
        phase[i].pz1 += 2.0 * phase[i].pz2 + 3.0 * phase[i].pz3 + 4.0 * phase[i].pz4;
        phase[i].px2 += 3.0 * phase[i].px3 + 6.0 * phase[i].px4;
        phase[i].py2 += 3.0 * phase[i].py3 + 6.0 * phase[i].py4;
        phase[i].pz2 += 3.0 * phase[i].pz3 + 6.0 * phase[i].pz4;
        phase[i].px3 += 4.0 * phase[i].px4;
        phase[i].py3 += 4.0 * phase[i].py4;
        phase[i].pz3 += 4.0 * phase[i].pz4;
    }

    // Periodic boundary conditions on the border cells
    int p1;
    for (i = 0; i < 2 * nfc2; i++)
    {
        p1 = fcellslist[bfc[0].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].x0);
            p1 = phase[p1].findex;
        }
        p1 = fcellslist[bfc[1].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].y0);
            p1 = phase[p1].findex;
        }
        p1 = fcellslist[bfc[2].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].z0);
            p1 = phase[p1].findex;
        }
    }

    // Forces at the predicted positions
    force(interact, phase, fcellslist, outfcells);

    // Gaussian thermostat multiplier (with proportional feedback)
    zeta = 0;
    znum = 0;
    zden = 0;
    if (constrain == 1)
    {
        for (i = 0; i < npart; i++)
        {
            znum += (phase[i].px0 * interact[i].fx) + (phase[i].py0 * interact[i].fy) + (phase[i].pz0 * interact[i].fz);
            zden += ((phase[i].px0 * phase[i].px0) + (phase[i].py0 * phase[i].py0) + (phase[i].pz0 * phase[i].pz0));
        }

        zeta = znum / zden;
        zeta += fback * ((zden - (3. * npart - freed) * temp) / ((3. * npart - freed) * temp));
    }

    // Corrector
    for (i = 0; i < npart; i++)
    {
        xcor = phase[i].x1 - phase[i].px0 * dt;
        ycor = phase[i].y1 - phase[i].py0 * dt;
        zcor = phase[i].z1 - phase[i].pz0 * dt;

        pxcor = phase[i].px1 - (interact[i].fx - (zeta * phase[i].px0)) * dt;
        pycor = phase[i].py1 - (interact[i].fy - (zeta * phase[i].py0)) * dt;
        pzcor = phase[i].pz1 - (interact[i].fz - (zeta * phase[i].pz0)) * dt;

        phase[i].x0 -= xcor * fo1;
        phase[i].y0 -= ycor * fo1;
        phase[i].z0 -= zcor * fo1;
        phase[i].x1 -= xcor;
        phase[i].y1 -= ycor;
        phase[i].z1 -= zcor;
        phase[i].x2 -= xcor * f21;
        phase[i].y2 -= ycor * f21;
        phase[i].z2 -= zcor * f21;
        phase[i].x3 -= xcor * f31;
        phase[i].y3 -= ycor * f31;
        phase[i].z3 -= zcor * f31;
        phase[i].x4 -= xcor * f41;
        phase[i].y4 -= ycor * f41;
        phase[i].z4 -= zcor * f41;

        phase[i].px0 -= pxcor * fo1;
        phase[i].py0 -= pycor * fo1;
        phase[i].pz0 -= pzcor * fo1;
        phase[i].px1 -= pxcor;
        phase[i].py1 -= pycor;
        phase[i].pz1 -= pzcor;
        phase[i].px2 -= pxcor * f21;
        phase[i].py2 -= pycor * f21;
        phase[i].pz2 -= pzcor * f21;
        phase[i].px3 -= pxcor * f31;
        phase[i].py3 -= pycor * f31;
        phase[i].pz3 -= pzcor * f31;
        phase[i].px4 -= pxcor * f41;
        phase[i].py4 -= pycor * f41;
        phase[i].pz4 -= pzcor * f41;
    }

    // Periodic boundary conditions on the border cells
    for (i = 0; i < 2 * nfc2; i++)
    {
        p1 = fcellslist[bfc[0].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].x0);
            p1 = phase[p1].findex;
        }
        p1 = fcellslist[bfc[1].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].y0);
            p1 = phase[p1].findex;
        }
        p1 = fcellslist[bfc[2].n[i]].flist;
        while (p1 != -1)
        {
            pbcs(&phase[p1].z0);
            p1 = phase[p1].findex;
        }
    }
}


/* ============================================================================
 *  Distribution functions
 * ========================================================================== */

// Build the linked cell lists used by Pair_function and Triplet_function.
// Ghost cells point to the list of the real cell they copy.
void Distrib_neighbourlist(Phase *phase, Dcellslist *dcellslist)
{
    int i, ix, iy, iz;
    int cx, cy, cz, index;

    for (i = 0; i < npart; i++)
    {
        phase[i].dindex = -1;
    }

    for (i = 0; i < tdc3; i++)
    {
        dcellslist[i].dlist = -1;
    }

    // Insert every particle at the head of the list of its cell
    for (i = 0; i < npart; i++)
    {
        cx = ((phase[i].x0 + (length / 2)) / ddcell) + dlayer;
        cy = ((phase[i].y0 + (length / 2)) / ddcell) + dlayer;
        cz = ((phase[i].z0 + (length / 2)) / ddcell) + dlayer;

        index = cx * tdc2 + cy * tdcell + cz;

        if (dcellslist[index].dlist == -1)
        {
            dcellslist[index].dlist = i;
        }
        else
        {
            phase[i].dindex = dcellslist[index].dlist;
            dcellslist[index].dlist = i;
        }
    }

    // Copy the list heads into the ghost cells
    for (cx = 0; cx < tdcell; cx++)
    {
        for (cy = 0; cy < tdcell; cy++)
        {
            for (cz = 0; cz < tdcell; cz++)
            {
                ix = cx;
                iy = cy;
                iz = cz;
                i = ix * tdc2 + iy * tdcell + iz;
                ddpbcs(&ix);
                ddpbcs(&iy);
                ddpbcs(&iz);
                index = ix * tdc2 + iy * tdcell + iz;
                dcellslist[i].dlist = dcellslist[index].dlist;
            }
        }
    }
}

// Histogram all triplets (p0, p1, p2) taken from cells (i, index1, index2) whose
// three side lengths are all below dthreshold.
// check = 10*d21 + d10 encodes which cells coincide:
//    0 : i == index1 == index2   -> p0 < p1 < p2
//   10 : i == index1 != index2   -> p0 < p1, any p2
//    1 : i != index1 == index2   -> any p1, p1 < p2
//   11 : all three cells differ  -> all combinations
void selectparticles(Phase *phase, Dcellslist *dcellslist, int *i, int *index1, int *index2, long long int *triad, int *check)
{
    int p0, p1, p2;
    double rx1, ry1, rz1;
    double rx2, ry2, rz2;
    double rx3, ry3, rz3;
    double d1, d2, d3;
    double x0, y0, z0;
    double x1, y1, z1;
    int i1, i2, i3;
    double r1, r2, r3;

    switch (*check)
    {
    case 11:
        p0 = dcellslist[*i].dlist;
        while (p0 != -1)
        {
            x0 = phase[p0].x0;
            y0 = phase[p0].y0;
            z0 = phase[p0].z0;

            p1 = dcellslist[*index1].dlist;
            while (p1 != -1)
            {
                x1 = phase[p1].x0 + dcellslist[*index1].x;
                y1 = phase[p1].y0 + dcellslist[*index1].y;
                z1 = phase[p1].z0 + dcellslist[*index1].z;

                rx1 = x1 - x0;
                ry1 = y1 - y0;
                rz1 = z1 - z0;
                r1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                if (r1 < dt2)
                {
                    p2 = dcellslist[*index2].dlist;
                    while (p2 != -1)
                    {
                        rx2 = (phase[p2].x0 + dcellslist[*index2].x) - x0;
                        ry2 = (phase[p2].y0 + dcellslist[*index2].y) - y0;
                        rz2 = (phase[p2].z0 + dcellslist[*index2].z) - z0;

                        r2 = (rx2 * rx2) + (ry2 * ry2) + (rz2 * rz2);

                        if (r2 < dt2)
                        {
                            rx3 = rx1 - rx2;
                            ry3 = ry1 - ry2;
                            rz3 = rz1 - rz2;

                            r3 = (rx3 * rx3) + (ry3 * ry3) + (rz3 * rz3);

                            if (r3 < dt2)
                            {
                                d1 = sqrt(r1);
                                d2 = sqrt(r2);
                                d3 = sqrt(r3);

                                // sort so that d1 <= d2 <= d3
                                if (d1 < d2)
                                {
                                    if (d3 < d1) swap(d1, d3);
                                }
                                else
                                {
                                    if (d2 < d3) swap(d1, d2);
                                    else swap(d1, d3);
                                }
                                if (d3 < d2) swap(d2, d3);

                                i1 = d1 / tsbin;
                                i2 = d2 / tsbin;
                                i3 = d3 / tsbin;

                                triad[((i3) * (i3 + 1) * (i3 + 2) / 6) + ((i2 + 1) * i2 / 2) + i1]++;
                            }
                        }
                        p2 = phase[p2].dindex;
                    }
                }
                p1 = phase[p1].dindex;
            }
            p0 = phase[p0].dindex;
        }
        break;

    case 10:
        p0 = dcellslist[*i].dlist;
        while (p0 != -1)
        {
            x0 = phase[p0].x0;
            y0 = phase[p0].y0;
            z0 = phase[p0].z0;

            p1 = dcellslist[*index1].dlist;
            while (p1 > p0 && p1 != -1)
            {
                x1 = phase[p1].x0 + dcellslist[*index1].x;
                y1 = phase[p1].y0 + dcellslist[*index1].y;
                z1 = phase[p1].z0 + dcellslist[*index1].z;

                rx1 = x1 - x0;
                ry1 = y1 - y0;
                rz1 = z1 - z0;
                r1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                if (r1 < dt2)
                {
                    p2 = dcellslist[*index2].dlist;
                    while (p2 != -1)
                    {
                        rx2 = (phase[p2].x0 + dcellslist[*index2].x) - x0;
                        ry2 = (phase[p2].y0 + dcellslist[*index2].y) - y0;
                        rz2 = (phase[p2].z0 + dcellslist[*index2].z) - z0;

                        r2 = (rx2 * rx2) + (ry2 * ry2) + (rz2 * rz2);

                        if (r2 < dt2)
                        {
                            rx3 = rx1 - rx2;
                            ry3 = ry1 - ry2;
                            rz3 = rz1 - rz2;

                            r3 = (rx3 * rx3) + (ry3 * ry3) + (rz3 * rz3);

                            if (r3 < dt2)
                            {
                                d1 = sqrt(r1);
                                d2 = sqrt(r2);
                                d3 = sqrt(r3);

                                // sort so that d1 <= d2 <= d3
                                if (d1 < d2)
                                {
                                    if (d3 < d1) swap(d1, d3);
                                }
                                else
                                {
                                    if (d2 < d3) swap(d1, d2);
                                    else swap(d1, d3);
                                }
                                if (d3 < d2) swap(d2, d3);

                                i1 = d1 / tsbin;
                                i2 = d2 / tsbin;
                                i3 = d3 / tsbin;

                                triad[((i3) * (i3 + 1) * (i3 + 2) / 6) + ((i2 + 1) * i2 / 2) + i1]++;
                            }
                        }
                        p2 = phase[p2].dindex;
                    }
                }
                p1 = phase[p1].dindex;
            }
            p0 = phase[p0].dindex;
        }
        break;

    case 1:
        p0 = dcellslist[*i].dlist;
        while (p0 != -1)
        {
            x0 = phase[p0].x0;
            y0 = phase[p0].y0;
            z0 = phase[p0].z0;

            p1 = dcellslist[*index1].dlist;
            while (p1 != -1)
            {
                x1 = phase[p1].x0 + dcellslist[*index1].x;
                y1 = phase[p1].y0 + dcellslist[*index1].y;
                z1 = phase[p1].z0 + dcellslist[*index1].z;

                rx1 = x1 - x0;
                ry1 = y1 - y0;
                rz1 = z1 - z0;
                r1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                if (r1 < dt2)
                {
                    p2 = dcellslist[*index2].dlist;
                    while (p2 > p1 && p2 != -1)
                    {
                        rx2 = (phase[p2].x0 + dcellslist[*index2].x) - x0;
                        ry2 = (phase[p2].y0 + dcellslist[*index2].y) - y0;
                        rz2 = (phase[p2].z0 + dcellslist[*index2].z) - z0;

                        r2 = (rx2 * rx2) + (ry2 * ry2) + (rz2 * rz2);

                        if (r2 < dt2)
                        {
                            rx3 = rx1 - rx2;
                            ry3 = ry1 - ry2;
                            rz3 = rz1 - rz2;

                            r3 = (rx3 * rx3) + (ry3 * ry3) + (rz3 * rz3);

                            if (r3 < dt2)
                            {
                                d1 = sqrt(r1);
                                d2 = sqrt(r2);
                                d3 = sqrt(r3);

                                // sort so that d1 <= d2 <= d3
                                if (d1 < d2)
                                {
                                    if (d3 < d1) swap(d1, d3);
                                }
                                else
                                {
                                    if (d2 < d3) swap(d1, d2);
                                    else swap(d1, d3);
                                }
                                if (d3 < d2) swap(d2, d3);

                                i1 = d1 / tsbin;
                                i2 = d2 / tsbin;
                                i3 = d3 / tsbin;

                                triad[((i3) * (i3 + 1) * (i3 + 2) / 6) + ((i2 + 1) * i2 / 2) + i1]++;
                            }
                        }
                        p2 = phase[p2].dindex;
                    }
                }
                p1 = phase[p1].dindex;
            }
            p0 = phase[p0].dindex;
        }
        break;

    case 0:
        p0 = dcellslist[*i].dlist;
        while (p0 != -1)
        {
            x0 = phase[p0].x0;
            y0 = phase[p0].y0;
            z0 = phase[p0].z0;

            p1 = dcellslist[*index1].dlist;
            while (p1 > p0 && p1 != -1)
            {
                x1 = phase[p1].x0 + dcellslist[*index1].x;
                y1 = phase[p1].y0 + dcellslist[*index1].y;
                z1 = phase[p1].z0 + dcellslist[*index1].z;

                rx1 = x1 - x0;
                ry1 = y1 - y0;
                rz1 = z1 - z0;
                r1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                if (r1 < dt2)
                {
                    p2 = dcellslist[*index2].dlist;
                    while (p2 > p1 && p2 != -1)
                    {
                        rx2 = (phase[p2].x0 + dcellslist[*index2].x) - x0;
                        ry2 = (phase[p2].y0 + dcellslist[*index2].y) - y0;
                        rz2 = (phase[p2].z0 + dcellslist[*index2].z) - z0;

                        r2 = (rx2 * rx2) + (ry2 * ry2) + (rz2 * rz2);

                        if (r2 < dt2)
                        {
                            rx3 = rx1 - rx2;
                            ry3 = ry1 - ry2;
                            rz3 = rz1 - rz2;

                            r3 = (rx3 * rx3) + (ry3 * ry3) + (rz3 * rz3);

                            if (r3 < dt2)
                            {
                                d1 = sqrt(r1);
                                d2 = sqrt(r2);
                                d3 = sqrt(r3);

                                // sort so that d1 <= d2 <= d3
                                if (d1 < d2)
                                {
                                    if (d3 < d1) swap(d1, d3);
                                }
                                else
                                {
                                    if (d2 < d3) swap(d1, d2);
                                    else swap(d1, d3);
                                }
                                if (d3 < d2) swap(d2, d3);

                                i1 = d1 / tsbin;
                                i2 = d2 / tsbin;
                                i3 = d3 / tsbin;

                                triad[((i3) * (i3 + 1) * (i3 + 2) / 6) + ((i2 + 1) * i2 / 2) + i1]++;
                            }
                        }
                        p2 = phase[p2].dindex;
                    }
                }
                p1 = phase[p1].dindex;
            }
            p0 = phase[p0].dindex;
        }
        break;
    }
}

// Add one sample to the triplet histogram (g3).
// Every unordered triplet of cells is visited once: for each real cell i,
// index1 runs over the half-shell (index1 >= i) of a spherically truncated
// stencil, and index2 runs over the cells with index2 >= index1 that lie within
// the stencil of both i and index1.
void Triplet_function(Phase *phase, Dcellslist *dcellslist, long long int *triad, int *ymax, int *zmax)
{
    int i;
    int nx1, ny1, nz1, nx2, ny2, nz2;
    int x1, y1, z1, x2, y2, z2;
    int cx, cy, cz;
    int index1, index2;
    int d10, d21;

    int indx1, indy1, indz1;
    int indx2, indy2, indz2;
    int diffx, diffy, diffz;

    Distrib_neighbourlist(phase, dcellslist);

    int check = 0;

    for (cx = dlayer; cx < tdcell - dlayer; cx++)
    {
        for (cy = dlayer; cy < tdcell - dlayer; cy++)
        {
            for (cz = dlayer; cz < tdcell - dlayer; cz++)
            {
                i = cx * tdc2 + cy * tdcell + cz;

                // second cell (index1): half-shell of the stencil of i
                for (nx1 = 0; nx1 <= dlayer; nx1++)
                {
                    indx1 = abs(nx1);
                    for (ny1 = -ymax[indx1]; ny1 <= ymax[indx1]; ny1++)
                    {
                        indy1 = abs(ny1);
                        indz1 = (dlayer + 1) * indx1 + indy1;
                        for (nz1 = -zmax[indz1]; nz1 <= zmax[indz1]; nz1++)
                        {
                            x1 = cx + nx1;
                            y1 = cy + ny1;
                            z1 = cz + nz1;

                            index1 = x1 * tdc2 + y1 * tdcell + z1;
                            if (index1 >= i)
                            {
                                if (index1 != i)
                                    d10 = 1;
                                else
                                    d10 = 0;

                                // third cell (index2): inside the stencils of both i and index1
                                for (nx2 = max(0, nx1 - dlayer); nx2 <= min(dlayer, nx1 + dlayer); nx2++)
                                {
                                    indx2 = abs(nx2);
                                    diffx = abs(nx1 - nx2);
                                    for (ny2 = max(-ymax[indx2], ny1 - ymax[diffx]); ny2 <= min(ymax[indx2], ny1 + ymax[diffx]); ny2++)
                                    {
                                        indy2 = abs(ny2);
                                        indz2 = (dlayer + 1) * indx2 + indy2;
                                        diffy = abs(ny1 - ny2);
                                        diffz = diffx * (dlayer + 1) + diffy;
                                        for (nz2 = max(-zmax[indz2], nz1 - zmax[diffz]); nz2 <= min(zmax[indz2], nz1 + zmax[diffz]); nz2++)
                                        {
                                            x2 = cx + nx2;
                                            y2 = cy + ny2;
                                            z2 = cz + nz2;

                                            index2 = x2 * tdc2 + y2 * tdcell + z2;
                                            if (index2 >= index1)
                                            {
                                                if (index2 != index1)
                                                    d21 = 1;
                                                else
                                                    d21 = 0;

                                                check = 10 * d21 + d10;
                                                selectparticles(phase, dcellslist, &i, &index1, &index2, triad, &check);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Add one sample to the pair histogram (g2). Each pair within dthreshold is
// counted once and adds 2 to its bin (one for each particle).
void Pair_function(Phase *phase, Dcellslist *dcellslist, long long int *dpair, int *ymax, int *zmax)
{
    int i;
    double rx1, ry1, rz1;
    double x0, y0, z0;
    double d1;

    int p0, p1;

    int nx1, ny1, nz1;
    int x1, y1, z1;

    int cx, cy, cz;
    int index1;

    int d;
    int indx1, indy1, indz1;
    int b;

    Distrib_neighbourlist(phase, dcellslist);

    for (cx = dlayer; cx <= upd; cx++)
    {
        for (cy = dlayer; cy <= upd; cy++)
        {
            for (cz = dlayer; cz <= upd; cz++)
            {
                i = cx * tdc2 + cy * tdcell + cz;

                for (nx1 = 0; nx1 <= dlayer; nx1++)
                {
                    indx1 = abs(nx1);
                    for (ny1 = -ymax[indx1]; ny1 <= ymax[indx1]; ny1++)
                    {
                        indy1 = abs(ny1);
                        indz1 = (dlayer + 1) * indx1 + indy1;
                        for (nz1 = -zmax[indz1]; nz1 <= zmax[indz1]; nz1++)
                        {
                            x1 = cx + nx1;
                            y1 = cy + ny1;
                            z1 = cz + nz1;

                            index1 = x1 * tdc2 + y1 * tdcell + z1;
                            if (index1 >= i)
                            {
                                if (index1 != i)
                                    d = 1;
                                else
                                    d = 0;

                                switch (d)
                                {
                                case 1:  // different cells: all pairs
                                    p0 = dcellslist[i].dlist;
                                    while (p0 != -1)
                                    {
                                        x0 = phase[p0].x0;
                                        y0 = phase[p0].y0;
                                        z0 = phase[p0].z0;

                                        p1 = dcellslist[index1].dlist;
                                        while (p1 != -1)
                                        {
                                            rx1 = phase[p1].x0 + dcellslist[index1].x - x0;
                                            ry1 = phase[p1].y0 + dcellslist[index1].y - y0;
                                            rz1 = phase[p1].z0 + dcellslist[index1].z - z0;

                                            d1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                                            if (d1 < dt2)
                                            {
                                                d1 = sqrt(d1);
                                                b = d1 / psbin;
                                                dpair[b] += 2;
                                            }
                                            p1 = phase[p1].dindex;
                                        }
                                        p0 = phase[p0].dindex;
                                    }
                                    break;

                                case 0:  // same cell: each pair once (p1 > p0)
                                    p0 = dcellslist[i].dlist;
                                    while (p0 != -1)
                                    {
                                        x0 = phase[p0].x0;
                                        y0 = phase[p0].y0;
                                        z0 = phase[p0].z0;

                                        p1 = dcellslist[index1].dlist;
                                        while (p1 > p0 && p1 != -1)
                                        {
                                            rx1 = phase[p1].x0 + dcellslist[index1].x - x0;
                                            ry1 = phase[p1].y0 + dcellslist[index1].y - y0;
                                            rz1 = phase[p1].z0 + dcellslist[index1].z - z0;

                                            d1 = (rx1 * rx1) + (ry1 * ry1) + (rz1 * rz1);

                                            if (d1 < dt2)
                                            {
                                                d1 = sqrt(d1);
                                                b = d1 / psbin;
                                                dpair[b] += 2;
                                            }
                                            p1 = phase[p1].dindex;
                                        }
                                        p0 = phase[p0].dindex;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


/* ============================================================================
 *  Diagnostics, sampling frequency and output
 * ========================================================================== */

// Append the total kinetic energy (sum of p^2) and the three components of the
// total momentum to the check file. Useful to spot bugs.
void check(Phase *phase, ofstream &chck)
{
    int i;
    double ksum = 0;
    double xsum = 0;
    double ysum = 0;
    double zsum = 0;

    for (i = 0; i < npart; i++)
    {
        ksum += phase[i].px0 * phase[i].px0 + phase[i].py0 * phase[i].py0 + phase[i].pz0 * phase[i].pz0;
        xsum += phase[i].px0;
        ysum += phase[i].py0;
        zsum += phase[i].pz0;
    }
    chck << ksum << "," << xsum << "," << ysum << "," << zsum << "," << endl;
}

// Measure the CPU cost of one time step, one g2 sample and one g3 sample, and
// set the sampling intervals so that sampling g2 and g3 each take roughly as
// much time as the dynamics:
//   dfreq[0] = cost(g2 sample) / cost(time step)
//   dfreq[1] = cost(g3 sample) / cost(time step)
void setfrequency(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, Dcellslist *dcellslist, long long int *dpair, long long int *triad, int *maxy, int *maxz, int *dfreq)
{
    double ptime = 0, ttime = 0, gtime = 0;
    int i, j;
    clock_t timegear;
    clock_t timepair;
    clock_t timetriad;

    int outcycle = 10;
    int incycle  = 1000;

    for (i = 0; i < outcycle; i++)
    {
        timegear = clock();
        for (j = 0; j < incycle; j++)
        {
            gear(phase, interact, fcellslist, bfc, outfcells);
        }
        timegear = clock() - timegear;
        gtime += timegear;

        timepair = clock();
        for (j = 0; j < incycle; j++)
        {
            Pair_function(phase, dcellslist, dpair, maxy, maxz);
        }
        timepair = clock() - timepair;
        ptime += timepair;

        timetriad = clock();
        Triplet_function(phase, dcellslist, triad, maxy, maxz);
        timetriad = clock() - timetriad;
        ttime += timetriad;
    }

    gtime /= (outcycle * incycle);
    ptime /= (outcycle * incycle);
    ttime /= outcycle;

    ptime /= gtime;
    ttime /= gtime;

    dfreq[0] = ptime;
    dfreq[1] = ttime;

    if (dfreq[0] == 0)
    {
        dfreq[0] = 1;
    }
}

// Save the current histograms, microstate and run report.
//   Pair.txt   : number of g2 samples, then the pbin g2 counts
//   Triad.txt  : number of g3 samples, then the ttbin g3 counts
//   Phase.txt  : x, y, z, px, py, pz of every particle
//   Report.txt : sample counts, time steps, CPU time and sampling intervals
void tempsave(long long int *dpair, long long int *triad, Phase *phase, int pc, int tc, int c, double timespent, int *dfreq)
{
    // Distribution functions (large files: written with fprintf)
    FILE *pairfile;
    FILE *tripletfile;

    pairfile    = fopen("Pair.txt", "w+");
    tripletfile = fopen("Triad.txt", "w+");

    int i;
    fprintf(pairfile, "%lld ", pc);
    for (i = 0; i < pbin; i++)
    {
        fprintf(pairfile, "%lld ", dpair[i]);
    }

    fprintf(tripletfile, "%lld ", tc);
    for (i = 0; i < ttbin; i++)
    {
        fprintf(tripletfile, "%lld ", triad[i]);
    }

    fclose(pairfile);
    fclose(tripletfile);

    // Microstate and report (small files: ofstream is fine here, but it is too
    // slow for g2/g3, which use fprintf or fwrite)
    ofstream tempphase("Phase.txt");
    ofstream report("Report.txt");

    for (i = 0; i < npart; i++)
    {
        tempphase << phase[i].x0 << "," << phase[i].y0 << "," << phase[i].z0 << "," << phase[i].px0 << "," << phase[i].py0 << "," << phase[i].pz0 << endl;
    }

    report << "Pair samples=" << pc << ", Triplet samples=" << tc << ", Timesteps=" << c << " , Simulation time (s) = " << timespent << " , pair frequency= " << dfreq[0] << " , triplet frequency = " << dfreq[1] << endl;

    tempphase.close();
    report.close();
}

// Write the simulation parameters to Simulation_information.txt
void savefile(long long int *dpair, long long int *triad, int *samp)
{
    ofstream info("Simulation_information.txt");
    info << "Number of particles = " << npart << " , density = " << dens << " , temperature = " << temp << " , Potential rcut= " << r_cut << endl;
    info << "Threshold = " << dthreshold << ", n. of bins (pair) = " << pbin << endl;
    info << "Threshold = " << dthreshold << ", n. of bins (triplet) = " << tbin << endl;
    info << "Total simulation time (s)= " << simtime << endl;
    info.close();
}


/* ============================================================================
 *  Simulation driver
 * ========================================================================== */

// Run the simulation:
//   1. 10000 equilibration steps to melt the initial lattice
//   2. setfrequency() to choose the g2 and g3 sampling intervals
//   3. reset the histograms
//   4. production run: integrate, sample g2 every dfreq[0] steps and g3 every
//      dfreq[1] steps, save every savetime CPU seconds, stop after simtime
void evolution(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, Dcellslist *dcellslist, long long int *dpair, long long int *triad, int *maxy, int *maxz, int *dfreq, int *samp, ofstream &chck)
{
    int i, j;

    // Equilibration: let the system decorrelate from the initial lattice
    force(interact, phase, fcellslist, outfcells);

    for (i = 0; i < 10000; i++)
    {
        gear(phase, interact, fcellslist, bfc, outfcells);
    }

    // Choose the sampling intervals for g2 and g3
    setfrequency(phase, interact, fcellslist, bfc, outfcells, dcellslist, dpair, triad, maxy, maxz, dfreq);

    clock_t start_time;
    clock_t time;

    i = 0;
    double timespent;
    int c    = 0;  // time steps
    int pc   = 0;  // g2 samples
    int tc   = 0;  // g3 samples
    int step = 1;  // index of the next save
    start_time = clock();

    // Reset the histograms (they were filled during setfrequency)
    for (j = 0; j < pbin; j++)
    {
        dpair[j] = 0;
    }
    for (j = 0; j < ttbin; j++)
    {
        triad[j] = 0;
    }

    // Production run
    while (i == 0)
    {
        gear(phase, interact, fcellslist, bfc, outfcells);

        if (c % dfreq[0] == 0)
        {
            Pair_function(phase, dcellslist, dpair, maxy, maxz);
            pc++;
        }
        if (c % dfreq[1] == 0)
        {
            Triplet_function(phase, dcellslist, triad, maxy, maxz);
            tc++;
        }

        c++;

        time = clock();
        timespent = (double)((time - start_time) / CLOCKS_PER_SEC);

        if (timespent >= savetime * step)
        {
            step++;
            tempsave(dpair, triad, phase, pc, tc, c, timespent, dfreq);
            check(phase, chck);
        }

        if (timespent > simtime)
        {
            i = 1;
        }
    }

    tempsave(dpair, triad, phase, pc, tc, c, timespent, dfreq);
}

// Set up the initial state and the static cell-list tables:
//   - setLattice    : particles on a bcc lattice
//   - randomMomenta : random momenta, zero total momentum, target temperature
//   - Bordercells   : periodic shifts, ghost cells, border cells and the
//                     distribution-function stencil. Periodic boundary
//                     conditions are then applied only to particles in border
//                     cells, which keeps their cost low.
void initialize(Phase *phase, Dcellslist *dcellslist, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, int *maxy, int *maxz)
{
    setLattice(phase);
    randomMomenta(phase);
    Bordercells(dcellslist, fcellslist, bfc, outfcells, maxy, maxz);
}


/* ============================================================================
 *  Main
 * ========================================================================== */

int main()
{
    Phase        *phase      = new Phase[npart];
    Fcellslist   *fcellslist = new Fcellslist[tfc3];
    Dcellslist   *dcellslist = new Dcellslist[tdc3];
    Borderfcells *bfc        = new Borderfcells[3];
    Interact     *interact   = new Interact[npart];
    Outerfcells  *outfcells  = new Outerfcells[ofc];

    int *maxy  = new int[dlayer + 1];
    int *maxz  = new int[(dlayer + 1) * (dlayer + 1)];
    int *dfreq = new int[2];
    int *samp  = new int[2];

    long long int *dpair = new long long int[pbin];
    long long int *triad = new long long int[ttbin];

    if (triad == nullptr)
    {
        cout << "fail" << endl;
    }
    else
    {
        cout << "success" << endl;
    }

    ofstream checkevol("checkevolution.txt");
    savefile(dpair, triad, samp);

    // Initial positions and momenta, and the cell-list tables
    initialize(phase, dcellslist, fcellslist, bfc, outfcells, maxy, maxz);

    evolution(phase, interact, fcellslist, bfc, outfcells, dcellslist, dpair, triad, maxy, maxz, dfreq, samp, checkevol);

    checkevol.close();

    delete[] phase;
    delete[] fcellslist;
    delete[] dcellslist;
    delete[] bfc;
    delete[] dpair;
    delete[] triad;
    return 0;
}
