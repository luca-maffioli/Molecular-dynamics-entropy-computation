# Two- and three-particle distribution functions of a WCA fluid

> **Reference.** This code was used to produce the results published in
>
> L. Maffioli, N. Clisby, F. Frascoli and B. D. Todd,
> *Computation of the equilibrium three-particle entropy for dense atomic fluids by
> molecular dynamics simulation*,
> J. Chem. Phys. **151**, 164102 (2019).
> [doi:10.1063/1.5124715](https://doi.org/10.1063/1.5124715)
>

```bibtex
@article{Maffioli2019,
  author  = {Maffioli, Luca and Clisby, Nathan and Frascoli, Federico and Todd, B. D.},
  title   = {Computation of the equilibrium three-particle entropy for dense atomic fluids by molecular dynamics simulation},
  journal = {The Journal of Chemical Physics},
  volume  = {151},
  number  = {16},
  pages   = {164102},
  year    = {2019},
  doi     = {10.1063/1.5124715}
}
```

This program is an equilibrium molecular dynamics (MD) simulation of a three-dimensional
fluid. It samples the **two-particle (g2)** and **three-particle (g3) distribution
functions** of the fluid.

- **System:** particles in a periodic cubic box, interacting through the
  Weeks–Chandler–Andersen (WCA) potential. This is the Lennard-Jones potential cut at
  r<sub>c</sub> = 2<sup>1/6</sup>σ.
- **Integrator:** 5-value Gear predictor–corrector.
- **Thermostat:** Gaussian isokinetic thermostat with a proportional feedback term.
  Isoenergetic dynamics are also available.
- **Neighbour search:** linked cell lists with ghost cells for the periodic boundaries.
  These are used both for the forces and for the distribution functions.
- **Run control:** the program runs for a fixed CPU-time budget and saves the histograms
  periodically, so a long run can be stopped at any time without losing data.

---

## Contents

1. [Building and running](#building-and-running)
2. [Parameters](#parameters)
3. [Program flow](#program-flow)
4. [Cell lists](#cell-lists)
5. [Storage of g3](#storage-of-g3)
6. [Data structures](#data-structures)
7. [Function reference](#function-reference)
8. [Output files](#output-files)
9. [Notes](#notes)

---

## Building and running

```bash
g++ -O3 -march=native main.cpp -o g3md
./g3md
```

The program reads no input: all parameters are compile-time constants at the top of
`main.cpp`. All output files are written to the working directory.

The g3 histogram has about 36 million entries, so the run needs **about 300 MB of
RAM**.

---

## Parameters

### Physical parameters

| Constant    | Default            | Meaning                                                              |
|-------------|--------------------|----------------------------------------------------------------------|
| `nside`     | 15                 | Lattice cells per side of the initial bcc lattice                    |
| `npart`     | 2·nside³ = 6750    | Number of particles                                                  |
| `dens`      | 1                  | Number density                                                       |
| `sigma`     | 1                  | Particle diameter σ                                                  |
| `epsilon`   | 1                  | LJ energy ε                                                          |
| `mass`      | 1                  | Particle mass. Unused: the equations of motion assume m = 1          |
| `temp`      | 1                  | Target temperature                                                   |
| `freed`     | 3                  | Degrees of freedom removed by the thermostat (= dimension)           |
| `fback`     | 5                  | Coefficient of the thermostat feedback term                          |
| `constrain` | 1                  | 1 = isokinetic (thermostatted), 0 = isoenergetic                     |

### Run-time parameters

| Constant   | Default       | Meaning                                         |
|------------|---------------|-------------------------------------------------|
| `simtime`  | 4.5 days (s)  | CPU time of the production run                  |
| `savetime` | 3600 s        | CPU time between two saves                      |
| `dt`       | 0.001         | Integration time step                           |
| `pbin`     | 600           | Number of bins of g2                            |
| `tbin`     | 600           | Number of bins per distance of g3               |

### Derived quantities

The values in brackets are for the default parameters.

| Constant     | Meaning                                                                                     |
|--------------|---------------------------------------------------------------------------------------------|
| `length`     | Box side, (npart/dens)<sup>1/3</sup> [18.90]                                                |
| `r_cut`      | WCA cutoff, 2<sup>1/6</sup>σ [1.122]                                                        |
| `nfcell`     | Force cells per side [16]                                                                   |
| `fdcell`     | Force cell size, ≥ r_cut [1.181]                                                            |
| `flayer`     | Layers of neighbour/ghost force cells [1]                                                   |
| `ndcell`     | Distribution cells per side, (npart/8)<sup>1/3</sup> [9]                                    |
| `ddcell`     | Distribution cell size [2.100]                                                              |
| `dlayer`     | Layers of neighbour/ghost distribution cells [4]                                            |
| `dthreshold` | Largest distance sampled in g2 and g3. It is the smaller of L/2 and `maxlayer`·`ddcell` [8.40] |
| `psbin`, `tsbin` | Bin widths of g2 and g3 [0.0140]                                                        |
| `ttbin`      | Number of stored g3 bins, tbin(tbin+1)(tbin+2)/6 [36 180 200]                               |

---

## Program flow

```
main
 ├─ savefile                    write Simulation_information.txt
 ├─ initialize
 │   ├─ setLattice              bcc lattice
 │   ├─ randomMomenta           random momenta, P = 0, T = temp
 │   └─ Bordercells             static cell-list tables
 └─ evolution
     ├─ force                   initial forces
     ├─ gear × 10000            equilibration (melts the lattice)
     ├─ setfrequency            choose the g2/g3 sampling intervals
     ├─ reset histograms
     └─ loop until simtime
         ├─ gear                one time step
         ├─ Pair_function       every dfreq[0] steps
         ├─ Triplet_function    every dfreq[1] steps
         └─ tempsave + check    every savetime seconds
```

### Equations of motion

With m = 1, the integrated equations are

```
dr_i/dt = p_i
dp_i/dt = F_i − α p_i

α = Σ F_i·p_i / Σ p_i²  +  fback · (Σ p_i² − (3N − freed)·T) / ((3N − freed)·T)
```

The first term of α is the Gaussian multiplier that keeps the kinetic energy constant.
The second term pulls the kinetic energy back to its target value and corrects the
numerical drift. With `constrain = 0`, α = 0 and the dynamics conserve the total energy.

The equations are integrated with a 5-value Gear predictor–corrector for first-order
equations. The corrector coefficients are 251/720, 1, 11/12, 1/3 and 1/24. The
derivatives are stored scaled: `x1 = dt·dx/dt`, `x2 = dt²/2·d²x/dt²`, …

---

## Cell lists

Two independent cell grids are used: one for the forces (cell size ≥ r_cut) and one for
the distribution functions (cell size `ddcell`).

**Ghost layers.** Each grid is padded with `flayer` (or `dlayer`) layers of ghost cells
on every side. A ghost cell is a periodic image of a real cell:

- its list head points to the particle list of that real cell, and
- it stores the shift (0 or ±`length` along each axis) that must be added to the
  positions of its particles.

With ghost cells, the inner loops never have to apply the minimum-image convention.

**Linked lists.** Each cell stores the index of its first particle (`flist`/`dlist`).
Each particle stores the index of the next particle in the same cell
(`findex`/`dindex`), and -1 marks the end of a list. Particles are inserted at the head,
so each list is in decreasing particle index. This ordering is used to visit each pair
inside a cell only once (`p1 > p0`).

**Half-shell traversal.** For each real cell `i`, only the neighbour cells `j ≥ i` are
visited, so each pair of cells is visited once. Newton's third law then gives the force
on both particles.

**Spherical stencil.** For the distribution functions, the neighbour stencil is
truncated to the cells that can contain particles closer than `dthreshold`:

- `ymax[|nx|]` is the largest |ny| for a given |nx|;
- `zmax[|nx|·(dlayer+1) + |ny|]` is the largest |nz| for given |nx| and |ny|.

**Periodic boundary conditions.** In one time step, only particles in the outermost real
cells can leave the box. `gear` therefore applies `pbcs` only to the particles in those
border cells.

---

## Storage of g3

g3(r<sub>1</sub>, r<sub>2</sub>, r<sub>3</sub>) is symmetric under any permutation of
its three arguments, so only the region i ≥ j ≥ k (bin indices) is stored. This region
is stored as a 1D array:

```
g3(i, j, k) = triad[ i(i+1)(i+2)/6 + (j+1)j/2 + k ],    i ≥ j ≥ k
```

For any other ordering, sort the three indices first. The histograms are raw counts:

- g2 adds 2 per pair, one for each particle of the pair;
- g3 adds 1 per unordered triplet.

To convert them into distribution functions, normalise by the number of samples, the
number of particles, the density and the bin volume.

---

## Data structures

| Type           | Fields                                                                                  | Description |
|----------------|------------------------------------------------------------------------------------------|-------------|
| `Phase`        | `x0…x4`, `y0…y4`, `z0…z4`, `px0…px4`, `py0…py4`, `pz0…pz4`, `findex`, `dindex`          | State of one particle: position, momentum and their scaled derivatives up to 4th order (Gear), plus the links to the next particle in the force and distribution cell lists |
| `Interact`     | `fx`, `fy`, `fz`                                                                         | Total force on one particle |
| `Fcellslist`   | `flist`, `x`, `y`, `z`                                                                   | Force cell: head of its particle list and its periodic shift |
| `Dcellslist`   | `dlist`, `x`, `y`, `z`                                                                   | Distribution cell: same as `Fcellslist` |
| `Borderfcells` | `n[2·nfc2]`                                                                              | Indices of the real force cells on the two faces normal to one axis. Three instances are used, for x, y and z |
| `Outerfcells`  | `n[2]`                                                                                   | `{ghost cell, real cell it copies}` |

---

## Function reference

### Periodic boundary conditions

#### `void pbcs(double *x)`
Wraps a coordinate back into [−L/2, L/2) by adding or subtracting `length` once. It
assumes the particle has moved by less than one box length.

#### `void fdpbcs(int *x)`
Maps a force-cell coordinate along one axis, including the ghost layers, onto the
coordinate of the corresponding real cell.

#### `void ddpbcs(int *x)`
Same as `fdpbcs`, for the distribution-function grid.

### Initialisation

#### `void rand0(void)`
Seeds the C random number generator with the current time.

#### `void setLattice(Phase *phase)`
Places the `npart` particles on a body-centred cubic lattice (2 particles per lattice
cell, `nside` cells per side). It then wraps them into the box and sets all position
derivatives `x1…x4`, `y1…y4`, `z1…z4` to zero.

#### `void randomMomenta(Phase *phase)`
Initialises the momenta:

1. draws each momentum component uniformly from [0, 1];
2. subtracts the mean, so the total momentum is zero;
3. rescales all momenta so that Σp² = 3·npart·temp;
4. sets all momentum derivatives `px1…pz4` to zero.

#### `void Bordercells(Dcellslist *dcellslist, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, int *y, int *z)`
Builds the static tables used by the cell lists. It is called once, at start-up.

| Argument     | Filled with |
|--------------|-------------|
| `fcellslist` | Periodic shift (`x`, `y`, `z`) of every force cell: 0 for real cells, ±`length` for ghost cells |
| `outfcells`  | For each of the `ofc` ghost force cells, its index and the index of the real cell it copies |
| `dcellslist` | Periodic shift of every distribution cell |
| `bfc`        | `bfc[0]`, `bfc[1]`, `bfc[2]`: real force cells on the faces normal to x, y and z |
| `y`          | `ymax` of the spherical stencil (size `dlayer+1`) |
| `z`          | `zmax` of the spherical stencil (size `(dlayer+1)²`) |

#### `void initialize(Phase *phase, Dcellslist *dcellslist, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, int *maxy, int *maxz)`
Calls `setLattice`, `randomMomenta` and `Bordercells` in that order.

### Forces and integration

#### `void force(Interact *interact, Phase *phase, Fcellslist *fcellslist, Outerfcells *outfcells)`
Computes the WCA force on every particle:

1. rebuilds the force cell lists from the current positions;
2. copies the real list heads into the ghost cells (`outfcells`);
3. for every real cell, visits the half-shell of neighbour cells:
   - within the same cell, each pair is taken once (`p1 > p0`);
   - in a neighbour cell, all pairs are taken, with the cell's periodic shift added
     to `p1`;
4. for every pair with r² < r<sub>c</sub>², computes

   F = 24ε/σ² · [2(σ²/r²)⁷ − (σ²/r²)⁴]

   and adds F·**r** to `p1` and −F·**r** to `p0`.

The result is written to `interact`.

#### `void gear(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells)`
Advances the system by one time step `dt`:

1. **Predictor:** Taylor-expands positions, momenta and their derivatives (Pascal
   triangle update).
2. **PBC:** applies `pbcs` to the particles in the border cells. The cell membership
   comes from the previous call to `force`.
3. **Forces:** calls `force` at the predicted positions.
4. **Thermostat:** if `constrain == 1`, computes the multiplier α (see
   [Equations of motion](#equations-of-motion)).
5. **Corrector:** computes the errors

   `xcor = x1 − p·dt` and `pcor = p1 − (F − αp)·dt`

   and corrects every derivative with the coefficients `fo1`, 1, `f21`, `f31`, `f41`.
6. **PBC:** applies `pbcs` again to the particles in the border cells.

### Distribution functions

#### `void Distrib_neighbourlist(Phase *phase, Dcellslist *dcellslist)`
Builds the distribution-function cell lists from the current positions. It then makes
every ghost cell point to the list of its real image. It is called at the start of
`Pair_function` and `Triplet_function`.

#### `void Pair_function(Phase *phase, Dcellslist *dcellslist, long long int *dpair, int *ymax, int *zmax)`
Adds one sample to the g2 histogram `dpair`:

1. rebuilds the distribution cell lists;
2. for each real cell `i`, visits the cells `index1 ≥ i` of the spherical stencil;
3. for every pair closer than `dthreshold`, adds 2 to bin `r / psbin`.

Within a cell, each pair is taken once (`p1 > p0`).

#### `void Triplet_function(Phase *phase, Dcellslist *dcellslist, long long int *triad, int *ymax, int *zmax)`
Adds one sample to the g3 histogram `triad`. It rebuilds the distribution cell lists,
then loops over triplets of cells `(i, index1, index2)`:

- `i` is a real cell;
- `index1 ≥ i` is in the half-shell of the stencil of `i`;
- `index2 ≥ index1` is within the stencil of both `i` and `index1`.

This visits each unordered triplet of cells once. For each triplet of cells, it records
which cells coincide (`d10`: `index1 ≠ i`, `d21`: `index2 ≠ index1`) and passes them to
`selectparticles` as `check = 10·d21 + d10`.

#### `void selectparticles(Phase *phase, Dcellslist *dcellslist, int *i, int *index1, int *index2, long long int *triad, int *check)`
Histograms the particle triplets for one triplet of cells. `p0` comes from cell `i`,
`p1` from `index1` and `p2` from `index2`. The value of `check` sets which particle
orderings are used, so that every unordered triplet of particles is counted exactly
once:

| `check` | Cells                        | Particle ordering  |
|---------|------------------------------|--------------------|
| 0       | i = index1 = index2          | p0 < p1 < p2       |
| 10      | i = index1 ≠ index2          | p0 < p1, any p2    |
| 1       | i ≠ index1 = index2          | any p1, p1 < p2    |
| 11      | all different                | all combinations   |

The three side lengths of the triangle are checked in turn against `dthreshold`, and the
loop stops at the first one that fails. For accepted triplets, the three distances are
sorted (d1 ≤ d2 ≤ d3) and binned with `tsbin`. The corresponding entry of `triad`
(see [Storage of g3](#storage-of-g3)) is then incremented.

### Sampling, diagnostics and output

#### `void setfrequency(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, Dcellslist *dcellslist, long long int *dpair, long long int *triad, int *maxy, int *maxz, int *dfreq)`
Measures, with `clock()`, the average CPU cost of:

- one `gear` step (10 × 1000 calls);
- one `Pair_function` call (10 × 1000 calls);
- one `Triplet_function` call (10 calls).

It then sets the sampling intervals

```
dfreq[0] = cost(Pair_function)    / cost(gear)
dfreq[1] = cost(Triplet_function) / cost(gear)
```

so that sampling g2 and sampling g3 each take roughly as much CPU time as the dynamics.
`dfreq[0]` is at least 1. The dynamics keep running during the measurement, and the
histograms filled here are discarded afterwards.

#### `void check(Phase *phase, ofstream &chck)`
Appends one line to `checkevolution.txt`:

```
Σp², ΣPx, ΣPy, ΣPz,
```

Σp² is twice the kinetic energy, and ΣPx, ΣPy, ΣPz are the components of the total
momentum. In the isokinetic case, Σp² should stay at its target value and the total
momentum at zero, which makes this a quick consistency check.

#### `void tempsave(long long int *dpair, long long int *triad, Phase *phase, int pc, int tc, int c, double timespent, int *dfreq)`
Overwrites the output files with the current state of the run: `Pair.txt`, `Triad.txt`,
`Phase.txt` and `Report.txt` (see [Output files](#output-files)).

| Argument    | Meaning                         |
|-------------|---------------------------------|
| `pc`, `tc`  | Number of g2 and g3 samples     |
| `c`         | Number of time steps            |
| `timespent` | CPU time of the production run  |

#### `void savefile(long long int *dpair, long long int *triad, int *samp)`
Writes the main simulation parameters to `Simulation_information.txt`. The arguments
are not used.

#### `void evolution(Phase *phase, Interact *interact, Fcellslist *fcellslist, Borderfcells *bfc, Outerfcells *outfcells, Dcellslist *dcellslist, long long int *dpair, long long int *triad, int *maxy, int *maxz, int *dfreq, int *samp, ofstream &chck)`
Runs the whole simulation:

1. computes the initial forces;
2. runs 10 000 equilibration steps;
3. calls `setfrequency`;
4. resets the histograms;
5. runs the production loop. At every step it calls `gear`, samples g2 every
   `dfreq[0]` steps and g3 every `dfreq[1]` steps, and calls `tempsave` and `check`
   every `savetime` CPU seconds. It stops after `simtime` CPU seconds and calls
   `tempsave` one last time.

#### `int main()`
Allocates all arrays and opens `checkevolution.txt`. It then calls `savefile`,
`initialize` and `evolution`.

---

## Output files

| File                         | Written by                | Content |
|------------------------------|---------------------------|---------|
| `Simulation_information.txt` | `savefile`                | Number of particles, density, temperature, cutoff, `dthreshold`, numbers of bins, `simtime` |
| `Pair.txt`                   | `tempsave`                | One line: the number of g2 samples, then the `pbin` counts of g2, space-separated |
| `Triad.txt`                  | `tempsave`                | One line: the number of g3 samples, then the `ttbin` counts of g3, space-separated (about 300 MB of text) |
| `Phase.txt`                  | `tempsave`                | One line per particle: `x,y,z,px,py,pz` |
| `Report.txt`                 | `tempsave`                | Numbers of samples, time steps, CPU time, sampling intervals |
| `checkevolution.txt`         | `check`                   | One line per save: Σp², total momentum |

The histogram files are overwritten at each save, so they always hold the cumulative
histograms of the run so far.

