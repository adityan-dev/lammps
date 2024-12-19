#include "pair_ntw.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
// #include "neighbor.h"
// #include "respa.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
// using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairNTW::PairNTW(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairNTW::~PairNTW()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rc);
    memory->destroy(rcsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(A);
    memory->destroy(Asq);
    memory->destroy(B);
    memory->destroy(C);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairNTW::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double r, rinv, rsq, r2inv, r6inv, forcelj, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < Asq[itype][jtype]) {
        r = sqrt(rsq);
        rinv = 1 / r;
        r2inv = 1.0 / rsq;
        r6inv = r2inv * r2inv * r2inv;

        if (rsq < rcsq[itype][jtype]) {
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          fpair = factor_lj * forcelj * r2inv;
        } else {
          forcelj = -3.0 * B[itype][jtype] * pow(r - A[itype][jtype], 2.0);
          fpair = factor_lj * forcelj * rinv;
        }

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          if (rsq < rcsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) + C[itype][jtype];
            evdwl *= factor_lj;
          } else {
            evdwl = B[itype][jtype] * pow(r - A[itype][jtype], 3.0);
            evdwl *= factor_lj;
          }
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNTW::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(cut, n, n, "pair:cut");
  memory->create(rc, n, n, "pair:rc");
  memory->create(rcsq, n, n, "pair:rcsq");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(A, n, n, "pair:A");
  memory->create(Asq, n, n, "pair:Asq");
  memory->create(B, n, n, "pair:B");
  memory->create(C, n, n, "pair:C");
  memory->create(lj1, n, n, "pair:lj1");
  memory->create(lj2, n, n, "pair:lj2");
  memory->create(lj3, n, n, "pair:lj3");
  memory->create(lj4, n, n, "pair:lj4");
  memory->create(offset, n, n, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNTW::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  gamma = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  // if (allocated) {
  //   int i, j;
  //   for (i = 1; i <= atom->ntypes; i++)
  //     for (j = i; j <= atom->ntypes; j++)
  //       if (setflag[i][j]) cut[i][j] = gamma * sigma[i][j];
  // }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNTW::coeff(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

  //std::cout << "Gamma = " << gamma << "\n";
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      // std::cout << "Sigma" << i << j << sigma_one << "\n";
      // std::cout << "Epsilon" << i << j << epsilon_one << "\n";
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNTW::init_one(int i, int j)
{

  double A0, B0, C0, g, g3, g6, g12, g15;

  g = gamma;
  g3 = g * g * g;
  g6 = g3 * g3;
  g12 = g6 * g6;
  g15 = g12 * g3;

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  if (i == j) {
    lj1[i][j] = 12.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj2[i][j] = 0.0;
    lj3[i][j] = epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj4[i][j] = 0.0;
  } else {
    lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
    lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  }

  A0 = sigma[i][j];
  C0 = epsilon[i][j];
  B0 = C0 / (A0 * A0 * A0);

  if (i == j) {
    A[i][j] = (15.0 / 13.0) * g * A0;
    B[i][j] = (-169.0 / g15) * B0;
    C[i][j] = (-5.0 / (13.0 * g12)) * C0;
  } else {
    A[i][j] = 3.0 * ((3.0 * g6 - 10.0) / (7.0 * g6 - 26.0)) * g * A0;
    B[i][j] = (2.0 / g15) * ((7.0 * g6 - 26.0) * (7.0 * g6 - 26.0) / (g6 - 2.0)) * B0;
    C[i][j] = (4.0 / g12) * ((3.0 * g6 - 2.0) * (g6 - 5.0) / (7.0 * g6 - 26.0)) * C0;
  }

  // std::cout << i << j << " " << A[i][j] << " " << B[i][j] << " " << C[i][j] << "\n";

  cut[i][j] = A[i][j];
  Asq[i][j] = A[i][j] * A[i][j];

  rc[i][j] = g * sigma[i][j];
  rcsq[i][j] = rc[i][j] * rc[i][j];

  // std::cout << i << j << " " << Asq[i][j] << " " << A[i][j] << "\n";

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  A[j][i] = A[i][j];
  Asq[j][i] = Asq[i][j];
  B[j][i] = B[i][j];
  C[j][i] = C[i][j];
  rc[j][i] = rc[i][j];
  rcsq[j][i] = rcsq[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNTW::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNTW::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNTW::write_restart_settings(FILE *fp)
{
  fwrite(&gamma, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNTW::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) { utils::sfread(FLERR, &gamma, sizeof(double), 1, fp, nullptr, error); }
  MPI_Bcast(&gamma, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairNTW::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairNTW::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j], A[i][j],
              B[i][j], C[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairNTW::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r, rinv, r2inv, r6inv, forcelj, philj;

  r = sqrt(rsq);
  rinv = 1.0 / r;
  r2inv = 1.0 / rsq;
  r6inv = r2inv * r2inv * r2inv;

  if (rsq < rcsq[itype][jtype]) {
    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
    fforce = factor_lj * forcelj * r2inv;
    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) + C[itype][jtype];
  } else if (rsq < Asq[itype][jtype]) {
    forcelj = -3.0 * B[itype][jtype] * pow(r - A[itype][jtype], 2.0);
    fforce = factor_lj * forcelj * rinv;
    //std::cout << itype << jtype << " " << factor_lj << "\n";
    philj = B[itype][jtype] * pow(r - A[itype][jtype], 3.0);
  } else {
    forcelj = 0.0;
    philj = 0.0;
  }
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void *PairNTW::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str, "switch_a") == 0) return (void *) A;
  if (strcmp(str, "switch_b") == 0) return (void *) B;
  if (strcmp(str, "switch_c") == 0) return (void *) C;
  if (strcmp(str, "cutoff") == 0) return (void *) rc;
  return nullptr;
}
