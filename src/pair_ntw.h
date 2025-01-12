#ifdef PAIR_CLASS
// clang-format off
PairStyle(ntw,PairNTW);
// clang-format on
#else

#ifndef PAIR_NTW_H
#define PAIR_NTW_H

#include "pair.h"

namespace LAMMPS_NS {

// Hello

class PairNTW : public Pair {
 public:
  PairNTW(class LAMMPS *);
  ~PairNTW() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  //void init_style() override;
  double init_one(int, int) override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;

  double single(int, int, int, int, double, double, double, double &) override;
  //void born_matrix(int, int, int, int, double, double, double, double &, double &) override;
  void *extract(const char *, int &) override;

  //void compute_inner() override;
  //void compute_middle() override;
  //void compute_outer(int, int) override;

 protected:
  double gamma;
  double **cut;
  double **rc, **rcsq;
  double **epsilon, **sigma;
  double **A, **B, **C, **Asq;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  /* double *cut_respa; */

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif    // PAIR_NTW_H
#endif
