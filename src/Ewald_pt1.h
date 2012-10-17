/******

   James McClain and Sebastian Wouters
   Doing the Ewald stuff for Waffle 0.1
   Hack-a-thon Oct 14 to Oct 16, 2012

*******/

#ifndef EWALD_PT1_H
#define EWALD_PT1_H

#include "PhfSolidDef.h"

class Ewald_pt1{

   public:
      //Ewald::Ewald takes in the parameter eta (see docs/EwaldTDR.tex) and info about the supercell.
      Ewald_pt1(double eta, FSolidModel & Solid);

      virtual ~Ewald_pt1();

      //Ewald::NN calculates and returns the nuclear-nuclear repulsion energy W per supercell (see docs/EwaldTDR.tex).
      double NN();

      //Ewald::eN takes in info of two gaussian primitives (ζa, lmax_A, vec(A), ζb , lmax_B , vec(B)) [cart co!]
      //With count defined as:
      //   int count = -1;
      //   for (int nx = lmax; nx >=0; nx--)
      //      for (int ny = lmax - nx; ny >=0; ny--)
      //         int nz = lmax - nx - ny;
      //         count++;
      //count_A is related to the unnormalized angular momentum basis function (x - A_x)^nx * (y - A_y)^ny * (z - A_z)^nz * exp(-eta * (vec{r} - vec{A})^2)
      //There are nTotal_A = (lmax+1)*(lmax+2)/2 such angular momentum basis functions. Result should hence be of size nTotal_A*nTotal_B (can be larger too)
      //The last term of Eq. (15) of docs/EwaldTDR.tex, including the -4 pi / Omega, sandwiched between unnormalized ang. mom. basisfunctions (count_A,count_B) is stored in result[count_A+nTotal_A*count_B]. Note that there should already be memory allocated for result!!
      void eN(double * result, double zeta_A, int lmax_A, double Ax, double Ay, double Az, double zeta_B, int lmax_B, double Bx, double By, double Bz);

   private:
      
      //! Number of nuclei in the supercell
      int nNuclei;

      //! List with charges of the nuclei [size: nNuclei]
      int * nZvals;

      //! Array with locations of nuclei [size: nNuclei*3]: nucleus i with charge nZvals[i] is located at (\vec{r}_i)_{j-co in cart co} = dZloc[i + nNuclei*j]
      double * dZloc;

      //! The supercell translation vectors defined in a right-handed cartesian coordinate system [size: 3*3]: (\vec{T}_i)_{j-coordinate in cart co} = dTvecs[i + 3*j]
      double * dTvecs;

      //! The supercell reciprocal vectors (right-handed cart co) [size: 3*3]: ( \vec{G}_i )_{j-co in cart co} = dGvecs[i + 3*j]
      double * dGvecs;

      //! The supercell volume
      double Volume;

      //! eta parameter of Ewald summation
      double eta;

      //! threshold for exp(-G^2/eta)
      double epsilon;

      //! Truncation parameter for the summation over the reciprocal lattice vectors in the LR Ewald part.
      int maxN;

      //Calculate the signed volume (the triple product without taking the absolute value) of the supercell, the dTvecs should be initialized though
      double CalcSignedVolume();

      //Calculate the reciprocal lattice vectors of the supercell in cart co, dTvecs and Volume should be OK then and dGvecs initialized.
      void CalcReciprocal();

      //if useG==true:
      //   Based on the metric of the skew reciprocal basis, an upper bound is determined for n1,n2,n3 (continuous veriables, |ni|<=upper bound) so that in the set (n1G1 + n2G2 + n3G3) certainly all vectors within the unit sphere are contained. Important for doing the summation over all G vectors in Ewald::eN.
      //else:
      //   Same thing but for the skew supercell translation vectors.
      double giveNmaxOverR(bool useG=true);

      //Work up the ladder for the factorized nested sine/cosine recursion relations
      void FillCosSinOneElectron(double * element_cos, double * element_sin, int lmax, double seed, double zeta, double Aco, double Bco, double Pco, double Gco);

};

#endif

