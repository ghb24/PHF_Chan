/******

   James McClain and Sebastian Wouters
   Doing the Ewald stuff for Waffle 0.1
   Hack-a-thon Oct 14 to Oct 16, 2012

*******/

#include <math.h>
#include "Ewald_pt1.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "Ewald_lapack.h"

Ewald_pt1::Ewald_pt1(double etaval){ //FIXME

   //FIXME: the following block should actually be passed/ filled in from passed data when Ewald is created.
   //For now: do NaCl lattice
   this->nNuclei = 2;
   this->nZvals = new int[this->nNuclei];
   this->nZvals[0] = 1;
   this->nZvals[1] = 1;
   this->dTvecs = new double[9];
   double a_pm = 4.0;
   this->dTvecs[0 + 3*0] = a_pm/2; //vec1 (c++:0) xco (c++:0)
   this->dTvecs[0 + 3*1] = a_pm/2; //vec1 (c++:0) yco (c++:1)
   this->dTvecs[0 + 3*2] = 0.0;
   this->dTvecs[1 + 3*0] = 0.0;
   this->dTvecs[1 + 3*1] = a_pm/2;
   this->dTvecs[1 + 3*2] = a_pm/2;
   this->dTvecs[2 + 3*0] = a_pm/2;
   this->dTvecs[2 + 3*1] = 0.0;
   this->dTvecs[2 + 3*2] = a_pm/2;
   this->dZloc = new double[this->nNuclei * 3];
   this->dZloc[0 + nNuclei*0] = this->dZloc[0 + 3*1] = this->dZloc[0 + 3*2] = 0.0;
   this->dZloc[1 + nNuclei*0] = a_pm/2;
   this->dZloc[1 + nNuclei*1] = this->dZloc[1 + 3*2] = 0.0;

   //This part is independent from any kind of input. It only depends on whether Tvecs is filled in correctly. Not that we use the Gi.Tj = 2 pi delta_ij convention.
   this->Volume = CalcSignedVolume();
   if ((this->Volume)<0.0)
      this->Volume = -(this->Volume);
   this->dGvecs = new double[9];
   CalcReciprocal();

   //Store eta and
   this->eta = etaval;
   this->epsilon = 1e-24; //Whenever exp(-G^2/eta ) is smaller than epsilon: don't include it in the summation
   double radius = sqrt( - eta * log(epsilon) );
   this->maxN = ceil(giveNmaxOverR()*radius)+0.1;

   std::cout << radius << "\t" << maxN << std::endl;

}

Ewald_pt1::~Ewald_pt1(){

   delete [] nZvals;
   delete [] dTvecs;
   delete [] dZloc;
   delete [] dGvecs;

}

double Ewald_pt1::CalcSignedVolume(){

   double result = dTvecs[0 + 3*0] * ( dTvecs[1 + 3*1] * dTvecs[2 + 3*2] - dTvecs[1 + 3*2] * dTvecs[2 + 3*1] )
                 + dTvecs[0 + 3*1] * ( dTvecs[1 + 3*2] * dTvecs[2 + 3*0] - dTvecs[1 + 3*0] * dTvecs[2 + 3*2] )
                 + dTvecs[0 + 3*2] * ( dTvecs[1 + 3*0] * dTvecs[2 + 3*1] - dTvecs[1 + 3*1] * dTvecs[2 + 3*0] );
   return result;

}

void Ewald_pt1::CalcReciprocal(){

   double prefact = 2 * M_PI / CalcSignedVolume();
   
   for (int cnt=0; cnt<=2; cnt++){
      for (int cnt2=0; cnt2<=2; cnt2++){
         dGvecs[cnt + 3*cnt2] = prefact * ( dTvecs[ (cnt+1)%3 + 3 * ((cnt2+1)%3) ] * dTvecs[ (cnt+2)%3 + 3 * ((cnt2+2)%3) ] -  dTvecs[ (cnt+2)%3 + 3 * ((cnt2+1)%3) ] * dTvecs[ (cnt+1)%3 + 3 * ((cnt2+2)%3) ] ) ;
      }
   }

   std::cout << "G-vecs" << std::endl;
   for (int cnt=0; cnt<3; cnt++)
      std::cout << dGvecs[cnt + 3*0] << "\t" << dGvecs[cnt + 3*1] << "\t" << dGvecs[cnt + 3*2] << std::endl;


}

double Ewald_pt1::NN(){

   int totalCharge = 0;
   int sumZsquared = 0;
   for (int cnt=0; cnt<nNuclei; cnt++){
      totalCharge += nZvals[cnt];
      sumZsquared += nZvals[cnt]*nZvals[cnt];
   }
   double result_pt1 = - (2 * totalCharge * totalCharge) * M_PI / Volume / eta
                     - 0.5 * sumZsquared * sqrt(eta/M_PI);

   double result_pt2 = 0.0;

   for (int cnt_n1 = 0; cnt_n1<=maxN; cnt_n1++){
      for (int cnt_n2 = -maxN; cnt_n2<=maxN; cnt_n2++){
         for (int cnt_n3 = -maxN; cnt_n3<=maxN; cnt_n3++){
            if (!((cnt_n1==0) && (cnt_n2==0) && (cnt_n3==0))){
               //This defines \vec{G}. Note that -\vec{G} should also be taken into account here, as the summation over cnt_n1 only does positive values.
               double Gx = cnt_n1*dGvecs[0 + 3*0] + cnt_n2*dGvecs[1 + 3*0] + cnt_n3*dGvecs[2 + 3*0];
               double Gy = cnt_n1*dGvecs[0 + 3*1] + cnt_n2*dGvecs[1 + 3*1] + cnt_n3*dGvecs[2 + 3*1];
               double Gz = cnt_n1*dGvecs[0 + 3*2] + cnt_n2*dGvecs[1 + 3*2] + cnt_n3*dGvecs[2 + 3*2];

               double Gsquared = Gx*Gx + Gy*Gy + Gz*Gz;
               double factor = exp( - Gsquared / eta );
               if (factor>=epsilon){

                  for (int cnt1=0; cnt1<nNuclei; cnt1++)
                     for (int cnt2=0; cnt2<nNuclei; cnt2++){

                        double seed = Gx * (dZloc[cnt1+nNuclei*0] - dZloc[cnt2+nNuclei*0]) + Gy * (dZloc[cnt1+nNuclei*1] - dZloc[cnt2+nNuclei*1])
                                    + Gz * (dZloc[cnt1+nNuclei*2] - dZloc[cnt2+nNuclei*2]);
                        result_pt2 += factor * cos( seed ) / Gsquared * nZvals[cnt1] * nZvals[cnt2];

                     }
               }
            }
         }
      }
   }

   result_pt2 *= 4 * M_PI / Volume;

   double result_pt3 = 0.0; //--> erfc part
   //erfc(10) = 10^(-45)
   int maxNvalForTSum = 3 + ceil(giveNmaxOverR(false)*20/sqrt(eta)) + 0.1;

   std::cout << maxNvalForTSum << std::endl;

   for (int N1=-maxNvalForTSum; N1<=maxNvalForTSum; N1++){
      for (int N2=-maxNvalForTSum; N2<=maxNvalForTSum; N2++){
         for (int N3=-maxNvalForTSum; N3<=maxNvalForTSum; N3++){
            double Tx = N1 * dTvecs[0 + 3*0] + N2 * dTvecs[1 + 3*0] + N3 * dTvecs[2 + 3*0];
            double Ty = N1 * dTvecs[0 + 3*1] + N2 * dTvecs[1 + 3*1] + N3 * dTvecs[2 + 3*1];
            double Tz = N1 * dTvecs[0 + 3*2] + N2 * dTvecs[1 + 3*2] + N3 * dTvecs[2 + 3*2];

            for (int cnt1=0; cnt1<nNuclei; cnt1++){
               for (int cnt2=0; cnt2<nNuclei; cnt2++){
                  double LengthVec = sqrt( (dZloc[cnt1+3*0] - dZloc[cnt2+3*0] + Tx)*(dZloc[cnt1+3*0] - dZloc[cnt2+3*0] + Tx)
                                         + (dZloc[cnt1+3*1] - dZloc[cnt2+3*1] + Ty)*(dZloc[cnt1+3*1] - dZloc[cnt2+3*1] + Ty)
                                         + (dZloc[cnt1+3*2] - dZloc[cnt2+3*2] + Tz)*(dZloc[cnt1+3*2] - dZloc[cnt2+3*2] + Tz) );
                  if (!((N1==0) && (N2==0) && (N3==0) && (cnt1==cnt2))){
                     result_pt3 += 0.5 * erfc( LengthVec * sqrt(eta) * 0.5 ) / LengthVec * nZvals[cnt1] * nZvals[cnt2];
                  }
               }
            }
         }
      }
   }

   std::cout << "NN\t" << result_pt1 + result_pt2 + result_pt3 << std::endl;

   return (result_pt1 + result_pt2 + result_pt3);

}

void Ewald_pt1::eN(double * result, double zeta_A, int lmax_A, double Ax, double Ay, double Az, double zeta_B, int lmax_B, double Bx, double By, double Bz){
 
   double zeta = zeta_A+zeta_B;

   double Px = (zeta_A * Ax + zeta_B * Bx)/zeta;
   double Py = (zeta_A * Ay + zeta_B * By)/zeta;
   double Pz = (zeta_A * Az + zeta_B * Bz)/zeta;

   double beta = zeta_A*zeta_B/zeta * ((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By) + (Az-Bz)*(Az-Bz));

   int lmax = std::max(lmax_A,lmax_B);
   double * elementx_cos = new double[(lmax+1)*(lmax+1)];
   double * elementx_sin = new double[(lmax+1)*(lmax+1)];
   double * elementy_cos = new double[(lmax+1)*(lmax+1)];
   double * elementy_sin = new double[(lmax+1)*(lmax+1)];
   double * elementz_cos = new double[(lmax+1)*(lmax+1)];
   double * elementz_sin = new double[(lmax+1)*(lmax+1)];

   //create two arrays cos/sin of size count_A times count_B in which the results will be copied.
   int dim_A = (lmax_A+1)*(lmax_A+2)/2;
   int dim_B = (lmax_B+1)*(lmax_B+2)/2;
   int dimension = dim_A*dim_B;
   int inc1 = 1;
   double * cos_term_rightsize = new double[dimension];
   double * sin_term_rightsize = new double[dimension];
   double * temp_rightsize = new double[dimension];
   for (int counter=0; counter<dimension; counter++)
      result[counter] = 0.0;

   for (int cnt_n1 = 0; cnt_n1<=maxN; cnt_n1++){
      for (int cnt_n2 = -maxN; cnt_n2<=maxN; cnt_n2++){
         for (int cnt_n3 = -maxN; cnt_n3<=maxN; cnt_n3++){
            if (!((cnt_n1==0) && (cnt_n2==0) && (cnt_n3==0))){
               //This defines \vec{G}. Note that -\vec{G} should also be taken into account here, as the summation over cnt_n1 only does positive values.
               double Gx = cnt_n1*dGvecs[0 + 3*0] + cnt_n2*dGvecs[1 + 3*0] + cnt_n3*dGvecs[2 + 3*0];
               double Gy = cnt_n1*dGvecs[0 + 3*1] + cnt_n2*dGvecs[1 + 3*1] + cnt_n3*dGvecs[2 + 3*1];
               double Gz = cnt_n1*dGvecs[0 + 3*2] + cnt_n2*dGvecs[1 + 3*2] + cnt_n3*dGvecs[2 + 3*2];
               
               double Gsquared = Gx*Gx + Gy*Gy + Gz*Gz;
               double factor = exp( - Gsquared / eta );
               if (factor>=epsilon){
               
                  FillCosSinOneElectron(elementx_cos, elementx_sin, lmax, Px*Gx, zeta, Ax, Bx, Px, Gx);
                  FillCosSinOneElectron(elementy_cos, elementy_sin, lmax, Py*Gy, zeta, Ay, By, Py, Gy);
                  FillCosSinOneElectron(elementz_cos, elementz_sin, lmax, Pz*Gz, zeta, Az, Bz, Pz, Gz);

                  //Combine the right elements into the arrays cos/sin of size count_A times count_B. Note that the factorization occured on the exp(i vec(G).vec(r)) level and that therefore the (cos_x + i sin_x)()() terms have to be worked out.
                  int count_A = -1;
                  for (int nx_A = lmax_A; nx_A>=0; nx_A--){
                     for (int ny_A = lmax_A - nx_A; ny_A>=0; ny_A--){
                        int nz_A = lmax_A - nx_A - ny_A;
                        count_A++;

                        int count_B = -1;
                        for (int nx_B = lmax_B; nx_B>=0; nx_B--){
                           for (int ny_B = lmax_B - nx_B; ny_B>=0; ny_B--){
                              int nz_B = lmax_B - nx_B - ny_B;
                              count_B++;
                              cos_term_rightsize[count_A + dim_A*count_B] = elementx_cos[nx_A+lmax*nx_B]*elementy_cos[ny_A+lmax*ny_B]*elementz_cos[nz_A+lmax*nz_B]
                                                                          - elementx_cos[nx_A+lmax*nx_B]*elementy_sin[ny_A+lmax*ny_B]*elementz_sin[nz_A+lmax*nz_B]
                                                                          - elementx_sin[nx_A+lmax*nx_B]*elementy_cos[ny_A+lmax*ny_B]*elementz_sin[nz_A+lmax*nz_B]
                                                                          - elementx_sin[nx_A+lmax*nx_B]*elementy_sin[ny_A+lmax*ny_B]*elementz_cos[nz_A+lmax*nz_B];
                              sin_term_rightsize[count_A+dim_A*count_B] = - elementx_sin[nx_A+lmax*nx_B]*elementy_sin[ny_A+lmax*ny_B]*elementz_sin[nz_A+lmax*nz_B]
                                                                          + elementx_sin[nx_A+lmax*nx_B]*elementy_cos[ny_A+lmax*ny_B]*elementz_cos[nz_A+lmax*nz_B]
                                                                          + elementx_cos[nx_A+lmax*nx_B]*elementy_sin[ny_A+lmax*ny_B]*elementz_cos[nz_A+lmax*nz_B]
                                                                          + elementx_cos[nx_A+lmax*nx_B]*elementy_cos[ny_A+lmax*ny_B]*elementz_sin[nz_A+lmax*nz_B];
                           }
                        }
                     }
                  }

                  //For the formule from the notes: 2cos(G(r-r_alpha)) = 2[ cos(Gr) cos(Gr_alpha) + sin(Gr) sin(Gr_alpha)  ]
                  //Make the matrix count_A times count_B containing sum_alpha Z_alpha  2[ cos(Gr) cos(Gr_alpha) + sin(Gr) sin(Gr_alpha)]

                  for (int cnter=0; cnter<dimension; cnter++)
                     temp_rightsize[cnter] = 0.0;

                  for (int cnt_nuclei = 0; cnt_nuclei < nNuclei; cnt_nuclei++){
                     double seed = Gx*dZloc[cnt_nuclei + nNuclei*0] + Gy*dZloc[cnt_nuclei + nNuclei*1] + Gz*dZloc[cnt_nuclei + nNuclei*2];
                     for (int cnter=0; cnter<dimension; cnter++)
                        temp_rightsize[cnter] += 2 * nZvals[cnt_nuclei] * ( cos(seed) * cos_term_rightsize[cnter] + sin(seed) * sin_term_rightsize[cnter]);
                  }

                  // Multiply with exp_min_alpha / G^2 and add it to the current sum_matrix
                  double prefactor = exp( - Gsquared * (0.25 / zeta + 1.0/eta) ) /Gsquared;
                  //daxpy: y = y + a*x
                  daxpy_(&dimension,&prefactor,temp_rightsize,&inc1,result,&inc1);
                  
               }
            }
         }
      }
   }

   //multiply the sum_matrix of size count_A times count_B with (pi/zeta)^(3/2) * exp(-beta) * (- 4 pi) / Omega
   double prefactor = - ( 4 * M_PI / (this->Volume) ) * pow(M_PI/zeta,1.5) * exp(-beta);
   dscal_(&dimension,&prefactor,result,&inc1);

   delete [] elementx_cos;
   delete [] elementx_sin;
   delete [] elementy_cos;
   delete [] elementy_sin;
   delete [] elementz_cos;
   delete [] elementz_sin;

   delete [] cos_term_rightsize;
   delete [] sin_term_rightsize;
   delete [] temp_rightsize;

   std::cout << "Whee at Ewald_pt1" << std::endl;

}


void Ewald_pt1::FillCosSinOneElectron(double * element_cos, double * element_sin, int lmax, double seed, double zeta, double Aco, double Bco, double Pco, double Gco){

   element_cos[0+lmax*0] = cos( seed );
   element_sin[0+lmax*0] = sin( seed );

   double PcoMinAco = Pco - Aco;
   double PcoMinBco = Pco - Bco;
   double halfOverZeta = 0.5/zeta;
   double Gfact = halfOverZeta * Gco;

   for (int cnt_l=1; cnt_l<=lmax; cnt_l++){
      for (int cnt_bound=0; cnt_bound<cnt_l; cnt_bound++){
         element_cos[cnt_l + lmax*cnt_bound] = PcoMinAco * element_cos[(cnt_l-1) + lmax*cnt_bound] - Gfact * element_sin[(cnt_l-1) + lmax*cnt_bound];
         element_cos[cnt_bound + lmax*cnt_l] = PcoMinBco * element_cos[cnt_bound + lmax*(cnt_l-1)] - Gfact * element_sin[cnt_bound + lmax*(cnt_l-1)];
         element_sin[cnt_l + lmax*cnt_bound] = PcoMinAco * element_sin[(cnt_l-1) + lmax*cnt_bound] + Gfact * element_cos[(cnt_l-1) + lmax*cnt_bound];
         element_sin[cnt_bound + lmax*cnt_l] = PcoMinBco * element_sin[cnt_bound + lmax*(cnt_l-1)] + Gfact * element_cos[cnt_bound + lmax*(cnt_l-1)];
         if ((cnt_l-2)>=0){
            element_cos[cnt_l + lmax*cnt_bound] += halfOverZeta * (cnt_l-1) * element_cos[(cnt_l-2) + lmax*cnt_bound];
            element_cos[cnt_bound + lmax*cnt_l] += halfOverZeta * (cnt_l-1) * element_cos[cnt_bound + lmax*(cnt_l-2)];
            element_sin[cnt_l + lmax*cnt_bound] += halfOverZeta * (cnt_l-1) * element_sin[(cnt_l-2) + lmax*cnt_bound];
            element_sin[cnt_bound + lmax*cnt_l] += halfOverZeta * (cnt_l-1) * element_sin[cnt_bound + lmax*(cnt_l-2)];
         }
         if ((cnt_bound-1)>=0){
            element_cos[cnt_l + lmax*cnt_bound] += halfOverZeta * cnt_bound * element_cos[(cnt_l-1) + lmax*(cnt_bound-1)];
            element_cos[cnt_bound + lmax*cnt_l] += halfOverZeta * cnt_bound * element_cos[(cnt_bound-1) + lmax*(cnt_l-1)];
            element_sin[cnt_l + lmax*cnt_bound] += halfOverZeta * cnt_bound * element_sin[(cnt_l-1) + lmax*(cnt_bound-1)];
            element_sin[cnt_bound + lmax*cnt_l] += halfOverZeta * cnt_bound * element_sin[(cnt_bound-1) + lmax*(cnt_l-1)];
         }
      }  
      element_cos[cnt_l + lmax*cnt_l] = PcoMinAco * element_cos[(cnt_l-1) + lmax*cnt_l] - Gfact * element_sin[(cnt_l-1) + lmax*cnt_l];
      element_sin[cnt_l + lmax*cnt_l] = PcoMinAco * element_sin[(cnt_l-1) + lmax*cnt_l] + Gfact * element_cos[(cnt_l-1) + lmax*cnt_l];
      element_cos[cnt_l + lmax*cnt_l] += halfOverZeta * cnt_l * element_cos[(cnt_l-1) + lmax*(cnt_l-1)];
      element_sin[cnt_l + lmax*cnt_l] += halfOverZeta * cnt_l * element_sin[(cnt_l-1) + lmax*(cnt_l-1)];
      if ((cnt_l-2)>=0){
         element_cos[cnt_l + lmax*cnt_l] = halfOverZeta * (cnt_l-1) * element_cos[(cnt_l-2) + lmax*cnt_l];
         element_sin[cnt_l + lmax*cnt_l] = halfOverZeta * (cnt_l-1) * element_sin[(cnt_l-2) + lmax*cnt_l];
      }
   }


}

double Ewald_pt1::giveNmaxOverR(bool useG){

   //standard useG = true
   //if false: use Tvecs

   /***
      (m1 G1 + m2 G2 + m3 G3) (n1 G1 + n2 G2 + n3 G3) = m^(dagger) G n
      x = G^{1/2} n --> m^(dagger) G n = y^dagger x
      if |xi| = r you move on the cube that binds the sphere you want to stay in
      so |n|_max = max_i sum_j | (G^{-1/2})_ij | radius --> worst case n values to stay within sphere      
   ***/

   double * metric = new double[3*3];
   for (int cnt=0; cnt<3; cnt++)
      for (int cnt2=cnt; cnt2<3; cnt2++){
         if (useG)
            metric[cnt+3*cnt2] = dGvecs[cnt+3*0]*dGvecs[cnt2+3*0]+dGvecs[cnt+3*1]*dGvecs[cnt2+3*1]+dGvecs[cnt+3*2]*dGvecs[cnt2+3*2];
         else
            metric[cnt+3*cnt2] = dTvecs[cnt+3*0]*dTvecs[cnt2+3*0]+dTvecs[cnt+3*1]*dTvecs[cnt2+3*1]+dTvecs[cnt+3*2]*dTvecs[cnt2+3*2];
         metric[cnt2+3*cnt] = metric[cnt+3*cnt2];
      }

   char jobz = 'V';
   char uplo = 'U';
   int matrixdim = 3;
   double * eigs = new double[matrixdim];
   int lwork = matrixdim*matrixdim;
   double * work = new double[lwork];
   int info = 0;

   dsyev_(&jobz,&uplo,&matrixdim,metric,&matrixdim,eigs,work,&lwork,&info);
   for (int cnt=0; cnt<3; cnt++){
      int inc = 1;
      double scal = pow(eigs[cnt],-0.25);
      dscal_(&matrixdim,&scal,metric+cnt*matrixdim,&inc);
   }

   char trans = 'T';
   char notrans = 'N';
   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&notrans,&trans,&matrixdim,&matrixdim,&matrixdim,&alpha,metric,&matrixdim,metric,&matrixdim,&beta,work,&matrixdim);

   double val1 = fabs(work[0 + 3*0]) + fabs(work[0+3*1]) + fabs(work[0+3*2]);
   double val2 = fabs(work[1 + 3*0]) + fabs(work[1+3*1]) + fabs(work[1+3*2]);
   double val3 = fabs(work[2 + 3*0]) + fabs(work[2+3*1]) + fabs(work[2+3*2]);

   delete [] work;
   delete [] eigs;
   delete [] metric;

   return std::max(std::max(val1,val2),val3);

}










