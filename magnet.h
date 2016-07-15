/********************************************************/
/* magnet.h -- Class definitions for magnet.            */
/* Supplement the unit cell (lattice and atomic basis)  */
/* from Crystal class with a derived neighbor list.     */
/* Or alternatively input neighbor information directly */
/* with no associated crystallographic information.     */
/********************************************************/
// Each atom in the simulation is indexed 1-NTOT,
// this list is associated with global arrays, gl***.
// 
#ifndef MAGNET_H
#define MAGNET_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "crystal.h"
#include "vector.h"
#include "KKR.h"

static const int DIMENSIONS=3; // spin dimension
                               // This can be varied.

using namespace std;

class Magnet: public Crystal {
public:
  Magnet();                               //Default constructor
  Magnet(const string m, const string c); //constructor from Crystal input file
  Magnet(const Magnet & ref);             //Copy constructor
  Magnet(int n1, int n2, int n3, const Magnet & ref);      //Copy constructor, replicated
  ~Magnet();                              //Deletes memory allocated 

  const Magnet &operator =(const Magnet &);  //Assignment operator
  bool operator ==(const Magnet &);          //Equality operator
  bool operator !=(const Magnet &);          //Inequality operator

  double Energy(int ith);   //Heisenberg energy of ith spin
  double totEnergy();       //Heisenberg energy of entire cell

  int getNumSpecies() const {return xtlData.getNumSpecies();}
  int getNumAtoms(const int i) const {return xtlData.getNumAtoms(i);}
  string getSpeciesName(int i) const { return xtlData.getSpeciesName(i); }
  double getIthSpin(int ith, double v[DIMENSIONS]) const {for (int i=0;i<DIMENSIONS;i++)v[i]=glSpins[ith][i];}
  double writeIthSpin(int ith, double v[DIMENSIONS]){for (int i=0; i<DIMENSIONS; i++) glSpins[ith][i]=v[i];}
  ///
  int getnTot() const {return nTot;}
  int getIJthShellSize(int ith, int jth) const {return shellSize[ith][jth];}
  int getIJKthNbr(int ith, int jth, int kth) const {return shellNbr[ith][jth][kth]; }
  double getIJKthJ(int ith, int jth, int kth) const {return J[ith][jth][kth]; }
  string getglName(int ith)const {return glName[ith];}
  double getglMoment(int ith) const{return glMoment[ith];}
   ///
  void  getglCoord(int ith, double * vec){ vec[0]=glCoord[ith][0]; vec[1]=glCoord[ith][1]; vec[2]=glCoord[ith][2]; } 
  double getRmax() const {return kkr.getRmax();}
  double setRmax() {rmax = kkr.getRmax();}
  int getIthLattice(int ith, double vecg[LATTICED]) const {xtlData.getIthLattice(ith,vecg);}
  int secondPass(int nTot, int nLines, int glSpeciesID[], int glAtomID[],
                 int* nShells, double** shellRadius2,
                 int** shellSize, int*** shellNbr, double rmax);
int getIthNShell(int i) const {return nShells[i];}
double setT(double Tin) {T=Tin;}
double getT() const {return T;}
double setr(double rin){r=rin;}
double getr()const{return r;} 
double settheta(double thetain){theta=thetain;}
double gettheta() const{return theta;}
double  setphi(double phiin){phi=phiin;}
double getphi() const{return phi;}
double setK(double Kin){K=Kin;}
double getK() const{return K;}
double setH(double h[3]) {H[0]=h[0]; H[1]=h[1]; H[2]=h[2];}
double getH(double h[3]) const {h[0]=H[0]; h[1]=H[1]; h[2]=H[2];} //%%% TYPO 
double setD(double d[3]) {D[0]=d[0]; D[1]=d[1]; D[2]=d[2];}
double getD(double d[3]) const {d[0]=D[0]; d[1]=D[1]; d[2]=D[2];} 

int getnMonte() const {return nMonte;}
int setnMonte(int i) {nMonte=i;}

private:

// Inherit crystallographic data.
Crystal xtlData;

// Reformat Crystal data:
// All nTot atoms are grouped into a global 1D list.
// ith global atom is indexed in Crystal as 
// glAtomID[i]'th member of glSpeciesID[i].
int nTot;
int*     glSpeciesID;
string*  glName;
int*     glAtomID;
 double** glCoord; //

// added 4/1/15 MPS local moment to be read from file.dat in crystal for each site, 1-nTot
double* glMoment;

double rmax;  // this is approx the cutoff for Heisenberg interactions

// Each site has a spin vector.
double** glSpins;          // ith = 1 to nTot

// canned procedure to convert Crystal to global variables.
int reformatCrystal();

// canned procedure to compute neighbor lists out to rmax
int makeNbrLists(double rmax);
int statusReport();

// Define neighbor shells and neighbor lists.
// If this data is input directly, then it is possible to have
// nShells=1 and shellRadius2=undefined;
                         // otherwise:
int*     nShells;        // # of shells around ith site
int**    shellSize;      // # of neighbors in each such shell, sitesXshells
double** shellRadius2;   // squared nbr distance of each shell, sitesXshells
int***   shellNbr;       // index of atoms in each shell, sitesXshellsXnbrs

// Define Heisenberg coupling constants over sitesXshellsXnbrs
double*** J;

// save temperature as one of the magnet parameters.
double T;
//save angle
double r;
double theta;
double phi;
// save external field
double H[3];
//save field direction
double D[3];
//save K
double K;
// save step count
int nMonte;

KKR kkr;

};

#endif //MAGNET_H
