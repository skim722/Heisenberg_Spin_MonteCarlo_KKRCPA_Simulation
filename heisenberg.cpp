#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <sstream>
#include "omp.h"
#include "heisenberg.h"

using namespace std;

// invoke as: heisenberg LofCell NofMC HScalar TinmeV
int main (int argc, char *argv[]) 
{
//=====================================[ Part 1 :Read Input]============================================//
clock_t start;
double duration_1, duration_2, duration_3, duration_4, duration_5, duration_6, duration_7,duration_8,duration_9, duration_10;
start= clock();
// Construct manget object from constructed spin map 
int ispin = 0;
double H[3];
double Hdirection[3];
double SpinI[3];
srand(time(NULL)); 
cout << " Automatically read in unit cell.\n";
string cfile = "file.dat";
string mfile= "mfile.dat";
Magnet  unitMag(cfile,mfile);
duration_1=(clock()-start)/(double)CLOCKS_PER_SEC;
cout<<"Time consumed for part 1 :"<<duration_1<<'\n'; 

//=====================================[ Part 2 :Replicate Unitecell]====================================//
start=clock();
cout<<"CHECK Point -3";
int n;
if (argc>0) n = atoi(argv[1]);  
cout<<"CHECK POINT 2";
cout << " Entered number of unit cells in each direction. " << n << endl;
cout << " Replicate supercell from unit cell.\n";
Magnet superMag(n,n,n,unitMag);
duration_2=(clock()-start)/(double)CLOCKS_PER_SEC;
cout<<"Time consumed for part 2 :"<<duration_2<<'\n';

//=====================================[ Part 3 :Initialization ]==========================================//
start=clock();
int nTot = n*n*n*unitMag.getnTot();
cout << " Total number of atoms is: " << nTot << endl;
int MCSteps;
if (argc>1) MCSteps = atoi(argv[2]);
cout << " Entered number of Monte Carlo steps: " << MCSteps << endl;
initialize(superMag,argc,argv);
double T = superMag.getT();
cout << " T " << T << endl;
double K = superMag.getK();
cout<< "K " <<K <<endl; 
superMag.getH(H);
superMag.getD(Hdirection);
cout<< " H "<< H[0]<<' ' <<H[1]<<' ' <<H[2]<<endl;
cout<< " H direction " << Hdirection[0]<< ' '<<Hdirection[1]<< ' ' <<Hdirection[2]<<endl;
print_coordinate("Coordinate.dat", superMag);
// Seek an optimal step size, delta_max.
// delta_max will depend on temperature.
// Then just use a fixed delta_max for the entire run.
duration_3=(clock()-start)/(double)CLOCKS_PER_SEC;
cout<<"Time consumed for part 3 :"<<duration_3<<'\n';

//=====================================[ Part 4 : Optimize Delta_max  ]========================================//
start=clock();
Magnet dumMag(superMag);
double delta_max;
delta_max = 0.10;
double edum;
for (int ii=0;ii<5000;ii++)
     {
          double ratio;
          ratio = oneMonteCarloStepPerSpin(dumMag,delta_max);
	  // ratio = oneMonteCarloStepPerSpin(superMag,delta_max);
           if     (ratio<0.30) {delta_max*=0.90;}
           else if(ratio<0.45) {delta_max*=0.99;}
           else if(ratio>0.70) {delta_max*=1.09;}
           else if(ratio>0.55) {delta_max*=1.01;}
          // limit stepsize at ultrahigh T
          // seems like a larger cutoff would be OK at moderate T...
          if (delta_max>0.5 ) {delta_max = 0.5 ;}
          for (int n=0;n<nTot;n++) edum += energyPerSpin(n, dumMag);
     }
cout << " delta Max " << delta_max << endl;
duration_4=(clock()-start)/(double)CLOCKS_PER_SEC;
cout<<"Time consumed for part 4 :"<<duration_4<<'\n';

//=====================================[ Part 5 :Thermalization ]=============================================//
start=clock();
for (int t=0;t<1;t++)
     {                            //loop for temperature control
          int thermSteps = int(0.2 * MCSteps);
          cout<<"Starting Spin Set at T="<<T<<'\n'<<endl;
          for (int i = 0; i < thermSteps; i++)
              {
                     oneMonteCarloStepPerSpin(superMag,delta_max);
                   //  cout << i;
              }
          cout << " end therm stage" << endl;
          duration_5=(clock()-start)/(double)CLOCKS_PER_SEC;
          cout<<"Time consumed for part 5 :"<<duration_5<<"\n";

//=====================================[ Part 6 :Monte Carlo Steps ]===========================================//
          start=clock();
          ostringstream h0,h1,h2;
          h0<<H[0]; h1<<H[1]; h2<<H[2];
          ostringstream Temp;
          Temp<<T;
          string name="output_";
          name=name+h0.str()+h1.str()+h2.str()+"_"+Temp.str();
          ofstream file(name.c_str());
          // total magnetization and energy, to be summed over replicas
          double mAv=0,mAv_d=0, m2Av = 0,m2Av_d=0, eAv = 0, e2Av = 0, MAv=0;
          double sum_mAv=0.,sum_mAv_d=0,sum_m2Av=0.0,sum_m2Av_d=0, sum_eAv=0, sum_e2Av=0, sum_MAv=0 ; 
          // species-specific magnetization, to be summed over replicas
          int nSp = superMag.getNumSpecies();
          assert(nSp>0);
          double mAv_Species[nSp],sum_mAv_Species[nSp];
          for (int iSp=0;iSp<nSp;iSp++)
              {
                   mAv_Species[iSp]=0.;
                   sum_mAv_Species[iSp]=0.;
              }
//=====================================[[[[[MONTE CARLO STEPS: START]]]]]========================================//
          for (int kMC = 0; kMC < MCSteps; kMC++) 
              {
                   start=clock();
                   for (int mmm=0;mmm<100;mmm++)
                         {
                              oneMonteCarloStepPerSpin(superMag,delta_max);
                         }
                   cout<<"Monte Carlo Step #: "<< kMC<<endl;
                   double t=(clock()-start)/(double)CLOCKS_PER_SEC;
                   cout<<"time "<<t<<endl;
                   // initialize species-dependent net magnetization for single snapshot
                   double siteSumM_Vector_Species[nSp][DIMENSIONS], siteSumM_Scalar_Species[nSp];
                   int natoms_Species[nSp];
                   for (int iSp=0;iSp<nSp;iSp++)
                         {
                             natoms_Species[iSp]=0;
                              for (int k=0;k<DIMENSIONS;k++) 
                                   for (int iSp=0;iSp<nSp;iSp++) 
                                         siteSumM_Vector_Species[iSp][k]=0.0;
                         }
                   // sum vector moments over sites
                   double spinI[DIMENSIONS];
                   int nTot=superMag.getnTot();
                   double K=superMag.getK();
                   int iSp;
                   for (int i=0;i<nTot;i++)
                        { 
                             string st1, st2;
                             st1 = superMag.getglName(i);
                             superMag.getIthSpin(i,spinI);      
                             assert(fabs(mynorm(spinI)-1.)<1.e-8);
                             double momentI = superMag.getglMoment(i); 
                             //match atom i with the correct species, sum the species-specific moment
                             int ismatch;
                             for (iSp=0;iSp<nSp;iSp++)
                                  {
                                       st2 = superMag.getSpeciesName(iSp);
                                       ismatch=0;
                                       if (strcmp(st1.c_str(),st2.c_str())==0)
                                            {
                                                  ismatch=1;
                                                  ++natoms_Species[iSp];
 	                                          for (int k=0;k<DIMENSIONS;k++) siteSumM_Vector_Species[iSp][k]+=momentI*spinI[k];
                                                  break;
                                            }
                                  }
                             assert(ismatch);
                        }
                   // normalize total energy, total magnetization, species magnetization
                   double ssum[3] = {0.,0.,0.};
                   assert(ssum[0]==0.&&ssum[1]==0.&&ssum[2]==0.);
                   for (int iSp=0;iSp<nSp;iSp++)
                       {
                            ssum[0]+=siteSumM_Vector_Species[iSp][0];
                            ssum[1]+=siteSumM_Vector_Species[iSp][1];
                            ssum[2]+=siteSumM_Vector_Species[iSp][2];
                            siteSumM_Scalar_Species[iSp] = mynorm(siteSumM_Vector_Species[iSp])
                                             /double(natoms_Species[iSp]);
                        }
                   double mtot = mynorm(ssum)/double(nTot);
                   double mtot_d=mydot(ssum,Hdirection)/double(nTot);
                   double etot = 0.;
                   for (int n=0;n<nTot;n++) etot += energyPerSpin(n, superMag);
                   edum = etot/double(nTot);
                   double sum_spin[3]={0.,0.,0.};
                   double sum_allspin[3]={0,0,0};
                   // sum over snapshots
                   sum_spin[0] += ssum[0]/double(nTot);
                   sum_spin[1] += ssum[1]/double(nTot);
                   sum_spin[2] += ssum[2]/double(nTot);  
                   sum_allspin[0]+=sum_spin[0];
                   sum_allspin[1]+=sum_spin[1];
                   sum_allspin[2]+=sum_spin[2];
                   sum_mAv += mtot; sum_m2Av += mtot * mtot;
                   sum_mAv_d += mtot_d; sum_m2Av_d += mtot_d*mtot_d;
                   sum_eAv += etot; sum_e2Av += etot * etot;
                   // sum over snapshots
                   for (int iSp=0;iSp<nSp;iSp++) sum_mAv_Species[iSp] += siteSumM_Scalar_Species[iSp];
                   double spinAvi[3];
                   spinAvi[0]=sum_allspin[0]/double(kMC+1);
                   spinAvi[1]=sum_allspin[1]/double(kMC+1);
                   spinAvi[2]=sum_allspin[2]/double(kMC+1);
                   mAv = sum_mAv/double(kMC+1);
                   m2Av = sum_m2Av/double(kMC+1);
                   mAv_d = sum_mAv_d / double(kMC+1);
                   for (int iSp=0;iSp<nSp;iSp++) mAv_Species[iSp]=sum_mAv_Species[iSp]/double(kMC+1);
                   eAv =sum_eAv/double(kMC+1);
                   e2Av =sum_e2Av/double(kMC+1);
                   if ((kMC+1)%100==0)
                        {
                              //PrintSpin(superMag); 
                              RotateMagnet(superMag,delta_max); 
                              //PrintSpin(superMag);  
                        }
                   if (kMC==0||(kMC+1)%500==0)
                        {
                              //PrintSpin(superMag);
                              oneMonteCarloGlobalSpinReversal(superMag); 
                              //PrintSpin(superMag);
                         }
                   file << kMC+1 <<"\t" <<"spin x y z\t"<<sum_spin[0]<<"\t"<<sum_spin[1]<<"\t"<<sum_spin[2]<<"\tspin Ave\t" <<spinAvi[0]<<"\t"<<spinAvi[1]<<"\t"<<spinAvi[2]<<"\tmtot_d\t"<< mtot_d <<"\t"<<mAv_d<<"\tmtot\t"<<mtot<<" \t"<<mAv<<"\tSpontaneous Energy\t"<<etot<<"\tAve_E\t"<<eAv/double(nTot)<<endl;
                  if((kMC+1)%2000==0)
                        print_magnetic("mfile",superMag,T,kMC+1,H[0],H[1],H[2]);          
 
            }
//=====================================[[[[[MONTE CARLO STEPS: END]]]]]========================================//
           cout << edum << " " << eAv << " wwe EEEE" << endl;
           file << "#nTot): " << nTot << "   H(Magnetic Field): " << H[0] << " " << H[1] << " " << H[2] ;
           file <<'\t'<<" MCSteps: "<< MCSteps <<"  "<<"T"<<'\t'<<  T <<'\t' << " <m> " <<'\t'<< mAv ;
           for (int iSp=0;iSp<nSp;iSp++) file <<"  <m_" << superMag.getSpeciesName(iSp) << "> " << mAv_Species[iSp] ;
           file <<"M"<<'\t'<<MAv<<'\t' <<" <e> = " <<'\t'<< eAv/double(nTot)<<'\t'<<" Cv "<<'\t'<<double(e2Av-eAv*eAv)/double(nTot)/T/T<< endl; 
           file.close();
           // repeat the last line, appending it to another file
           ofstream afile("heisenberg.data",ios::app);
           afile << " " <<K<<"    \t"<< nTot << "                     " << H[0] << " " << H[1] << " " << H[2] ;
           afile <<'\t'<<"          "<< MCSteps <<"  "<<" "<<'\t'<<  T <<'\t' << "     " <<'\t'<< mAv_d<<'\t'<<"    "<<mAv<<'\t';
           for (int iSp=0;iSp<nSp;iSp++) afile <<"     " << iSp << "  " << mAv_Species[iSp] ;
           afile <<" "<<'\t'<<MAv<<'\t' <<"       " <<'\t'<< eAv/double(nTot)<<'\t'<<"    "<<'\t'<<double(e2Av-eAv*eAv)/double(nTot)/T/T<< endl; 
           afile.close();
     }
 duration_6=(clock()-start)/(double)CLOCKS_PER_SEC;
 cout<<"Time consumed for part 6 :"<<duration_6<<'\n';
//======================================DONE and record time consumption========================================//
 ofstream file("time_consumption.out");
 file<<"Part 1 : Get data and construct magnet object \n"<<duration_1<<'\n';
 file<<"Part 2 : Replicate supercell from unit cell \n"<<duration_2<<'\n';
 file<<"Part 3 : Initialize spin as a random fluctuation \n"<<duration_3<<'\n';
 file<<"Part 4 : Copy magnet object to dummag and optimize deltamax \n" <<duration_4<<'\n';
 file<<"Part 5 : Multiply m_Fe and m_B to each spins \n"<<duration_5<<'\n';
 file<<"Part 6 : All Monte Carlo Steps for MCStep sizes \n"<<duration_6<<endl;
 file.close();
}

