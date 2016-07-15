#include <cmath>
#include <cstdlib>
#include <iostream>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random/variate_generator.hpp>
#include <fstream>
#include <float.h>
#include <string>
//#include <limits>
#include <assert.h>
#include <sstream>
#include <string> 
#include "magnet.h"
#include "crystal.h"
#include "vector.h"
#define PI 3.14159265
using namespace std;

//const long double expmin = logl(numeric_limits<long double>::min());
//const long double expmax = logl(numeric_limits<long double>::max());
//Toots
inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}


inline double mydot(double v[], double w[])
{
 double s = 0.;
  for (int i=0;i<DIMENSIONS;i++) s += v[i]*w[i];
   return s;
   }

inline double mynorm(double v[])
{
 double s=0.0;
  for (int i=0;i<DIMENSIONS;i++) s +=v[i]*v[i];
  return sqrt(s);
}

inline double std_rand();
int normalize(double a[3],double n[3]);
void cross(const double a[3], const double b[3], double c[3]);
float rand_gauss (void);
void rotation(const double u[3], double theta, double s[3],double rot_s[3]);
double angle( double v[], double u[]);
void conv_spherical(double xyz[3],double spherical[3]);
void conv_cartesian(double spherical[3],double xyz[3]);
//Metropolis Monte Carlo 
void initialize_spinmap_prevrun(Magnet &,int argc,char* argv[] );
void initialize_spinmap(Magnet &);
void initialize_spinmap_random(Magnet &);
void initialize(Magnet&,int argc, char* argv[]);
bool MetropolisStep(Magnet& max, double delta_max);
void try_move(double * try_spin, double delta_max);
void try_move_sphere(double * try_spin);
void print_magnetic(string filename, Magnet& mag, double T,int MCSteps,double h0, double h1, double h2 );
int findnRhomb(double rmax, double lattice[3][3], int nRhomb[3]);
//Calculate Observable
double magnetizationPerSpin(Magnet&);
double magnetization(Magnet&);
double energyPerSpin(int i, Magnet& max);
double dEnergydSpin(int i, Magnet& max);
double acceptanceRatio;
int steps = 0;                  // steps so far



void oneMonteCarloGlobalSpinReversal(Magnet& max){
     double spini[3];
     double T=max.getT();
     int nTot=max.getnTot();
     int i=int(nTot*std_rand());
     double e0=0,e1=0;
     for(int i=0; i<nTot;i++){e0+=energyPerSpin(i,max);}
     Magnet dummy(max);

     for (int i=0;i<nTot;i++){
       dummy.getIthSpin(i,spini);
       spini[0]*=-1;
       spini[1]*=-1;
       spini[2]*=-1;
       dummy.writeIthSpin(i,spini);  
   
      }

     for (int i=0;i<nTot;i++) e1+=energyPerSpin(i, dummy);

     long double ratexp=e1-e0;
     if (ratexp<0) {
        double SpinTemp[3];
        for (int j=0;j<nTot;j++){dummy.getIthSpin(j,SpinTemp); max.writeIthSpin(j,SpinTemp);}     
     cout<<"Let's make it upside down!!"<<endl;//accept reversed movement
     }
     else  cout<<"Not reverse!"<<endl;
        //return to the previous status
}


void print_magnetic(string filename, Magnet& mag, double T, int MCSteps,double h0, double h1, double h2 ){
  ostringstream oss,mc,H0,H1,H2;
  string str;
  oss<<T; mc<<MCSteps;H0<<h0;H1<<h1;H2<<h2;
  filename=filename+"_"+oss.str()+"_"+H0.str()+H1.str()+H2.str()+"_"+mc.str()+".data";
  ofstream file(filename.c_str());
  int nTot = mag.getnTot();
  double spin[3];

  cout<<"Target T="<< T << endl;
      for (int j = 0; j < nTot; j++)
        {
          mag.getIthSpin(j,spin);
          str=mag.getglName(j);
          file<<str<<"\t";
          for (int k=0;k<DIMENSIONS;k++)
            {
              file<<spin[k]<<"\t";
            }
          file<<'\n';
        }
  file.close();

}

void print_coordinate(string filename, Magnet & mag){
  double glCoord[3];
  string str; 
  int nTot=mag.getnTot();
  ofstream file(filename.c_str());
  
  for (int i=0;i<nTot;i++)
   {
      str=mag.getglName(i);
      mag.getglCoord(i,glCoord);
      file<<str<<"\t";
      for (int k=0; k<3;k++)
       {  
           file<<glCoord[k]<<'\t';
       }       
       file<<'\n';
    }
   file.close();      
}


void initialize(Magnet& max, int argc, char* argv[]) {
           //(1) Setting up parameters

           double r,theta,phi; //theta and phi in degree
           //double u[3],a[3];
           //double easy[3]={0.0,0.0,1.0};

           if (argc>2) r = atof(argv[3]);
           if (argc>3) theta = atof(argv[4]);
           if (argc>4) phi= atof(argv[5]);
          
           double H[3],a[3];
           H[0] = r*sin(theta*PI/180)*cos(phi*PI/180);
           H[1] = r*sin(theta*PI/180)*sin(phi*PI/180);
           H[2] = r*cos(theta*PI/180);
           cout << " Entered magnetic field H: " << H[0] << " " << H[1] << " " << H[2] << endl;
           double r1,r2;
           r2=abs(sqrt(H[0]*H[0]+H[1]*H[1]+H[2]*H[2]));
           if (abs(r)<0.0001) {
                 r1=1;
                 a[0]=0; a[1]=0; a[2]=1;
                }
           else{
                r1 = 1/r2;
                a[0]=H[0]*r1; a[1]=H[1]*r1; a[2]=H[2]*r1;  
               }
           cout<< "Normalized field direction a: " <<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;      
           cout<<"Okay, I am going to set the filed.!";
           max.setr(r);
           max.settheta(theta);
           max.setphi(phi);
           max.setH(H);
           max.setD(a);
           double T;
           if (argc>5) T = atof(argv[6]);
           max.setT(T);
           cout << " Entered the highest temperature T: " << max.getT() << endl;
           double K;
           if (argc>6) K = atof(argv[7]);
           max.setK(K);
           cout<<" Entered K value: " <<max.getK()<<endl;
           double geth[3];
           max.getH(geth);
           cout<<"I will get the field  H="<<geth[0]<<" "<<geth[1]<<" "<<geth[2]<<endl;
           int number;
           cout<<"CHOOSE INITIALIZE OPTION? \n 1. Start from random spinmap \n 2. Start from minimum energy spinmap \n 3. Start from previous run\n Enter number and go for it!"<<endl;
           if (argc>7) number = atoi(argv[8]);
           cout<<number<<endl;
           if (number==1)       initialize_spinmap_random(max);  
           else if (number==2)  initialize_spinmap(max);
           else if (number==3)  initialize_spinmap_prevrun(max,argc,argv); 
           else initialize_spinmap(max);
}         


void initialize_spinmap(Magnet& max){
           double r=max.getr();
           cout<<"##########INITIALIZE SPIN MAP! #########################"<<endl; 
           double easy[3]={0.0,0.0,1.0};
           double u[3],a[3];
           double T=max.getT();
           max.getD(a);
           //cross product btw easy and H direction
           cross(easy,a,u);
           cout<<easy[0]<<" "<<easy[1]<<" "<<easy[2]<<endl;
           cout<<a[0]<<" " <<a[1]<<" "<<a[2]<<endl;
           cout<<"u_norm"<<u[0]<<" " <<u[1]<<" " <<u[2]<<endl;;
 
           // Initialize spin map 
           int nTot = max.getnTot();
           int i=1,n=1;
           double myangle;
           myangle=angle(easy,a);
           cout<<"####MYANGLE######"<<myangle<<endl;  
           double spini[3],rot_spini[3];
           double hard[3]={1/sqrt(3),1/sqrt(3),1/sqrt(3)};
           for (int j = 0; j < nTot; j++)
           { max.writeIthSpin(j,easy);}
           double stoping=50;           
                 
         //while(stoping >5.0){
         while(r!=0&&stoping>0.01){              
              i= int(nTot*std_rand()) ;
              double e0 = energyPerSpin( i,  max);
              cout<<"e0: " <<e0<<endl;
              //Rotate initial spin as myangle 
              max.getIthSpin(i, spini);
              rotation(u,myangle,spini,rot_spini);
              double r3=1./abs(sqrt(rot_spini[0]*rot_spini[0]+rot_spini[1]*rot_spini[1]+rot_spini[2]*rot_spini[2]));
              rot_spini[0]*=r3; rot_spini[1]*=r3; rot_spini[2]*=r3;//normalize function in vector.h doen't work somehow  
              cout<<"rotation axis" <<u[0]<<'\t'<<u[1]<<'\t'<<u[2]<<endl;
              cout<<"spin i" <<spini[0]<<'\t'<<spini[1]<<'\t'<<spini[2];
              cout<<"rot_spin"<<rot_spini[0]<<'\t'<<rot_spini[1]<<'\t'<<rot_spini[1]<<endl;
              
              for (int j=0;j<nTot;j++){max.writeIthSpin(j,rot_spini);}
              double e1=  energyPerSpin( i,  max);
              cout<<"e1 : "<<e1<<endl;       
  
              if (e1 < e0) {
                           myangle=-1.0*myangle/2.0; //search for the opposite angle
               }
              else{
                          for(int j=0;j<nTot;j++){
                                max.writeIthSpin(j,spini);} //return to the unmodified sataus
                             myangle=myangle/2.0;
                  }
                 stoping = abs(myangle); 
                 n++;

          // ############# DEBUGING ######################################################
		 max.getIthSpin(i,spini);             
            cout<< "############SPEEPEST DESCENT STEP : "<<n<<"##################"<<endl;
            cout<< "SPIN DIRECTION  : "<<spini[0]<<"\t"<<spini[1]<<"\t"<<spini[2]<<endl;
            cout <<"e0  :"<<e0 <<"   e1    :"<<e1<<endl;
            cout<<"abs(e0-e1)"<<abs(e0-e1)<<endl;
            cout<< "ANGLE  :"<<myangle<<endl;
            cout<<"OPTIMIZE STEP NUMBER   : " <<n<<endl;
            cout<<"################################################"<<endl;             
           }
           double H[3];
           max.getH(H);
	 //#######################################################
           print_magnetic("initial",max,T,0,H[0],H[1],H[2]);
           
}       


void initialize_spinmap_random(Magnet& max) {

	   double T=max.getT();
	   int nTot = max.getnTot();
	   double spin[DIMENSIONS];
	   for (int j = 0; j < nTot; j++)
	       {
		// choose a random vector, uniformly distributed inside unit sphere
		// then normalize to a unit vector uniformly distributed over 4pi
		double s = 10.;
		while (s>1.0 || s==0.)
		      {
		       for (int i=0;i<DIMENSIONS;i++) spin[i] = 2.*std_rand()-1.0;
		       s = mynorm(spin);
		      }
		for (int k=0;k<DIMENSIONS;k++) spin[k] /= s;
			max.writeIthSpin(j,spin);
          	       }

	             steps = 0;
		     max.setnMonte(steps);
                     double H_field[3]; max.getH(H_field); 
		     print_magnetic("initial",max,T,0,H_field[1],H_field[2],H_field[3]);
}


void initialize_spinmap_prevrun(Magnet & max,int argc,char * argv[]){
    cout<<"######################################################"<<endl;
    cout<<"Initialize Spin from Previous Run!! "<<endl;
    cout<<"#####################################################"<<endl;
    double T=max.getT();
    double H_temp[3];
    max.getH(H_temp);
    int MCSteps;
    cout<< "Step numbers to start from : ";
    if (argc>8) MCSteps = atoi(argv[9]);
    cout<<MCSteps;
    ostringstream oss,mc,H0,H1,H2;
    oss<<T; mc<<MCSteps; H0<<H_temp[0]; H1<<H_temp[1]; H2<<H_temp[2];
    string filename="mfile";
    filename=filename+"_"+oss.str()+"_"+H0.str()+H1.str()+H2.str()+"_"+mc.str()+".data"; 
    ifstream file(filename.c_str());
    
    int nTot=max.getnTot();
    string line;
    string temp_spin; 
    double spin,SpinI[3]; 
    if (file.is_open()){ 
         for(int j;j<nTot;j++)
         {
               getline(file,line);
               for (int i =0; i<4;i++)
                    {
                       size_t index=line.find("\t");
                       temp_spin=line.substr(0,index);
                       if (i>0)
                          {
                                SpinI[i-1]=atof(temp_spin.c_str());
                                cout<<SpinI[i-1]<<"\n";
                          }
                      line=line.substr(index+1);
                      max.writeIthSpin(j,SpinI);  
          } 
     }
    file.close();
    cout<<"\n\n\nSTARTING FROM "<<MCSteps<<"th TIME STEPS\n\n\n"<<endl;
    double H[3]; max.getH(H);
    print_magnetic("initial",max,T,0,H[0],H[1],H[2]);
  }
}

	// looks OK 7/20/14 MPS
	//Random Movement through the Sphere
void try_move(double * try_spin, double delta_max)
{
  double delta [DIMENSIONS];
    // choose a random vector, uniformly distribute inside unit sphere
  double s = 10.;
  while (s>1.0 || s==0.)
  {
	   for (int i=0;i<DIMENSIONS;i++) delta[i] = 2.*std_rand()-1.0;
	   s = mynorm(delta);
   }

	    for (int k=0;k<DIMENSIONS;k++) delta[k] *= delta_max/s; // dividing by s overweights large angular moves...
	    for (int k=0;k<DIMENSIONS;k++) try_spin[k]=try_spin[k]+delta[k];
	    s=mynorm(try_spin);
	    for (int k=0;k<DIMENSIONS;k++) try_spin[k]/=s;

}


bool MetropolisStep (Magnet& max, double delta_max) {
  int nTot=max.getnTot(); 
  double T=max.getT();
  double SpinI[3];
  int i= int(nTot*std_rand()) ;
//cout << i << " random site " << endl;
//Compute energy of initial configuration
  double e0 = dEnergydSpin( i,  max);
//Random move to spin i,j
  double try_spin[3];
// save old spin
  max.getIthSpin(i, SpinI);    
  for (int k=0;k<3;k++){try_spin[k]=SpinI[k];}
// Move spin direction and overwrite our spin database
  try_move(try_spin, delta_max);
  max.writeIthSpin(i,try_spin);
//Compute energy of trial configuration
  double e1=  dEnergydSpin( i,  max);
// ratio of Boltzmann factors
  long double ratexp =  -(e1-e0)/(T+1e-8);
// cout<<"DEBUGGING e0:  "<<e0 <<"e1  "<<e1<<endl; 
  if (std_rand() < expl(ratexp)) {
    return true; // leave as the modified dataset
    } 
  else{
	       max.writeIthSpin(i,SpinI); //return to the unmodified sataus
		return false;
    }
	    
}
        
bool RotateMagnet(Magnet& max, double delta_max)
{
     int nTot=max.getnTot();
     double T=max.getT();
     double SpinI[3],SpinJ[3],try_spin[3],SpinI_sph[3],SpinJ_sph[3],try_spin_sph[3];
     double theta,phi;
     int i=int(nTot*std_rand());
     double e0=0,e1=0; 
     for(int j;j<nTot;j++) e0+=energyPerSpin(j,max);//total energy before rotate
        
     max.getIthSpin(i,SpinI);
     for (int k=0;k<3;k++) try_spin[k]=SpinI[k];
     try_move(try_spin,delta_max);         
     conv_spherical(SpinI,SpinI_sph);
     conv_spherical(try_spin,try_spin_sph);
     theta=try_spin_sph[1]-SpinI_sph[1];
     phi=try_spin_sph[2]-SpinI_sph[2]; //get angle
     Magnet dummy(max);
     for( int j=0;j<nTot;j++)
     {
          dummy.getIthSpin(j,SpinJ);
          conv_spherical(SpinJ,SpinJ_sph);
          SpinJ_sph[1]+=theta;
          SpinJ_sph[2]+=phi;
          conv_cartesian(SpinJ_sph,SpinJ); 
          dummy.writeIthSpin(j,SpinJ);
     } //rotate all spins
     for (int j=0;j<nTot;j++) e1+=energyPerSpin(j,dummy); //total energy after rotate
          long double ratexp=-(e1-e0)/(T+1e-8);
          if (std_rand()<expl(ratexp))
           {
              double SpinTemp[3];
              for(int j=0;j<nTot;j++) {dummy.getIthSpin(j,SpinTemp); max.writeIthSpin(j,SpinTemp);}
              cout<<"ROTATE MAGNET!"<<endl;
              return true; //leave as the rotated dataset                 
           }
          else{
                  cout<<"DO NOT ROTATE MAGNET"<<endl;
                  return false;//return to the unmodified dataset  
              } 
} 


double oneMonteCarloStepPerSpin_2 (Magnet& max, const double delta_max) {
     int accepts = 0;
     for (int i = 0; i < max.getnTot(); i++)
     if (MetropolisStep(max,delta_max))
             ++accepts;
     acceptanceRatio = double(accepts)/double(max.getnTot());
     ++steps;
     return acceptanceRatio;
}

double oneMonteCarloStepPerSpin (Magnet& max, const double delta_max) {
     double T=max.getT();
     double SpinI[3],spinJ[3];
     int nTot=max.getnTot();
     double e0,e1;     
     double H[3];
     max.getH(H); 
     double K=max.getK();       
     int accepts=0;

     for (int ii = 0; ii < nTot; ii++)
       {
           bool Metropolis;
           int i=int(nTot*std_rand());
           //CORRECTED
           max.getIthSpin(i,SpinI);
           string element_i=max.getglName(i);
           string element_j; 
           //
           int nShellI = max.getIthNShell(i);
           int shellSizeI[nShellI];
           for (int j=0;j<nShellI;j++)
             shellSizeI[j] = max.getIJthShellSize(i,j);
           //CORRECTED
           double T=max.getT(); 
           //
           double energyI=0.0;
           double energyI_2=0.0;
           //**//
          //Random move to spin i j
          double try_spin[3];
          for (int k=0;k<3;k++){try_spin[k]=SpinI[k];}
          try_move(try_spin, delta_max);
          //max.writeIthSpin(i,try_spin); 
          //Compute energy of trial configuration
           for (int j=0;j<nShellI;j++)
             {
              for (int k=0;k<shellSizeI[j];k++)
                  {
                   int m = max.getIJKthNbr(i,j,k);
                   double J = max.getIJKthJ(i,j,k);
                   
                   //CORRECTED:.scaling factor;
                   max.getIthSpin(m,spinJ);
                   element_j=max.getglName(j);
                   
                   if (element_j=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_j=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;
                   if (element_i=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_i=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;  

                   // 
                   double nij[3];
                   nij[2]=spinJ[2]-SpinI[2];
                   nij[1]=spinJ[1]-SpinI[1];
                   nij[0]=spinJ[0]-SpinI[0];
                   double div = sqrt(mydot(nij,nij))+1e-10;
                   nij[0]=nij[0]/div;
                   nij[1]=nij[1]/div;
                   nij[2]=nij[2]/div;

                   double nij_2[3];
                   nij_2[2]=spinJ[2]-try_spin[2];
                   nij_2[1]=spinJ[1]-try_spin[1];
                   nij_2[0]=spinJ[0]-try_spin[0];
                   double div_2 = sqrt(mydot(nij,nij))+1e-10;
                   nij_2[0]=nij_2[0]/div_2;
                   nij_2[1]=nij_2[1]/div_2;
                   nij_2[2]=nij_2[2]/div_2;
                   // energyI += J * mydot(SpinI,spinJ) + K * mydot(SpinI,nij)*mydot(spinJ,nij);
                   // energyI_2 += J * mydot(try_spin,spinJ) + K * mydot(try_spin,nij_2)*mydot(spinJ,nij_2);
                   // CORRECTED
                   energyI += J * mydot(SpinI,spinJ);
                   energyI_2 += J * mydot(try_spin,spinJ);                
                  }
             }
           e0=-energyI+mydot(H,SpinI);  
           e1=-energyI_2+mydot(H,try_spin);
 // Use KKR sign convention for Heisenberg J
 //  // Return the energy difference between having
 //   // spinI at site i versus having (0,0,0) at site i.
 //    return -energyI + mydot(H,spinI);
// ratio of Boltzmann factors
          long double ratexp =  -(e1-e0)/(T+1e-8);       
          if (std_rand() < expl(ratexp))
                  {  
                     max.writeIthSpin(i,try_spin);
                     Metropolis= true; //  modified dataset
                     ++accepts;   
                  }
              else{
                     //max.writeIthSpin(i,SpinI); //return to the unmodified sataus
                     Metropolis=false;
                  }
        }       
     acceptanceRatio = double(accepts)/double(nTot);
     ++steps; 
     return acceptanceRatio;
}

double magnetizationPerSpin(Magnet& max)
{
    double spinI[3];
    int nTot=max.getnTot();
	  
    double ssum[DIMENSIONS];
    for (int k=0;k<DIMENSIONS;k++) ssum[k]=0.0;
	  for (int i=0;i<nTot;i++)
	     {
		max.getIthSpin(i,spinI);
		for (int k=0;k<DIMENSIONS;k++)
		   ssum[k] +=  spinI[k];
	     }

    return mynorm(ssum)/double(nTot);
}

double magnetization(Magnet& max)
{
     double spinI[3];
     int nTot=max.getnTot();
	  
     double ssum[DIMENSIONS];
     for (int k=0;k<DIMENSIONS;k++) ssum[k]=0.0;
	  for (int i=0;i<nTot;i++)
	     {
		max.getIthSpin(i,spinI);
		for (int k=0;k<DIMENSIONS;k++)
		   ssum[k] +=  spinI[k];
	     }

     return mynorm(ssum);
}

double energyPerSpin(int i, Magnet& max)
{
	 // energy for i'th spin, 
	 // Heisenberg terms are split 50/50 between site i and neighbor sites
	 // Sum this function over all sites i to get the total energy.
          string element_i,element_j;
	  double spinI[3], spinJ[3];
         //CORRECTED
          element_i=max.getglName(i);
          double T=max.getT();
          //; 
	  int nShellI = max.getIthNShell(i);
	  int shellSizeI[nShellI];
	  for (int j=0;j<nShellI;j++)
	      shellSizeI[j] = max.getIJthShellSize(i,j);

	  max.getIthSpin(i,spinI);
	  
	 double energyI = 0.;
	 for (int j=0;j<nShellI;j++)
	     {
	      for (int k=0;k<shellSizeI[j];k++)
		  {
		   int m = max.getIJKthNbr(i,j,k);
		   double J = max.getIJKthJ(i,j,k);
		   max.getIthSpin(m,spinJ);
		   //CORRECTED
	           element_j=max.getglName(j);

                   if (element_j=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_j=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;
                   if (element_i=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_i=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;

                   // 
                   //
		   double K;
		   K=max.getK();           

		   double nij[3];
		   double coordI[3],coordM[3];
		   max.getglCoord(i,coordI);
		   max.getglCoord(m,coordM);
		   
		  nij[2]=coordM[2]-coordI[2];
		   nij[1]=coordM[1]-coordI[1];
		   nij[0]=coordM[0]-coordI[0];
		   double div = sqrt(mydot(nij,nij));
		   nij[0]=nij[0]/div;
		   nij[1]=nij[1]/div;
		   nij[2]=nij[2]/div;


//		   energyI += J * mydot(spinI,spinJ) + K * mydot(spinI,nij)*mydot(spinJ,nij);
                   energyI += J * mydot(spinI,spinJ); 
		  }
	     }

	 double H[3];
	 max.getH(H);

	 // Use KKR sign convention for Heisenberg J
	 return -energyI/2. + mydot(H,spinI);
	 
}

double dEnergydSpin(int i, Magnet& max)
{
	 // Change in energy relative to setting the i'th spin to zero.
	 // Use this function to compute the energy difference for altering a single spin.
	 double H[3];
         double T=max.getT();
	 max.getH(H);
	 double spinI[3], spinJ[3];
         //CORRECTED
         string element_i,element_j;
         element_i=max.getglName(i);
         //  
	 int nShellI = max.getIthNShell(i);
	 int shellSizeI[nShellI];
	 for (int j=0;j<nShellI;j++)
	     shellSizeI[j] = max.getIJthShellSize(i,j);
	/*
	 int** shellNbrI;
	 shellNbrI = new int* [nShellI];
	 for (int j=0;j<nShellI;j++)
	     {
	      shellNbrI[j] = new int [shellSizeI[j]];
	      for (int k=0;k<shellSizeI[j];k++)
		  shellNbrI[j][k] = max.getIJKthNbr(i,j,k);
	     }
	*/
	 max.getIthSpin(i,spinI);
	 double energyI = 0.;
	 for (int j=0;j<nShellI;j++)
	     {
	      for (int k=0;k<shellSizeI[j];k++)
		  {
		   int m = max.getIJKthNbr(i,j,k);
		   double J = max.getIJKthJ(i,j,k);

                    //CORRECTED
		   max.getIthSpin(m,spinJ);
                   element_j=max.getglName(j);


                   if (element_j=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_j=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;
                   if (element_i=="Co") J=J*(1.0 - 0.000923628358514078*T);
                   else if (element_i=="Pt") J=J*(1.0 - 0.008880410000955151*T);
                   else cout<<"Wrong global Name!"<<endl;
                   //
		   double K;
		   K=max.getK();          
                   double nij[3];
                   nij[2]=spinJ[2]-spinI[2];
                   nij[1]=spinJ[1]-spinI[1];
                   nij[0]=spinJ[0]-spinI[0];
         	   double div = sqrt(mydot(nij,nij))+1e-10;
                   nij[0]=nij[0]/div; 
	           nij[1]=nij[1]/div;
         	   nij[2]=nij[2]/div;
                   
                   if (0.0==K) energyI += J * mydot(spinI,spinJ);  
                   else energyI += J * mydot(spinI,spinJ) + K * mydot(spinI,nij)*mydot(spinJ,nij);

                  }     
             }
 // Use KKR sign convention for Heisenberg J
 // Return the energy difference between having
 // spinI at site i versus having (0,0,0) at site i.
 return -energyI + mydot(H,spinI);
 
}

void PrintSpin(Magnet & max)
{
     double nTot,T,K;
     double H[3],SpinI[3];
     nTot=max.getnTot();
     T=max.getT(); 
     max.getH(H);
     K=max.getK();  
     cout<<"======SPONTANEOUS SPIN SNAPSHOT========"<<endl;
     for(int i=0;i<nTot;i++)
          {
              max.getIthSpin(i,SpinI); 
              cout<<SpinI[0]<<",\t"<<SpinI[1]<<",\t"<<SpinI[2]<<"|\t\t\t";
          }
     cout<<"======================================="<<endl;   
     cout<<"T: "<<T<<'\t'<<"K: "<<K<<'\t'<<"nTot: "<<nTot<<endl;
     cout<<"H: "<<H[0]<<'\t'<<H[1]<<'\t'<<H[2]<<endl;
}
