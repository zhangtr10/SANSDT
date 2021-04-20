#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <random>
#include <complex>
#include <iostream>
//#include <boost/math/distributions/triangular.hpp>
#define pow2(x) ((x)*(x))
#define PI 3.141592653589793238462643383

using Eigen::MatrixXd;
using namespace std;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using namespace Eigen;
//typedef complex<double> dcomp;

double triangulardist(double, double, double,double);
double gaussian(double sigma, double mean);
double ran2( void );

void read_input_file(char*);
long int idum ;
//initialization
int nlayer,   nNeutrons, j;
double H, d, wg, AL, *L, dd, *H_n, *wg_n;

/*
const int j=25;
const int nlayer = 9;
int nNeutrons=10000;
const double d =380, H = 360, wg = 170, L = 75;*/

/*double dd=13000;// for intensity calculation
double H_n[ nlayer ] ={5,5,5,30,100,100,100,15,75};
double wg_n[nlayer] = {330,280,185,180,175,170,165,145,170} ;*/

int hh=100;

complex<double> I = complex<double>( 0.0 , 1.0 ) ;
complex<double> F_x(int, double, double);


complex<double> F_x(int n, double x, double y){
    if(n==0){
      return y-x;
    }else{
      return -d/(2*PI*n*I)*exp(-2*PI/d*n*I*y)+d/(2*PI*n*I)*exp(-2*PI/d*n*I*x);
    }
    
}

int main(int argc, char* argv[]){
  if (argc < 2) { // We expect 3 arguments: the program name, the source path and the destination path
        std::cerr << "Usage: " << argv[0] << "SOURCE DESTINATION" << std::endl;
        return 1;
    }
	/*dcomp i;
	I = -1;
	I = sqrt(I); */

  read_input_file(argv[1]);
  printf("finish reading");fflush( stdout ) ;
  double g=2*PI/d, ww=d-wg, f=ww/d, rho_p = 0, rho_al = 0.0003877, rho_s=0.0003523, rho_air=0;
  //different layers, FFTW 
/////////////////////////////////////////////////////////////////////////////////////////////
  complex<double> rho_n[nlayer][4*j+1];
  
  int count_n=0;
  for(int m=0; m< nlayer; m++){
    count_n += H_n[m];

    if(count_n<= AL){
      for(int n=-2*j; n<=2*j; n++){

        rho_n[m][n+2*j]=1/d * (rho_al *F_x(n,-d/2, -wg_n[m]/2+L[m]) + rho_p *F_x(n,-wg_n[m]/2+L[m], wg_n[m]/2-L[m]) + rho_al *F_x(n, wg_n[m]/2-L[m],d/2)); 
    
      }
    }
    else if(count_n > H){   
      for(int n=-2*j; n<=2*j; n++){

        rho_n[m][n+2*j]=1/d * (rho_al * F_x(n,-wg_n[m]/2, wg_n[m]/2) + rho_s *F_x(n,-d/2, -wg_n[m]/2) + rho_s* F_x(n, wg_n[m]/2,d/2));
      }
    } else{
      for(int n=-2*j; n<=2*j; n++){

    
        rho_n[m][n+2*j]=1/d * ( (rho_al *F_x(n,-d/2, -wg_n[m]/2+L[m]) + rho_p *F_x(n,-wg_n[m]/2+L[m], wg_n[m]/2-L[m]) + rho_al *F_x(n, wg_n[m]/2-L[m],d/2)) +((-1)*rho_al + rho_s)*F_x(n,-d/2, -wg_n[m]/2)+((-1)*rho_al + rho_s)*F_x(n, wg_n[m]/2,d/2));
      }   
    }
   
  }

//angle set up
/////////////////////////////////////////////////////////////////////////////////////////////

    //Determining incident angles and wavelengths and their distributions
   // *Incident Angle, 90 = perpendicular to sample surface, normal distribution*)
    double avgsita = 90 * PI / 180, deltasita = 0.18* PI / 180;
   // std::default_random_engine generator;
    //std::normal_distribution<double> distribution(avgtheta,deltatheta/8);
    
    double *sita;
    sita = ( double* ) calloc(  nNeutrons, sizeof( double ) );

    for(int i=0; i< nNeutrons; i++){
      sita[i] = gaussian(deltasita,avgsita);
      

    }
    
//(*Wavelength, triangular distribution*)
    double avggamma = 0.6, deltagamma = avggamma*0.138;
    double *gamma;
    gamma = ( double* ) calloc(  nNeutrons, sizeof( double ) );
    double R;
    
    for(int i=0; i< nNeutrons; i++){
      R = ran2();
    	gamma[i] = triangulardist(avggamma - deltagamma , avggamma, avggamma + deltagamma, R) ;
    }

    
//(*Azimuthal Angle, 0-360 deg (only use 0-180 the mirror for second half of image), uniform distribution*)
    double *thi;
    thi=( double* ) calloc(  nNeutrons, sizeof( double ) );
    double thimin=0, thimax=PI;
    
    for(int i=0; i< nNeutrons; i++){
    	thi[i] = ran2()*thimax;
    }
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
/////////////////////main loop/////////////////////////
///////////////////////////////////////////////////////
    double kx, ky, kz;
    
    complex<double> **Intensin;
    Intensin = (complex<double>**)calloc( nNeutrons , sizeof( complex<double>*) ) ;
    for(int n=0; n<nNeutrons; n++){
      Intensin[n] = (complex<double>*)calloc( 2*j + 1 , sizeof( complex<double>) ) ;
    }

    double **xparr, **yparr;
    xparr = (double**)calloc( nNeutrons , sizeof( double*) ) ;
    yparr = (double**)calloc( nNeutrons , sizeof( double*) ) ;
    for(int n=0; n<nNeutrons; n++){
      xparr[n] = (double*)calloc( 2*j + 1 , sizeof( double) ) ;
      yparr[n] = (double*)calloc( 2*j + 1 , sizeof( double) ) ;
    }

    double INnm[hh*2+1][hh*2+1];

    for(int n=-hh; n<=hh; n++){

       for(int m=-hh; m<=hh; m++){
         INnm[n+hh][m+hh]=0;

    }
  }


  //#pragma omp parallel for
    for(int i = 0; i < nNeutrons; i ++){

    	kx=2*PI/gamma[i]*cos(thi[i])*cos(sita[i]);
    	ky=2*PI/gamma[i]*sin(thi[i])*cos(sita[i]);
    	kz=2*PI/gamma[i]*sin(sita[i]);

      MatrixXcd M_n( 2*j+1, 2*j+1);
  
      complex<double> pz_n[nlayer][2*j + 1];
      complex<double> pz0n[2*j + 1];
      complex<double> pzfn[2*j + 1];

      complex<double> b_nm[nlayer][2*j + 1][2*j + 1];

      for (int h=0; h<nlayer; h++){

        for(int n=-j; n<=j; n++){
          for(int m=-j; m<=j; m++){
            if(n==m){
              M_n(n+j, m+j)= -pow2( ky + n * g )+ pow2(ky) + pow2(kz) -4*PI*rho_n[h][2*j];
            }else{
              M_n(n+j, m+j)= -4 * PI * rho_n[h][ n-m + 2*j ];

            }
          }

        }
        ComplexEigenSolver<MatrixXcd> es(M_n);

        for(int n=-j; n<=j; n++){

          pz_n[h][n+j] = sqrt(es.eigenvalues()[ n + j].real());
          //printf("kzn %lf\n", kzn[ n+j]);

          pzfn[n+j] = sqrt( pow2(ky) + pow2(kz) - 4 * PI * rho_s -pow2(ky + n * g));
          pz0n[n+j] = sqrt( pow2(ky) + pow2(kz) -pow2(ky + n * g));
          //printf("ptn p0n %lf %lf\n", ptn[ n+j],p0n[ n+j] );
          for(int m=-j; m<=j; m++){

            b_nm[h][n+j][m+j]=es.eigenvectors().col(n + j)[m + j].real();
  


          }

        }



      }

       //cout << "Here is matrix, M:" << endl << M1 << endl << endl;
       

       //cout << "The eigenvalues of M are:" << endl << es1.eigenvalues() << endl;
       //cout << "The matrix of eigenvectors, is:" << endl << es1.eigenvectors() << endl << endl;

 /////////////////////////////////////////////////////////////////////////////////////////////    
       // initialize the r_n and t_n at different unknown number , in different euqations 
     /* complex<double> ****r_n, ****t_n;
       r_n = (complex<double>****)calloc( nlayer + 1 , sizeof( complex<double>*) ) ;
       t_n = (complex<double>****)calloc( nlayer + 1 , sizeof( complex<double>*) ) ;

       for(int n=0; n<nlayer+1; n++){
         r_n[n] = (complex<double>***)calloc( 2*(nlayer + 1) , sizeof( complex<double>) ) ;
         t_n[n] = (complex<double>***)calloc( 2*(nlayer + 1) , sizeof( complex<double>) ) ;
         for(int m=0; m<2*(nlayer + 1); m++){
           r_n[n][m] = (complex<double>**)calloc( 2*j + 1 , sizeof( complex<double>) ) ;
           t_n[n][m] = (complex<double>**)calloc( 2*j + 1 , sizeof( complex<double>) ) ;

           for(int h=0; j<2*j+1; h++){
             r_n[n][m][h]= (complex<double>*)calloc( 2*j + 1 , sizeof( complex<double>) ) ;
             t_n[n][m][h]=(complex<double>*)calloc( 2*j + 1 , sizeof( complex<double>) ) ;
             for(int a=0; a<2*j+1;a++){
                r_n[n][m][h][a]=0;
                t_n[n][m][h][a]=0;
             }

           }
         }
       }*/
/////////////////////////////////////////////////////////////////////////////////////////////




        //cout << "Here is matrix, bnm:" << endl << bnm << endl << endl;


        MatrixXcd arr((2*j+1)*(2*nlayer+2),(2*j+1)*(2*nlayer+2));
        MatrixXcd bb((2*j+1)*(2*nlayer+2),1);
        for(int m=0; m< (2*j+1)*(2*nlayer+2); m++){
          bb(m,0)=0;
          for(int n=0; n< (2*j+1)*(2*nlayer+2); n++){
            arr(m,n)=0;
          }
        }

        //cout << "Here is matrix, arr:" << endl << arr << endl << endl;
       /*
        int cc=0;
        for(int m=-j; m<=j; m++){
          r_n[0][0][m+j][m+j]=1;
          r_n[0][1][m+j][m+j]=pz0n[m+j];

          t_n[ nlayer ][2*nlayer][m+j][m+j]=-1;
          t_n[ nlayer ][2*nlayer+1][m+j][m+j]=pzfn[m+j]

        }
       for(int h=1; h<nlayer+1; h++){

         //for(int e=0; e<2*(nlayer + 1); e++){
         

           for(int m=-j; m<=j; m++){

 
             for(int n=-j; n<=j;n++){

                r_n[h][cc][m][n]=-exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j];
                r_n[h][cc+1][m][n]=-exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];
                r_n[h][cc+2][m][n]=b_nm[h-1][n+j][m+j];
                r_n[h][cc+3][m][n]=b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];

                t_n[h-1][cc][m][n] = -b_nm[h-1][n+j][m+j];
                t_n[h-1][cc+1][m][n] = b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];
                t_n[h-1][cc+2][m][n] = exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j];
                t_n[h-1][cc+3][m][n] = -exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];

             }

           //}
         }
          cc = h+1;

       }
       */
/////////////////////////////////////////////////////////////////////////////////////////////
// start to allocate the unknowns on the array
/////////////////////////////////////////////////////////////////////////////////////////////
          int cc=0;
          for(int m=-j; m<=j; m++){
            arr(m+j,m+j)=1;
            arr( m + 3*j+1 ,m + j) =pz0n[m+j];

           arr( m + (2*j+1)*(2*nlayer+2) -3*j-2, m + (2*j+1)*(2*nlayer+2) -j-1)=-1;
           arr( m + (2*j+1)*(2*nlayer+2) -j-1, m + (2*j+1)*(2*nlayer+2) -j-1) =pzfn[m+j];

          }
          for(int h=1; h<nlayer+1; h++){
            //for(int e=0; e<2*(nlayer + 1); e++ ){
              for(int m=-j; m<=j; m++){
                for(int n=-j; n<=j; n++){ 
                  arr( cc*2*j + cc + m+j, h*2*j+h+j+n)=-exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j];
                  arr( (cc+1)*2*j + cc+1 + m+j, h*2*j+h+j+n)=-exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];
                  arr( (cc+2)*2*j + cc+2 + m+j, h*2*j+h+j+n)=b_nm[h-1][n+j][m+j];
                  arr( (cc+3)*2*j + cc+3 + m+j, h*2*j+h+j+n)=b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];


                  arr( cc*2*j + cc + m+j,  (2*nlayer+1)*j+nlayer+ h*2*j+h+n)=-b_nm[h-1][n+j][m+j];
                  arr( (cc+1)*2*j + cc + 1 + m+j,  (2*nlayer+1)*j+nlayer+ h*2*j+h+n)=b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];
                  arr( (cc+2)*2*j + cc + 2 + m+j,  (2*nlayer+1)*j+nlayer+ h*2*j+h+n)=exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j];
                  arr( (cc+3)*2*j + cc + 3 + m+j,  (2*nlayer+1)*j+nlayer+ h*2*j+h+n)= -exp(I * pz_n[h-1][n+j] * H_n[h-1]) * b_nm[h-1][n+j][m+j] * pz_n[h-1][n+j];

                }
 
              }
            cc = cc+2;   
           // }
          }

           bb(j, 0) = -1;
           bb(3*j+1,0) = kz;

           //or colPivHouseholderQr(): high accuracy; householderQr():faster
           //cout << "Here is matrix, arr:" << endl << arr << endl << endl;
           MatrixXcd solv = arr.householderQr().solve(bb);
           //cout << "The solution is:\n" << solv << endl;

           MatrixXcd RFin(2*j+1,1);
           MatrixXcd TFin(2*j+1,1);
           MatrixXcd Rfin(2*j+1,1);
           MatrixXcd Tfin(2*j+1,1);

           MatrixXcd thtin(2*j+1,1);
           MatrixXcd phitin(2*j+1,1);

 
           for(int n=-j; n<=j; n++){

              RFin(n+j, 0) = solv(n+j,0);
              TFin(n+j, 0) = solv( n + (2*j+1)*(2*nlayer+2) -j-1,0);
              thtin(n+j, 0)= atan( sqrt( pow2(kx) + pow2(ky + n * g))/pzfn[n+j] );

              if(kx!=0){
                phitin(n+j, 0)=atan((ky + n * g)/kx);
              }else{
                phitin(n+j, 0)=PI/2;
              }

              if(pz0n[n+j].imag()!=0){
                Rfin(n+j, 0)=0;
              }else{
                Rfin(n+j, 0)=RFin(n+j, 0);
              }
              
              if(pzfn[n+j].imag()!=0){
                Tfin(n+j,0)=0;
              }else{
                Tfin(n+j,0)= TFin(n+j, 0);
              }


           }
            //(*calculate intensity*)
           for(int n=-j; n<=j; n++){
              Intensin[i][n+j] = pzfn[n+j]/kz * pow2( abs( Tfin(n+j,0) ) );

              double xx = (dd*tan(thtin(n+j, 0))*cos(phitin(n+j, 0))).real();
              int xp = floor(xx/5 +0.5);
              int xpneg = -floor(xx/5 +0.5);
              
              double yy = (dd* tan(thtin(n+j, 0)) * sin(phitin(n+j, 0))).real();
              int yp = floor(yy/5 +0.5);
              int ypneg = -floor(yy/5 +0.5);

              // not quite sure with the syntax in mathematica
              if(abs(xp)>hh||abs(yp)>hh){
                continue;
              }
              
              else if(abs(xp)==0){
                INnm[xp + hh][yp + hh] += (Intensin[i][n+j]).real();

              }else{
                INnm[xp + hh][yp + hh] += (Intensin[i][n+j]).real();
                INnm[xpneg + hh][yp + hh] += (Intensin[i][n+j]).real();
                if(INnm[xp + hh][yp + hh]!=INnm[xpneg + hh][yp + hh]){
                  printf("%lf %lf\n", INnm[xp + 100][yp + 100], INnm[xpneg + 100][yp + 100]);
                }
                
              }
               ///////////////////////////////////////////////////////////
              xparr[i][n+j] = xp;
              yparr[i][n+j] = yp;
           }
       	//}
           //(*organizing reflectance and transmission coefficients*)
       }
      FILE *fid;
      fid = fopen(argv[2], "w");
      double highNumx=0;
      double highNumy=0;
      double minNumx=xparr[0][0];
      double minNumy=yparr[0][0];

      for(int n=0; n< nNeutrons; n++){
         for(int m=-j; m<= j; m++){
           if (xparr[n][m+j] > highNumx)
              highNumx = xparr[n][m+j];
           
           if (xparr[n][m+j] < minNumx)
              minNumx = xparr[n][m+j];
           
           if (yparr[n][m+j] > highNumy)
              highNumy = yparr[n][m+j];
           
           if (yparr[n][m+j] < minNumy)
              minNumy = yparr[n][m+j];

         }

      }
      printf("xpmin, xpmax, ypmin, ypmax%lf %lf %lf %lf\n",minNumx, highNumx, minNumy, highNumy);


     for(int n=-hh; n<=hh; n++){

       for(int m=-hh; m<=hh; m++){

         if(m==hh){
          fprintf (fid, "{%d, %d, %lf}\n",n, m, INnm[n+hh][m+hh] );

         }else{
          fprintf (fid, "{%d, %d, %lf} ",n, m, INnm[n+hh][m+hh] );
         }
         

       }
     }

    return 0;

}

double triangulardist(double  Min, double  Mode, double  Max, double R){
    
          //between 0.0 and 1.0 gaussian
    
    //    Triangular
    if ( R == (( Mode -  Min) / ( Max -  Min)))
    {
        return  Mode;
    }
    else if ( R < (( Mode -  Min) / ( Max -  Min)))
    {
        return  Min + sqrt( R * ( Max -  Min) * ( Mode -  Min));
    }
    else
    {
        return  Max - sqrt((1 -  R) * ( Max -  Min) * ( Max -  Mode));
    }


}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

double ran2 (void)
{
  int j;
  long int k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  extern long idum;
  
  if (idum <= 0) {
    if (-(idum) < 1) idum = 1;
    else idum = -(idum);
    idum2 = (idum);
    for (j = NTAB + 7; j >=0; j--) {
      k = (idum) / IQ1;
      idum = IA1 * (idum - k * IQ1) - k * IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  
  k = (idum) / IQ1;
  idum = IA1 * (idum - k * IQ1) - k * IR1;
  if (idum < 0) idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
}

double gaussian ( double sigma, double mean){

    double r=2.0, v1, v2,l;
    while (r > 1.0){
   	v1=double(ran2()*2-1);
	v2=double(ran2()*2-1);
	r=v1*v1+v2*v2;
    }
    l=v1*sqrt(-2*log(r)/r);
    l=mean+sigma*l;
    return l;

}
void read_input_file(char* file){
  FILE *inp ;
  int i;
  double d1 ;

  char tt[80] ;

  inp = fopen( file , "r" ) ;


  fscanf( inp , "%d" , &nlayer ) ;
  fgets( tt , 80 , inp ) ;
  printf("nalyer: %d \n" , nlayer ) ;fflush( stdout ) ;

  fscanf( inp , "%d" , &nNeutrons ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &j ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &H ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &d ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &wg ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &AL ) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &dd ) ;
  fgets( tt , 80 , inp ) ;

  printf("dd: %lf \n" , dd ) ;fflush( stdout ) ;
  
   H_n= ( double* ) calloc( nlayer , sizeof( double ) ) ;

  for ( i=0 ; i<nlayer ; i++ ){
    fscanf( inp , "%lf" , &H_n[i] ) ; 
    printf("height: %lf \n" , H_n[i] ) ;fflush( stdout ) ;
  }
  fgets( tt , 80 , inp ) ;


  wg_n= ( double* ) calloc( nlayer , sizeof( double ) ) ;
  for ( i=0 ; i<nlayer ; i++ ){
    fscanf( inp , "%lf" , &wg_n[i] ) ; 
    printf("width: %lf\n" ,wg_n[i] ) ;fflush( stdout ) ;
  }
  printf("finish_pre_pre" ) ;fflush( stdout ) ;
  fgets( tt , 80 , inp ) ;

  L= ( double* ) calloc( nlayer , sizeof( double ) ) ;
  printf("finish_pre" ) ;fflush( stdout ) ;
  for ( i=0 ; i<nlayer ; i++ ){
    fscanf( inp , "%lf" , &L[i] ) ; 
    printf("al thickness: %lf\n" ,L[i] ) ;fflush( stdout ) ;
  }
  fgets( tt , 80 , inp ) ;
  printf("finish" ) ;fflush( stdout ) ;
}

/*

double KroneckerDelta(int x, int y){

	if (x==y){
		return 1;
	}else{
		return 0
	}
}

double triangulardist(double  Min,double  Mode, double  Max){
    double  R=0.0;
    
     R = r.NextDouble();       //between 0.0 and 1.0 gaussian
    //    Triangular
    if ( R == (( Mode -  Min) / ( Max -  Min)))
    {
        return  Mode;
    }
    else if ( R < (( Mode -  Min) / ( Max -  Min)))
    {
        return  Min + Math.Sqrt( R * ( Max -  Min) * ( Mode -  Min));
    else
    {
        return  Max - Math.Sqrt((1 -  R) * ( Max -  Min) * ( Max -  Mode));
    }
}
*/


















