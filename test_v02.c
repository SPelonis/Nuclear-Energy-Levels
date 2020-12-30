#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXREP 100
 
double x_t[10000], y_t[10000], z_t[10000], E_t[10000], value[10000];    
    
/* Hear we define the initial values of the parameters below */
double y = 0;
double z = 0.1;
double E = -50;

double charge = 8;
double A = 16;
int l_max = 2;

double G(double E_i, int l_i);

double f(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double f1(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double f2(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double f3(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double f4(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double g(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double g1(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double g2(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double g3(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double g4(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i);
double R(double A);


int main()
{
int flag = 0;

double E1, E2, x;
double E1_;
double y1_, y2_, y_;
double E_keep = 0;

int i = 1;
int j = 0;
int l = 0; 

x = R(A) + 3; /*This is x_start*/

/*E1 and E2 are the energy limits that we will use.*/
printf(" Give me the value of E1:");
scanf("%lf", &E1);

printf(" Give me the value of E2:");
scanf("%lf", &E2);

printf("\n Initial Parameters: ");
printf("\n x = %lf", x);
printf("\t y = %lf", y);
printf("\t z = %lf\n", z); 

while(l<=l_max){
    
printf("\n l = %d \n",l);

    while(E1<E2){
    flag = 0;
    E1_ = E1 + 1; /*The values E1 and E1_ are needed to run every loop.*/
    E_keep = E1;
    printf("\n Energy: %lf MeV", E1);
    y1_ = G(E1,l);
    y2_ = G(E1_,l);
        if(y1_*y2_ > 0)
        {
            printf("\n The method cannot be applied.");
        
        }
        else
        {
            printf("\n     i         E[MeV]              y\n");
            printf("\n------------------------------------\n");
 
            do
            {
                E = 0.5 * (E1 + E1_);
                y_ = G(E,l);
                printf("%6d %14.8f %14.8f \n", i, E, y_);
 
                if (y_ == 0){
            
                    flag = 1;
                    break;
                }
                else if (y_*y1_ > 0){
                    E1 = E;
                }
                else{
                    E1_ = E;
                }
                y1_ = G(E1,l);
 
                if (i >= MAXREP){
            
                    printf("\n Max number of iterations exited.\n ");
                    flag = 2;
                    break;
                }
 
                i++;
 
            }while (flag == 0);
        }
        /*Here we set a high threshold for the Eigenvalues.*/
        if (flag == 2 && fabs(y_) < 5){                    
            printf("\n The Energy Eigenvalue is E = %.14f\n", E);
            flag = 0;
        }
        else if (flag == 2 && fabs(y_) > 5){
        printf("\n The desired accuracy was not achieved.\n");
        }

    E1 = E_keep + 1;
    i = 1;
    j++;
    }
E1 = E1 - j; 
j = 0;
l++;
}
return 0;

} 

double G(double E_i, int l_i)
{
double x = R(A) + 3;
double a,b;
int N = 100;
int i = 0;
int j = 0;
x_t[0] = x; 
y_t[0] = y;
z_t[0] = z;

double dt;
dt = (x - 0.1)/N;

int flag = 0;

E_t[0] = E_i;
   
    
while(i < N){
    
    y_t[i+1] = y_t[i] + ( f1(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) + 2*f2(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) +2*f3(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) + f4(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) )/6;    
    
    z_t[i+1] = z_t[i] + ( g1(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) + 2*g2(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) +2*g3(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) + g4(x_t[i],y_t[i],z_t[i],dt,E_i,l_i,R(A),charge) )/6; 
    
    x_t[i+1] = x_t[i] - dt;
    i++;
    }
    
return y_t[i]; 
}


/*They need to be changed according to the second ODE. */
/*This is the 1st differential*/
double f(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i) 
{
    return z_i;
}
/*This is the 2nd differential*/
double g(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
if(x_i>R_i){
    return    ( 20.75*l_i*(l_i+1)/pow(x_i,2) + 1.44*(charge_i)/x_i - 43/(pow(2.72,(x_i - R_i)/0.6)+1) - Q)*y_i;
}
if(x_i<=R_i){
    return    ( 20.75*l_i*(l_i+1)/pow(x_i,2) + (1.44/2)*( 3 - pow(x_i,2)/pow(R_i,2) )*(charge_i)/R_i - 43/(pow(2.72,(x_i - R_i)/0.6)+1) - Q)*y_i;
}
}


double f1(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*f(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i);
}
double g1(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*g(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i);
}


double f2(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*f(x_i + 0.5*dt_i,y_i + 0.5*f1(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + 0.5*g1(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}
double g2(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*g(x_i + 0.5*dt_i,y_i + 0.5*f1(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + 0.5*g1(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}


double f3(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*f(x_i + 0.5*dt_i,y_i + 0.5*f2(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + 0.5*g2(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}
double g3(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*g(x_i + 0.5*dt_i,y_i + 0.5*f2(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + 0.5*g2(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}


double f4(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*f(x_i + dt_i,y_i + f3(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + g3(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}
double g4(double x_i,double y_i,double z_i,double dt_i, double E_i, double l_i, double R_i, double charge_i)
{
    return dt_i*g(x_i + dt_i,y_i + f3(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),z_i + g3(x_i,y_i,z_i,dt_i,E_i,l_i,R_i,charge_i),dt_i,E_i,l_i,R_i,charge_i);
}

double R(double A)
{
    return 1.2 * pow(A,0.333);
}



























