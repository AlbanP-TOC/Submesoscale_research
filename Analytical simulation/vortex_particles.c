
#include "view.h"
#include "run.h"
#include "tracer-particles.h"
#include "stokes-particles.h"
#include "write_class_particles.h"
#include "scatter2.h"
#include "poisson.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define LEVEL 9 
#define alph 1.25643// Lamb-Oseen vortex constant
#define grp_n 2 //Number of interial particles groups

double rc = 42000.;
double gmma = 92600.;

// Particle characteristic time and color
double tparticles[grp_n] = {100000.,1000000.};
double color_particles[grp_n][3] = {{0.004, 0.796, 0.882},{0.011, 0.564, 0.631}};//,{0.011, 0.404, 0.45}} ;

//Stokes drift parameters
double lmda = 100.;
double a = 2;
double angle = 0.;

//Simulation time parameters
double Tc = 10000.;
double time_sim = 10;
double time_div = 300.;

//particle objects (physics and data export)
Particles inertial[grp_n];
Particles flow;
TextWriter lpwriters[grp_n+1];

// We do not use the values of these fields:
face vector mu[];
scalar rho[], psi[], omega[];
vector u[];
double nu = 1e-6;
double pos;
double re = 2000.; 

void iso_contour (scalar s, double isoval){
  vertex scalar vs[];
  scalar il[];
  boundary ({s}); // Just in case ... 
  foreach_vertex()
    vs[] = interpolate (s, x , y) - isoval;
  boundary ({vs});
  fractions (vs, il);
  boundary ({il});
  draw_vof("il", lw = 3);
}
bool isOut(double x, double y, double eps_) {
  bool isOut = false;
  if ((x + eps_) > X0 + L0 || (x - eps_) < X0 || 
      (y + eps_) > Y0 + L0 || (y - eps_) < Y0)
    isOut = true;
  return isOut;
}



double lamb_oseen (double r, double gma) {
  //return the tangent velocity at distance r from the vortex center (gma = +/- gmma depending on the vortex rotation)
  return gma / (2. * pi * r) * (1 - exp(-sq(r) / sq(rc)));
}

event double_lamb_oseen (i = 0) {
  //define the velocity field induced by four Lamb-Oseen vortices
  foreach () {
    double r0 = sqrt( sq(x - pos) + sq(y - pos) );
    double ut0 = lamb_oseen (r0, -gmma);
    double th0 = atan2 (y - pos, x - pos);
    double r1 = sqrt( sq(x + pos) + sq(y - pos) );
    double ut1 = lamb_oseen (r1, gmma);
    double th1 = atan2 (y - pos, x + pos);
    double r2 = sqrt( sq(x + pos) + sq(y + pos) );
    double ut2 = lamb_oseen (r2, -gmma);
    double th2 = atan2 (y + pos, x + pos);
    double r3 = sqrt( sq(x - pos) + sq(y + pos) );
    double ut3 = lamb_oseen (r3, gmma);
    double th3 = atan2 (y + pos, x - pos);
    
    u.x[] = - ut0 * sin(th0) - ut1 * sin(th1) - ut2 * sin(th2) - ut3 * sin(th3) ;
    u.y[] =   ut0 * cos(th0) + ut1 * cos(th1) + ut2 * cos(th2) + ut3 * cos(th3);
  }
}

event compute_vorticity (i = 0) {
  //define the vorticity field
  double coef = gmma / (pi * sq(rc));
  foreach () {
    double r0 ,r1 ,r2, r3;
    r0 = sqrt( sq(x - pos) + sq(y - pos) );
    r1 = sqrt( sq(x + pos) + sq(y - pos) );
    r2 = sqrt( sq(x + pos) + sq(y + pos) );
    r3 = sqrt( sq(x - pos) + sq(y + pos) );
    omega[] =  coef * (-exp (-sq(r0) / sq(rc))+exp (-sq(r1) / sq(rc))-exp (-sq(r2) / sq(rc))+exp (-sq(r3) / sq(rc)));
  }
}

int main(int argc, char *argv[]) {
  //arguments: gmma, rc, lmda, time_sim
  if (argc > 1)
    gmma = atoi(argv[1]);
    
  if (argc > 2)
    rc = atoi(argv[2]);
    
  if (argc > 3)
    lmda = atoi(argv[3]);
    
  if (argc > 4)
    time_sim = atoi(argv[4]);
    
    
  //Characteristic time corresponding to a lagrangian revolution of the Lamb-Oseen vortex
  Tc = (4.*alph*sq(pi* rc)/((1.-exp(-alph))*24000));
  //Total domain corresponds to 4 vortex with one vortex domain being 20rc*20rc
  L0 = 20*2*rc;
  X0 = Y0 = -L0/2;
  N = 64;

  periodic (left);
  periodic (top);
  DT = Tc / time_div;
  run();
}

event vis (i = 0) {
  foreach_face ()
    mu.x[] = 1e-6 ;
}
event init (t = 0) {
  pos = L0/4.;
  int n_part = 64;
  double dl = L0/(2.*n_part);
  
  flow = new_tracer_particles (0);
  
  for (int j = 1; j <= n_part; j++) {
    for (int k = 1; k <= n_part; k++){
      particle p = {
                .x = (j-0.5) * dl ,
                .y = (k-0.5) * dl ,
              };
      add_particle (p, flow);
    }
  } 
  
  particle_boundary (flow);
  
  lpwriters[0] = init_particle_files("w",flow,"tracer","tracer"); 
  
  for (int l = 0; l < grp_n; l++)
  {
    inertial[l] = new_inertial_particles (0);
  
    for (int j = 1; j <= n_part; j++) {
      for (int k = 1; k <= n_part; k++){
        particle p = {
                  .x = (j-0.5) * dl ,
                  .y = (k-0.5) * dl ,
                  .u2.z = tparticles[l], //Set relaxation timescale 
                  .inertia = true
                };
        add_particle (p, inertial[l]);
      }
    }
    
    particle_boundary(inertial[l]);
    
    char particle_name[80];
    sprintf(particle_name,"inertial_t%d",(int) tparticles[l]);
    
    lpwriters[l+1] = init_particle_files("w",inertial[l],particle_name,particle_name);
  }
  

}
event kinematic (t += Tc/time_div) {
  printf("time = %g, progress = %g %\n", t, 100*t/(time_sim*Tc));
  particles_kinematics(lpwriters[0],t);
  for (int l = 0; l < grp_n; l++)
  {
    particles_kinematics(lpwriters[l+1],t);
  }
}

event set_dtmax (i++, last) {dt = dtnext(Tc / time_div);}

event set_stokes_drift (i = 1) {
  foreach () {
    boundary ({u});
    double g = 9.81;
    double k = 2. * pi / lmda;
    double om = sqrt(g*k);
    double u0 = gmma/(2*pi*rc);
    u.x[] += 0.*om* k * sq(a) + 0.0 * u0 * noise();
    u.y[] += 1.*om* k * sq(a) + 0.0 * u0 * noise();
    }
}


event text_field (i = 1) {

  char name[80];
  sprintf (name, "omega");
  FILE * fp = fopen (name, "w");
  
  output_field ({omega}, fp, pow(2, LEVEL), linear = true);
  fclose (fp);

  sprintf (name, "ux");
  fp = fopen (name, "w");
  
  output_field ({u.x}, fp, pow(2, LEVEL), linear = true);
  fclose (fp);

  sprintf (name, "uy");
  fp = fopen (name, "w");
  
  output_field ({u.y}, fp, pow(2, LEVEL), linear = true);
  fclose (fp);
  
}

event bviewer (t += Tc/time_div) {
   
 
   poisson(psi, omega);
 
   stats statpsi;
   statpsi = statsf (psi);
 
   if (true) {
     view (width = 1000, height = 1000);
     squares ("omega", map = blue_white_red, min = -1, max = 1);
     /*for (int l = 0; l < grp_n; l++)
       {
         scatter (inertial[l], s = l+1 , pc = color_particles[l]);
       }*/
     scatter (inertial[0], s = 3 , pc = {0.004, 0.796, 0.882});
     scatter (inertial[1], s = 3 , pc = {0.011, 0.564, 0.631});
     //scatter (inertial[2], s = 3 , pc = {0.011, 0.404, 0.45});
     scatter (flow , s = 1);
     box();
     save ("parts.mp4");
     save ("parts.png");
     if (i == 0)
       save ("parts0.png");
     
   clear();
   squares ("omega", map = jet);
   /*for (int l = 0; l < grp_n; l++)
       {
         scatter (inertial[l], s = l+1 , pc = color_particles[l]);
       }*/
   scatter (inertial[0], s = 3 , pc = {0.004, 0.796, 0.882});
   scatter (inertial[1], s = 3 , pc = {0.011, 0.564, 0.631});
   //scatter (inertial[2], s = 3 , pc = {0.011, 0.404, 0.45});
   scatter (flow , s = 1);
   box();
   save ("omega.mp4");
   save ("ux.png");
   if (i == 0)
     save ("omega0.png");
   
   for (double iv = statpsi.min; iv <= statpsi.max; iv += (statpsi.max - statpsi.min) / 10.) 
     iso_contour(psi, iv);
     
     // scatter (flow , s = 1);
   if (i == 0)
     save("stream0.png");
     //save("stream.png");
     
    // cells();
     // box();
     // save ("cells.mp4");
   }
 }
event stop (t = Tc * time_sim);
 
 
