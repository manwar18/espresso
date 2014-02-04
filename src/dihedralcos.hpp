// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef DIHEDRALCOS_H
#define DIHEDRALCOS_H
/** \file dihedral.h Routines to calculate the dihedral energy or/and
 *  and force for a particle quadruple.  Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!  \ref forces.c
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

#define ANGLE_NOT_DEFINED -100

/// set dihedral parameters
int dihedralcos_set_params(int bond_type, double bend, double bend1, double bend2);

/** Calculates the dihedral angle between particle quadruple p1, p2,
p3 and p4. The dihedral angle is the angle between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). 
Vectors a, b and c are the bond vectors between consequtive particles.
If the a,b or b,c are parallel the dihedral angle is not defined in which
case the routine returns phi=-1. Calling functions should check for that
(Written by: Arijit Maitra) */


inline void calc_dihedralcos_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
				  double a[3], double b[3], double c[3], 
				  double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
				  double *cosphi, double *phi)
{
    int i;
    get_mi_vector(a, p2->r.p, p1->r.p);
    get_mi_vector(b, p3->r.p, p2->r.p);
    get_mi_vector(c, p4->r.p, p3->r.p);

    /* calculate vector product a X b and b X c */
    vector_product(a, b, aXb);
    vector_product(b, c, bXc);

    /* calculate the unit vectors */
    *l_aXb = sqrt(sqrlen(aXb));
    *l_bXc = sqrt(sqrlen(bXc));

    /* catch case of undefined dihedral angle */
    if ( *l_aXb <= TINY_LENGTH_VALUE || *l_bXc <= TINY_LENGTH_VALUE ) { *phi = -1.0; *cosphi = 0; return;}
    for (i=0;i<3;i++) {
        aXb[i] /= *l_aXb;
        bXc[i] /= *l_bXc;
    }

    *cosphi = scalar(aXb, bXc);

    if ( fabs(fabs(*cosphi)-1)  < TINY_SIN_VALUE  ) *cosphi = dround(*cosphi);
        /* Calculate dihedral angle */
        *phi = acos(*cosphi);
        if( scalar(aXb, c) < 0.0 ) {*phi = (2.0*PI) - *phi;}
 }


inline int calc_dihedralcos_force(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force2[3],
				 double force1[2], double force3[2])
{
    int i;
    /* vectors for dihedral angle calculation */
    double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
    double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];

    /* dihedral angle, cosine of the dihedral angle */
    double phi, cosphi ;
  
    /* force factors sinmphi_sinphi */
    double fac, f1[3], f4[3];

    /* dihedral angle */
    calc_dihedralcos_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
    /* dihedral angle not defined - force zero */
    if ( phi == -1.0 ) { 
        for(i=0;i<3;i++) { force1[i] = 0.0; force2[i] = 0.0; force3[i] = 0.0; }
        return 0;
    }

    /* calculate force components (directions) */
    for(i=0;i<3;i++)  {
        f1[i] = (v23Xv34[i] - cosphi*v12Xv23[i])/l_v12Xv23;
        f4[i] = (v12Xv23[i] - cosphi*v23Xv34[i])/l_v23Xv34;
    }

    vector_product(v23, f1, v23Xf1);
    vector_product(v23, f4, v23Xf4);
    vector_product(v34, f4, v34Xf4);
    vector_product(v12, f1, v12Xf1);

    if(fabs(sin(phi)) < TINY_SIN_VALUE) {
        /* First term of Maclaurin Series*///Just first derivative of force from sinphi(I multiplied odd terms of cos with -1 due to sign convension(-1)^n for cos(n*phi))
        fac  = -iaparams->p.dihedralcos.bend*cos(phi);
        fac +=  4*iaparams->p.dihedralcos.bend1*cos(2*phi);
        fac += -9*iaparams->p.dihedralcos.bend2*cos(3*phi);
        fac *= -0.5;
        fac /=  cos(phi);
    }

    else { //(I multiplied odd terms of cos with -1 due to sign convension, (-1)^n for cos(n*phi))
        fac  = -iaparams->p.dihedralcos.bend*sin(phi);
        fac +=  2*iaparams->p.dihedralcos.bend1*sin(2*phi);
        fac += -3*iaparams->p.dihedralcos.bend2*sin(3*phi);
        fac *= -0.5;
        fac /=  sin(phi);
    }
    /* store dihedral forces */
    for(i=0;i<3;i++) {
        force1[i] = fac*v23Xf1[i];
        force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]); 
        force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
    }
    return 0;
}


inline int dihedralcos_energy(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
			     Bonded_ia_parameters *iaparams, double *_energy) 
{
    /* vectors for dihedral calculations. */
    double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
    /* dihedral angle, cosine of the dihedral angle */
    double phi, cosphi;
    /* energy factors */
    double fac;
    calc_dihedralcos_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
    //(I multiplied odd terms of cos with -1 due to sign convension, (-1)^n for cos(n*phi))
    fac  = iaparams->p.dihedralcos.bend * (1 + cos(phi));
    fac += iaparams->p.dihedralcos.bend1 * (1 - cos(2*phi));
    fac += iaparams->p.dihedralcos.bend2  * (1 + cos(3*phi));
    fac *= 0.5;
    *_energy = fac;
    return 0;
}

#endif
