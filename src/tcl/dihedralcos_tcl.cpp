/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file dihedral_tcl.c
 *
 *  Parser for the dihedral potential
 */

#include "dihedralcos_tcl.hpp"
#include "dihedralcos.hpp"



/// parse parameters for the dihedral potential
int tclcommand_inter_parse_dihedralcos(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, bend1, bend2;

  if (argc < 4 ) {
    Tcl_AppendResult(interp, "dihedralcos needs 3 parameters: "
		     "<bend> <bend1> <bend2>", (char *) NULL);
    return (ES_ERROR);
  }
  if ( !ARG_IS_D(1, bend) || !ARG_IS_D(2, bend1) || !ARG_IS_D(3, bend2) ) {
    Tcl_AppendResult(interp, "dihedralcos needs 3 parameters of types DOUBLE DOUBLE DOUBLE: "
		     "<bend> <bend1> <bend2> ", (char *) NULL);
    return ES_ERROR;
  }
  
  CHECK_VALUE(dihedralcos_set_params(bond_type, bend, bend1, bend2), "bond type must be nonnegative");
}



int tclprint_to_result_dihedralcosIA(Tcl_Interp *interp,
				  Bonded_ia_parameters *params)
{
  char buffer[TCL_DOUBLE_SPACE];
  sprintf(buffer, "%e", (double)(params->p.dihedralcos.bend));
  Tcl_AppendResult(interp, "dihedralcos ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.dihedralcos.bend1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, params->p.dihedralcos.bend2, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  return (TCL_OK);

}
