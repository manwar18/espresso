#include <mpi.h>
#include <string.h>
#include "layered.h"
#include "global.h"
#include "communication.h"
#include "debug.h"
#include "ghosts.h"
#include "forces.h"
#include "pressure.h"
#include "energy.h"
#include "constraint.h"
#include "domain_decomposition.h"

/* Organization: Layers only in one direction.
   ghost_bottom
   c1
   c2
   c3
   .
   .
   .
   cn
   ghost_top

   First, all nodes send downwards, then upwards. Within these actions,
   first the odd nodes send. For even n_nodes, this algorithm is straight
   forward: first the odd nodes send, the even receive, then vice versa.
   For odd n_nodes, we have
   1) 1->2 3->4 5->1
   2) 2->3 4->5
   So in the first round node 5 has to wait for node 1 to
   complete the send and get ready to receive. In other words,
   what physically happens is:
   1) 1->2 3->4 5->*
   2) *->1 2->3 4->*
   3) *->5
   This means that one pending receive has to be done in addition
   provided that all other transmissions can happen in parallel.

*/

/** wether we are the lowest node */
#define LAYERED_BOTTOM 1
/** wether we are the highest node */
#define LAYERED_TOP    2
/** same as PERIODIC(2) */
#define LAYERED_PERIODIC 4
#define LAYERED_BTM_MASK (LAYERED_BOTTOM|LAYERED_PERIODIC)
#define LAYERED_TOP_MASK (LAYERED_TOP|LAYERED_PERIODIC)
/** node has a neighbor above (modulo n_nodes) */
#define LAYERED_TOP_NEIGHBOR ((layered_flags & LAYERED_TOP_MASK) != LAYERED_TOP)
/** node has a neighbor below (modulo n_nodes) */
#define LAYERED_BTM_NEIGHBOR ((layered_flags & LAYERED_BTM_MASK) != LAYERED_BOTTOM)

int layered_flags = 0;
int n_layers = -1, determine_n_layers = 1, btm, top;
double layer_h = 0, layer_h_i = 0;

static void layered_get_mi_vector(double res[3], double a[3], double b[3])
{
  int i;

  for(i=0;i<2;i++) {
    res[i] = a[i] - b[i];
#ifdef PARTIAL_PERIODIC
    if (PERIODIC(i))
#endif
      res[i] -= dround(res[i]*box_l_i[i])*box_l[i];
  }
  res[2] = a[2] - b[2];
}

Cell *layered_position_to_cell(double pos[3])
{
  int cpos = (int)((pos[2] - my_left[2])*layer_h_i) + 1;
  if (cpos < 1) {
    if (!LAYERED_BTM_NEIGHBOR)
      cpos = 1;
    else
      return NULL;
  }
  else if (cpos > n_layers) {
    /* not periodic, but at top */
    if (!LAYERED_TOP_NEIGHBOR)
      cpos = n_layers;
    else
      return NULL;
  }
  return &(cells[cpos]);
}

void layered_topology_release()
{
  CELL_TRACE(fprintf(stderr,"%d: layered_topology_release:\n", this_node));
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void layered_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int even_odd;
  int c, n;

  if (n_nodes > 1) {
    /* more than one node => no local transfers */

    /* how many communications to do: up even/odd, down even/odd */
    n = 4;
    /* one communication missing if not periodic but on border */
    if (!LAYERED_TOP_NEIGHBOR)
      n -= 2;
    if (!LAYERED_BTM_NEIGHBOR)
      n -= 2;

    prepare_comm(comm, data_parts, n);

    /* always sending/receiving 1 cell per time step */
    for(c = 0; c < n; c++) {
      comm->comm[c].part_lists = malloc(sizeof(ParticleList *));
      comm->comm[c].n_part_lists = 1;
      comm->comm[c].mpi_comm = MPI_COMM_WORLD;
    }

    c = 0;

    CELL_TRACE(fprintf(stderr, "%d: ghostrec new comm\n", this_node));
    /* downwards */
    for (even_odd = 0; even_odd < 2; even_odd++) {
      /* send */
      if (this_node % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
	comm->comm[c].type = GHOST_SEND;
	/* round 1 use prefetched data */
	if (c == 1)
	  comm->comm[c].type |= GHOST_PREFETCH;
	comm->comm[c].node = btm;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[0];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send force to %d btmg\n", this_node, btm));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[1];

	  /* if periodic and bottom or top, send shifted */
	  comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	  if (((layered_flags & LAYERED_BTM_MASK) == LAYERED_BTM_MASK) &&
	      (data_parts & GHOSTTRANS_POSITION)) {
	    comm->data_parts |= GHOSTTRANS_POSSHFTD;
	    comm->comm[c].shift[2] = box_l[2];
	  }
	  else
	    comm->comm[c].shift[2] = 0;
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send to %d shift %f btml\n", this_node, btm, comm->comm[c].shift[2]));
	}
	c++;
      }
      /* recv. Note we test r_node as we always have to test the sender
	 as for odd n_nodes maybe we send AND receive. */
      if (top % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
	comm->comm[c].type = GHOST_RECV;
	/* round 0 prefetch */
	if (c == 0)
	  comm->comm[c].type |= GHOST_PREFETCH;
	comm->comm[c].node = top;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[n_layers];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec get force from %d topl\n", this_node, top));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec recv from %d topg\n", this_node, btm));
	}
	c++;
      }
    }

    /* upwards */
    for (even_odd = 0; even_odd < 2; even_odd++) {
      /* send */
      if (this_node % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
	comm->comm[c].type = GHOST_SEND;
	/* round 1 use prefetched data from round 0.
	   But this time there may already have been two transfers downwards */
	if (c % 2 == 1)
	  comm->comm[c].type |= GHOST_PREFETCH;
	comm->comm[c].node = top;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send force to %d topg\n", this_node, top));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[n_layers];

	  /* if periodic and bottom or top, send shifted */
	  comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	  if (((layered_flags & LAYERED_TOP_MASK) == LAYERED_TOP_MASK) &&
	      (data_parts & GHOSTTRANS_POSITION)) {
	    comm->data_parts |= GHOSTTRANS_POSSHFTD;
	    comm->comm[c].shift[2] = -box_l[2];
	  }
	  else
	    comm->comm[c].shift[2] = 0;
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec send to %d shift %f topl\n", this_node, top, comm->comm[c].shift[2]));
	}
	c++;
      }
      /* recv. Note we test r_node as we always have to test the sender
	 as for odd n_nodes maybe we send AND receive. */
      if (btm % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
	comm->comm[c].type = GHOST_RECV;
	/* round 0 prefetch. But this time there may already have been two transfers downwards */
	if (c % 2 == 0)
	  comm->comm[c].type |= GHOST_PREFETCH;
	comm->comm[c].node = btm;
	if (data_parts == GHOSTTRANS_FORCE) {
	  comm->comm[c].part_lists[0] = &cells[1];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec get force from %d btml\n", this_node, btm));
	}
	else {
	  comm->comm[c].part_lists[0] = &cells[0];
	  CELL_TRACE(fprintf(stderr, "%d: ghostrec recv from %d btmg\n", this_node, btm));
	}
	c++;
      }
    }
  }
  else {
    /* one node => local transfers, either 2 (up and down, periodic) or zero*/

    n = (layered_flags & LAYERED_PERIODIC) ? 2 : 0;

    prepare_comm(comm, data_parts, n);

    if (n != 0) {
      /* two cells: from and to */
      for(c = 0; c < n; c++) {
	comm->comm[c].part_lists = malloc(2*sizeof(ParticleList *));
	comm->comm[c].n_part_lists = 2;
	comm->comm[c].mpi_comm = MPI_COMM_WORLD;
	comm->comm[c].node = this_node;
      }

      c = 0;

      /* downwards */
      comm->comm[c].type = GHOST_LOCL;
      if (data_parts == GHOSTTRANS_FORCE) {
	comm->comm[c].part_lists[0] = &cells[0];
	comm->comm[c].part_lists[1] = &cells[n_layers];
      }
      else {
	comm->comm[c].part_lists[0] = &cells[1];
	comm->comm[c].part_lists[1] = &cells[n_layers + 1];
	/* here it is periodic */
	if (data_parts & GHOSTTRANS_POSITION)
	  comm->data_parts |= GHOSTTRANS_POSSHFTD;
	comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	comm->comm[c].shift[2] = box_l[2];
      }
      c++;

      /* upwards */
      comm->comm[c].type = GHOST_LOCL;
      if (data_parts == GHOSTTRANS_FORCE) {
	comm->comm[c].part_lists[0] = &cells[n_layers + 1];
	comm->comm[c].part_lists[1] = &cells[1];
      }
      else {
	comm->comm[c].part_lists[0] = &cells[n_layers];
	comm->comm[c].part_lists[1] = &cells[0];
	/* here it is periodic */
	if (data_parts & GHOSTTRANS_POSITION)
	  comm->data_parts |= GHOSTTRANS_POSSHFTD;
	comm->comm[c].shift[0] = comm->comm[c].shift[1] = 0;
	comm->comm[c].shift[2] = -box_l[2];
      }
    }
  }
}

void layered_topology_init(CellPList *old)
{
  Particle *part;
  int c, p, np;

  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init, %d old particle lists\n", this_node, old->n));

  if (cell_structure.type != CELL_STRUCTURE_LAYERED) {
    cell_structure.type = CELL_STRUCTURE_LAYERED;
    cell_structure.position_to_node = map_position_node_array;
    cell_structure.position_to_cell = layered_position_to_cell;
  }

  /* check node grid. All we can do is 1x1xn. */
  if (node_grid[0] != 1 || node_grid[1] != 1) {
    fprintf(stderr, "ERROR: selected node grid is not suitable for layered cell structure (needs 1x1xn grid)\n");
    errexit();
  }

  if (determine_n_layers) {
    if (max_range > 0) {
      n_layers = (int)floor(local_box_l[2]/max_range);
      if (n_layers < 1) {
	fprintf(stderr, "ERROR: layered: maximal interaction range %f larger than local box length %f\n", max_range, local_box_l[2]);
	errexit();
      }
      if (n_layers > max_num_cells - 2)
	n_layers = max_num_cells - 2;
    }
    else
      n_layers = 1;
  }

  top = (this_node == n_nodes - 1) ? 0 : this_node + 1;
  btm = (this_node == 0) ? n_nodes - 1 : this_node - 1;

  layer_h = local_box_l[2]/(double)(n_layers);
  layer_h_i = 1/layer_h;

  if (layer_h < max_range) {
    fprintf(stderr, "ERROR: layered: maximal interaction range %f larger than layer height %f\n", max_range, layer_h);
    errexit();
  }

  /* check wether node is top and/or bottom */
  layered_flags = 0;
  if (this_node == 0)
    layered_flags |= LAYERED_BOTTOM;
  if (this_node == n_nodes - 1)
    layered_flags |= LAYERED_TOP;

  if (PERIODIC(2))
    layered_flags |= LAYERED_PERIODIC;

  /* allocate cells and mark them */
  realloc_cells(n_layers + 2);
  realloc_cellplist(&local_cells, local_cells.n = n_layers);
  for (c = 0; c < n_layers; c++)
    local_cells.cell[c] = &cells[c + 1];
  realloc_cellplist(&ghost_cells, ghost_cells.n = 2);
  ghost_cells.cell[0] = &cells[0];
  ghost_cells.cell[1] = &cells[n_layers + 1];

  /* create communicators */
  layered_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
  layered_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++) {
      Cell *nc = layered_position_to_cell(part[p].r.p);
      /* particle does not belong to this node. Just stow away
	 somewhere for the moment */
      if (nc == NULL)
	nc = local_cells.cell[0];
      append_unindexed_particle(nc, &part[p]);
    }
  }
  for (c = 1; c <= n_layers; c++)
    update_local_particles(&cells[c]);

  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init done\n", this_node));
}
 
static void layered_append_particles(ParticleList *pl)
{
  int p;

  for(p = 0; p < pl->n; p++) {
    if (LAYERED_BTM_NEIGHBOR && pl->part[p].r.p[2] < my_left[2])
      continue;
    if (LAYERED_TOP_NEIGHBOR && pl->part[p].r.p[2] > my_right[2])
      continue;

    move_indexed_particle(layered_position_to_cell(pl->part[p].r.p), pl, p);
  }
}

void layered_exchange_and_sort_particles(int global_flag)
{
  Particle *part;
  Cell *nc, *oc;
  int c, p, flag, redo;
  ParticleList send_buf_dn, send_buf_up;
  ParticleList recv_buf;

  init_particlelist(&send_buf_dn);
  init_particlelist(&send_buf_up);

  init_particlelist(&recv_buf);

  for (;;) {
    /* sort local particles */
    for (c = 1; c <= n_layers; c++) {
      oc = &cells[c];

      for (p = 0; p < oc->n; p++) {
	part = &oc->part[p];
	fold_position(part->r.p, part->l.i);
	nc = layered_position_to_cell(part->r.p);

	if (LAYERED_BTM_NEIGHBOR && part->r.p[2] < my_left[2])
	  move_indexed_particle(&send_buf_dn, oc, p);
	else if (LAYERED_TOP_NEIGHBOR && part->r.p[2] > my_right[2])
	  move_indexed_particle(&send_buf_up, oc, p);
	else {
	  if (nc != oc) {
	    move_indexed_particle(nc, oc, p);
	    if (p < oc->n) p--;
	  }
	}
      }
    }

    /* transfer */
    if (n_nodes > 1) {
      if (this_node % 2 == 0) {
	/* send down */
	if (LAYERED_BTM_NEIGHBOR)
	  send_particles(&send_buf_dn, btm);
	if (LAYERED_TOP_NEIGHBOR)
	  recv_particles(&recv_buf, top);
	/* send up */
	if (LAYERED_BTM_NEIGHBOR)
	  send_particles(&send_buf_up, top);
	if (LAYERED_TOP_NEIGHBOR)
	  recv_particles(&recv_buf, btm);
      }
      else {
	if (LAYERED_TOP_NEIGHBOR)
	  recv_particles(&recv_buf, top);
	if (LAYERED_BTM_NEIGHBOR)
	  send_particles(&send_buf_dn, btm);
	if (LAYERED_BTM_NEIGHBOR)
	  recv_particles(&recv_buf, btm);
	if (LAYERED_TOP_NEIGHBOR)
	send_particles(&send_buf_up, top);
      }
    }
#ifdef ADDITIONAL_CHECKS
    else {
      if (recv_buf.n != 0 || send_buf_dn.n != 0 || send_buf_up.n != 0) {
	fprintf(stderr, "1 node but transfer buffers not empty\n");
	errexit();
      }
    }
#endif
    layered_append_particles(&recv_buf);

    /* handshake redo */
    flag = (recv_buf.n != 0);

    if (global_flag == LAYERED_FULL_EXCHANGE) {
      MPI_Allreduce(&flag, &redo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (!redo)
	break;
      CELL_TRACE(fprintf(stderr, "%d: another exchange round\n", this_node));
    }
    else {
      if (flag) {
	fprintf(stderr,"%d: layered_exchange_and_sort_particles: particle moved to far\n", this_node);
	errexit();
      }
      break;
    }
  }

  realloc_particlelist(&recv_buf, 0);
}

/** nonbonded and bonded force calculation using the verlet list */
void layered_calculate_ia()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    if (LAYERED_BTM_NEIGHBOR) {
      cellb = &cells[c-1];
      pb    = cellb->part;
      npb   = cellb->n;
    }
    else {
      npb   = 0;
      pb    = NULL;
    }

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      add_bonded_force(p1);
#ifdef CONSTRAINTS
      add_constraints_forces(p1);
#endif

      /* cell itself and bonded / constraints */
      CELL_TRACE(fprintf(stderr, "%d: inter local\n", this_node));
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_force(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      CELL_TRACE(fprintf(stderr, "%d: inter bottom\n", this_node));
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_force(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
}

void layered_calculate_energies()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    if (LAYERED_BTM_NEIGHBOR) {
      cellb = &cells[c-1];
      pb    = cellb->part;
      npb   = cellb->n;
    }
    else {
      npb   = 0;
      pb    = NULL;
    }

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      add_kinetic_energy(p1);

      add_bonded_energy(p1);
#ifdef CONSTRAINTS
      add_constraints_energy(p1);
#endif

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_energy(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_energy(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
}

void layered_calculate_virials()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];
 
  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    if (LAYERED_BTM_NEIGHBOR) {
      cellb = &cells[c-1];
      pb    = cellb->part;
      npb   = cellb->n;
    }
    else {
      npb   = 0;
      pb    = NULL;
    }

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      add_kinetic_virials(p1);

      add_bonded_virials(p1);

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_virials(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
	dist2 = sqrlen(d);
	add_non_bonded_pair_virials(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
}
