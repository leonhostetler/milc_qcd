/************************ setup.c ****************************/
/* MIMD version 7 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//  $Log: setup.c,v $
//  Revision 1.24  2012/11/24 00:20:09  detar
//  Add support for future HYPISQ action.
//
//  Revision 1.23  2011/11/29 22:30:03  detar
//  Support new fermion links scheme.
//
//  Revision 1.22  2009/05/31 00:50:16  detar
//  Change naik_term_mass to naik_term_epsilon
//
//  Revision 1.21  2009/04/14 17:09:24  detar
//  AB: Calculate structure of multi_x array in terms of different Naik corrections
//  AB: Preset Naik corrections in eps_naik[] controlled by HISQ_NAIT_ADJUSTABLE
//
//  Revision 1.20  2008/03/28 15:20:22  detar
//  Bazavov changes 3/08.  Debugging and Naik term mass.
//
//  Revision 1.19  2007/12/14 04:51:45  detar
//  Add HISQ code.
//
//  Revision 1.18  2007/11/09 16:06:28  detar
//  Pull FN link calculation out of inverter.
//
//  Revision 1.17  2007/10/09 19:51:02  detar
//  Add ferm_links_t and ks_action_paths structures and pass them as params
//
//  Revision 1.16  2007/05/21 04:23:32  detar
//  Add Precision selection for fermion force in QOP and QDP
//
//  Revision 1.15  2007/03/27 20:48:54  detar
//  Fix n_pseudo check.
//
//  Revision 1.14  2006/12/13 18:41:40  detar
//  Add precision arg to mat_invert_uml and mat_invert_cg
//
//  Revision 1.13  2006/12/09 14:10:57  detar
//  Add mixed precision capability for QDP and QOP inverters
//
//  Revision 1.12  2006/11/16 06:07:50  detar
//  Moved optimization selections to Makefile
//  Put rational function parameters in preamble of input file
//  Recreated fiducial output files and remade errtol files
//  Remade rationals.sample files
//  Removed obsolete files
//
//  Revision 1.11  2006/11/07 02:30:58  detar
//  Fix some omissions to complete the previous update.
//
//  Revision 1.10  2006/11/04 23:35:14  detar
//  Add separate CG control for MD, FA, GR
//  Add nrestart parameter
//  Remove some debugging lines from sample output files and errtol files
//
//  Revision 1.9  2006/10/29 02:35:54  detar
//  Abandon rationals.h and load parameters from file instead.
//
//  Revision 1.8  2006/10/09 03:44:04  detar
//  Move fermion_force_fn to generic_ks/fermion_force_fn_multi.c
//  and path_transport_field to generic/path_transport.c
//  Change ks_multicg selection method to a set_opts call.
//
//  Revision 1.7  2006/10/02 04:13:50  detar
//  Distinguish inverter residuals for molecular dynamics and action
//
//  Revision 1.6  2006/09/19 03:06:58  detar
//  Upgrade for concurrent EOS calculations
//
//  Revision 1.5  2006/08/24 04:29:02  detar
//  Remove or protect unused variable declarations.
//  Remove unwanted Makefile
//
//  Revision 1.4  2006/08/22 19:32:22  detar
//  Upgrade to QDP and train error tolerance file
//
//  Revision 1.3  2006/08/13 15:16:06  detar
//  Train error files.  Trivia.
//
//  Revision 1.2  2006/08/13 04:02:32  detar
//  Switch to function pointer for selecting multicg inverter species
//
//  Revision 1.1  2006/08/09 04:22:19  detar
//  Adding Doug's code to repository
//
//  Revision 1.3  2005/11/10 16:58:44  detar
//  Experimenting with tags and versions
//
//  Revision 1.2  2005/11/10 16:55:09  detar
//  Add CVS version and log tags to file
//
//
*/
/* MIMD version 7 */
#define IF_OK if(status==0)
#define _POSIX_C_SOURCE 200112L // for gethostname

#include "ks_imp_includes.h"
#include "quark_action.h"
#include "lattice_qdp.h"
#include "../include/fermion_links.h"
#define SU3_MAT_OP_NO_STORAGE
#include "../include/su3_mat_op.h"
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

/* Each node has a params structure for passing simulation parameters */
#include "params.h"

/* Forward declarations */
static int initial_set(void);
static void third_neighbor(int, int, int, int, int *, int, int *, int *, int *, int *);
static void make_3n_gathers(void);

int
setup(void)
{
  int prompt;

  
  /* print banner, get volume, seed */
  prompt = initial_set();
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();
  /* initialize the node random number generator */
  initialize_prn( &node_prn, iseed, volume+mynode() );
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();

  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  /* set up 3rd nearest neighbor pointers and comlink structures
     code for this routine is below  */
  make_3n_gathers();
  /* set up K-S phase vectors, boundary conditions */
  phaseset();
  
  return( prompt );
}

static int n_naiks;
static double eps_naik[MAX_NAIK];

/* SETUP ROUTINES */
static int 
initial_set(void)
{
  int prompt=0,status,i,tmporder;
  Real current_naik_epsilon;

  /* On node zero, read lattice size, seed, and send to others */
  if(mynode()==0){
    /* print banner */
    printf("SU3 with improved KS action\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("Rational function hybrid Monte Carlo algorithm\n");
    printf("MIMD version %s\n",MILC_CODE_VERSION);
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");

    /* Print list of options selected */
    if(mynode()==0)printf("Options selected...\n");
    show_generic_opts();
    show_generic_ks_opts();
    show_generic_ks_md_opts();
#ifdef INT_ALG
    if(mynode()==0)printf("INT_ALG=%s\n",ks_int_alg_opt_chr());
#endif
#if FERM_ACTION == HISQ
    show_su3_mat_opts();
    show_hisq_links_opts();
    show_hisq_force_opts();
#elif FERM_ACTION == HYPISQ
    show_su3_mat_opts();
    show_hypisq_links_opts();
    show_hypisq_force_opts();
#endif

    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &param.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &param.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &param.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &param.nt );
#ifdef FIX_NODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "node_geometry", 
			   param.node_geometry, 4);
#ifdef FIX_IONODE_GEOM
    IF_OK status += get_vi(stdin, prompt, "ionode_geometry", 
			   param.ionode_geometry, 4);
#endif
#endif
    {
      int iseed;
      IF_OK {
	status += get_i(stdin, prompt,"iseed", &iseed );
	param.iseed = iseed;  /* implicitly convert int to unsigned */
      }
    }
    /* Number of pseudofermions */
    IF_OK status += get_i(stdin, prompt,"n_pseudo", &param.n_pseudo );
    if(param.n_pseudo > MAX_N_PSEUDO){
      printf("Error:  Too many pseudofermion fields.  Recompile. Current max is %d\n"
	     ,MAX_N_PSEUDO);
      terminate(1);
    }
    /* get name of file containing rational function parameters */
    IF_OK status += get_s(stdin, prompt, "load_rhmc_params", 
			  param.rparamfile);
    /* beta, quark masses */
    IF_OK status += get_f(stdin, prompt,"beta", &param.beta );

    IF_OK status += get_i(stdin, prompt,"n_dyn_masses", &param.n_dyn_masses );
    IF_OK status += get_vf(stdin, prompt, "dyn_mass", param.dyn_mass, param.n_dyn_masses);
    IF_OK status += get_vi(stdin, prompt, "dyn_flavors", param.dyn_flavors, param.n_dyn_masses);

    IF_OK status += get_f(stdin, prompt,"u0", &param.u0 );

    /* Use antiperiodic boundary conditions */
    for(int k = 0; k < 4; k++)
      param.coord_origin[k] = 0;
    param.time_bc = 0;

    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  fflush(stdout);
  
  if( param.stopflag != 0 )
    normal_exit(0);
  
  nx        = param.nx;
  ny        = param.ny;
  nz        = param.nz;
  nt        = param.nt;
#ifdef FIX_NODE_GEOM
  for(i = 0; i < 4; i++)
    node_geometry[i] = param.node_geometry[i];
#ifdef FIX_IONODE_GEOM
  for(i = 0; i < 4; i++)
    ionode_geometry[i] = param.ionode_geometry[i];
#endif
#endif
  iseed     = param.iseed;
  n_pseudo  = param.n_pseudo;
  strcpy(rparamfile,param.rparamfile);
  
  number_of_nodes = numnodes();
  volume=(size_t)nx*ny*nz*nt;
  total_iters=0;
#ifdef HISQ_SVD_COUNTER
  hisq_svd_counter = 0;
#endif
      
#ifdef HISQ_FORCE_FILTER_COUNTER
  hisq_force_filter_counter = 0;
#endif

  /* Load rational function parameters */
  rparam = load_rhmc_params(rparamfile, n_pseudo);  
  if(rparam == NULL)terminate(1);

  /* Determine the maximum rational fcn order */
  max_rat_order = 0;
  for(i = 0; i < n_pseudo; i++){
    if(rparam[i].MD.order > max_rat_order)max_rat_order = rparam[i].MD.order;
    if(rparam[i].GR.order > max_rat_order)max_rat_order = rparam[i].GR.order;
    if(rparam[i].FA.order > max_rat_order)max_rat_order = rparam[i].FA.order;
  }
  if(mynode()==0)printf("Maximum rational func order is %d\n",max_rat_order);
  fflush(stdout);

  /* Determine the number of different Naik masses
     and fill in n_orders_naik and n_pseudo_naik        */
  current_naik_epsilon = rparam[0].naik_term_epsilon;
  tmporder = 0;
  n_naiks = 0;
  n_order_naik_total = 0;
  for( i=0; i<n_pseudo; i++ ) {
    if( rparam[i].naik_term_epsilon != current_naik_epsilon ) {
      if( tmporder > 0 ) {
        n_orders_naik[n_naiks] = tmporder;
	eps_naik[n_naiks] = current_naik_epsilon;
        current_naik_epsilon = rparam[i].naik_term_epsilon;
        n_naiks++;
        n_order_naik_total += tmporder;
        tmporder = 0;
      }
    }
    tmporder += rparam[i].MD.order;
    n_pseudo_naik[n_naiks]++;
  }
  if( tmporder > 0 ) {
    n_orders_naik[n_naiks] = tmporder;
    eps_naik[n_naiks] = current_naik_epsilon;
    n_order_naik_total += tmporder;
    n_naiks++;
  }
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  // calculate epsilon corrections for different Naik terms
  if( 0 != eps_naik[0] ) {
    if(mynode()==0)printf("IN HISQ AND HYPISQ ACTIONS FIRST SET OF PSEUDO FERMION FIELDS SHOULD HAVE EPSILON CORRECTION TO NAIK TERM ZERO.\n");
    fflush(stdout);
    terminate(1);
  }
#endif
  if(mynode()==0)printf("Naik term correction structure of multi_x:\n");
  if(mynode()==0)printf("n_naiks %d\n",n_naiks);
  for( i=0; i<n_naiks; i++ ) {
    if(mynode()==0)printf("n_pseudo_naik[%d]=%d\n", i, n_pseudo_naik[i]);
    if(mynode()==0)printf("n_orders_naik[%d]=%d\n", i, n_orders_naik[i]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
    if(mynode()==0)printf("eps_naik[%d]=%f\n", i, eps_naik[i]);
#endif
  }
  fflush(stdout);

  if(mynode()==0)printf("n_order_naik_total %d\n",n_order_naik_total);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  if( n_naiks+1 > MAX_NAIK ) {
    if(mynode()==0)printf("MAX_NAIK=%d < n_naiks+1=%d\n", MAX_NAIK, n_naiks+1 );
    if(mynode()==0)printf("Increase MAX_NAIK\n");
    fflush(stdout);
    terminate(1);
  }
#else /* non HISQ */
  if( n_naiks>1 ) {
    if(mynode()==0)printf("FOR ACTIONS OTHER THAN HISQ AND HYPISQ EPSILON CORRECTION IS NOT USED.\n");
    if(mynode()==0)printf("ONLY ONE SET OF X LINKS IS USED.\n");
    if(mynode()==0)printf("SET ALL naik_mass TO 0 IN RATIONAL FUNCTION FILE.\n");
    fflush(stdout);
    terminate(1);
  }
#endif /* HISQ */

  beta = param.beta;
  
  n_dyn_masses = param.n_dyn_masses;
  for(i = 0; i < n_dyn_masses; i++){
    dyn_mass[i] = param.dyn_mass[i];
    dyn_flavors[i] = param.dyn_flavors[i];
  }
  u0 = param.u0;

  return(prompt);
}

/* read in parameters and coupling constants	*/
int
readin(int prompt)
{
  /* read in parameters for su3 monte carlo	*/
  /* argument "prompt" is 1 if prompts are to be given for input	*/
  
  int status;
  int i;
  
  /* On node zero, read parameters and send to all other nodes */
  if(mynode()==0) {
    
    printf("\n\n");
    status=0;
    
    /* warms, trajecs */
    IF_OK status += get_i(stdin, prompt,"warms", &param.warms );
    IF_OK status += get_i(stdin, prompt,"trajecs", &param.trajecs );
    
    /* trajectories between propagator measurements */
    IF_OK status += 
      get_i(stdin, prompt,"traj_between_meas", &param.propinterval );
    
    /* microcanonical time step */
    IF_OK status += 
      get_f(stdin, prompt,"microcanonical_time_step", &param.epsilon );
    
    /*microcanonical steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &param.steps );
    
    /* Data for each pseudofermion */

    for(i = 0; i < param.n_pseudo; i++){
      Real tmp[3]; int itmp[3];

      /* Residuals for multicg solves */
      IF_OK status += get_vf(stdin, prompt,"cgresid_md_fa_gr", tmp, 3 );
      /* rsqmin is r**2 in conjugate gradient */
      IF_OK {
	param.rsqmin_md[i] = tmp[0]*tmp[0];
	param.rsqmin_fa[i] = tmp[1]*tmp[1];
	param.rsqmin_gr[i] = tmp[2]*tmp[2];
      }

      /* Max CG iterations for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "max_multicg_md_fa_gr", itmp, 3);
      IF_OK {
	param.niter_md[i] = itmp[0];
	param.niter_fa[i] = itmp[1];
	param.niter_gr[i] = itmp[2];
      }

      /* Precision for multicg solves */
      IF_OK status += get_vi(stdin, prompt, "cgprec_md_fa_gr", itmp, 3);
      IF_OK {
	param.prec_md[i] = itmp[0];
	param.prec_fa[i] = itmp[1];
	param.prec_gr[i] = itmp[2];
      }
    }

    /* Max restarts for cleanup solves */
    IF_OK param.nrestart = 5;
    
    /* Precision for fermion force calculation */
    IF_OK status = get_i(stdin, prompt, "prec_ff", &param.prec_ff);

    /*------------------------------------------------------------*/
    /* Chiral condensate and related quantities                   */
    /*------------------------------------------------------------*/

    /* Eigenpairs for HMC not supported at present */
    param.eigen_param.Nvecs = 0;  
    
    IF_OK status += get_i(stdin, prompt, "number_of_pbp_masses",
			  &param.num_pbp_masses);
    if(param.num_pbp_masses > MAX_MASS_PBP){
      printf("Number of masses exceeds dimension %d\n",MAX_MASS_PBP);
      status++;
    }
    IF_OK if(param.num_pbp_masses > 0){
      IF_OK status += get_i(stdin, prompt, "max_cg_prop",
			    &param.qic_pbp[0].max);
      IF_OK status += get_i(stdin, prompt, "max_cg_prop_restarts",
			    &param.qic_pbp[0].nrestart);
      IF_OK status += get_i(stdin, prompt, "npbp_reps", &param.npbp_reps );
      IF_OK status += get_i(stdin, prompt, "prec_pbp", &param.prec_pbp);
      IF_OK for(i = 0; i < param.num_pbp_masses; i++){
	IF_OK status += get_f(stdin, prompt, "mass", &param.ksp_pbp[i].mass);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	IF_OK status += get_f(stdin, prompt, "naik_term_epsilon", 
			      &param.ksp_pbp[i].naik_term_epsilon);
#endif
	param.qic_pbp[i].min = 0;
	param.qic_pbp[i].start_flag = 0;
	param.qic_pbp[i].nsrc = 1;
	param.qic_pbp[i].max = param.qic_pbp[0].max;
	param.qic_pbp[i].nrestart = param.qic_pbp[0].nrestart;
	param.qic_pbp[i].prec = param.prec_pbp;
	IF_OK status += get_f(stdin, prompt, "error_for_propagator", &param.qic_pbp[i].resid);
	IF_OK status += get_f(stdin, prompt, "rel_error_for_propagator", &param.qic_pbp[i].relresid );
      }
    }

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(param.startflag),
					 param.startfile );
    
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
				       param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				 param.stringLFN );
    
    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */
  
    /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));
  
  if( param.stopflag != 0 )return param.stopflag;
  
  warms = param.warms;
  trajecs = param.trajecs;
  steps = param.steps;
  propinterval = param.propinterval;
  niter = param.niter;
  nrestart = param.nrestart;
  for(i = 0; i< n_pseudo; i++){
    niter_md[i] = param.niter_md[i];
    niter_fa[i] = param.niter_fa[i];
    niter_gr[i] = param.niter_gr[i];

    rsqmin_md[i] = param.rsqmin_md[i];
    rsqmin_fa[i] = param.rsqmin_fa[i];
    rsqmin_gr[i] = param.rsqmin_gr[i];

    prec_md[i] = param.prec_md[i];
    prec_fa[i] = param.prec_fa[i];
    prec_gr[i] = param.prec_gr[i];
  }
  prec_ff = param.prec_ff;
  rsqprop = param.rsqprop;
  epsilon = param.epsilon;
  n_pseudo = param.n_pseudo;
  startflag = param.startflag;
  saveflag = param.saveflag;
  strcpy(startfile,param.startfile);
  strcpy(savefile,param.savefile);
  strcpy(stringLFN, param.stringLFN);

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  int isite, idir;
  site *s;
  FORALLSITES(isite,s) {
    for( idir=XUP;idir<=TUP;idir++ ) {
      lattice[isite].on_step_Y[idir] = 0;
      lattice[isite].on_step_W[idir] = 0;
      lattice[isite].on_step_V[idir] = 0;
      /* zero out information from previous time step
         if fresh lattice; keep everything for continuation */
      if( startflag != CONTINUE ){
        lattice[isite].phase_Y_previous[idir] = 0.0;
        lattice[isite].phase_Y[idir] = 0.0;
        lattice[isite].Vdet[idir] = 0.0;
        clear_su3mat( &(lattice[isite].Wlink[idir]) );
        clear_su3mat( &(lattice[isite].Wlink_previous[idir]) );
      }
    }
  }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */

#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  /* Add PBP quantities to the eps_naik table of unique Naik epsilon
     coefficients .  Also build the hash table for mapping a mass term to
     its Naik epsilon index */

  /* Contribution from the chiral condensate epsilons */
  for(i = 0; i < param.num_pbp_masses; i++){
    param.ksp_pbp[i].naik_term_epsilon_index = 
      fill_eps_naik(eps_naik, &n_naiks, 
		    param.ksp_pbp[i].naik_term_epsilon);
  }

#endif

  /* Do whatever is needed to get lattice */
  if( startflag == CONTINUE ){
    rephase( OFF );
  }
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile );
  /* if a lattice was read in, put in KS phases and AP boundary condition */

  phases_in = OFF;
  rephase( ON );
  
  /* Copy gauge links from site structure to field */
  /* (This is a transitional step.  As we upgrade the code, we will 
     read/construct the lattice directly in the field and drop the
     site structure */

  /* Set options for fermion links */
#ifdef DM_DU0
  /* We want to calculate both the links and their u0 derivatives */
  fermion_links_want_du0(1);
#endif
  
#ifdef DBLSTORE_FN
  /* We want to double-store the links for optimization */
  fermion_links_want_back(1);
#endif
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )& defined(DM_DEPS)
  /* We want to calculate both the links and their Naik eps
     derivatives (HISQ only) */
  fermion_links_want_deps(1);
#endif
  
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fn_links = create_fermion_links_from_site(MILC_PRECISION, n_naiks, eps_naik);
#else
  fn_links = create_fermion_links_from_site(MILC_PRECISION, 0, NULL);
#endif
    
  /* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();
  
  return(0);
}

/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
static void 
make_3n_gathers(void)
{
  int i;
  
  for(i=XUP; i<=TUP; i++) {
    make_gather(third_neighbor, &i, WANT_INVERSE,
		ALLOW_EVEN_ODD, SWITCH_PARITY);
  }
  
  /* Sort into the order we want for nearest neighbor gathers,
     so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  
  sort_eight_gathers(X3UP);
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

static void 
third_neighbor(int x, int y, int z, int t, int *dirpt, int FB,
	       int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
	"forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+3)%nx; break;
  case XDOWN: *xp = (x+4*nx-3)%nx; break;
  case YUP: *yp = (y+3)%ny; break;
  case YDOWN: *yp = (y+4*ny-3)%ny; break;
  case ZUP: *zp = (z+3)%nz; break;
  case ZDOWN: *zp = (z+4*nz-3)%nz; break;
  case TUP: *tp = (t+3)%nt; break;
  case TDOWN: *tp = (t+4*nt-3)%nt; break;
  default: printf("third_neighb: bad direction\n"); exit(1);
  }
}
