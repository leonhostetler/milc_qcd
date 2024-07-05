/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This file is for lattice generation with the RHMC algorithm */

#define CONTROL
#include "ks_spectrum_includes.h"
#include <string.h>


EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int
main( int argc, char **argv )
{
  int prompt; 
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  g_sync();

  /* set up */
  prompt = setup();
  
  /* loop over input sets */
  while( readin(prompt) == 0) {
    
    // Create and fill fermion vector for testing
    srand(1);
    su3_vector *src = create_v_field();
    int i, j;
    site *s;
    FORALLSITESDOMAIN(i,s) {
      for(j=0; j<3; j++) {
        src[i].c[j].real = (float)rand()/(float)(RAND_MAX);
        src[i].c[j].imag = (float)rand()/(float)(RAND_MAX);
      }
    }

    // Print random fermion vector
    if(mynode()==0){
      printf("\nRandom fermion vector:\n");
      FORALLSITESDOMAIN(i,s) {
        for(j=0; j<3; j++) {
          printf("%f%+fi\n", src[i].c[j].real, src[i].c[j].imag);
        }
      }
      fflush(stdout);
    }

    // Create solution vector
    su3_vector *dest = create_v_field();
    // Get fermion links
    imp_ferm_links_t **fn;
    fn = get_fm_links(fn_links);
    // Compute Dslash*src
    dslash_fn_field(src, dest, EVENANDODD, fn[0]);
    FORALLSITESDOMAIN(i,s) { // Factor of 1/2
      for(j=0; j<3; j++) {
          dest[i].c[j].real = dest[i].c[j].real/2.0;
          dest[i].c[j].imag = dest[i].c[j].imag/2.0;
        }
    }
    FORALLSITESDOMAIN(i,s) {
      scalar_mult_add_su3_vector(dest+i, src+i, -0.92, dest+i); // diagonal/mass term
    }
    // Print solution vector
    if(mynode()==0){
      printf("\nSolution vector:\n");
      FORALLSITESDOMAIN(i,s) {
        for(j=0; j<3; j++) {
          printf("%f%+fi\n", dest[i].c[j].real, dest[i].c[j].imag);
        }
      }
      fflush(stdout);
    }
    // Clean up
    node0_printf("DONE\n"); fflush(stdout);
    destroy_v_field(src);
    destroy_v_field(dest);

    /* Destroy fermion links (created in readin() */
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
  }
  free_lattice();
  normal_exit(0);
  return 0;
}

