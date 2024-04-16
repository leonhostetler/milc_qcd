/************************* control.c *******************************/
/* MIMD version 7 */

#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */


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
    su3_vector *src = create_v_field();
    int i, j;
    site *s;
    FORALLSITESDOMAIN(i,s) {
      for(j=0; j<3; j++) {
        src[i].c[j].real = 1.0;
        src[i].c[j].imag = 0.0;
      }
    }

    // Create solution vector
    su3_vector *dest = create_v_field();

    // Compute Dslash*src
    dslash_field(src, dest, EVENANDODD);
    FORALLSITESDOMAIN(i,s) { // Factor of 0.5
      for(j=0; j<3; j++) {
          dest[i].c[j].real = 0.5*dest[i].c[j].real;
          dest[i].c[j].imag = 0.5*dest[i].c[j].imag;
        }
    }
    FORALLSITESDOMAIN(i,s) {
      scalar_mult_add_su3_vector(dest+i, src+i, -0.92, dest+i); // diagonal/mass term
    }

    // Print solution vector
    // Note that dest contains the odd sites first and then the even sites
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
  }
  free_lattice();  
  normal_exit(0);
  return 0;
}

