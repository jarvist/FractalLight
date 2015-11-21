/* Cavity Sim
 * MSci Physics Project QOLS07, Imperial College ; Supervised by Prof Geoff New
 * Jarvist Frost & Benjamin Hall 2005-2006
 * Rewritten and made slightly less horrible, Jarvist Moore Frost 2015
 */

#include "cavity_sim.c" 

main ()
{
    setup_fftw(); // Build plans for this size FFT
    //setup_cavity_conjugate(); // Construct conjugate planes + etc. Not
    //working currently?
    setup_cavity_direct(); // Construct FOCUS and propogation
	 
	make_filter (5); // Make aperture polygon with this size
	fprintf (stderr, "Npolygon: %d M: %f Focal: %f\n", n, M,FOCAL);

//   input_ap_picture(); //Lena
	generate_initial_intensity ();
	 
	 for (passes = 0; passes < 10000; passes++)
	  {
        simulate_laser_partA();
        sprintf(name,"%.5d.pnm",framecount);
               output_ap_picture(name);	     
        simulate_laser_partB();
	  }
	fprintf (stderr, "Reset\n");

    destroy_fftw();
}
