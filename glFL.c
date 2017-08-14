/* Cavity Sim
 * MSci Physics Project QOLS07, Imperial College ; Supervised by Prof Geoff New
 * Jarvist Frost & Benjamin Hall 2005-2006
 * Rewritten and made slightly less horrible, Jarvist Moore Frost 2015
 */
/* Has been conjoined with: (gl) Video Feedback by Jarvist Moore Frost
 *  This code written 2004-2010.
 *  GLUT skeleton found at http://stackoverflow.com/questions/503816/linux-fastest-way-to-draw
*/
/* Install GLUT development libraries, on Ubuntu/Deb something like this should work:
   #sudo apt-get install libglut3-dev
*/

#include "cavity_sim.c" 
#include <GL/glut.h>

#include <stdio.h>
#include <math.h>
#include <limits.h>

#define X_RES N/2 
#define Y_RES N/2

//#define MAKEMOVIE // if you want lots of .pbm outputfiles :^)

//original on: http://stackoverflow.com/questions/503816/linux-fastest-way-to-draw
void renderScene() {    
    char picfile[50];

    simulate_laser_partA();
    curpic_ap_picture();

#ifdef MAKEMOVIE
        sprintf(picfile,"pic-%06D.pbm",framecount); // framecount global variable updated elsewhere
        output_ap_picture(picfile);
#endif 

    passes++;

    // render the texture here
    glEnable (GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glTexImage2D (
        GL_TEXTURE_2D,
        0,
        GL_RGB,
        X_RES,
        Y_RES,
        0,
        GL_RGB,
        GL_UNSIGNED_BYTE,
        &curpic[0][0][0]
    );

    glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0, -1.0);
        glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0, -1.0);
        glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0,  1.0);
        glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0,  1.0);
    glEnd();

    glFlush();
    glutSwapBuffers();

    simulate_laser_partB();
}

int main(int argc, char **argv) 
{
    setup_fftw(); // Build plans for this size FFT
    //setup_cavity_conjugate(); // Construct conjugate planes + etc. Not
    //working currently?
    setup_cavity_direct(); // Construct FOCUS and propogation
	 
	make_filter (6); // Make aperture polygon with this size
	fprintf (stderr, "Npolygon: %d M: %f Focal: %f\n", n, M,FOCAL);

//   input_ap_picture(); //Lena
	generate_initial_intensity ();
	
    // GLUT setup
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

    glutInitWindowPosition(100, 100);
    glutInitWindowSize(X_RES,Y_RES);
    glutCreateWindow("VIDEO FRACTALS");

    glutDisplayFunc(renderScene);

    //glutTimerFunc(100, renderScene, 0);
	// max FPS
    glutIdleFunc(renderScene);

    glutMainLoop();

    destroy_fftw();
    return 0;
}
