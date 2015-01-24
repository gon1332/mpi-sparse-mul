#include "test_suite.h"

/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

#define PARALLEL //SERIAL
#define NTASKS 4

int main(int argc, char *argv[])
{

#ifdef SERIAL
    /* Serial code. */
    test_suite();
#else
    /* Parallel code. */
    parallel_tests(NTASKS, argc, argv);
#endif

    return 0;
}
