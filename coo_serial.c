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

#define TEST

int main(void)
{
#ifdef TEST
    test_suite();
#endif

    return 0;
}
