#ifndef TEST_SUITE_H_JGP7CH2X
#define TEST_SUITE_H_JGP7CH2X

extern void test_suite(void);
extern void parallel_tests(int numtasks, int argc, char *argv[]);
int readMtx( char filename[], int **Ii, int **Ji, double **vali, int *size);

#endif /* TEST_SUITE_H_JGP7CH2X */

