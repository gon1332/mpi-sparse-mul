#include "test_suite.h"

#define TEST

int main(void)
{
#ifdef TEST
    test_suite();
#endif

    return 0;
}
