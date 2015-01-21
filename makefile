CC = clang
CFLAGS = -std=c99 -Wall -Wextra -ggdb3

main: coo_sparse.h test_suite.h
	$(CC) $(CFLAGS) coo_sparse.c test_suite.c coo_serial.c -o serial

docs: docs_conf
	doxygen docs_conf

clean:
	rm -rf *~ serial
