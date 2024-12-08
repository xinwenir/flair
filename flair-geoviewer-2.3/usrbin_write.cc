#include <stdio.h>
#include <stdio.h>
#include <assert.h>

#include "os.h"
#include "usrbin.h"

using namespace std;

/** main */
int main(int, char *av[])
{
	Usrbin usrbin(av[1], 1);
	usrbin.save(av[2]);
} // main
