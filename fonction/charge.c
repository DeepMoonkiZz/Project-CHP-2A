#include <stdio.h>
#include <stdlib.h>

void charge(int me, int n, int nb_proc, int *iBeg, int *iEnd)
{
    int ss_int;
    int reste;
    
    ss_int = n/nb_proc;
    reste = n%nb_proc;

    if (reste == 0) {
        *iBeg = me * ss_int;
        *iEnd = (me + 1) * ss_int - 1;
    }

    else if (me < reste) {
        *iBeg = me * (ss_int + 1);
        *iEnd = (me + 1) * (ss_int + 1) - 1;   
    }

    else if (me >= reste) {
        *iBeg = reste + me * ss_int;
        *iEnd = *iBeg + ss_int - 1;   
    }

    else {
        printf("error in charge");
    }
}