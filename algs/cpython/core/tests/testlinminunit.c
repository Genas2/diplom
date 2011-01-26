

#include <stdafx.h>
#include <stdio.h>
#include "testlinminunit.h"


/*$ Declarations $*/


/*$ Body $*/


ae_bool testlinmin(ae_bool silent, ae_state *_state)
{
    ae_bool waserrors;
    ae_bool result;


    waserrors = ae_false;
    if( !silent )
    {
        printf("TESTING LINMIN\n");
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*$ End $*/
