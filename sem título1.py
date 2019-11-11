

/******************************************************************************

                            Online C Compiler.
                Code, Compile, Run and Debug C program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <stdio.h>
#include <string.h>
char* faltante(){
    char *a[55] = {"1","2","3","5","111","2221","2223","4445","6422","7445", "8855"};
    static int i = -1;
    i++;
    i=i%55;
    return a[i];
}
int main()
{
    char *a=0;
    
    char msg[42];

    int i = 0, cnt = 0, faltantes = 1;

    while((a = faltante()) != NULL && cnt < 42){
        
        if(strlen(a) + cnt == 42 || (strlen(a) + cnt < 42 && 2*strlen(a) + cnt  > 42)){
            
            cnt = cnt + sprintf(msg+cnt, "%s", a);
            break;
            
        }else if(strlen(a) + cnt < 42){

            cnt = cnt + sprintf(msg+cnt, "%s-", a);
            //printf("%s %d\n", msg, cnt);

        }
        
    }
    
    
    if(a == NULL){
            
        printf("Chegou aque\n");
        msg[cnt-1] = '\0';

    }
        
    
    printf("%s %d\n", msg, strlen(msg));

    
    
    return 0;
}
