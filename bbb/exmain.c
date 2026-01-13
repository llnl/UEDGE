/**
 * @file exmain.c
 *
 * Purpose: Allow the trapping of control-C to support Basis-like debug
 *          mode in the python version.
 *
 * $Id: exmain.c,v 1.6 2021/03/23 18:43:14 meyer8 Exp $
 *
 */

#ifdef FORTHON
#include <Python.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <setjmp.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAS_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif



static struct sigaction act,oact;
static sigjmp_buf ev;


/* 
    Handler for SIGINT signal
*/
void int_handler(int signum) {
   (void)signum;
   char mymyline[200],*ret;
   sigset_t block_mask;
   printf("\nType \"cont\" to continue exmain(), \"abort\" (not compatible with openmp) or \"stop\" (with openmp) to return to Python prompt \n");
   printf("or a single line to be evaluated by Python.\n");
#pragma omp master
    {
int condition;
condition=1;
   while(1){
#ifdef HAS_READLINE
       ret = readline("Debug>>> ");
       if(ret == (char *)NULL)return;
       add_history(ret); 
       strncpy(mymyline,ret,sizeof(mymyline)-1); 
       free(ret); 
#else
       printf("Debug>>> ");
       ret = fgets(mymyline,150,stdin);
       if(ret == (char *)NULL)return;
#endif
       if(strncmp(mymyline,"cont",4) == 0){
           return;
       } else if (strncmp(mymyline,"abort",5) == 0) {
          PyRun_SimpleString("bbb.exmain_aborted = True");
          siglongjmp(ev,1);
       } else if (strncmp(mymyline,"stop",4) == 0) {
	  PyRun_SimpleString("print(\"Stopping exmain ... Please wait...\")");
          PyRun_SimpleString("bbb.exmain_aborted = True");
          return;
       } else if (strncmp(mymyline,"exit",4) == 0) {
          PyRun_SimpleString("bbb.exmain_aborted = True");
          siglongjmp(ev,1);
       } else {
          PyRun_SimpleString(mymyline);
          /* matplotlib seems to unset the hander
             so it is set again just in case */
          sigfillset(&block_mask);
          act.sa_handler = int_handler;
          act.sa_mask = block_mask;
          act.sa_flags = 0;
          sigaction(SIGINT,&act,NULL);
       }
   }
   
}
}
#endif



/* FORTHON is defined by the Python build. This exmain does nothing when 
   compiled for the basis version of the code, it just drops through to
   the Fortran routine.  */

#if defined(FC_FUNC)
void FC_FUNC(exmain_f, EXMAIN_F)(void); 
#else
void exmain_f_(void); 
#endif

#if defined(FC_FUNC)
void FC_FUNC(exmain, EXMAIN)() {
#else
void exmain_() {
#endif
#ifdef FORTHON
   sigset_t block_mask;
   int ival;
#pragma omp master
{
   ival = sigsetjmp(ev,1);
   if(ival != 0){
       sigaction(SIGINT,&oact,NULL);
       return;
   }
   }
   

/* setup to catch SIGINT and save the previous handler to be restored
   on return */
#pragma omp master
{
   sigfillset(&block_mask);
   act.sa_handler = int_handler;
   act.sa_mask = block_mask;
   act.sa_flags = 0;
   sigaction(SIGINT,&act,&oact);

   PyRun_SimpleString("from uedge import bbb");
   PyRun_SimpleString("bbb.exmain_aborted = False");
   }
   

#endif


/*  now call the Fortran version of exmain  */

#if defined(FC_FUNC)
   FC_FUNC(exmain_f, EXMAIN_F)();
#else
   exmain_f_(); 
#endif
#ifdef FORTHON
   sigaction(SIGINT,&oact,NULL);
#endif
}


