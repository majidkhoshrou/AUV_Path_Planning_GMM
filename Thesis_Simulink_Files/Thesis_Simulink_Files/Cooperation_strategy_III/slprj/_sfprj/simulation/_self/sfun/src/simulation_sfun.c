/* Include files */

#include "simulation_sfun.h"
#include "simulation_sfun_debug_macros.h"
#include "c1_simulation.h"
#include "c2_simulation.h"
#include "c3_simulation.h"
#include "c4_simulation.h"
#include "c5_simulation.h"
#include "c6_simulation.h"
#include "c7_simulation.h"
#include "c8_simulation.h"
#include "c9_simulation.h"
#include "c10_simulation.h"
#include "c11_simulation.h"
#include "c12_simulation.h"
#include "c13_simulation.h"
#include "c14_simulation.h"
#include "c15_simulation.h"
#include "c16_simulation.h"
#include "c17_simulation.h"
#include "c18_simulation.h"
#include "c19_simulation.h"
#include "c20_simulation.h"
#include "c21_simulation.h"
#include "c22_simulation.h"
#include "c23_simulation.h"
#include "c24_simulation.h"
#include "c25_simulation.h"
#include "c26_simulation.h"
#include "c27_simulation.h"
#include "c28_simulation.h"
#include "c29_simulation.h"
#include "c30_simulation.h"
#include "c31_simulation.h"
#include "c32_simulation.h"
#include "c34_simulation.h"
#include "c35_simulation.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _simulationMachineNumber_;
real_T _sfTime_;

/* Function Declarations */

/* Function Definitions */
void simulation_initializer(void)
{
}

void simulation_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_simulation_method_dispatcher(SimStruct *simstructPtr, unsigned
  int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==11) {
    c11_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==12) {
    c12_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==13) {
    c13_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==16) {
    c16_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==18) {
    c18_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==19) {
    c19_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==20) {
    c20_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==21) {
    c21_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==22) {
    c22_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==23) {
    c23_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==24) {
    c24_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==25) {
    c25_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==26) {
    c26_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==27) {
    c27_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==28) {
    c28_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==29) {
    c29_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==30) {
    c30_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==31) {
    c31_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==32) {
    c32_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==34) {
    c34_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==35) {
    c35_simulation_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_simulation_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2926814656U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(321913351U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1291019589U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(146373733U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4206893689U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3644054092U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4023821809U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1738556633U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_simulation_get_check_sum(mxArray *plhs[]);
          sf_c1_simulation_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_simulation_get_check_sum(mxArray *plhs[]);
          sf_c2_simulation_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_simulation_get_check_sum(mxArray *plhs[]);
          sf_c3_simulation_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_simulation_get_check_sum(mxArray *plhs[]);
          sf_c4_simulation_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_simulation_get_check_sum(mxArray *plhs[]);
          sf_c5_simulation_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_simulation_get_check_sum(mxArray *plhs[]);
          sf_c6_simulation_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_simulation_get_check_sum(mxArray *plhs[]);
          sf_c7_simulation_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_simulation_get_check_sum(mxArray *plhs[]);
          sf_c8_simulation_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_simulation_get_check_sum(mxArray *plhs[]);
          sf_c9_simulation_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_simulation_get_check_sum(mxArray *plhs[]);
          sf_c10_simulation_get_check_sum(plhs);
          break;
        }

       case 11:
        {
          extern void sf_c11_simulation_get_check_sum(mxArray *plhs[]);
          sf_c11_simulation_get_check_sum(plhs);
          break;
        }

       case 12:
        {
          extern void sf_c12_simulation_get_check_sum(mxArray *plhs[]);
          sf_c12_simulation_get_check_sum(plhs);
          break;
        }

       case 13:
        {
          extern void sf_c13_simulation_get_check_sum(mxArray *plhs[]);
          sf_c13_simulation_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_simulation_get_check_sum(mxArray *plhs[]);
          sf_c14_simulation_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_simulation_get_check_sum(mxArray *plhs[]);
          sf_c15_simulation_get_check_sum(plhs);
          break;
        }

       case 16:
        {
          extern void sf_c16_simulation_get_check_sum(mxArray *plhs[]);
          sf_c16_simulation_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_simulation_get_check_sum(mxArray *plhs[]);
          sf_c17_simulation_get_check_sum(plhs);
          break;
        }

       case 18:
        {
          extern void sf_c18_simulation_get_check_sum(mxArray *plhs[]);
          sf_c18_simulation_get_check_sum(plhs);
          break;
        }

       case 19:
        {
          extern void sf_c19_simulation_get_check_sum(mxArray *plhs[]);
          sf_c19_simulation_get_check_sum(plhs);
          break;
        }

       case 20:
        {
          extern void sf_c20_simulation_get_check_sum(mxArray *plhs[]);
          sf_c20_simulation_get_check_sum(plhs);
          break;
        }

       case 21:
        {
          extern void sf_c21_simulation_get_check_sum(mxArray *plhs[]);
          sf_c21_simulation_get_check_sum(plhs);
          break;
        }

       case 22:
        {
          extern void sf_c22_simulation_get_check_sum(mxArray *plhs[]);
          sf_c22_simulation_get_check_sum(plhs);
          break;
        }

       case 23:
        {
          extern void sf_c23_simulation_get_check_sum(mxArray *plhs[]);
          sf_c23_simulation_get_check_sum(plhs);
          break;
        }

       case 24:
        {
          extern void sf_c24_simulation_get_check_sum(mxArray *plhs[]);
          sf_c24_simulation_get_check_sum(plhs);
          break;
        }

       case 25:
        {
          extern void sf_c25_simulation_get_check_sum(mxArray *plhs[]);
          sf_c25_simulation_get_check_sum(plhs);
          break;
        }

       case 26:
        {
          extern void sf_c26_simulation_get_check_sum(mxArray *plhs[]);
          sf_c26_simulation_get_check_sum(plhs);
          break;
        }

       case 27:
        {
          extern void sf_c27_simulation_get_check_sum(mxArray *plhs[]);
          sf_c27_simulation_get_check_sum(plhs);
          break;
        }

       case 28:
        {
          extern void sf_c28_simulation_get_check_sum(mxArray *plhs[]);
          sf_c28_simulation_get_check_sum(plhs);
          break;
        }

       case 29:
        {
          extern void sf_c29_simulation_get_check_sum(mxArray *plhs[]);
          sf_c29_simulation_get_check_sum(plhs);
          break;
        }

       case 30:
        {
          extern void sf_c30_simulation_get_check_sum(mxArray *plhs[]);
          sf_c30_simulation_get_check_sum(plhs);
          break;
        }

       case 31:
        {
          extern void sf_c31_simulation_get_check_sum(mxArray *plhs[]);
          sf_c31_simulation_get_check_sum(plhs);
          break;
        }

       case 32:
        {
          extern void sf_c32_simulation_get_check_sum(mxArray *plhs[]);
          sf_c32_simulation_get_check_sum(plhs);
          break;
        }

       case 34:
        {
          extern void sf_c34_simulation_get_check_sum(mxArray *plhs[]);
          sf_c34_simulation_get_check_sum(plhs);
          break;
        }

       case 35:
        {
          extern void sf_c35_simulation_get_check_sum(mxArray *plhs[]);
          sf_c35_simulation_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(814460797U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(400623215U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1072597456U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1176453921U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3476620980U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1482397629U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2907754226U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4108439520U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_simulation_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "BFFgMwtUirVpzJyOPqQVOF") == 0) {
          extern mxArray *sf_c1_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c1_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "L1U9KWlF16TwgOJsxCOdKC") == 0) {
          extern mxArray *sf_c2_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c2_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "OeR1lMmOHDhnKxpINKMQ5G") == 0) {
          extern mxArray *sf_c3_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c3_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "MtA2Skv1UMVFHJdWgTiSX") == 0) {
          extern mxArray *sf_c4_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c4_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "BFFgMwtUirVpzJyOPqQVOF") == 0) {
          extern mxArray *sf_c5_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c5_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "L1U9KWlF16TwgOJsxCOdKC") == 0) {
          extern mxArray *sf_c6_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c6_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "OeR1lMmOHDhnKxpINKMQ5G") == 0) {
          extern mxArray *sf_c7_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c7_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "MtA2Skv1UMVFHJdWgTiSX") == 0) {
          extern mxArray *sf_c8_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c8_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "Lu03muUfFEvEfeDapiDzm") == 0) {
          extern mxArray *sf_c9_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c9_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "hFv827iKDSdAJrds9GSJlF") == 0) {
          extern mxArray *sf_c10_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c10_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 11:
      {
        if (strcmp(aiChksum, "1i0keOOofqyCVWTwJwtZXH") == 0) {
          extern mxArray *sf_c11_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c11_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 12:
      {
        if (strcmp(aiChksum, "lqOmPKQg6fA36l4NJ2CljF") == 0) {
          extern mxArray *sf_c12_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c12_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 13:
      {
        if (strcmp(aiChksum, "y8e0cPHhL8wowOVjmeuxWB") == 0) {
          extern mxArray *sf_c13_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c13_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "uG2trjyup9dMsOszMva5p") == 0) {
          extern mxArray *sf_c14_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c14_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "pBMiY0EJcj4eioiEh38iGC") == 0) {
          extern mxArray *sf_c15_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c15_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 16:
      {
        if (strcmp(aiChksum, "BFFgMwtUirVpzJyOPqQVOF") == 0) {
          extern mxArray *sf_c16_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c16_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "Lu03muUfFEvEfeDapiDzm") == 0) {
          extern mxArray *sf_c17_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c17_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 18:
      {
        if (strcmp(aiChksum, "hFv827iKDSdAJrds9GSJlF") == 0) {
          extern mxArray *sf_c18_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c18_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 19:
      {
        if (strcmp(aiChksum, "L1U9KWlF16TwgOJsxCOdKC") == 0) {
          extern mxArray *sf_c19_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c19_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 20:
      {
        if (strcmp(aiChksum, "OeR1lMmOHDhnKxpINKMQ5G") == 0) {
          extern mxArray *sf_c20_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c20_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 21:
      {
        if (strcmp(aiChksum, "1i0keOOofqyCVWTwJwtZXH") == 0) {
          extern mxArray *sf_c21_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c21_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 22:
      {
        if (strcmp(aiChksum, "MtA2Skv1UMVFHJdWgTiSX") == 0) {
          extern mxArray *sf_c22_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c22_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 23:
      {
        if (strcmp(aiChksum, "Lu03muUfFEvEfeDapiDzm") == 0) {
          extern mxArray *sf_c23_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c23_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 24:
      {
        if (strcmp(aiChksum, "lqOmPKQg6fA36l4NJ2CljF") == 0) {
          extern mxArray *sf_c24_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c24_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 25:
      {
        if (strcmp(aiChksum, "hFv827iKDSdAJrds9GSJlF") == 0) {
          extern mxArray *sf_c25_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c25_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 26:
      {
        if (strcmp(aiChksum, "1i0keOOofqyCVWTwJwtZXH") == 0) {
          extern mxArray *sf_c26_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c26_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 27:
      {
        if (strcmp(aiChksum, "lqOmPKQg6fA36l4NJ2CljF") == 0) {
          extern mxArray *sf_c27_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c27_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 28:
      {
        if (strcmp(aiChksum, "y8e0cPHhL8wowOVjmeuxWB") == 0) {
          extern mxArray *sf_c28_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c28_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 29:
      {
        if (strcmp(aiChksum, "uG2trjyup9dMsOszMva5p") == 0) {
          extern mxArray *sf_c29_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c29_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 30:
      {
        if (strcmp(aiChksum, "apEnm6j64xX2lsw7rrqo6D") == 0) {
          extern mxArray *sf_c30_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c30_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 31:
      {
        if (strcmp(aiChksum, "y8e0cPHhL8wowOVjmeuxWB") == 0) {
          extern mxArray *sf_c31_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c31_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 32:
      {
        if (strcmp(aiChksum, "pBMiY0EJcj4eioiEh38iGC") == 0) {
          extern mxArray *sf_c32_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c32_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 34:
      {
        if (strcmp(aiChksum, "uG2trjyup9dMsOszMva5p") == 0) {
          extern mxArray *sf_c34_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c34_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 35:
      {
        if (strcmp(aiChksum, "pBMiY0EJcj4eioiEh38iGC") == 0) {
          extern mxArray *sf_c35_simulation_get_autoinheritance_info(void);
          plhs[0] = sf_c35_simulation_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_simulation_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray *sf_c1_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray *sf_c10_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 11:
      {
        extern const mxArray *sf_c11_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c11_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 12:
      {
        extern const mxArray *sf_c12_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c12_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 13:
      {
        extern const mxArray *sf_c13_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray *sf_c14_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 16:
      {
        extern const mxArray *sf_c16_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c16_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray *sf_c17_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 18:
      {
        extern const mxArray *sf_c18_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c18_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 19:
      {
        extern const mxArray *sf_c19_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c19_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 20:
      {
        extern const mxArray *sf_c20_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c20_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 21:
      {
        extern const mxArray *sf_c21_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c21_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 22:
      {
        extern const mxArray *sf_c22_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c22_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 23:
      {
        extern const mxArray *sf_c23_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c23_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 24:
      {
        extern const mxArray *sf_c24_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c24_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 25:
      {
        extern const mxArray *sf_c25_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c25_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 26:
      {
        extern const mxArray *sf_c26_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c26_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 27:
      {
        extern const mxArray *sf_c27_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c27_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 28:
      {
        extern const mxArray *sf_c28_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c28_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 29:
      {
        extern const mxArray *sf_c29_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c29_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 30:
      {
        extern const mxArray *sf_c30_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c30_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 31:
      {
        extern const mxArray *sf_c31_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c31_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 32:
      {
        extern const mxArray *sf_c32_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c32_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 34:
      {
        extern const mxArray *sf_c34_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c34_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 35:
      {
        extern const mxArray *sf_c35_simulation_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c35_simulation_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_simulation_third_party_uses_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c1_simulation_third_party_uses_info(void);
          plhs[0] = sf_c1_simulation_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c2_simulation_third_party_uses_info(void);
          plhs[0] = sf_c2_simulation_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c3_simulation_third_party_uses_info(void);
          plhs[0] = sf_c3_simulation_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c4_simulation_third_party_uses_info(void);
          plhs[0] = sf_c4_simulation_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c5_simulation_third_party_uses_info(void);
          plhs[0] = sf_c5_simulation_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c6_simulation_third_party_uses_info(void);
          plhs[0] = sf_c6_simulation_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c7_simulation_third_party_uses_info(void);
          plhs[0] = sf_c7_simulation_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c8_simulation_third_party_uses_info(void);
          plhs[0] = sf_c8_simulation_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c9_simulation_third_party_uses_info(void);
          plhs[0] = sf_c9_simulation_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c10_simulation_third_party_uses_info(void);
          plhs[0] = sf_c10_simulation_third_party_uses_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c11_simulation_third_party_uses_info(void);
          plhs[0] = sf_c11_simulation_third_party_uses_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c12_simulation_third_party_uses_info(void);
          plhs[0] = sf_c12_simulation_third_party_uses_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c13_simulation_third_party_uses_info(void);
          plhs[0] = sf_c13_simulation_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c14_simulation_third_party_uses_info(void);
          plhs[0] = sf_c14_simulation_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c15_simulation_third_party_uses_info(void);
          plhs[0] = sf_c15_simulation_third_party_uses_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c16_simulation_third_party_uses_info(void);
          plhs[0] = sf_c16_simulation_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c17_simulation_third_party_uses_info(void);
          plhs[0] = sf_c17_simulation_third_party_uses_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c18_simulation_third_party_uses_info(void);
          plhs[0] = sf_c18_simulation_third_party_uses_info();
          break;
        }
      }

     case 19:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c19_simulation_third_party_uses_info(void);
          plhs[0] = sf_c19_simulation_third_party_uses_info();
          break;
        }
      }

     case 20:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c20_simulation_third_party_uses_info(void);
          plhs[0] = sf_c20_simulation_third_party_uses_info();
          break;
        }
      }

     case 21:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c21_simulation_third_party_uses_info(void);
          plhs[0] = sf_c21_simulation_third_party_uses_info();
          break;
        }
      }

     case 22:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c22_simulation_third_party_uses_info(void);
          plhs[0] = sf_c22_simulation_third_party_uses_info();
          break;
        }
      }

     case 23:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c23_simulation_third_party_uses_info(void);
          plhs[0] = sf_c23_simulation_third_party_uses_info();
          break;
        }
      }

     case 24:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c24_simulation_third_party_uses_info(void);
          plhs[0] = sf_c24_simulation_third_party_uses_info();
          break;
        }
      }

     case 25:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c25_simulation_third_party_uses_info(void);
          plhs[0] = sf_c25_simulation_third_party_uses_info();
          break;
        }
      }

     case 26:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c26_simulation_third_party_uses_info(void);
          plhs[0] = sf_c26_simulation_third_party_uses_info();
          break;
        }
      }

     case 27:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c27_simulation_third_party_uses_info(void);
          plhs[0] = sf_c27_simulation_third_party_uses_info();
          break;
        }
      }

     case 28:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c28_simulation_third_party_uses_info(void);
          plhs[0] = sf_c28_simulation_third_party_uses_info();
          break;
        }
      }

     case 29:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c29_simulation_third_party_uses_info(void);
          plhs[0] = sf_c29_simulation_third_party_uses_info();
          break;
        }
      }

     case 30:
      {
        if (strcmp(tpChksum, "Zdcsl9NFBPhOh43b6o46NC") == 0) {
          extern mxArray *sf_c30_simulation_third_party_uses_info(void);
          plhs[0] = sf_c30_simulation_third_party_uses_info();
          break;
        }
      }

     case 31:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c31_simulation_third_party_uses_info(void);
          plhs[0] = sf_c31_simulation_third_party_uses_info();
          break;
        }
      }

     case 32:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c32_simulation_third_party_uses_info(void);
          plhs[0] = sf_c32_simulation_third_party_uses_info();
          break;
        }
      }

     case 34:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c34_simulation_third_party_uses_info(void);
          plhs[0] = sf_c34_simulation_third_party_uses_info();
          break;
        }
      }

     case 35:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c35_simulation_third_party_uses_info(void);
          plhs[0] = sf_c35_simulation_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_simulation_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c1_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c2_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c3_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c4_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c5_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c6_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c7_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c8_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c9_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c10_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c11_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c11_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c12_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c12_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c13_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c13_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c14_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c14_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c15_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "efPVx35z9ooeRynC39kn0E") == 0) {
          extern mxArray *sf_c16_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c16_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c17_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c17_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c18_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c18_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 19:
      {
        if (strcmp(tpChksum, "fBk9zgNY5ZyqKmoxPyrvRF") == 0) {
          extern mxArray *sf_c19_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c19_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 20:
      {
        if (strcmp(tpChksum, "KbHOOpjeqG6w7b0b79p3g") == 0) {
          extern mxArray *sf_c20_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c20_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 21:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c21_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c21_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 22:
      {
        if (strcmp(tpChksum, "3tpUe8B5p7tcdlXcRFHEED") == 0) {
          extern mxArray *sf_c22_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c22_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 23:
      {
        if (strcmp(tpChksum, "MqMLod0FmHsYrBUgfZ425F") == 0) {
          extern mxArray *sf_c23_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c23_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 24:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c24_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c24_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 25:
      {
        if (strcmp(tpChksum, "CIyXD3hc1iCmjdmHrS8xRD") == 0) {
          extern mxArray *sf_c25_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c25_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 26:
      {
        if (strcmp(tpChksum, "MVWvJkFnvz9zAFQJSl6x4E") == 0) {
          extern mxArray *sf_c26_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c26_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 27:
      {
        if (strcmp(tpChksum, "K4ZOrZMb09CExf3z87c9rF") == 0) {
          extern mxArray *sf_c27_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c27_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 28:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c28_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c28_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 29:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c29_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c29_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 30:
      {
        if (strcmp(tpChksum, "Zdcsl9NFBPhOh43b6o46NC") == 0) {
          extern mxArray *sf_c30_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c30_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 31:
      {
        if (strcmp(tpChksum, "QTzmt1qx9WtHjQwVoCpiiG") == 0) {
          extern mxArray *sf_c31_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c31_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 32:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c32_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c32_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 34:
      {
        if (strcmp(tpChksum, "VRvN14rWB8XykPyyzWPveD") == 0) {
          extern mxArray *sf_c34_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c34_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     case 35:
      {
        if (strcmp(tpChksum, "enD2oeLM1OTy9hAeoXHiqB") == 0) {
          extern mxArray *sf_c35_simulation_updateBuildInfo_args_info(void);
          plhs[0] = sf_c35_simulation_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void simulation_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _simulationMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "simulation","sfun",0,34,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_simulationMachineNumber_,
    0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_simulationMachineNumber_,0);
}

void simulation_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_simulation_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("simulation",
      "simulation");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_simulation_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
