/* Include files */

#include "fleet_model_sfun.h"
#include "fleet_model_sfun_debug_macros.h"
#include "c1_fleet_model.h"
#include "c2_fleet_model.h"
#include "c3_fleet_model.h"
#include "c4_fleet_model.h"
#include "c5_fleet_model.h"
#include "c6_fleet_model.h"
#include "c7_fleet_model.h"
#include "c8_fleet_model.h"
#include "c9_fleet_model.h"
#include "c10_fleet_model.h"
#include "c15_fleet_model.h"
#include "c16_fleet_model.h"
#include "c17_fleet_model.h"
#include "c18_fleet_model.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _fleet_modelMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void fleet_model_initializer(void)
{
}

void fleet_model_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_fleet_model_method_dispatcher(SimStruct *simstructPtr, unsigned
  int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==16) {
    c16_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==18) {
    c18_fleet_model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_fleet_model_process_check_sum_call( int nlhs, mxArray * plhs[],
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
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(718980213U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(784913590U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(151588553U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2169697620U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1612149201U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1926047018U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3224725308U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1935695481U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c1_fleet_model_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c2_fleet_model_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c3_fleet_model_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c4_fleet_model_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c5_fleet_model_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c6_fleet_model_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c7_fleet_model_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c8_fleet_model_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c9_fleet_model_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c10_fleet_model_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c15_fleet_model_get_check_sum(plhs);
          break;
        }

       case 16:
        {
          extern void sf_c16_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c16_fleet_model_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c17_fleet_model_get_check_sum(plhs);
          break;
        }

       case 18:
        {
          extern void sf_c18_fleet_model_get_check_sum(mxArray *plhs[]);
          sf_c18_fleet_model_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3031367619U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4001028638U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3978939492U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(838979348U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2615792668U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3983241017U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3371106623U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(496645936U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_fleet_model_autoinheritance_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
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
        if (strcmp(aiChksum, "hZTv82GXTzMy3biWPxXCXB") == 0) {
          extern mxArray *sf_c1_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c1_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "8rXTa16BDUf0o5sJR5OP0E") == 0) {
          extern mxArray *sf_c2_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c2_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "hZTv82GXTzMy3biWPxXCXB") == 0) {
          extern mxArray *sf_c3_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c3_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "P1zv5GYmaV1kzj4dTv8F5") == 0) {
          extern mxArray *sf_c4_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c4_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "hZTv82GXTzMy3biWPxXCXB") == 0) {
          extern mxArray *sf_c5_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c5_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "UDLF52EBWLu1Wb0fimMzmG") == 0) {
          extern mxArray *sf_c6_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c6_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "rHcddR4JbH1lTrBnW2LjqE") == 0) {
          extern mxArray *sf_c7_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c7_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "h87v3JDLyLFW7kUTWLsYHF") == 0) {
          extern mxArray *sf_c8_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c8_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "edJUxuZYWgfG5xVpVe3C2F") == 0) {
          extern mxArray *sf_c9_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c9_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "EDY3OpvD6VM6zUkc9jKpWG") == 0) {
          extern mxArray *sf_c10_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c10_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "YJeCjyeZTsBeGhl2ylgi5C") == 0) {
          extern mxArray *sf_c15_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c15_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 16:
      {
        if (strcmp(aiChksum, "YJeCjyeZTsBeGhl2ylgi5C") == 0) {
          extern mxArray *sf_c16_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c16_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "YJeCjyeZTsBeGhl2ylgi5C") == 0) {
          extern mxArray *sf_c17_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c17_fleet_model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 18:
      {
        if (strcmp(aiChksum, "YJeCjyeZTsBeGhl2ylgi5C") == 0) {
          extern mxArray *sf_c18_fleet_model_get_autoinheritance_info(void);
          plhs[0] = sf_c18_fleet_model_get_autoinheritance_info();
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

unsigned int sf_fleet_model_get_eml_resolved_functions_info( int nlhs, mxArray *
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
        extern const mxArray *sf_c1_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray *sf_c10_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 16:
      {
        extern const mxArray *sf_c16_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c16_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray *sf_c17_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_fleet_model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 18:
      {
        extern const mxArray *sf_c18_fleet_model_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c18_fleet_model_get_eml_resolved_functions_info();
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

unsigned int sf_fleet_model_third_party_uses_info( int nlhs, mxArray * plhs[],
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
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c1_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c1_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "F3TfhItzVW0lv1EOAcoQgH") == 0) {
          extern mxArray *sf_c2_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c2_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c3_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c3_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "mpf6zSEj4cY344ZZFyKf0") == 0) {
          extern mxArray *sf_c4_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c4_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c5_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c5_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "qgGzdzCaUx3abOIwPKyTvD") == 0) {
          extern mxArray *sf_c6_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c6_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "s5oTCkGhFJJJWqSmpnCnaG") == 0) {
          extern mxArray *sf_c7_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c7_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "EdQ92DG9wZth2dYgK5SWB") == 0) {
          extern mxArray *sf_c8_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c8_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "uGDt08L789reduRVk42TvD") == 0) {
          extern mxArray *sf_c9_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c9_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "klD9sVlWYN5tTPFDO8XMIE") == 0) {
          extern mxArray *sf_c10_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c10_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c15_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c15_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c16_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c16_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c17_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c17_fleet_model_third_party_uses_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c18_fleet_model_third_party_uses_info(void);
          plhs[0] = sf_c18_fleet_model_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_fleet_model_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
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
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c1_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "F3TfhItzVW0lv1EOAcoQgH") == 0) {
          extern mxArray *sf_c2_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c3_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "mpf6zSEj4cY344ZZFyKf0") == 0) {
          extern mxArray *sf_c4_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "yqZjO6dl88XxrAoqvZWeyG") == 0) {
          extern mxArray *sf_c5_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "qgGzdzCaUx3abOIwPKyTvD") == 0) {
          extern mxArray *sf_c6_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "s5oTCkGhFJJJWqSmpnCnaG") == 0) {
          extern mxArray *sf_c7_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "EdQ92DG9wZth2dYgK5SWB") == 0) {
          extern mxArray *sf_c8_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "uGDt08L789reduRVk42TvD") == 0) {
          extern mxArray *sf_c9_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "klD9sVlWYN5tTPFDO8XMIE") == 0) {
          extern mxArray *sf_c10_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c15_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c16_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c16_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c17_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c17_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "zqiLc8QnCfsh2OPxEZZFF") == 0) {
          extern mxArray *sf_c18_fleet_model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c18_fleet_model_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void fleet_model_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _fleet_modelMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "fleet_model","sfun",0,14,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_fleet_modelMachineNumber_,
    0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_fleet_modelMachineNumber_,
    0);
}

void fleet_model_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_fleet_model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("fleet_model",
      "fleet_model");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_fleet_model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
