/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c35_simulation.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c35_debug_family_names[31] = { "l", "covll", "Pl", "xh",
  "yh", "xlh", "ylh", "dxhat", "dyhat", "rhat", "H", "R", "S", "K", "lambda",
  "cov11", "nargin", "nargout", "state_predict", "Sigma_predict", "z", "covl",
  "statel", "nsonar", "j", "i", "kl", "Sigma", "state_est", "state_i", "cov_i" };

/* Function Declarations */
static void initialize_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance);
static void enable_c35_simulation(SFc35_simulationInstanceStruct *chartInstance);
static void disable_c35_simulation(SFc35_simulationInstanceStruct *chartInstance);
static void c35_update_debugger_state_c35_simulation
  (SFc35_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c35_simulation
  (SFc35_simulationInstanceStruct *chartInstance);
static void set_sim_state_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_st);
static void finalize_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance);
static void sf_c35_simulation(SFc35_simulationInstanceStruct *chartInstance);
static void c35_chartstep_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc35_simulation(SFc35_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c35_machineNumber, uint32_T
  c35_chartNumber);
static const mxArray *c35_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static void c35_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_cov_i, const char_T *c35_identifier, real_T c35_y[4]);
static void c35_b_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[4]);
static void c35_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_b_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static void c35_c_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_state_i, const char_T *c35_identifier, real_T c35_y[2]);
static void c35_d_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[2]);
static void c35_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_c_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static void c35_e_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_state_est, const char_T *c35_identifier, real_T c35_y[9]);
static void c35_f_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[9]);
static void c35_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_d_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static void c35_g_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_Sigma, const char_T *c35_identifier, real_T c35_y[81]);
static void c35_h_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[81]);
static void c35_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_e_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static real_T c35_i_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId);
static void c35_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_f_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static void c35_j_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[4]);
static void c35_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static const mxArray *c35_g_sf_marshallOut(void *chartInstanceVoid, real_T
  c35_inData_data[17], int32_T c35_inData_sizes[2]);
static void c35_k_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T
  c35_y_data[17], int32_T c35_y_sizes[2]);
static void c35_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, real_T c35_outData_data[17],
  int32_T c35_outData_sizes[2]);
static void c35_info_helper(const mxArray **c35_info);
static const mxArray *c35_emlrt_marshallOut(char * c35_u);
static const mxArray *c35_b_emlrt_marshallOut(uint32_T c35_u);
static void c35_b_info_helper(const mxArray **c35_info);
static void c35_c_info_helper(const mxArray **c35_info);
static real_T c35_mpower(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_a);
static void c35_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_eml_error(SFc35_simulationInstanceStruct *chartInstance);
static void c35_b_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static boolean_T c35_eml_use_refblas(SFc35_simulationInstanceStruct
  *chartInstance);
static void c35_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance, int32_T
  c35_k, real_T c35_A_data[17], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[9]);
static void c35_check_forloop_overflow_error(SFc35_simulationInstanceStruct
  *chartInstance, boolean_T c35_overflow);
static void c35_c_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_mldivide(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_A[81], real_T c35_B[81], real_T c35_Y[81]);
static void c35_realmin(SFc35_simulationInstanceStruct *chartInstance);
static void c35_eps(SFc35_simulationInstanceStruct *chartInstance);
static void c35_eml_matlab_zgetrf(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_b_A[81], int32_T c35_ipiv[9], int32_T *c35_info);
static int32_T c35_eml_ixamax(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_n, real_T c35_x[81], int32_T c35_ix0);
static void c35_eml_xger(SFc35_simulationInstanceStruct *chartInstance, int32_T
  c35_m, int32_T c35_n, real_T c35_alpha1, int32_T c35_ix0, int32_T c35_iy0,
  real_T c35_A[81], int32_T c35_ia0, real_T c35_b_A[81]);
static void c35_eml_warning(SFc35_simulationInstanceStruct *chartInstance);
static void c35_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_A[81], real_T c35_B[81], real_T c35_b_B[81]);
static void c35_below_threshold(SFc35_simulationInstanceStruct *chartInstance);
static void c35_b_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81], real_T c35_b_B[81]);
static void c35_d_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_b_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B_data[17], int32_T c35_B_sizes[1], real_T c35_C
  [9], real_T c35_b_C[9]);
static void c35_b_below_threshold(SFc35_simulationInstanceStruct *chartInstance);
static void c35_e_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_tanh(SFc35_simulationInstanceStruct *chartInstance, real_T
                     c35_x[9], real_T c35_b_x[9]);
static void c35_f_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_g_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance);
static void c35_c_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_k, real_T c35_A_data[153], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[81], real_T c35_b_C[81]);
static const mxArray *c35_h_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData);
static int32_T c35_l_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId);
static void c35_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData);
static uint8_T c35_m_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_b_is_active_c35_simulation, const char_T
  *c35_identifier);
static uint8_T c35_n_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId);
static void c35_b_eml_matlab_zgetrf(SFc35_simulationInstanceStruct
  *chartInstance, real_T c35_A[81], int32_T c35_ipiv[9], int32_T *c35_info);
static void c35_b_eml_xger(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_m, int32_T c35_n, real_T c35_alpha1, int32_T c35_ix0, int32_T
  c35_iy0, real_T c35_A[81], int32_T c35_ia0);
static void c35_c_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81]);
static void c35_d_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81]);
static void c35_d_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B_data[17], int32_T c35_B_sizes[1], real_T c35_C
  [9]);
static void c35_b_tanh(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_x[9]);
static void c35_e_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_k, real_T c35_A_data[153], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[81]);
static void init_dsm_address_info(SFc35_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c35_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c35_is_active_c35_simulation = 0U;
}

static void initialize_params_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance)
{
  real_T c35_d0;
  real_T c35_d1;
  real_T c35_dv0[81];
  int32_T c35_i0;
  sf_set_error_prefix_string(
    "Error evaluating data 'nsonar' in the parent workspace.\n");
  sf_mex_import_named("nsonar", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c35_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c35_nsonar = c35_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'i' in the parent workspace.\n");
  sf_mex_import_named("i", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c35_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c35_i = c35_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'kl' in the parent workspace.\n");
  sf_mex_import_named("kl", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      c35_dv0, 0, 0, 0U, 1, 0U, 2, 9, 9);
  for (c35_i0 = 0; c35_i0 < 81; c35_i0++) {
    chartInstance->c35_kl[c35_i0] = c35_dv0[c35_i0];
  }

  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c35_simulation(SFc35_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c35_simulation(SFc35_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c35_update_debugger_state_c35_simulation
  (SFc35_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c35_simulation
  (SFc35_simulationInstanceStruct *chartInstance)
{
  const mxArray *c35_st;
  const mxArray *c35_y = NULL;
  int32_T c35_i1;
  real_T c35_u[81];
  const mxArray *c35_b_y = NULL;
  int32_T c35_i2;
  real_T c35_b_u[4];
  const mxArray *c35_c_y = NULL;
  int32_T c35_i3;
  real_T c35_c_u[9];
  const mxArray *c35_d_y = NULL;
  int32_T c35_i4;
  real_T c35_d_u[2];
  const mxArray *c35_e_y = NULL;
  uint8_T c35_hoistedGlobal;
  uint8_T c35_e_u;
  const mxArray *c35_f_y = NULL;
  real_T (*c35_state_i)[2];
  real_T (*c35_state_est)[9];
  real_T (*c35_cov_i)[4];
  real_T (*c35_Sigma)[81];
  c35_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c35_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c35_state_est = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c35_Sigma = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 1);
  c35_st = NULL;
  c35_st = NULL;
  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_createcellarray(5), FALSE);
  for (c35_i1 = 0; c35_i1 < 81; c35_i1++) {
    c35_u[c35_i1] = (*c35_Sigma)[c35_i1];
  }

  c35_b_y = NULL;
  sf_mex_assign(&c35_b_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 2, 9, 9),
                FALSE);
  sf_mex_setcell(c35_y, 0, c35_b_y);
  for (c35_i2 = 0; c35_i2 < 4; c35_i2++) {
    c35_b_u[c35_i2] = (*c35_cov_i)[c35_i2];
  }

  c35_c_y = NULL;
  sf_mex_assign(&c35_c_y, sf_mex_create("y", c35_b_u, 0, 0U, 1U, 0U, 1, 4),
                FALSE);
  sf_mex_setcell(c35_y, 1, c35_c_y);
  for (c35_i3 = 0; c35_i3 < 9; c35_i3++) {
    c35_c_u[c35_i3] = (*c35_state_est)[c35_i3];
  }

  c35_d_y = NULL;
  sf_mex_assign(&c35_d_y, sf_mex_create("y", c35_c_u, 0, 0U, 1U, 0U, 1, 9),
                FALSE);
  sf_mex_setcell(c35_y, 2, c35_d_y);
  for (c35_i4 = 0; c35_i4 < 2; c35_i4++) {
    c35_d_u[c35_i4] = (*c35_state_i)[c35_i4];
  }

  c35_e_y = NULL;
  sf_mex_assign(&c35_e_y, sf_mex_create("y", c35_d_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c35_y, 3, c35_e_y);
  c35_hoistedGlobal = chartInstance->c35_is_active_c35_simulation;
  c35_e_u = c35_hoistedGlobal;
  c35_f_y = NULL;
  sf_mex_assign(&c35_f_y, sf_mex_create("y", &c35_e_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c35_y, 4, c35_f_y);
  sf_mex_assign(&c35_st, c35_y, FALSE);
  return c35_st;
}

static void set_sim_state_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_st)
{
  const mxArray *c35_u;
  real_T c35_dv1[81];
  int32_T c35_i5;
  real_T c35_dv2[4];
  int32_T c35_i6;
  real_T c35_dv3[9];
  int32_T c35_i7;
  real_T c35_dv4[2];
  int32_T c35_i8;
  real_T (*c35_Sigma)[81];
  real_T (*c35_cov_i)[4];
  real_T (*c35_state_est)[9];
  real_T (*c35_state_i)[2];
  c35_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c35_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c35_state_est = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c35_Sigma = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c35_doneDoubleBufferReInit = TRUE;
  c35_u = sf_mex_dup(c35_st);
  c35_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c35_u, 0)),
    "Sigma", c35_dv1);
  for (c35_i5 = 0; c35_i5 < 81; c35_i5++) {
    (*c35_Sigma)[c35_i5] = c35_dv1[c35_i5];
  }

  c35_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c35_u, 1)),
                       "cov_i", c35_dv2);
  for (c35_i6 = 0; c35_i6 < 4; c35_i6++) {
    (*c35_cov_i)[c35_i6] = c35_dv2[c35_i6];
  }

  c35_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c35_u, 2)),
    "state_est", c35_dv3);
  for (c35_i7 = 0; c35_i7 < 9; c35_i7++) {
    (*c35_state_est)[c35_i7] = c35_dv3[c35_i7];
  }

  c35_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c35_u, 3)),
    "state_i", c35_dv4);
  for (c35_i8 = 0; c35_i8 < 2; c35_i8++) {
    (*c35_state_i)[c35_i8] = c35_dv4[c35_i8];
  }

  chartInstance->c35_is_active_c35_simulation = c35_m_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c35_u, 4)),
     "is_active_c35_simulation");
  sf_mex_destroy(&c35_u);
  c35_update_debugger_state_c35_simulation(chartInstance);
  sf_mex_destroy(&c35_st);
}

static void finalize_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c35_simulation(SFc35_simulationInstanceStruct *chartInstance)
{
  int32_T c35_i9;
  int32_T c35_i10;
  int32_T c35_i11;
  int32_T c35_i12;
  int32_T c35_i13;
  int32_T c35_i14;
  int32_T c35_i15;
  int32_T c35_i16;
  int32_T c35_i17;
  real_T *c35_z;
  real_T *c35_j;
  real_T (*c35_cov_i)[4];
  real_T (*c35_state_i)[2];
  real_T (*c35_statel)[2];
  real_T (*c35_covl)[4];
  real_T (*c35_state_est)[9];
  real_T (*c35_Sigma_predict)[81];
  real_T (*c35_state_predict)[9];
  real_T (*c35_Sigma)[81];
  c35_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c35_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c35_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c35_statel = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
  c35_covl = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c35_state_est = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c35_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c35_Sigma_predict = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 1);
  c35_state_predict = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  c35_Sigma = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 33U, chartInstance->c35_sfEvent);
  for (c35_i9 = 0; c35_i9 < 81; c35_i9++) {
    _SFD_DATA_RANGE_CHECK((*c35_Sigma)[c35_i9], 0U);
  }

  for (c35_i10 = 0; c35_i10 < 9; c35_i10++) {
    _SFD_DATA_RANGE_CHECK((*c35_state_predict)[c35_i10], 1U);
  }

  for (c35_i11 = 0; c35_i11 < 81; c35_i11++) {
    _SFD_DATA_RANGE_CHECK((*c35_Sigma_predict)[c35_i11], 2U);
  }

  _SFD_DATA_RANGE_CHECK(*c35_z, 3U);
  for (c35_i12 = 0; c35_i12 < 9; c35_i12++) {
    _SFD_DATA_RANGE_CHECK((*c35_state_est)[c35_i12], 4U);
  }

  for (c35_i13 = 0; c35_i13 < 4; c35_i13++) {
    _SFD_DATA_RANGE_CHECK((*c35_covl)[c35_i13], 5U);
  }

  for (c35_i14 = 0; c35_i14 < 2; c35_i14++) {
    _SFD_DATA_RANGE_CHECK((*c35_statel)[c35_i14], 6U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c35_nsonar, 7U);
  _SFD_DATA_RANGE_CHECK(*c35_j, 8U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c35_i, 9U);
  for (c35_i15 = 0; c35_i15 < 2; c35_i15++) {
    _SFD_DATA_RANGE_CHECK((*c35_state_i)[c35_i15], 10U);
  }

  for (c35_i16 = 0; c35_i16 < 4; c35_i16++) {
    _SFD_DATA_RANGE_CHECK((*c35_cov_i)[c35_i16], 11U);
  }

  for (c35_i17 = 0; c35_i17 < 81; c35_i17++) {
    _SFD_DATA_RANGE_CHECK(chartInstance->c35_kl[c35_i17], 12U);
  }

  chartInstance->c35_sfEvent = CALL_EVENT;
  c35_chartstep_c35_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c35_chartstep_c35_simulation(SFc35_simulationInstanceStruct
  *chartInstance)
{
  real_T c35_hoistedGlobal;
  real_T c35_b_hoistedGlobal;
  real_T c35_c_hoistedGlobal;
  real_T c35_d_hoistedGlobal;
  int32_T c35_i18;
  real_T c35_state_predict[9];
  int32_T c35_i19;
  real_T c35_Sigma_predict[81];
  real_T c35_z;
  int32_T c35_i20;
  real_T c35_covl[4];
  int32_T c35_i21;
  real_T c35_statel[2];
  real_T c35_b_nsonar;
  real_T c35_j;
  real_T c35_b_i;
  int32_T c35_i22;
  real_T c35_b_kl[81];
  uint32_T c35_debug_family_var_map[31];
  real_T c35_l;
  real_T c35_covll[4];
  real_T c35_Pl[4];
  real_T c35_xh;
  real_T c35_yh;
  real_T c35_xlh;
  real_T c35_ylh;
  real_T c35_dxhat;
  real_T c35_dyhat;
  real_T c35_rhat;
  int32_T c35_H_sizes[2];
  real_T c35_H_data[17];
  real_T c35_R;
  real_T c35_S;
  real_T c35_K[9];
  real_T c35_lambda;
  real_T c35_cov11[4];
  real_T c35_nargin = 9.0;
  real_T c35_nargout = 4.0;
  real_T c35_Sigma[81];
  real_T c35_state_est[9];
  real_T c35_state_i[2];
  real_T c35_cov_i[4];
  int32_T c35_i23;
  int32_T c35_i24;
  const mxArray *c35_y = NULL;
  int32_T c35_i25;
  real_T c35_b;
  real_T c35_b_y;
  int32_T c35_i26;
  int32_T c35_i27;
  real_T c35_x[4];
  int32_T c35_k;
  int32_T c35_b_k;
  real_T c35_b_b;
  real_T c35_d2;
  real_T c35_c_b;
  real_T c35_d3;
  real_T c35_d_b;
  real_T c35_d4;
  real_T c35_e_b;
  real_T c35_d5;
  int32_T c35_i28;
  int32_T c35_i29;
  real_T c35_f_b;
  real_T c35_c_y;
  real_T c35_g_b;
  real_T c35_d_y;
  real_T c35_b_x;
  real_T c35_c_x;
  real_T c35_A;
  real_T c35_B;
  real_T c35_d_x;
  real_T c35_e_y;
  real_T c35_e_x;
  real_T c35_f_y;
  real_T c35_g_y;
  real_T c35_b_A;
  real_T c35_b_B;
  real_T c35_f_x;
  real_T c35_h_y;
  real_T c35_g_x;
  real_T c35_i_y;
  real_T c35_j_y;
  real_T c35_h_b;
  real_T c35_k_y;
  real_T c35_c_A;
  real_T c35_c_B;
  real_T c35_h_x;
  real_T c35_l_y;
  real_T c35_i_x;
  real_T c35_m_y;
  real_T c35_n_y;
  real_T c35_d_A;
  real_T c35_d_B;
  real_T c35_j_x;
  real_T c35_o_y;
  real_T c35_k_x;
  real_T c35_p_y;
  real_T c35_q_y;
  real_T c35_i_b;
  real_T c35_r_y;
  int32_T c35_iv0[2];
  int32_T c35_iv1[2];
  int32_T c35_loop_ub;
  int32_T c35_i30;
  int32_T c35_b_loop_ub;
  int32_T c35_i31;
  real_T c35_j_b;
  int32_T c35_a_sizes[2];
  int32_T c35_a;
  int32_T c35_b_a;
  int32_T c35_c_loop_ub;
  int32_T c35_i32;
  real_T c35_a_data[17];
  int32_T c35_i33;
  real_T c35_k_b[81];
  boolean_T c35_innerDimOk;
  int32_T c35_i34;
  static char_T c35_cv0[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  char_T c35_u[21];
  const mxArray *c35_s_y = NULL;
  int32_T c35_c_k;
  int32_T c35_b_a_sizes[2];
  int32_T c35_c_a;
  int32_T c35_d_a;
  int32_T c35_d_loop_ub;
  int32_T c35_i35;
  real_T c35_b_a_data[17];
  int32_T c35_i36;
  real_T c35_l_b[81];
  real_T c35_t_y[9];
  int32_T c35_b_sizes;
  int32_T c35_e_loop_ub;
  int32_T c35_i37;
  real_T c35_b_data[17];
  boolean_T c35_b_innerDimOk;
  int32_T c35_i38;
  char_T c35_b_u[21];
  const mxArray *c35_u_y = NULL;
  real_T c35_v_y;
  int32_T c35_d_k;
  int32_T c35_e_k;
  int32_T c35_i39;
  real_T c35_c_kl[81];
  int32_T c35_i40;
  real_T c35_b_Sigma_predict[81];
  int32_T c35_f_loop_ub;
  int32_T c35_i41;
  boolean_T c35_c_innerDimOk;
  int32_T c35_i42;
  char_T c35_c_u[21];
  const mxArray *c35_w_y = NULL;
  int32_T c35_i43;
  real_T c35_x_y[9];
  int32_T c35_i44;
  real_T c35_m_b[81];
  int32_T c35_b_b_sizes;
  int32_T c35_g_loop_ub;
  int32_T c35_i45;
  real_T c35_b_b_data[17];
  real_T c35_e_B;
  real_T c35_y_y;
  real_T c35_ab_y;
  int32_T c35_i46;
  int32_T c35_i47;
  int32_T c35_i48;
  int32_T c35_i49;
  int32_T c35_i50;
  real_T c35_C[9];
  int32_T c35_i51;
  int32_T c35_i52;
  int32_T c35_i53;
  int32_T c35_i54;
  int32_T c35_i55;
  int32_T c35_i56;
  int32_T c35_i57;
  real_T c35_n_b;
  int32_T c35_i58;
  int32_T c35_i59;
  int32_T c35_i60;
  int32_T c35_e_a;
  int32_T c35_f_a;
  int32_T c35_h_loop_ub;
  int32_T c35_i61;
  int32_T c35_y_sizes[2];
  int32_T c35_i62;
  int32_T c35_i_loop_ub;
  int32_T c35_i63;
  real_T c35_y_data[153];
  int32_T c35_i64;
  boolean_T c35_d_innerDimOk;
  int32_T c35_i65;
  char_T c35_d_u[21];
  const mxArray *c35_bb_y = NULL;
  int32_T c35_f_k;
  int32_T c35_i66;
  real_T c35_cb_y[81];
  int32_T c35_b_y_sizes[2];
  int32_T c35_db_y;
  int32_T c35_eb_y;
  int32_T c35_j_loop_ub;
  int32_T c35_i67;
  real_T c35_b_y_data[153];
  int32_T c35_i68;
  real_T c35_o_b[81];
  int32_T c35_i69;
  int32_T c35_i70;
  int32_T c35_i71;
  int32_T c35_i72;
  int32_T c35_i73;
  int32_T c35_i74;
  int32_T c35_i75;
  int32_T c35_i76;
  int32_T c35_i77;
  int32_T c35_i78;
  int32_T c35_i79;
  real_T *c35_b_z;
  real_T *c35_b_j;
  real_T (*c35_b_Sigma)[81];
  real_T (*c35_b_state_est)[9];
  real_T (*c35_b_state_i)[2];
  real_T (*c35_b_cov_i)[4];
  real_T (*c35_b_statel)[2];
  real_T (*c35_b_covl)[4];
  real_T (*c35_c_Sigma_predict)[81];
  real_T (*c35_b_state_predict)[9];
  boolean_T guard1 = FALSE;
  c35_b_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c35_b_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c35_b_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c35_b_statel = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
  c35_b_covl = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
  c35_b_state_est = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c35_b_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c35_c_Sigma_predict = (real_T (*)[81])ssGetInputPortSignal(chartInstance->S, 1);
  c35_b_state_predict = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  c35_b_Sigma = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 33U, chartInstance->c35_sfEvent);
  c35_hoistedGlobal = *c35_b_z;
  c35_b_hoistedGlobal = chartInstance->c35_nsonar;
  c35_c_hoistedGlobal = *c35_b_j;
  c35_d_hoistedGlobal = chartInstance->c35_i;
  for (c35_i18 = 0; c35_i18 < 9; c35_i18++) {
    c35_state_predict[c35_i18] = (*c35_b_state_predict)[c35_i18];
  }

  for (c35_i19 = 0; c35_i19 < 81; c35_i19++) {
    c35_Sigma_predict[c35_i19] = (*c35_c_Sigma_predict)[c35_i19];
  }

  c35_z = c35_hoistedGlobal;
  for (c35_i20 = 0; c35_i20 < 4; c35_i20++) {
    c35_covl[c35_i20] = (*c35_b_covl)[c35_i20];
  }

  for (c35_i21 = 0; c35_i21 < 2; c35_i21++) {
    c35_statel[c35_i21] = (*c35_b_statel)[c35_i21];
  }

  c35_b_nsonar = c35_b_hoistedGlobal;
  c35_j = c35_c_hoistedGlobal;
  c35_b_i = c35_d_hoistedGlobal;
  for (c35_i22 = 0; c35_i22 < 81; c35_i22++) {
    c35_b_kl[c35_i22] = chartInstance->c35_kl[c35_i22];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 31U, 31U, c35_debug_family_names,
    c35_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_l, 0U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_covll, 1U, c35_sf_marshallOut,
    c35_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_Pl, 2U, c35_f_sf_marshallOut,
    c35_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_xh, 3U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_yh, 4U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_xlh, 5U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_ylh, 6U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_dxhat, 7U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_dyhat, 8U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_rhat, 9U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_DYN_IMPORTABLE(c35_H_data, (const int32_T *)
    &c35_H_sizes, NULL, 0, 10, (void *)c35_g_sf_marshallOut, (void *)
    c35_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_R, 11U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_S, 12U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_K, 13U, c35_c_sf_marshallOut,
    c35_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_lambda, 14U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_cov11, 15U, c35_f_sf_marshallOut,
    c35_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_nargin, 16U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_nargout, 17U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c35_state_predict, 18U, c35_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c35_Sigma_predict, 19U, c35_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c35_z, 20U, c35_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c35_covl, 21U, c35_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c35_statel, 22U, c35_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_b_nsonar, 23U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c35_j, 24U, c35_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c35_b_i, 25U, c35_e_sf_marshallOut,
    c35_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_b_kl, 26U, c35_d_sf_marshallOut,
    c35_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_Sigma, 27U, c35_d_sf_marshallOut,
    c35_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_state_est, 28U, c35_c_sf_marshallOut,
    c35_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_state_i, 29U, c35_b_sf_marshallOut,
    c35_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c35_cov_i, 30U, c35_sf_marshallOut,
    c35_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 5);
  for (c35_i23 = 0; c35_i23 < 9; c35_i23++) {
    c35_state_est[c35_i23] = c35_state_predict[c35_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 7);
  for (c35_i24 = 0; c35_i24 < 81; c35_i24++) {
    c35_Sigma[c35_i24] = c35_Sigma_predict[c35_i24];
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 9);
  if (CV_EML_IF(0, 1, 0, c35_j == c35_b_i)) {
    _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 10);
    c35_l = 0.0;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 11);
    if (CV_EML_IF(0, 1, 1, c35_j > c35_b_i)) {
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 12);
      c35_l = c35_j - 1.0;
    } else {
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 14);
      c35_l = c35_j;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 16);
  if (c35_l <= 5.0) {
  } else {
    c35_y = NULL;
    sf_mex_assign(&c35_y, sf_mex_create("y", "Assertion failed.", 15, 0U, 0U, 0U,
      2, 1, strlen("Assertion failed.")), FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, c35_y);
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 17);
  CV_EML_COND(0, 1, 0, c35_l > 5.0);
  if (CV_EML_COND(0, 1, 1, c35_l < 1.0)) {
    CV_EML_MCDC(0, 1, 0, TRUE);
    CV_EML_IF(0, 1, 2, TRUE);
    _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 18);
    c35_l = 1.0;
  } else {
    CV_EML_MCDC(0, 1, 0, FALSE);
    CV_EML_IF(0, 1, 2, FALSE);
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 20);
  for (c35_i25 = 0; c35_i25 < 4; c35_i25++) {
    c35_covll[c35_i25] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 23);
  guard1 = FALSE;
  if (CV_EML_COND(0, 1, 2, c35_b_i != c35_j)) {
    c35_b = c35_b_nsonar;
    c35_b_y = 20.0 * c35_b;
    if (CV_EML_COND(0, 1, 3, c35_z > c35_b_y)) {
      CV_EML_MCDC(0, 1, 1, TRUE);
      CV_EML_IF(0, 1, 3, TRUE);
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 26);
      for (c35_i26 = 0; c35_i26 < 4; c35_i26++) {
        c35_covll[c35_i26] = c35_covl[c35_i26];
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 27);
      for (c35_i27 = 0; c35_i27 < 4; c35_i27++) {
        c35_x[c35_i27] = c35_covll[c35_i27];
      }

      for (c35_k = 1; c35_k < 5; c35_k++) {
        c35_b_k = c35_k - 1;
        c35_Pl[c35_b_k] = c35_x[c35_b_k];
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 28);
      c35_b_b = c35_l;
      c35_d2 = 2.0 * c35_b_b;
      c35_c_b = c35_l;
      c35_d3 = 2.0 * c35_c_b;
      if (c35_d2 == c35_d3) {
      } else {
        _SFD_RUNTIME_ERROR_MSGID(
          "EMLRT:runTime:RepeatedExprWithDifferentResultsInColonIndexing");
      }

      c35_d_b = c35_l;
      c35_d4 = 2.0 * c35_d_b;
      c35_e_b = c35_l;
      c35_d5 = 2.0 * c35_e_b;
      if (c35_d4 == c35_d5) {
      } else {
        _SFD_RUNTIME_ERROR_MSGID(
          "EMLRT:runTime:RepeatedExprWithDifferentResultsInColonIndexing");
      }

      for (c35_i28 = 0; c35_i28 < 2; c35_i28++) {
        for (c35_i29 = 0; c35_i29 < 2; c35_i29++) {
          c35_Sigma_predict[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
            "Sigma_predict", (int32_T)_SFD_INTEGER_CHECK("2*l+2:2*l+3", c35_d2 +
                               (2.0 + (real_T)c35_i29)), 1, 9, 1, 0) + 9 *
                             ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
            "Sigma_predict", (int32_T)_SFD_INTEGER_CHECK("2*l+2:2*l+3", c35_d4 +
                                (2.0 + (real_T)c35_i28)), 1, 9, 2, 0) - 1)) - 1]
            = c35_Pl[c35_i29 + (c35_i28 << 1)];
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 32);
      c35_xh = c35_state_predict[0];
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 33);
      c35_yh = c35_state_predict[1];
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 34);
      c35_xlh = c35_statel[0];
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 35);
      c35_ylh = c35_statel[1];
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 36);
      c35_f_b = c35_l;
      c35_c_y = 2.0 * c35_f_b;
      c35_state_predict[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
        "state_predict", (int32_T)_SFD_INTEGER_CHECK("2*l+2", c35_c_y + 2.0), 1,
        9, 1, 0) - 1] = c35_xlh;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 37);
      c35_g_b = c35_l;
      c35_d_y = 2.0 * c35_g_b;
      c35_state_predict[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
        "state_predict", (int32_T)_SFD_INTEGER_CHECK("2*l+3", c35_d_y + 3.0), 1,
        9, 1, 0) - 1] = c35_ylh;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 40);
      c35_dxhat = c35_xh - c35_xlh;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 41);
      c35_dyhat = c35_yh - c35_ylh;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 43);
      c35_b_x = c35_mpower(chartInstance, c35_dxhat) + c35_mpower(chartInstance,
        c35_dyhat);
      c35_rhat = c35_b_x;
      if (c35_rhat < 0.0) {
        c35_eml_error(chartInstance);
      }

      c35_c_x = c35_rhat;
      c35_rhat = c35_c_x;
      c35_rhat = muDoubleScalarSqrt(c35_rhat);
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 44);
      c35_A = c35_dxhat;
      c35_B = c35_rhat;
      c35_d_x = c35_A;
      c35_e_y = c35_B;
      c35_e_x = c35_d_x;
      c35_f_y = c35_e_y;
      c35_g_y = c35_e_x / c35_f_y;
      c35_b_A = c35_dyhat;
      c35_b_B = c35_rhat;
      c35_f_x = c35_b_A;
      c35_h_y = c35_b_B;
      c35_g_x = c35_f_x;
      c35_i_y = c35_h_y;
      c35_j_y = c35_g_x / c35_i_y;
      c35_h_b = c35_l;
      c35_k_y = 2.0 * c35_h_b;
      c35_c_A = -c35_dxhat;
      c35_c_B = c35_rhat;
      c35_h_x = c35_c_A;
      c35_l_y = c35_c_B;
      c35_i_x = c35_h_x;
      c35_m_y = c35_l_y;
      c35_n_y = c35_i_x / c35_m_y;
      c35_d_A = -c35_dyhat;
      c35_d_B = c35_rhat;
      c35_j_x = c35_d_A;
      c35_o_y = c35_d_B;
      c35_k_x = c35_j_x;
      c35_p_y = c35_o_y;
      c35_q_y = c35_k_x / c35_p_y;
      c35_i_b = c35_l;
      c35_r_y = 2.0 * c35_i_b;
      c35_iv0[0] = 1;
      c35_iv0[1] = (int32_T)_SFD_INTEGER_CHECK("2*l-2", c35_k_y - 2.0);
      c35_iv1[0] = 1;
      c35_iv1[1] = (int32_T)_SFD_INTEGER_CHECK("4-(2*l-2)",
        _SFD_NON_NEGATIVE_CHECK("4-(2*l-2)", 4.0 - (c35_r_y - 2.0)));
      c35_H_sizes[0] = 1;
      c35_H_sizes[1] = (c35_iv0[1] + c35_iv1[1]) + 5;
      c35_H_data[0] = c35_g_y;
      c35_H_data[c35_H_sizes[0]] = c35_j_y;
      c35_H_data[c35_H_sizes[0] << 1] = 0.0;
      c35_loop_ub = c35_iv0[1] - 1;
      for (c35_i30 = 0; c35_i30 <= c35_loop_ub; c35_i30++) {
        c35_H_data[c35_H_sizes[0] * (c35_i30 + 3)] = 0.0;
      }

      c35_H_data[c35_H_sizes[0] * (3 + c35_iv0[1])] = c35_n_y;
      c35_H_data[c35_H_sizes[0] * (c35_iv0[1] + 4)] = c35_q_y;
      c35_b_loop_ub = c35_iv1[1] - 1;
      for (c35_i31 = 0; c35_i31 <= c35_b_loop_ub; c35_i31++) {
        c35_H_data[c35_H_sizes[0] * ((c35_i31 + c35_iv0[1]) + 5)] = 0.0;
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 47);
      c35_j_b = c35_b_nsonar;
      c35_R = 100.0 * c35_j_b;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 50);
      c35_a_sizes[0] = 1;
      c35_a_sizes[1] = c35_H_sizes[1];
      c35_a = c35_a_sizes[0];
      c35_b_a = c35_a_sizes[1];
      c35_c_loop_ub = c35_H_sizes[0] * c35_H_sizes[1] - 1;
      for (c35_i32 = 0; c35_i32 <= c35_c_loop_ub; c35_i32++) {
        c35_a_data[c35_i32] = c35_H_data[c35_i32];
      }

      for (c35_i33 = 0; c35_i33 < 81; c35_i33++) {
        c35_k_b[c35_i33] = c35_Sigma_predict[c35_i33];
      }

      c35_innerDimOk = ((real_T)c35_a_sizes[1] == 9.0);
      if (!c35_innerDimOk) {
        for (c35_i34 = 0; c35_i34 < 21; c35_i34++) {
          c35_u[c35_i34] = c35_cv0[c35_i34];
        }

        c35_s_y = NULL;
        sf_mex_assign(&c35_s_y, sf_mex_create("y", c35_u, 10, 0U, 1U, 0U, 2, 1,
          21), FALSE);
        sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U,
          1U, 14, c35_s_y));
      }

      c35_c_k = (int32_T)(real_T)c35_a_sizes[1];
      c35_b_eml_scalar_eg(chartInstance);
      c35_b_eml_scalar_eg(chartInstance);
      c35_b_a_sizes[0] = 1;
      c35_b_a_sizes[1] = c35_a_sizes[1];
      c35_c_a = c35_b_a_sizes[0];
      c35_d_a = c35_b_a_sizes[1];
      c35_d_loop_ub = c35_a_sizes[0] * c35_a_sizes[1] - 1;
      for (c35_i35 = 0; c35_i35 <= c35_d_loop_ub; c35_i35++) {
        c35_b_a_data[c35_i35] = c35_a_data[c35_i35];
      }

      for (c35_i36 = 0; c35_i36 < 81; c35_i36++) {
        c35_l_b[c35_i36] = c35_k_b[c35_i36];
      }

      c35_eml_xgemm(chartInstance, c35_c_k, c35_b_a_data, c35_b_a_sizes, c35_l_b,
                    c35_c_k, c35_t_y);
      c35_b_sizes = c35_H_sizes[1];
      c35_e_loop_ub = c35_H_sizes[1] - 1;
      for (c35_i37 = 0; c35_i37 <= c35_e_loop_ub; c35_i37++) {
        c35_b_data[c35_i37] = c35_H_data[c35_H_sizes[0] * c35_i37];
      }

      c35_b_innerDimOk = (9.0 == (real_T)c35_b_sizes);
      if (!c35_b_innerDimOk) {
        for (c35_i38 = 0; c35_i38 < 21; c35_i38++) {
          c35_b_u[c35_i38] = c35_cv0[c35_i38];
        }

        c35_u_y = NULL;
        sf_mex_assign(&c35_u_y, sf_mex_create("y", c35_b_u, 10, 0U, 1U, 0U, 2, 1,
          21), FALSE);
        sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U,
          1U, 14, c35_u_y));
      }

      c35_c_eml_scalar_eg(chartInstance);
      c35_c_eml_scalar_eg(chartInstance);
      c35_v_y = 0.0;
      for (c35_d_k = 1; c35_d_k < 10; c35_d_k++) {
        c35_e_k = c35_d_k;
        c35_v_y += c35_t_y[c35_e_k - 1] * c35_b_data[_SFD_EML_ARRAY_BOUNDS_CHECK
          ("", c35_e_k, 1, c35_b_sizes, 1, 0) - 1];
      }

      c35_S = c35_v_y + c35_R;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 52);
      for (c35_i39 = 0; c35_i39 < 81; c35_i39++) {
        c35_c_kl[c35_i39] = c35_b_kl[c35_i39];
      }

      for (c35_i40 = 0; c35_i40 < 81; c35_i40++) {
        c35_b_Sigma_predict[c35_i40] = c35_Sigma_predict[c35_i40];
      }

      c35_mldivide(chartInstance, c35_c_kl, c35_b_Sigma_predict, c35_k_b);
      c35_b_sizes = c35_H_sizes[1];
      c35_f_loop_ub = c35_H_sizes[1] - 1;
      for (c35_i41 = 0; c35_i41 <= c35_f_loop_ub; c35_i41++) {
        c35_b_data[c35_i41] = c35_H_data[c35_H_sizes[0] * c35_i41];
      }

      c35_c_innerDimOk = (9.0 == (real_T)c35_b_sizes);
      if (!c35_c_innerDimOk) {
        for (c35_i42 = 0; c35_i42 < 21; c35_i42++) {
          c35_c_u[c35_i42] = c35_cv0[c35_i42];
        }

        c35_w_y = NULL;
        sf_mex_assign(&c35_w_y, sf_mex_create("y", c35_c_u, 10, 0U, 1U, 0U, 2, 1,
          21), FALSE);
        sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U,
          1U, 14, c35_w_y));
      }

      c35_d_eml_scalar_eg(chartInstance);
      c35_d_eml_scalar_eg(chartInstance);
      for (c35_i43 = 0; c35_i43 < 9; c35_i43++) {
        c35_x_y[c35_i43] = 0.0;
      }

      for (c35_i44 = 0; c35_i44 < 81; c35_i44++) {
        c35_m_b[c35_i44] = c35_k_b[c35_i44];
      }

      c35_b_b_sizes = c35_b_sizes;
      c35_g_loop_ub = c35_b_sizes - 1;
      for (c35_i45 = 0; c35_i45 <= c35_g_loop_ub; c35_i45++) {
        c35_b_b_data[c35_i45] = c35_b_data[c35_i45];
      }

      c35_d_eml_xgemm(chartInstance, c35_m_b, c35_b_b_data, *(int32_T (*)[1])&
                      c35_b_b_sizes, c35_x_y);
      c35_e_B = c35_S;
      c35_y_y = c35_e_B;
      c35_ab_y = c35_y_y;
      for (c35_i46 = 0; c35_i46 < 9; c35_i46++) {
        c35_x_y[c35_i46] /= c35_ab_y;
      }

      for (c35_i47 = 0; c35_i47 < 81; c35_i47++) {
        c35_k_b[c35_i47] = c35_b_kl[c35_i47];
      }

      c35_b_tanh(chartInstance, c35_x_y);
      c35_f_eml_scalar_eg(chartInstance);
      c35_f_eml_scalar_eg(chartInstance);
      for (c35_i48 = 0; c35_i48 < 9; c35_i48++) {
        c35_K[c35_i48] = 0.0;
      }

      for (c35_i49 = 0; c35_i49 < 9; c35_i49++) {
        c35_K[c35_i49] = 0.0;
      }

      for (c35_i50 = 0; c35_i50 < 9; c35_i50++) {
        c35_C[c35_i50] = c35_K[c35_i50];
      }

      for (c35_i51 = 0; c35_i51 < 9; c35_i51++) {
        c35_K[c35_i51] = c35_C[c35_i51];
      }

      for (c35_i52 = 0; c35_i52 < 9; c35_i52++) {
        c35_C[c35_i52] = c35_K[c35_i52];
      }

      for (c35_i53 = 0; c35_i53 < 9; c35_i53++) {
        c35_K[c35_i53] = c35_C[c35_i53];
      }

      for (c35_i54 = 0; c35_i54 < 9; c35_i54++) {
        c35_K[c35_i54] = 0.0;
        c35_i55 = 0;
        for (c35_i56 = 0; c35_i56 < 9; c35_i56++) {
          c35_K[c35_i54] += c35_k_b[c35_i55 + c35_i54] * c35_x_y[c35_i56];
          c35_i55 += 9;
        }
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 54);
      c35_lambda = c35_z - c35_rhat;
      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 56);
      for (c35_i57 = 0; c35_i57 < 9; c35_i57++) {
        c35_x_y[c35_i57] = c35_K[c35_i57];
      }

      c35_n_b = c35_lambda;
      for (c35_i58 = 0; c35_i58 < 9; c35_i58++) {
        c35_x_y[c35_i58] *= c35_n_b;
      }

      for (c35_i59 = 0; c35_i59 < 9; c35_i59++) {
        c35_state_est[c35_i59] = c35_state_predict[c35_i59] + c35_x_y[c35_i59];
      }

      _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 58);
      for (c35_i60 = 0; c35_i60 < 9; c35_i60++) {
        c35_x_y[c35_i60] = c35_K[c35_i60];
      }

      c35_a_sizes[0] = 1;
      c35_a_sizes[1] = c35_H_sizes[1];
      c35_e_a = c35_a_sizes[0];
      c35_f_a = c35_a_sizes[1];
      c35_h_loop_ub = c35_H_sizes[0] * c35_H_sizes[1] - 1;
      for (c35_i61 = 0; c35_i61 <= c35_h_loop_ub; c35_i61++) {
        c35_a_data[c35_i61] = c35_H_data[c35_i61];
      }

      c35_y_sizes[0] = 9;
      c35_y_sizes[1] = c35_a_sizes[1];
      for (c35_i62 = 0; c35_i62 < 9; c35_i62++) {
        c35_i_loop_ub = c35_a_sizes[1] - 1;
        for (c35_i63 = 0; c35_i63 <= c35_i_loop_ub; c35_i63++) {
          c35_y_data[c35_i62 + c35_y_sizes[0] * c35_i63] = c35_x_y[c35_i62] *
            c35_a_data[c35_a_sizes[0] * c35_i63];
        }
      }

      for (c35_i64 = 0; c35_i64 < 81; c35_i64++) {
        c35_k_b[c35_i64] = c35_Sigma_predict[c35_i64];
      }

      c35_d_innerDimOk = ((real_T)c35_y_sizes[1] == 9.0);
      if (!c35_d_innerDimOk) {
        for (c35_i65 = 0; c35_i65 < 21; c35_i65++) {
          c35_d_u[c35_i65] = c35_cv0[c35_i65];
        }

        c35_bb_y = NULL;
        sf_mex_assign(&c35_bb_y, sf_mex_create("y", c35_d_u, 10, 0U, 1U, 0U, 2,
          1, 21), FALSE);
        sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U,
          1U, 14, c35_bb_y));
      }

      c35_f_k = (int32_T)(real_T)c35_y_sizes[1];
      c35_g_eml_scalar_eg(chartInstance);
      c35_g_eml_scalar_eg(chartInstance);
      for (c35_i66 = 0; c35_i66 < 81; c35_i66++) {
        c35_cb_y[c35_i66] = 0.0;
      }

      c35_b_y_sizes[0] = 9;
      c35_b_y_sizes[1] = c35_y_sizes[1];
      c35_db_y = c35_b_y_sizes[0];
      c35_eb_y = c35_b_y_sizes[1];
      c35_j_loop_ub = c35_y_sizes[0] * c35_y_sizes[1] - 1;
      for (c35_i67 = 0; c35_i67 <= c35_j_loop_ub; c35_i67++) {
        c35_b_y_data[c35_i67] = c35_y_data[c35_i67];
      }

      for (c35_i68 = 0; c35_i68 < 81; c35_i68++) {
        c35_o_b[c35_i68] = c35_k_b[c35_i68];
      }

      c35_e_eml_xgemm(chartInstance, c35_f_k, c35_b_y_data, c35_b_y_sizes,
                      c35_o_b, c35_f_k, c35_cb_y);
      for (c35_i69 = 0; c35_i69 < 81; c35_i69++) {
        c35_Sigma[c35_i69] = c35_Sigma_predict[c35_i69] - c35_cb_y[c35_i69];
      }
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    CV_EML_MCDC(0, 1, 1, FALSE);
    CV_EML_IF(0, 1, 3, FALSE);
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 62);
  for (c35_i70 = 0; c35_i70 < 2; c35_i70++) {
    c35_state_i[c35_i70] = c35_state_est[c35_i70];
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 63);
  c35_i71 = 0;
  c35_i72 = 0;
  for (c35_i73 = 0; c35_i73 < 2; c35_i73++) {
    for (c35_i74 = 0; c35_i74 < 2; c35_i74++) {
      c35_cov11[c35_i74 + c35_i71] = c35_Sigma[c35_i74 + c35_i72];
    }

    c35_i71 += 2;
    c35_i72 += 9;
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, 64);
  for (c35_i75 = 0; c35_i75 < 4; c35_i75++) {
    c35_cov_i[c35_i75] = c35_cov11[c35_i75];
  }

  _SFD_EML_CALL(0U, chartInstance->c35_sfEvent, -64);
  _SFD_SYMBOL_SCOPE_POP();
  for (c35_i76 = 0; c35_i76 < 81; c35_i76++) {
    (*c35_b_Sigma)[c35_i76] = c35_Sigma[c35_i76];
  }

  for (c35_i77 = 0; c35_i77 < 9; c35_i77++) {
    (*c35_b_state_est)[c35_i77] = c35_state_est[c35_i77];
  }

  for (c35_i78 = 0; c35_i78 < 2; c35_i78++) {
    (*c35_b_state_i)[c35_i78] = c35_state_i[c35_i78];
  }

  for (c35_i79 = 0; c35_i79 < 4; c35_i79++) {
    (*c35_b_cov_i)[c35_i79] = c35_cov_i[c35_i79];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 33U, chartInstance->c35_sfEvent);
}

static void initSimStructsc35_simulation(SFc35_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c35_machineNumber, uint32_T
  c35_chartNumber)
{
}

static const mxArray *c35_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_i80;
  real_T c35_b_inData[4];
  int32_T c35_i81;
  real_T c35_u[4];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  for (c35_i80 = 0; c35_i80 < 4; c35_i80++) {
    c35_b_inData[c35_i80] = (*(real_T (*)[4])c35_inData)[c35_i80];
  }

  for (c35_i81 = 0; c35_i81 < 4; c35_i81++) {
    c35_u[c35_i81] = c35_b_inData[c35_i81];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_cov_i, const char_T *c35_identifier, real_T c35_y[4])
{
  emlrtMsgIdentifier c35_thisId;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_cov_i), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_cov_i);
}

static void c35_b_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[4])
{
  real_T c35_dv5[4];
  int32_T c35_i82;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), c35_dv5, 1, 0, 0U, 1, 0U, 1, 4);
  for (c35_i82 = 0; c35_i82 < 4; c35_i82++) {
    c35_y[c35_i82] = c35_dv5[c35_i82];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_cov_i;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y[4];
  int32_T c35_i83;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_cov_i = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_cov_i), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_cov_i);
  for (c35_i83 = 0; c35_i83 < 4; c35_i83++) {
    (*(real_T (*)[4])c35_outData)[c35_i83] = c35_y[c35_i83];
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_b_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_i84;
  real_T c35_b_inData[2];
  int32_T c35_i85;
  real_T c35_u[2];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  for (c35_i84 = 0; c35_i84 < 2; c35_i84++) {
    c35_b_inData[c35_i84] = (*(real_T (*)[2])c35_inData)[c35_i84];
  }

  for (c35_i85 = 0; c35_i85 < 2; c35_i85++) {
    c35_u[c35_i85] = c35_b_inData[c35_i85];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_c_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_state_i, const char_T *c35_identifier, real_T c35_y[2])
{
  emlrtMsgIdentifier c35_thisId;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_state_i), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_state_i);
}

static void c35_d_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[2])
{
  real_T c35_dv6[2];
  int32_T c35_i86;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), c35_dv6, 1, 0, 0U, 1, 0U, 1, 2);
  for (c35_i86 = 0; c35_i86 < 2; c35_i86++) {
    c35_y[c35_i86] = c35_dv6[c35_i86];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_state_i;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y[2];
  int32_T c35_i87;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_state_i = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_state_i), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_state_i);
  for (c35_i87 = 0; c35_i87 < 2; c35_i87++) {
    (*(real_T (*)[2])c35_outData)[c35_i87] = c35_y[c35_i87];
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_c_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_i88;
  real_T c35_b_inData[9];
  int32_T c35_i89;
  real_T c35_u[9];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  for (c35_i88 = 0; c35_i88 < 9; c35_i88++) {
    c35_b_inData[c35_i88] = (*(real_T (*)[9])c35_inData)[c35_i88];
  }

  for (c35_i89 = 0; c35_i89 < 9; c35_i89++) {
    c35_u[c35_i89] = c35_b_inData[c35_i89];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 1, 9), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_e_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_state_est, const char_T *c35_identifier, real_T c35_y[9])
{
  emlrtMsgIdentifier c35_thisId;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_state_est), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_state_est);
}

static void c35_f_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[9])
{
  real_T c35_dv7[9];
  int32_T c35_i90;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), c35_dv7, 1, 0, 0U, 1, 0U, 1, 9);
  for (c35_i90 = 0; c35_i90 < 9; c35_i90++) {
    c35_y[c35_i90] = c35_dv7[c35_i90];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_state_est;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y[9];
  int32_T c35_i91;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_state_est = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_state_est), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_state_est);
  for (c35_i91 = 0; c35_i91 < 9; c35_i91++) {
    (*(real_T (*)[9])c35_outData)[c35_i91] = c35_y[c35_i91];
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_d_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_i92;
  int32_T c35_i93;
  int32_T c35_i94;
  real_T c35_b_inData[81];
  int32_T c35_i95;
  int32_T c35_i96;
  int32_T c35_i97;
  real_T c35_u[81];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  c35_i92 = 0;
  for (c35_i93 = 0; c35_i93 < 9; c35_i93++) {
    for (c35_i94 = 0; c35_i94 < 9; c35_i94++) {
      c35_b_inData[c35_i94 + c35_i92] = (*(real_T (*)[81])c35_inData)[c35_i94 +
        c35_i92];
    }

    c35_i92 += 9;
  }

  c35_i95 = 0;
  for (c35_i96 = 0; c35_i96 < 9; c35_i96++) {
    for (c35_i97 = 0; c35_i97 < 9; c35_i97++) {
      c35_u[c35_i97 + c35_i95] = c35_b_inData[c35_i97 + c35_i95];
    }

    c35_i95 += 9;
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 2, 9, 9), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_g_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_Sigma, const char_T *c35_identifier, real_T c35_y[81])
{
  emlrtMsgIdentifier c35_thisId;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_Sigma), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_Sigma);
}

static void c35_h_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[81])
{
  real_T c35_dv8[81];
  int32_T c35_i98;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), c35_dv8, 1, 0, 0U, 1, 0U, 2, 9,
                9);
  for (c35_i98 = 0; c35_i98 < 81; c35_i98++) {
    c35_y[c35_i98] = c35_dv8[c35_i98];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_Sigma;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y[81];
  int32_T c35_i99;
  int32_T c35_i100;
  int32_T c35_i101;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_Sigma = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_Sigma), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_Sigma);
  c35_i99 = 0;
  for (c35_i100 = 0; c35_i100 < 9; c35_i100++) {
    for (c35_i101 = 0; c35_i101 < 9; c35_i101++) {
      (*(real_T (*)[81])c35_outData)[c35_i101 + c35_i99] = c35_y[c35_i101 +
        c35_i99];
    }

    c35_i99 += 9;
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_e_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  real_T c35_u;
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  c35_u = *(real_T *)c35_inData;
  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", &c35_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static real_T c35_i_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId)
{
  real_T c35_y;
  real_T c35_d6;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), &c35_d6, 1, 0, 0U, 0, 0U, 0);
  c35_y = c35_d6;
  sf_mex_destroy(&c35_u);
  return c35_y;
}

static void c35_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_b_i;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_b_i = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_y = c35_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_b_i), &c35_thisId);
  sf_mex_destroy(&c35_b_i);
  *(real_T *)c35_outData = c35_y;
  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_f_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_i102;
  int32_T c35_i103;
  int32_T c35_i104;
  real_T c35_b_inData[4];
  int32_T c35_i105;
  int32_T c35_i106;
  int32_T c35_i107;
  real_T c35_u[4];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  c35_i102 = 0;
  for (c35_i103 = 0; c35_i103 < 2; c35_i103++) {
    for (c35_i104 = 0; c35_i104 < 2; c35_i104++) {
      c35_b_inData[c35_i104 + c35_i102] = (*(real_T (*)[4])c35_inData)[c35_i104
        + c35_i102];
    }

    c35_i102 += 2;
  }

  c35_i105 = 0;
  for (c35_i106 = 0; c35_i106 < 2; c35_i106++) {
    for (c35_i107 = 0; c35_i107 < 2; c35_i107++) {
      c35_u[c35_i107 + c35_i105] = c35_b_inData[c35_i107 + c35_i105];
    }

    c35_i105 += 2;
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_j_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T c35_y[4])
{
  real_T c35_dv9[4];
  int32_T c35_i108;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), c35_dv9, 1, 0, 0U, 1, 0U, 2, 2,
                2);
  for (c35_i108 = 0; c35_i108 < 4; c35_i108++) {
    c35_y[c35_i108] = c35_dv9[c35_i108];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_cov11;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  real_T c35_y[4];
  int32_T c35_i109;
  int32_T c35_i110;
  int32_T c35_i111;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_cov11 = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_cov11), &c35_thisId,
    c35_y);
  sf_mex_destroy(&c35_cov11);
  c35_i109 = 0;
  for (c35_i110 = 0; c35_i110 < 2; c35_i110++) {
    for (c35_i111 = 0; c35_i111 < 2; c35_i111++) {
      (*(real_T (*)[4])c35_outData)[c35_i111 + c35_i109] = c35_y[c35_i111 +
        c35_i109];
    }

    c35_i109 += 2;
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

static const mxArray *c35_g_sf_marshallOut(void *chartInstanceVoid, real_T
  c35_inData_data[17], int32_T c35_inData_sizes[2])
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_b_inData_sizes[2];
  int32_T c35_loop_ub;
  int32_T c35_i112;
  real_T c35_b_inData_data[17];
  int32_T c35_u_sizes[2];
  int32_T c35_b_loop_ub;
  int32_T c35_i113;
  real_T c35_u_data[17];
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  c35_b_inData_sizes[0] = 1;
  c35_b_inData_sizes[1] = c35_inData_sizes[1];
  c35_loop_ub = c35_inData_sizes[1] - 1;
  for (c35_i112 = 0; c35_i112 <= c35_loop_ub; c35_i112++) {
    c35_b_inData_data[c35_b_inData_sizes[0] * c35_i112] =
      c35_inData_data[c35_inData_sizes[0] * c35_i112];
  }

  c35_u_sizes[0] = 1;
  c35_u_sizes[1] = c35_b_inData_sizes[1];
  c35_b_loop_ub = c35_b_inData_sizes[1] - 1;
  for (c35_i113 = 0; c35_i113 <= c35_b_loop_ub; c35_i113++) {
    c35_u_data[c35_u_sizes[0] * c35_i113] =
      c35_b_inData_data[c35_b_inData_sizes[0] * c35_i113];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u_data, 0, 0U, 1U, 0U, 2,
    c35_u_sizes[0], c35_u_sizes[1]), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static void c35_k_emlrt_marshallIn(SFc35_simulationInstanceStruct *chartInstance,
  const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId, real_T
  c35_y_data[17], int32_T c35_y_sizes[2])
{
  int32_T c35_i114;
  uint32_T c35_uv0[2];
  int32_T c35_i115;
  static boolean_T c35_bv0[2] = { FALSE, TRUE };

  boolean_T c35_bv1[2];
  int32_T c35_tmp_sizes[2];
  real_T c35_tmp_data[17];
  int32_T c35_y;
  int32_T c35_b_y;
  int32_T c35_loop_ub;
  int32_T c35_i116;
  for (c35_i114 = 0; c35_i114 < 2; c35_i114++) {
    c35_uv0[c35_i114] = 1U + ((uint32_T)c35_i114 << 4);
  }

  for (c35_i115 = 0; c35_i115 < 2; c35_i115++) {
    c35_bv1[c35_i115] = c35_bv0[c35_i115];
  }

  sf_mex_import_vs(c35_parentId, sf_mex_dup(c35_u), c35_tmp_data, 1, 0, 0U, 1,
                   0U, 2, c35_bv1, c35_uv0, c35_tmp_sizes);
  c35_y_sizes[0] = 1;
  c35_y_sizes[1] = c35_tmp_sizes[1];
  c35_y = c35_y_sizes[0];
  c35_b_y = c35_y_sizes[1];
  c35_loop_ub = c35_tmp_sizes[0] * c35_tmp_sizes[1] - 1;
  for (c35_i116 = 0; c35_i116 <= c35_loop_ub; c35_i116++) {
    c35_y_data[c35_i116] = c35_tmp_data[c35_i116];
  }

  sf_mex_destroy(&c35_u);
}

static void c35_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, real_T c35_outData_data[17],
  int32_T c35_outData_sizes[2])
{
  const mxArray *c35_H;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  int32_T c35_y_sizes[2];
  real_T c35_y_data[17];
  int32_T c35_loop_ub;
  int32_T c35_i117;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_H = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_H), &c35_thisId,
    c35_y_data, c35_y_sizes);
  sf_mex_destroy(&c35_H);
  c35_outData_sizes[0] = 1;
  c35_outData_sizes[1] = c35_y_sizes[1];
  c35_loop_ub = c35_y_sizes[1] - 1;
  for (c35_i117 = 0; c35_i117 <= c35_loop_ub; c35_i117++) {
    c35_outData_data[c35_outData_sizes[0] * c35_i117] = c35_y_data[c35_y_sizes[0]
      * c35_i117];
  }

  sf_mex_destroy(&c35_mxArrayInData);
}

const mxArray *sf_c35_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c35_nameCaptureInfo = NULL;
  c35_nameCaptureInfo = NULL;
  sf_mex_assign(&c35_nameCaptureInfo, sf_mex_createstruct("structure", 2, 187, 1),
                FALSE);
  c35_info_helper(&c35_nameCaptureInfo);
  c35_b_info_helper(&c35_nameCaptureInfo);
  c35_c_info_helper(&c35_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c35_nameCaptureInfo);
  return c35_nameCaptureInfo;
}

static void c35_info_helper(const mxArray **c35_info)
{
  const mxArray *c35_rhs0 = NULL;
  const mxArray *c35_lhs0 = NULL;
  const mxArray *c35_rhs1 = NULL;
  const mxArray *c35_lhs1 = NULL;
  const mxArray *c35_rhs2 = NULL;
  const mxArray *c35_lhs2 = NULL;
  const mxArray *c35_rhs3 = NULL;
  const mxArray *c35_lhs3 = NULL;
  const mxArray *c35_rhs4 = NULL;
  const mxArray *c35_lhs4 = NULL;
  const mxArray *c35_rhs5 = NULL;
  const mxArray *c35_lhs5 = NULL;
  const mxArray *c35_rhs6 = NULL;
  const mxArray *c35_lhs6 = NULL;
  const mxArray *c35_rhs7 = NULL;
  const mxArray *c35_lhs7 = NULL;
  const mxArray *c35_rhs8 = NULL;
  const mxArray *c35_lhs8 = NULL;
  const mxArray *c35_rhs9 = NULL;
  const mxArray *c35_lhs9 = NULL;
  const mxArray *c35_rhs10 = NULL;
  const mxArray *c35_lhs10 = NULL;
  const mxArray *c35_rhs11 = NULL;
  const mxArray *c35_lhs11 = NULL;
  const mxArray *c35_rhs12 = NULL;
  const mxArray *c35_lhs12 = NULL;
  const mxArray *c35_rhs13 = NULL;
  const mxArray *c35_lhs13 = NULL;
  const mxArray *c35_rhs14 = NULL;
  const mxArray *c35_lhs14 = NULL;
  const mxArray *c35_rhs15 = NULL;
  const mxArray *c35_lhs15 = NULL;
  const mxArray *c35_rhs16 = NULL;
  const mxArray *c35_lhs16 = NULL;
  const mxArray *c35_rhs17 = NULL;
  const mxArray *c35_lhs17 = NULL;
  const mxArray *c35_rhs18 = NULL;
  const mxArray *c35_lhs18 = NULL;
  const mxArray *c35_rhs19 = NULL;
  const mxArray *c35_lhs19 = NULL;
  const mxArray *c35_rhs20 = NULL;
  const mxArray *c35_lhs20 = NULL;
  const mxArray *c35_rhs21 = NULL;
  const mxArray *c35_lhs21 = NULL;
  const mxArray *c35_rhs22 = NULL;
  const mxArray *c35_lhs22 = NULL;
  const mxArray *c35_rhs23 = NULL;
  const mxArray *c35_lhs23 = NULL;
  const mxArray *c35_rhs24 = NULL;
  const mxArray *c35_lhs24 = NULL;
  const mxArray *c35_rhs25 = NULL;
  const mxArray *c35_lhs25 = NULL;
  const mxArray *c35_rhs26 = NULL;
  const mxArray *c35_lhs26 = NULL;
  const mxArray *c35_rhs27 = NULL;
  const mxArray *c35_lhs27 = NULL;
  const mxArray *c35_rhs28 = NULL;
  const mxArray *c35_lhs28 = NULL;
  const mxArray *c35_rhs29 = NULL;
  const mxArray *c35_lhs29 = NULL;
  const mxArray *c35_rhs30 = NULL;
  const mxArray *c35_lhs30 = NULL;
  const mxArray *c35_rhs31 = NULL;
  const mxArray *c35_lhs31 = NULL;
  const mxArray *c35_rhs32 = NULL;
  const mxArray *c35_lhs32 = NULL;
  const mxArray *c35_rhs33 = NULL;
  const mxArray *c35_lhs33 = NULL;
  const mxArray *c35_rhs34 = NULL;
  const mxArray *c35_lhs34 = NULL;
  const mxArray *c35_rhs35 = NULL;
  const mxArray *c35_lhs35 = NULL;
  const mxArray *c35_rhs36 = NULL;
  const mxArray *c35_lhs36 = NULL;
  const mxArray *c35_rhs37 = NULL;
  const mxArray *c35_lhs37 = NULL;
  const mxArray *c35_rhs38 = NULL;
  const mxArray *c35_lhs38 = NULL;
  const mxArray *c35_rhs39 = NULL;
  const mxArray *c35_lhs39 = NULL;
  const mxArray *c35_rhs40 = NULL;
  const mxArray *c35_lhs40 = NULL;
  const mxArray *c35_rhs41 = NULL;
  const mxArray *c35_lhs41 = NULL;
  const mxArray *c35_rhs42 = NULL;
  const mxArray *c35_lhs42 = NULL;
  const mxArray *c35_rhs43 = NULL;
  const mxArray *c35_lhs43 = NULL;
  const mxArray *c35_rhs44 = NULL;
  const mxArray *c35_lhs44 = NULL;
  const mxArray *c35_rhs45 = NULL;
  const mxArray *c35_lhs45 = NULL;
  const mxArray *c35_rhs46 = NULL;
  const mxArray *c35_lhs46 = NULL;
  const mxArray *c35_rhs47 = NULL;
  const mxArray *c35_lhs47 = NULL;
  const mxArray *c35_rhs48 = NULL;
  const mxArray *c35_lhs48 = NULL;
  const mxArray *c35_rhs49 = NULL;
  const mxArray *c35_lhs49 = NULL;
  const mxArray *c35_rhs50 = NULL;
  const mxArray *c35_lhs50 = NULL;
  const mxArray *c35_rhs51 = NULL;
  const mxArray *c35_lhs51 = NULL;
  const mxArray *c35_rhs52 = NULL;
  const mxArray *c35_lhs52 = NULL;
  const mxArray *c35_rhs53 = NULL;
  const mxArray *c35_lhs53 = NULL;
  const mxArray *c35_rhs54 = NULL;
  const mxArray *c35_lhs54 = NULL;
  const mxArray *c35_rhs55 = NULL;
  const mxArray *c35_lhs55 = NULL;
  const mxArray *c35_rhs56 = NULL;
  const mxArray *c35_lhs56 = NULL;
  const mxArray *c35_rhs57 = NULL;
  const mxArray *c35_lhs57 = NULL;
  const mxArray *c35_rhs58 = NULL;
  const mxArray *c35_lhs58 = NULL;
  const mxArray *c35_rhs59 = NULL;
  const mxArray *c35_lhs59 = NULL;
  const mxArray *c35_rhs60 = NULL;
  const mxArray *c35_lhs60 = NULL;
  const mxArray *c35_rhs61 = NULL;
  const mxArray *c35_lhs61 = NULL;
  const mxArray *c35_rhs62 = NULL;
  const mxArray *c35_lhs62 = NULL;
  const mxArray *c35_rhs63 = NULL;
  const mxArray *c35_lhs63 = NULL;
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c35_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c35_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("reshape"), "name", "name", 2);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822368U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c35_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 3);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c35_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size"),
                  "context", "context", 4);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 4);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c35_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!varargin_nempty"),
                  "context", "context", 5);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 5);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c35_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size"),
                  "context", "context", 6);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 6);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1368186630U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c35_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 7);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 7);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c35_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 8);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("isinf"), "name", "name", 8);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717456U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c35_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 9);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c35_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 10);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 10);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822382U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c35_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 11);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name", 11);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c35_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 12);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmin"), "name", "name", 12);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c35_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 13);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 13);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c35_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 14);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 14);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c35_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 15);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmin"), "name", "name", 15);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c35_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!numel_for_size"),
                  "context", "context", 16);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name", 16);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c35_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 17);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c35_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 18);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name", 18);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c35_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size"),
                  "context", "context", 19);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_prod"), "name",
                  "name", 19);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 19);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c35_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 20);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c35_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 21);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c35_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 22);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name", 22);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c35_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 23);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 23);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c35_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 24);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c35_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 25);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c35_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 26);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c35_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 27);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mpower"), "name", "name", 27);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c35_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 28);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c35_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("ismatrix"), "name", "name",
                  29);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c35_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("power"), "name", "name", 30);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c35_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 31);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c35_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 32);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 32);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c35_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 33);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 33);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c35_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 34);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("floor"), "name", "name", 34);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 34);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c35_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 35);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c35_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 36);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822326U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c35_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 37);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 37);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c35_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 38);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name", 38);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c35_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 39);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("sqrt"), "name", "name", 39);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c35_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_error"), "name", "name",
                  40);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1343833958U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c35_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 41);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822338U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c35_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 42);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mrdivide"), "name", "name",
                  42);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c35_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 43);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("rdivide"), "name", "name",
                  43);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 43);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c35_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 44);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c35_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 45);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c35_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_div"), "name", "name",
                  46);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c35_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 47);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c35_rhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 48);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c35_rhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 49);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  49);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c35_rhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 50);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c35_rhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 51);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 51);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c35_rhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 52);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 52);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c35_rhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 53);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c35_rhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 54);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c35_rhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 55);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 55);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c35_rhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 56);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c35_rhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 57);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c35_rhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 58);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 58);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 58);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c35_rhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 59);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 59);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 59);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c35_rhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 60);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 60);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c35_rhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 61);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c35_rhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 62);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c35_rhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xdotu"), "name", "name",
                  63);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c35_rhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c35_rhs0);
  sf_mex_destroy(&c35_lhs0);
  sf_mex_destroy(&c35_rhs1);
  sf_mex_destroy(&c35_lhs1);
  sf_mex_destroy(&c35_rhs2);
  sf_mex_destroy(&c35_lhs2);
  sf_mex_destroy(&c35_rhs3);
  sf_mex_destroy(&c35_lhs3);
  sf_mex_destroy(&c35_rhs4);
  sf_mex_destroy(&c35_lhs4);
  sf_mex_destroy(&c35_rhs5);
  sf_mex_destroy(&c35_lhs5);
  sf_mex_destroy(&c35_rhs6);
  sf_mex_destroy(&c35_lhs6);
  sf_mex_destroy(&c35_rhs7);
  sf_mex_destroy(&c35_lhs7);
  sf_mex_destroy(&c35_rhs8);
  sf_mex_destroy(&c35_lhs8);
  sf_mex_destroy(&c35_rhs9);
  sf_mex_destroy(&c35_lhs9);
  sf_mex_destroy(&c35_rhs10);
  sf_mex_destroy(&c35_lhs10);
  sf_mex_destroy(&c35_rhs11);
  sf_mex_destroy(&c35_lhs11);
  sf_mex_destroy(&c35_rhs12);
  sf_mex_destroy(&c35_lhs12);
  sf_mex_destroy(&c35_rhs13);
  sf_mex_destroy(&c35_lhs13);
  sf_mex_destroy(&c35_rhs14);
  sf_mex_destroy(&c35_lhs14);
  sf_mex_destroy(&c35_rhs15);
  sf_mex_destroy(&c35_lhs15);
  sf_mex_destroy(&c35_rhs16);
  sf_mex_destroy(&c35_lhs16);
  sf_mex_destroy(&c35_rhs17);
  sf_mex_destroy(&c35_lhs17);
  sf_mex_destroy(&c35_rhs18);
  sf_mex_destroy(&c35_lhs18);
  sf_mex_destroy(&c35_rhs19);
  sf_mex_destroy(&c35_lhs19);
  sf_mex_destroy(&c35_rhs20);
  sf_mex_destroy(&c35_lhs20);
  sf_mex_destroy(&c35_rhs21);
  sf_mex_destroy(&c35_lhs21);
  sf_mex_destroy(&c35_rhs22);
  sf_mex_destroy(&c35_lhs22);
  sf_mex_destroy(&c35_rhs23);
  sf_mex_destroy(&c35_lhs23);
  sf_mex_destroy(&c35_rhs24);
  sf_mex_destroy(&c35_lhs24);
  sf_mex_destroy(&c35_rhs25);
  sf_mex_destroy(&c35_lhs25);
  sf_mex_destroy(&c35_rhs26);
  sf_mex_destroy(&c35_lhs26);
  sf_mex_destroy(&c35_rhs27);
  sf_mex_destroy(&c35_lhs27);
  sf_mex_destroy(&c35_rhs28);
  sf_mex_destroy(&c35_lhs28);
  sf_mex_destroy(&c35_rhs29);
  sf_mex_destroy(&c35_lhs29);
  sf_mex_destroy(&c35_rhs30);
  sf_mex_destroy(&c35_lhs30);
  sf_mex_destroy(&c35_rhs31);
  sf_mex_destroy(&c35_lhs31);
  sf_mex_destroy(&c35_rhs32);
  sf_mex_destroy(&c35_lhs32);
  sf_mex_destroy(&c35_rhs33);
  sf_mex_destroy(&c35_lhs33);
  sf_mex_destroy(&c35_rhs34);
  sf_mex_destroy(&c35_lhs34);
  sf_mex_destroy(&c35_rhs35);
  sf_mex_destroy(&c35_lhs35);
  sf_mex_destroy(&c35_rhs36);
  sf_mex_destroy(&c35_lhs36);
  sf_mex_destroy(&c35_rhs37);
  sf_mex_destroy(&c35_lhs37);
  sf_mex_destroy(&c35_rhs38);
  sf_mex_destroy(&c35_lhs38);
  sf_mex_destroy(&c35_rhs39);
  sf_mex_destroy(&c35_lhs39);
  sf_mex_destroy(&c35_rhs40);
  sf_mex_destroy(&c35_lhs40);
  sf_mex_destroy(&c35_rhs41);
  sf_mex_destroy(&c35_lhs41);
  sf_mex_destroy(&c35_rhs42);
  sf_mex_destroy(&c35_lhs42);
  sf_mex_destroy(&c35_rhs43);
  sf_mex_destroy(&c35_lhs43);
  sf_mex_destroy(&c35_rhs44);
  sf_mex_destroy(&c35_lhs44);
  sf_mex_destroy(&c35_rhs45);
  sf_mex_destroy(&c35_lhs45);
  sf_mex_destroy(&c35_rhs46);
  sf_mex_destroy(&c35_lhs46);
  sf_mex_destroy(&c35_rhs47);
  sf_mex_destroy(&c35_lhs47);
  sf_mex_destroy(&c35_rhs48);
  sf_mex_destroy(&c35_lhs48);
  sf_mex_destroy(&c35_rhs49);
  sf_mex_destroy(&c35_lhs49);
  sf_mex_destroy(&c35_rhs50);
  sf_mex_destroy(&c35_lhs50);
  sf_mex_destroy(&c35_rhs51);
  sf_mex_destroy(&c35_lhs51);
  sf_mex_destroy(&c35_rhs52);
  sf_mex_destroy(&c35_lhs52);
  sf_mex_destroy(&c35_rhs53);
  sf_mex_destroy(&c35_lhs53);
  sf_mex_destroy(&c35_rhs54);
  sf_mex_destroy(&c35_lhs54);
  sf_mex_destroy(&c35_rhs55);
  sf_mex_destroy(&c35_lhs55);
  sf_mex_destroy(&c35_rhs56);
  sf_mex_destroy(&c35_lhs56);
  sf_mex_destroy(&c35_rhs57);
  sf_mex_destroy(&c35_lhs57);
  sf_mex_destroy(&c35_rhs58);
  sf_mex_destroy(&c35_lhs58);
  sf_mex_destroy(&c35_rhs59);
  sf_mex_destroy(&c35_lhs59);
  sf_mex_destroy(&c35_rhs60);
  sf_mex_destroy(&c35_lhs60);
  sf_mex_destroy(&c35_rhs61);
  sf_mex_destroy(&c35_lhs61);
  sf_mex_destroy(&c35_rhs62);
  sf_mex_destroy(&c35_lhs62);
  sf_mex_destroy(&c35_rhs63);
  sf_mex_destroy(&c35_lhs63);
}

static const mxArray *c35_emlrt_marshallOut(char * c35_u)
{
  const mxArray *c35_y = NULL;
  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c35_u)), FALSE);
  return c35_y;
}

static const mxArray *c35_b_emlrt_marshallOut(uint32_T c35_u)
{
  const mxArray *c35_y = NULL;
  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", &c35_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c35_y;
}

static void c35_b_info_helper(const mxArray **c35_info)
{
  const mxArray *c35_rhs64 = NULL;
  const mxArray *c35_lhs64 = NULL;
  const mxArray *c35_rhs65 = NULL;
  const mxArray *c35_lhs65 = NULL;
  const mxArray *c35_rhs66 = NULL;
  const mxArray *c35_lhs66 = NULL;
  const mxArray *c35_rhs67 = NULL;
  const mxArray *c35_lhs67 = NULL;
  const mxArray *c35_rhs68 = NULL;
  const mxArray *c35_lhs68 = NULL;
  const mxArray *c35_rhs69 = NULL;
  const mxArray *c35_lhs69 = NULL;
  const mxArray *c35_rhs70 = NULL;
  const mxArray *c35_lhs70 = NULL;
  const mxArray *c35_rhs71 = NULL;
  const mxArray *c35_lhs71 = NULL;
  const mxArray *c35_rhs72 = NULL;
  const mxArray *c35_lhs72 = NULL;
  const mxArray *c35_rhs73 = NULL;
  const mxArray *c35_lhs73 = NULL;
  const mxArray *c35_rhs74 = NULL;
  const mxArray *c35_lhs74 = NULL;
  const mxArray *c35_rhs75 = NULL;
  const mxArray *c35_lhs75 = NULL;
  const mxArray *c35_rhs76 = NULL;
  const mxArray *c35_lhs76 = NULL;
  const mxArray *c35_rhs77 = NULL;
  const mxArray *c35_lhs77 = NULL;
  const mxArray *c35_rhs78 = NULL;
  const mxArray *c35_lhs78 = NULL;
  const mxArray *c35_rhs79 = NULL;
  const mxArray *c35_lhs79 = NULL;
  const mxArray *c35_rhs80 = NULL;
  const mxArray *c35_lhs80 = NULL;
  const mxArray *c35_rhs81 = NULL;
  const mxArray *c35_lhs81 = NULL;
  const mxArray *c35_rhs82 = NULL;
  const mxArray *c35_lhs82 = NULL;
  const mxArray *c35_rhs83 = NULL;
  const mxArray *c35_lhs83 = NULL;
  const mxArray *c35_rhs84 = NULL;
  const mxArray *c35_lhs84 = NULL;
  const mxArray *c35_rhs85 = NULL;
  const mxArray *c35_lhs85 = NULL;
  const mxArray *c35_rhs86 = NULL;
  const mxArray *c35_lhs86 = NULL;
  const mxArray *c35_rhs87 = NULL;
  const mxArray *c35_lhs87 = NULL;
  const mxArray *c35_rhs88 = NULL;
  const mxArray *c35_lhs88 = NULL;
  const mxArray *c35_rhs89 = NULL;
  const mxArray *c35_lhs89 = NULL;
  const mxArray *c35_rhs90 = NULL;
  const mxArray *c35_lhs90 = NULL;
  const mxArray *c35_rhs91 = NULL;
  const mxArray *c35_lhs91 = NULL;
  const mxArray *c35_rhs92 = NULL;
  const mxArray *c35_lhs92 = NULL;
  const mxArray *c35_rhs93 = NULL;
  const mxArray *c35_lhs93 = NULL;
  const mxArray *c35_rhs94 = NULL;
  const mxArray *c35_lhs94 = NULL;
  const mxArray *c35_rhs95 = NULL;
  const mxArray *c35_lhs95 = NULL;
  const mxArray *c35_rhs96 = NULL;
  const mxArray *c35_lhs96 = NULL;
  const mxArray *c35_rhs97 = NULL;
  const mxArray *c35_lhs97 = NULL;
  const mxArray *c35_rhs98 = NULL;
  const mxArray *c35_lhs98 = NULL;
  const mxArray *c35_rhs99 = NULL;
  const mxArray *c35_lhs99 = NULL;
  const mxArray *c35_rhs100 = NULL;
  const mxArray *c35_lhs100 = NULL;
  const mxArray *c35_rhs101 = NULL;
  const mxArray *c35_lhs101 = NULL;
  const mxArray *c35_rhs102 = NULL;
  const mxArray *c35_lhs102 = NULL;
  const mxArray *c35_rhs103 = NULL;
  const mxArray *c35_lhs103 = NULL;
  const mxArray *c35_rhs104 = NULL;
  const mxArray *c35_lhs104 = NULL;
  const mxArray *c35_rhs105 = NULL;
  const mxArray *c35_lhs105 = NULL;
  const mxArray *c35_rhs106 = NULL;
  const mxArray *c35_lhs106 = NULL;
  const mxArray *c35_rhs107 = NULL;
  const mxArray *c35_lhs107 = NULL;
  const mxArray *c35_rhs108 = NULL;
  const mxArray *c35_lhs108 = NULL;
  const mxArray *c35_rhs109 = NULL;
  const mxArray *c35_lhs109 = NULL;
  const mxArray *c35_rhs110 = NULL;
  const mxArray *c35_lhs110 = NULL;
  const mxArray *c35_rhs111 = NULL;
  const mxArray *c35_lhs111 = NULL;
  const mxArray *c35_rhs112 = NULL;
  const mxArray *c35_lhs112 = NULL;
  const mxArray *c35_rhs113 = NULL;
  const mxArray *c35_lhs113 = NULL;
  const mxArray *c35_rhs114 = NULL;
  const mxArray *c35_lhs114 = NULL;
  const mxArray *c35_rhs115 = NULL;
  const mxArray *c35_lhs115 = NULL;
  const mxArray *c35_rhs116 = NULL;
  const mxArray *c35_lhs116 = NULL;
  const mxArray *c35_rhs117 = NULL;
  const mxArray *c35_lhs117 = NULL;
  const mxArray *c35_rhs118 = NULL;
  const mxArray *c35_lhs118 = NULL;
  const mxArray *c35_rhs119 = NULL;
  const mxArray *c35_lhs119 = NULL;
  const mxArray *c35_rhs120 = NULL;
  const mxArray *c35_lhs120 = NULL;
  const mxArray *c35_rhs121 = NULL;
  const mxArray *c35_lhs121 = NULL;
  const mxArray *c35_rhs122 = NULL;
  const mxArray *c35_lhs122 = NULL;
  const mxArray *c35_rhs123 = NULL;
  const mxArray *c35_lhs123 = NULL;
  const mxArray *c35_rhs124 = NULL;
  const mxArray *c35_lhs124 = NULL;
  const mxArray *c35_rhs125 = NULL;
  const mxArray *c35_lhs125 = NULL;
  const mxArray *c35_rhs126 = NULL;
  const mxArray *c35_lhs126 = NULL;
  const mxArray *c35_rhs127 = NULL;
  const mxArray *c35_lhs127 = NULL;
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 64);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 64);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c35_rhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 65);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xdot"), "name", "name",
                  65);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 65);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c35_rhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "context",
                  "context", 66);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 66);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c35_rhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 67);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c35_rhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xdot"), "name",
                  "name", 68);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080372U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c35_rhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xdotx"), "name",
                  "name", 69);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c35_rhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 70);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 70);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 70);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c35_rhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 71);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c35_rhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 72);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 72);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c35_rhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 73);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 73);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c35_rhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 74);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 74);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c35_rhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 75);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c35_rhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context", 76);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mldivide"), "name", "name",
                  76);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 76);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c35_rhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 77);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_lusolve"), "name",
                  "name", 77);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 77);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1309454796U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c35_rhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "context",
                  "context", 78);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 78);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c35_rhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 79);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 79);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c35_rhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 80);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  80);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822406U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c35_rhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 81);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822410U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c35_rhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 82);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 82);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1302692594U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c35_rhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("realmin"), "name", "name",
                  83);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 83);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1307654842U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c35_rhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 84);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_realmin"), "name",
                  "name", 84);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1307654844U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c35_rhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 85);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 85);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c35_rhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 86);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eps"), "name", "name", 86);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 86);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c35_rhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 87);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 87);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822382U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c35_rhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 88);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_eps"), "name", "name",
                  88);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 88);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c35_rhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 89);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 89);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731596U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c35_rhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 90);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("min"), "name", "name", 90);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 90);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 90);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c35_rhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 91);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 91);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c35_rhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 92);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 92);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 92);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c35_rhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 93);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 93);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 93);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c35_rhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 94);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 94);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c35_rhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 95);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 95);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 95);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 95);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c35_rhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 96);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 96);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 96);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c35_rhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("colon"), "name", "name", 97);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 97);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1366165842U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c35_rhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 98);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("colon"), "name", "name", 98);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 98);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1366165842U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c35_rhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 99);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 99);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c35_rhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 100);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 100);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 100);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c35_rhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs100), "rhs",
                  "rhs", 100);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs100), "lhs",
                  "lhs", 100);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 101);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("floor"), "name", "name", 101);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 101);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c35_rhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs101), "rhs",
                  "rhs", 101);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs101), "lhs",
                  "lhs", 101);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 102);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmin"), "name", "name",
                  102);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 102);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c35_rhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs102), "rhs",
                  "rhs", 102);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs102), "lhs",
                  "lhs", 102);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 103);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name",
                  103);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 103);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c35_rhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs103), "rhs",
                  "rhs", 103);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs103), "lhs",
                  "lhs", 103);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 104);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmin"), "name", "name",
                  104);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 104);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c35_rhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs104), "rhs",
                  "rhs", 104);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs104), "lhs",
                  "lhs", 104);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 105);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name",
                  105);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 105);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c35_rhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs105), "rhs",
                  "rhs", 105);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs105), "lhs",
                  "lhs", 105);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 106);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_isa_uint"), "name",
                  "name", 106);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 106);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 106);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822384U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c35_rhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs106), "rhs",
                  "rhs", 106);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs106), "lhs",
                  "lhs", 106);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 107);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 107);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174180U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c35_rhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs107), "rhs",
                  "rhs", 107);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs107), "lhs",
                  "lhs", 107);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 108);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c35_rhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs108), "rhs",
                  "rhs", 108);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs108), "lhs",
                  "lhs", 108);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 109);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 109);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c35_rhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs109), "rhs",
                  "rhs", 109);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs109), "lhs",
                  "lhs", 109);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 110);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name",
                  110);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 110);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c35_rhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs110), "rhs",
                  "rhs", 110);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs110), "lhs",
                  "lhs", 110);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 111);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_isa_uint"), "name",
                  "name", 111);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 111);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 111);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822384U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c35_rhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs111), "rhs",
                  "rhs", 111);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs111), "lhs",
                  "lhs", 111);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 112);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 112);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 112);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c35_rhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs112), "rhs",
                  "rhs", 112);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs112), "lhs",
                  "lhs", 112);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 113);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 113);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c35_rhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs113), "rhs",
                  "rhs", 113);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs113), "lhs",
                  "lhs", 113);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 114);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 114);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c35_rhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs114), "rhs",
                  "rhs", 114);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs114), "lhs",
                  "lhs", 114);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 115);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 115);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c35_rhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs115), "rhs",
                  "rhs", 115);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs115), "lhs",
                  "lhs", 115);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 116);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 116);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c35_rhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs116), "rhs",
                  "rhs", 116);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs116), "lhs",
                  "lhs", 116);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 117);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 117);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c35_rhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs117), "rhs",
                  "rhs", 117);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs117), "lhs",
                  "lhs", 117);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 118);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 118);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 118);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c35_rhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs118), "rhs",
                  "rhs", 118);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs118), "lhs",
                  "lhs", 118);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 119);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 119);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 119);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c35_rhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs119), "rhs",
                  "rhs", 119);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs119), "lhs",
                  "lhs", 119);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 120);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 120);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 120);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c35_rhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs120), "rhs",
                  "rhs", 120);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs120), "lhs",
                  "lhs", 120);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 121);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  121);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 121);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c35_rhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs121), "rhs",
                  "rhs", 121);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs121), "lhs",
                  "lhs", 121);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 122);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 122);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c35_rhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs122), "rhs",
                  "rhs", 122);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs122), "lhs",
                  "lhs", 122);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m!below_threshold"),
                  "context", "context", 123);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("length"), "name", "name",
                  123);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1303149806U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c35_rhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs123), "rhs",
                  "rhs", 123);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs123), "lhs",
                  "lhs", 123);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 124);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 124);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c35_rhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs124), "rhs",
                  "rhs", 124);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs124), "lhs",
                  "lhs", 124);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"),
                  "context", "context", 125);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 125);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c35_rhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs125), "rhs",
                  "rhs", 125);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs125), "lhs",
                  "lhs", 125);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"),
                  "context", "context", 126);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_ixamax"), "name",
                  "name", 126);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080370U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c35_rhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs126), "rhs",
                  "rhs", 126);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs126), "lhs",
                  "lhs", 126);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 127);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 127);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c35_rhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs127), "rhs",
                  "rhs", 127);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs127), "lhs",
                  "lhs", 127);
  sf_mex_destroy(&c35_rhs64);
  sf_mex_destroy(&c35_lhs64);
  sf_mex_destroy(&c35_rhs65);
  sf_mex_destroy(&c35_lhs65);
  sf_mex_destroy(&c35_rhs66);
  sf_mex_destroy(&c35_lhs66);
  sf_mex_destroy(&c35_rhs67);
  sf_mex_destroy(&c35_lhs67);
  sf_mex_destroy(&c35_rhs68);
  sf_mex_destroy(&c35_lhs68);
  sf_mex_destroy(&c35_rhs69);
  sf_mex_destroy(&c35_lhs69);
  sf_mex_destroy(&c35_rhs70);
  sf_mex_destroy(&c35_lhs70);
  sf_mex_destroy(&c35_rhs71);
  sf_mex_destroy(&c35_lhs71);
  sf_mex_destroy(&c35_rhs72);
  sf_mex_destroy(&c35_lhs72);
  sf_mex_destroy(&c35_rhs73);
  sf_mex_destroy(&c35_lhs73);
  sf_mex_destroy(&c35_rhs74);
  sf_mex_destroy(&c35_lhs74);
  sf_mex_destroy(&c35_rhs75);
  sf_mex_destroy(&c35_lhs75);
  sf_mex_destroy(&c35_rhs76);
  sf_mex_destroy(&c35_lhs76);
  sf_mex_destroy(&c35_rhs77);
  sf_mex_destroy(&c35_lhs77);
  sf_mex_destroy(&c35_rhs78);
  sf_mex_destroy(&c35_lhs78);
  sf_mex_destroy(&c35_rhs79);
  sf_mex_destroy(&c35_lhs79);
  sf_mex_destroy(&c35_rhs80);
  sf_mex_destroy(&c35_lhs80);
  sf_mex_destroy(&c35_rhs81);
  sf_mex_destroy(&c35_lhs81);
  sf_mex_destroy(&c35_rhs82);
  sf_mex_destroy(&c35_lhs82);
  sf_mex_destroy(&c35_rhs83);
  sf_mex_destroy(&c35_lhs83);
  sf_mex_destroy(&c35_rhs84);
  sf_mex_destroy(&c35_lhs84);
  sf_mex_destroy(&c35_rhs85);
  sf_mex_destroy(&c35_lhs85);
  sf_mex_destroy(&c35_rhs86);
  sf_mex_destroy(&c35_lhs86);
  sf_mex_destroy(&c35_rhs87);
  sf_mex_destroy(&c35_lhs87);
  sf_mex_destroy(&c35_rhs88);
  sf_mex_destroy(&c35_lhs88);
  sf_mex_destroy(&c35_rhs89);
  sf_mex_destroy(&c35_lhs89);
  sf_mex_destroy(&c35_rhs90);
  sf_mex_destroy(&c35_lhs90);
  sf_mex_destroy(&c35_rhs91);
  sf_mex_destroy(&c35_lhs91);
  sf_mex_destroy(&c35_rhs92);
  sf_mex_destroy(&c35_lhs92);
  sf_mex_destroy(&c35_rhs93);
  sf_mex_destroy(&c35_lhs93);
  sf_mex_destroy(&c35_rhs94);
  sf_mex_destroy(&c35_lhs94);
  sf_mex_destroy(&c35_rhs95);
  sf_mex_destroy(&c35_lhs95);
  sf_mex_destroy(&c35_rhs96);
  sf_mex_destroy(&c35_lhs96);
  sf_mex_destroy(&c35_rhs97);
  sf_mex_destroy(&c35_lhs97);
  sf_mex_destroy(&c35_rhs98);
  sf_mex_destroy(&c35_lhs98);
  sf_mex_destroy(&c35_rhs99);
  sf_mex_destroy(&c35_lhs99);
  sf_mex_destroy(&c35_rhs100);
  sf_mex_destroy(&c35_lhs100);
  sf_mex_destroy(&c35_rhs101);
  sf_mex_destroy(&c35_lhs101);
  sf_mex_destroy(&c35_rhs102);
  sf_mex_destroy(&c35_lhs102);
  sf_mex_destroy(&c35_rhs103);
  sf_mex_destroy(&c35_lhs103);
  sf_mex_destroy(&c35_rhs104);
  sf_mex_destroy(&c35_lhs104);
  sf_mex_destroy(&c35_rhs105);
  sf_mex_destroy(&c35_lhs105);
  sf_mex_destroy(&c35_rhs106);
  sf_mex_destroy(&c35_lhs106);
  sf_mex_destroy(&c35_rhs107);
  sf_mex_destroy(&c35_lhs107);
  sf_mex_destroy(&c35_rhs108);
  sf_mex_destroy(&c35_lhs108);
  sf_mex_destroy(&c35_rhs109);
  sf_mex_destroy(&c35_lhs109);
  sf_mex_destroy(&c35_rhs110);
  sf_mex_destroy(&c35_lhs110);
  sf_mex_destroy(&c35_rhs111);
  sf_mex_destroy(&c35_lhs111);
  sf_mex_destroy(&c35_rhs112);
  sf_mex_destroy(&c35_lhs112);
  sf_mex_destroy(&c35_rhs113);
  sf_mex_destroy(&c35_lhs113);
  sf_mex_destroy(&c35_rhs114);
  sf_mex_destroy(&c35_lhs114);
  sf_mex_destroy(&c35_rhs115);
  sf_mex_destroy(&c35_lhs115);
  sf_mex_destroy(&c35_rhs116);
  sf_mex_destroy(&c35_lhs116);
  sf_mex_destroy(&c35_rhs117);
  sf_mex_destroy(&c35_lhs117);
  sf_mex_destroy(&c35_rhs118);
  sf_mex_destroy(&c35_lhs118);
  sf_mex_destroy(&c35_rhs119);
  sf_mex_destroy(&c35_lhs119);
  sf_mex_destroy(&c35_rhs120);
  sf_mex_destroy(&c35_lhs120);
  sf_mex_destroy(&c35_rhs121);
  sf_mex_destroy(&c35_lhs121);
  sf_mex_destroy(&c35_rhs122);
  sf_mex_destroy(&c35_lhs122);
  sf_mex_destroy(&c35_rhs123);
  sf_mex_destroy(&c35_lhs123);
  sf_mex_destroy(&c35_rhs124);
  sf_mex_destroy(&c35_lhs124);
  sf_mex_destroy(&c35_rhs125);
  sf_mex_destroy(&c35_lhs125);
  sf_mex_destroy(&c35_rhs126);
  sf_mex_destroy(&c35_lhs126);
  sf_mex_destroy(&c35_rhs127);
  sf_mex_destroy(&c35_lhs127);
}

static void c35_c_info_helper(const mxArray **c35_info)
{
  const mxArray *c35_rhs128 = NULL;
  const mxArray *c35_lhs128 = NULL;
  const mxArray *c35_rhs129 = NULL;
  const mxArray *c35_lhs129 = NULL;
  const mxArray *c35_rhs130 = NULL;
  const mxArray *c35_lhs130 = NULL;
  const mxArray *c35_rhs131 = NULL;
  const mxArray *c35_lhs131 = NULL;
  const mxArray *c35_rhs132 = NULL;
  const mxArray *c35_lhs132 = NULL;
  const mxArray *c35_rhs133 = NULL;
  const mxArray *c35_lhs133 = NULL;
  const mxArray *c35_rhs134 = NULL;
  const mxArray *c35_lhs134 = NULL;
  const mxArray *c35_rhs135 = NULL;
  const mxArray *c35_lhs135 = NULL;
  const mxArray *c35_rhs136 = NULL;
  const mxArray *c35_lhs136 = NULL;
  const mxArray *c35_rhs137 = NULL;
  const mxArray *c35_lhs137 = NULL;
  const mxArray *c35_rhs138 = NULL;
  const mxArray *c35_lhs138 = NULL;
  const mxArray *c35_rhs139 = NULL;
  const mxArray *c35_lhs139 = NULL;
  const mxArray *c35_rhs140 = NULL;
  const mxArray *c35_lhs140 = NULL;
  const mxArray *c35_rhs141 = NULL;
  const mxArray *c35_lhs141 = NULL;
  const mxArray *c35_rhs142 = NULL;
  const mxArray *c35_lhs142 = NULL;
  const mxArray *c35_rhs143 = NULL;
  const mxArray *c35_lhs143 = NULL;
  const mxArray *c35_rhs144 = NULL;
  const mxArray *c35_lhs144 = NULL;
  const mxArray *c35_rhs145 = NULL;
  const mxArray *c35_lhs145 = NULL;
  const mxArray *c35_rhs146 = NULL;
  const mxArray *c35_lhs146 = NULL;
  const mxArray *c35_rhs147 = NULL;
  const mxArray *c35_lhs147 = NULL;
  const mxArray *c35_rhs148 = NULL;
  const mxArray *c35_lhs148 = NULL;
  const mxArray *c35_rhs149 = NULL;
  const mxArray *c35_lhs149 = NULL;
  const mxArray *c35_rhs150 = NULL;
  const mxArray *c35_lhs150 = NULL;
  const mxArray *c35_rhs151 = NULL;
  const mxArray *c35_lhs151 = NULL;
  const mxArray *c35_rhs152 = NULL;
  const mxArray *c35_lhs152 = NULL;
  const mxArray *c35_rhs153 = NULL;
  const mxArray *c35_lhs153 = NULL;
  const mxArray *c35_rhs154 = NULL;
  const mxArray *c35_lhs154 = NULL;
  const mxArray *c35_rhs155 = NULL;
  const mxArray *c35_lhs155 = NULL;
  const mxArray *c35_rhs156 = NULL;
  const mxArray *c35_lhs156 = NULL;
  const mxArray *c35_rhs157 = NULL;
  const mxArray *c35_lhs157 = NULL;
  const mxArray *c35_rhs158 = NULL;
  const mxArray *c35_lhs158 = NULL;
  const mxArray *c35_rhs159 = NULL;
  const mxArray *c35_lhs159 = NULL;
  const mxArray *c35_rhs160 = NULL;
  const mxArray *c35_lhs160 = NULL;
  const mxArray *c35_rhs161 = NULL;
  const mxArray *c35_lhs161 = NULL;
  const mxArray *c35_rhs162 = NULL;
  const mxArray *c35_lhs162 = NULL;
  const mxArray *c35_rhs163 = NULL;
  const mxArray *c35_lhs163 = NULL;
  const mxArray *c35_rhs164 = NULL;
  const mxArray *c35_lhs164 = NULL;
  const mxArray *c35_rhs165 = NULL;
  const mxArray *c35_lhs165 = NULL;
  const mxArray *c35_rhs166 = NULL;
  const mxArray *c35_lhs166 = NULL;
  const mxArray *c35_rhs167 = NULL;
  const mxArray *c35_lhs167 = NULL;
  const mxArray *c35_rhs168 = NULL;
  const mxArray *c35_lhs168 = NULL;
  const mxArray *c35_rhs169 = NULL;
  const mxArray *c35_lhs169 = NULL;
  const mxArray *c35_rhs170 = NULL;
  const mxArray *c35_lhs170 = NULL;
  const mxArray *c35_rhs171 = NULL;
  const mxArray *c35_lhs171 = NULL;
  const mxArray *c35_rhs172 = NULL;
  const mxArray *c35_lhs172 = NULL;
  const mxArray *c35_rhs173 = NULL;
  const mxArray *c35_lhs173 = NULL;
  const mxArray *c35_rhs174 = NULL;
  const mxArray *c35_lhs174 = NULL;
  const mxArray *c35_rhs175 = NULL;
  const mxArray *c35_lhs175 = NULL;
  const mxArray *c35_rhs176 = NULL;
  const mxArray *c35_lhs176 = NULL;
  const mxArray *c35_rhs177 = NULL;
  const mxArray *c35_lhs177 = NULL;
  const mxArray *c35_rhs178 = NULL;
  const mxArray *c35_lhs178 = NULL;
  const mxArray *c35_rhs179 = NULL;
  const mxArray *c35_lhs179 = NULL;
  const mxArray *c35_rhs180 = NULL;
  const mxArray *c35_lhs180 = NULL;
  const mxArray *c35_rhs181 = NULL;
  const mxArray *c35_lhs181 = NULL;
  const mxArray *c35_rhs182 = NULL;
  const mxArray *c35_lhs182 = NULL;
  const mxArray *c35_rhs183 = NULL;
  const mxArray *c35_lhs183 = NULL;
  const mxArray *c35_rhs184 = NULL;
  const mxArray *c35_lhs184 = NULL;
  const mxArray *c35_rhs185 = NULL;
  const mxArray *c35_lhs185 = NULL;
  const mxArray *c35_rhs186 = NULL;
  const mxArray *c35_lhs186 = NULL;
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 128);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  128);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822306U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c35_rhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs128), "rhs",
                  "rhs", 128);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs128), "lhs",
                  "lhs", 128);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 129);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("abs"), "name", "name", 129);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 129);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c35_rhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs129), "rhs",
                  "rhs", 129);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs129), "lhs",
                  "lhs", 129);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 130);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 130);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c35_rhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs130), "rhs",
                  "rhs", 130);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs130), "lhs",
                  "lhs", 130);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 131);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 131);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c35_rhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs131), "rhs",
                  "rhs", 131);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs131), "lhs",
                  "lhs", 131);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 132);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 132);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c35_rhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs132), "rhs",
                  "rhs", 132);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs132), "lhs",
                  "lhs", 132);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m"),
                  "context", "context", 133);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 133);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 133);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c35_rhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs133), "rhs",
                  "rhs", 133);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs133), "lhs",
                  "lhs", 133);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 134);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xswap"), "name", "name",
                  134);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717474U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c35_rhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs134), "rhs",
                  "rhs", 134);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs134), "lhs",
                  "lhs", 134);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 135);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c35_rhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs135), "rhs",
                  "rhs", 135);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs135), "lhs",
                  "lhs", 135);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m"),
                  "context", "context", 136);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 136);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c35_rhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs136), "rhs",
                  "rhs", 136);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs136), "lhs",
                  "lhs", 136);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m"),
                  "context", "context", 137);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xswap"), "name",
                  "name", 137);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080386U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c35_rhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs137), "rhs",
                  "rhs", 137);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs137), "lhs",
                  "lhs", 137);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 138);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 138);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 138);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c35_rhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs138), "rhs",
                  "rhs", 138);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs138), "lhs",
                  "lhs", 138);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 139);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("abs"), "name", "name", 139);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 139);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 139);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c35_rhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs139), "rhs",
                  "rhs", 139);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs139), "lhs",
                  "lhs", 139);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 140);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 140);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 140);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 140);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c35_rhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs140), "rhs",
                  "rhs", 140);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs140), "lhs",
                  "lhs", 140);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 141);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 141);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 141);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 141);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c35_rhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs141), "rhs",
                  "rhs", 141);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs141), "lhs",
                  "lhs", 141);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 142);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 142);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 142);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c35_rhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs142), "rhs",
                  "rhs", 142);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs142), "lhs",
                  "lhs", 142);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"),
                  "context", "context", 143);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 143);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 143);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c35_rhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs143), "rhs",
                  "rhs", 143);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs143), "lhs",
                  "lhs", 143);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 144);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_div"), "name", "name",
                  144);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c35_rhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs144), "rhs",
                  "rhs", 144);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs144), "lhs",
                  "lhs", 144);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 145);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  145);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 145);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717472U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c35_rhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs145), "rhs",
                  "rhs", 145);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs145), "lhs",
                  "lhs", 145);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 146);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 146);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c35_rhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs146), "rhs",
                  "rhs", 146);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs146), "lhs",
                  "lhs", 146);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 147);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xger"), "name", "name",
                  147);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c35_rhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs147), "rhs",
                  "rhs", 147);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs147), "lhs",
                  "lhs", 147);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m"), "context",
                  "context", 148);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 148);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c35_rhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs148), "rhs",
                  "rhs", 148);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs148), "lhs",
                  "lhs", 148);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 149);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmax"), "name", "name",
                  149);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 149);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c35_rhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs149), "rhs",
                  "rhs", 149);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs149), "lhs",
                  "lhs", 149);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 150);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("min"), "name", "name", 150);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1311258918U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c35_rhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs150), "rhs",
                  "rhs", 150);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs150), "lhs",
                  "lhs", 150);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 151);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 151);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 151);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c35_rhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs151), "rhs",
                  "rhs", 151);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs151), "lhs",
                  "lhs", 151);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 152);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 152);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c35_rhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs152), "rhs",
                  "rhs", 152);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs152), "lhs",
                  "lhs", 152);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 153);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 153);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c35_rhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs153), "rhs",
                  "rhs", 153);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs153), "lhs",
                  "lhs", 153);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 154);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 154);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c35_rhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs154), "rhs",
                  "rhs", 154);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs154), "lhs",
                  "lhs", 154);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold"),
                  "context", "context", 155);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name",
                  155);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c35_rhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs155), "rhs",
                  "rhs", 155);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs155), "lhs",
                  "lhs", 155);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m"),
                  "context", "context", 156);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 156);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c35_rhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs156), "rhs",
                  "rhs", 156);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs156), "lhs",
                  "lhs", 156);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xger"), "name",
                  "name", 157);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080376U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c35_rhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs157), "rhs",
                  "rhs", 157);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs157), "lhs",
                  "lhs", 157);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xgerx"), "name",
                  "name", 158);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c35_rhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs158), "rhs",
                  "rhs", 158);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs158), "lhs",
                  "lhs", 158);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 159);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 159);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c35_rhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs159), "rhs",
                  "rhs", 159);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs159), "lhs",
                  "lhs", 159);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 160);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("abs"), "name", "name", 160);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 160);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 160);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c35_rhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs160), "rhs",
                  "rhs", 160);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs160), "lhs",
                  "lhs", 160);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 161);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 161);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 161);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c35_rhs161, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs161, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs161), "rhs",
                  "rhs", 161);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs161), "lhs",
                  "lhs", 161);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 162);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 162);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c35_rhs162, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs162, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs162), "rhs",
                  "rhs", 162);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs162), "lhs",
                  "lhs", 162);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 163);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 163);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c35_rhs163, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs163, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs163), "rhs",
                  "rhs", 163);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs163), "lhs",
                  "lhs", 163);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m"),
                  "context", "context", 164);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 164);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 164);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 164);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c35_rhs164, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs164, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs164), "rhs",
                  "rhs", 164);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs164), "lhs",
                  "lhs", 164);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 165);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_warning"), "name",
                  "name", 165);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 165);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822402U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c35_rhs165, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs165, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs165), "rhs",
                  "rhs", 165);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs165), "lhs",
                  "lhs", 165);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 166);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 166);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 166);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c35_rhs166, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs166, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs166), "rhs",
                  "rhs", 166);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs166), "lhs",
                  "lhs", 166);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 167);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 167);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 167);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c35_rhs167, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs167, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs167), "rhs",
                  "rhs", 167);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs167), "lhs",
                  "lhs", 167);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN"),
                  "context", "context", 168);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  168);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 168);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080378U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c35_rhs168, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs168, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs168), "rhs",
                  "rhs", 168);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs168), "lhs",
                  "lhs", 168);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 169);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 169);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c35_rhs169, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs169, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs169), "rhs",
                  "rhs", 169);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs169), "lhs",
                  "lhs", 169);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m!below_threshold"),
                  "context", "context", 170);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name",
                  170);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 170);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c35_rhs170, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs170, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs170), "rhs",
                  "rhs", 170);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs170), "lhs",
                  "lhs", 170);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 171);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 171);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c35_rhs171, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs171, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs171), "rhs",
                  "rhs", 171);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs171), "lhs",
                  "lhs", 171);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 172);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 172);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 172);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c35_rhs172, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs172, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs172), "rhs",
                  "rhs", 172);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs172), "lhs",
                  "lhs", 172);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"),
                  "context", "context", 173);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_refblas_xtrsm"), "name",
                  "name", 173);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1360285952U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c35_rhs173, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs173, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs173), "rhs",
                  "rhs", 173);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs173), "lhs",
                  "lhs", 173);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 174);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 174);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 174);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c35_rhs174, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs174, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs174), "rhs",
                  "rhs", 174);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs174), "lhs",
                  "lhs", 174);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 175);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 175);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c35_rhs175, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs175, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs175), "rhs",
                  "rhs", 175);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs175), "lhs",
                  "lhs", 175);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 176);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 176);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c35_rhs176, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs176, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs176), "rhs",
                  "rhs", 176);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs176), "lhs",
                  "lhs", 176);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 177);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 177);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c35_rhs177, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs177, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs177), "rhs",
                  "rhs", 177);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs177), "lhs",
                  "lhs", 177);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 178);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 178);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 178);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c35_rhs178, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs178, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs178), "rhs",
                  "rhs", 178);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs178), "lhs",
                  "lhs", 178);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 179);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 179);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 179);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c35_rhs179, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs179, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs179), "rhs",
                  "rhs", 179);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs179), "lhs",
                  "lhs", 179);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 180);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 180);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c35_rhs180, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs180, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs180), "rhs",
                  "rhs", 180);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs180), "lhs",
                  "lhs", 180);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 181);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("intmin"), "name", "name",
                  181);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 181);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 181);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c35_rhs181, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs181, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs181), "rhs",
                  "rhs", 181);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs181), "lhs",
                  "lhs", 181);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m"),
                  "context", "context", 182);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_div"), "name", "name",
                  182);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 182);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 182);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c35_rhs182, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs182, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs182), "rhs",
                  "rhs", 182);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs182), "lhs",
                  "lhs", 182);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 183);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("mtimes"), "name", "name",
                  183);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 183);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c35_rhs183, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs183, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs183), "rhs",
                  "rhs", 183);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs183), "lhs",
                  "lhs", 183);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "context", "context", 184);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 184);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1326731922U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c35_rhs184, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs184, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs184), "rhs",
                  "rhs", 184);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs184), "lhs",
                  "lhs", 184);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(""), "context", "context",
                  185);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("tanh"), "name", "name", 185);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "resolved",
                  "resolved", 185);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1343833988U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c35_rhs185, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs185, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs185), "rhs",
                  "rhs", 185);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs185), "lhs",
                  "lhs", 185);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "context",
                  "context", 186);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("eml_scalar_tanh"), "name",
                  "name", 186);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c35_info, c35_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tanh.m"),
                  "resolved", "resolved", 186);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(1286822340U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c35_info, c35_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c35_rhs186, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c35_lhs186, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_rhs186), "rhs",
                  "rhs", 186);
  sf_mex_addfield(*c35_info, sf_mex_duplicatearraysafe(&c35_lhs186), "lhs",
                  "lhs", 186);
  sf_mex_destroy(&c35_rhs128);
  sf_mex_destroy(&c35_lhs128);
  sf_mex_destroy(&c35_rhs129);
  sf_mex_destroy(&c35_lhs129);
  sf_mex_destroy(&c35_rhs130);
  sf_mex_destroy(&c35_lhs130);
  sf_mex_destroy(&c35_rhs131);
  sf_mex_destroy(&c35_lhs131);
  sf_mex_destroy(&c35_rhs132);
  sf_mex_destroy(&c35_lhs132);
  sf_mex_destroy(&c35_rhs133);
  sf_mex_destroy(&c35_lhs133);
  sf_mex_destroy(&c35_rhs134);
  sf_mex_destroy(&c35_lhs134);
  sf_mex_destroy(&c35_rhs135);
  sf_mex_destroy(&c35_lhs135);
  sf_mex_destroy(&c35_rhs136);
  sf_mex_destroy(&c35_lhs136);
  sf_mex_destroy(&c35_rhs137);
  sf_mex_destroy(&c35_lhs137);
  sf_mex_destroy(&c35_rhs138);
  sf_mex_destroy(&c35_lhs138);
  sf_mex_destroy(&c35_rhs139);
  sf_mex_destroy(&c35_lhs139);
  sf_mex_destroy(&c35_rhs140);
  sf_mex_destroy(&c35_lhs140);
  sf_mex_destroy(&c35_rhs141);
  sf_mex_destroy(&c35_lhs141);
  sf_mex_destroy(&c35_rhs142);
  sf_mex_destroy(&c35_lhs142);
  sf_mex_destroy(&c35_rhs143);
  sf_mex_destroy(&c35_lhs143);
  sf_mex_destroy(&c35_rhs144);
  sf_mex_destroy(&c35_lhs144);
  sf_mex_destroy(&c35_rhs145);
  sf_mex_destroy(&c35_lhs145);
  sf_mex_destroy(&c35_rhs146);
  sf_mex_destroy(&c35_lhs146);
  sf_mex_destroy(&c35_rhs147);
  sf_mex_destroy(&c35_lhs147);
  sf_mex_destroy(&c35_rhs148);
  sf_mex_destroy(&c35_lhs148);
  sf_mex_destroy(&c35_rhs149);
  sf_mex_destroy(&c35_lhs149);
  sf_mex_destroy(&c35_rhs150);
  sf_mex_destroy(&c35_lhs150);
  sf_mex_destroy(&c35_rhs151);
  sf_mex_destroy(&c35_lhs151);
  sf_mex_destroy(&c35_rhs152);
  sf_mex_destroy(&c35_lhs152);
  sf_mex_destroy(&c35_rhs153);
  sf_mex_destroy(&c35_lhs153);
  sf_mex_destroy(&c35_rhs154);
  sf_mex_destroy(&c35_lhs154);
  sf_mex_destroy(&c35_rhs155);
  sf_mex_destroy(&c35_lhs155);
  sf_mex_destroy(&c35_rhs156);
  sf_mex_destroy(&c35_lhs156);
  sf_mex_destroy(&c35_rhs157);
  sf_mex_destroy(&c35_lhs157);
  sf_mex_destroy(&c35_rhs158);
  sf_mex_destroy(&c35_lhs158);
  sf_mex_destroy(&c35_rhs159);
  sf_mex_destroy(&c35_lhs159);
  sf_mex_destroy(&c35_rhs160);
  sf_mex_destroy(&c35_lhs160);
  sf_mex_destroy(&c35_rhs161);
  sf_mex_destroy(&c35_lhs161);
  sf_mex_destroy(&c35_rhs162);
  sf_mex_destroy(&c35_lhs162);
  sf_mex_destroy(&c35_rhs163);
  sf_mex_destroy(&c35_lhs163);
  sf_mex_destroy(&c35_rhs164);
  sf_mex_destroy(&c35_lhs164);
  sf_mex_destroy(&c35_rhs165);
  sf_mex_destroy(&c35_lhs165);
  sf_mex_destroy(&c35_rhs166);
  sf_mex_destroy(&c35_lhs166);
  sf_mex_destroy(&c35_rhs167);
  sf_mex_destroy(&c35_lhs167);
  sf_mex_destroy(&c35_rhs168);
  sf_mex_destroy(&c35_lhs168);
  sf_mex_destroy(&c35_rhs169);
  sf_mex_destroy(&c35_lhs169);
  sf_mex_destroy(&c35_rhs170);
  sf_mex_destroy(&c35_lhs170);
  sf_mex_destroy(&c35_rhs171);
  sf_mex_destroy(&c35_lhs171);
  sf_mex_destroy(&c35_rhs172);
  sf_mex_destroy(&c35_lhs172);
  sf_mex_destroy(&c35_rhs173);
  sf_mex_destroy(&c35_lhs173);
  sf_mex_destroy(&c35_rhs174);
  sf_mex_destroy(&c35_lhs174);
  sf_mex_destroy(&c35_rhs175);
  sf_mex_destroy(&c35_lhs175);
  sf_mex_destroy(&c35_rhs176);
  sf_mex_destroy(&c35_lhs176);
  sf_mex_destroy(&c35_rhs177);
  sf_mex_destroy(&c35_lhs177);
  sf_mex_destroy(&c35_rhs178);
  sf_mex_destroy(&c35_lhs178);
  sf_mex_destroy(&c35_rhs179);
  sf_mex_destroy(&c35_lhs179);
  sf_mex_destroy(&c35_rhs180);
  sf_mex_destroy(&c35_lhs180);
  sf_mex_destroy(&c35_rhs181);
  sf_mex_destroy(&c35_lhs181);
  sf_mex_destroy(&c35_rhs182);
  sf_mex_destroy(&c35_lhs182);
  sf_mex_destroy(&c35_rhs183);
  sf_mex_destroy(&c35_lhs183);
  sf_mex_destroy(&c35_rhs184);
  sf_mex_destroy(&c35_lhs184);
  sf_mex_destroy(&c35_rhs185);
  sf_mex_destroy(&c35_lhs185);
  sf_mex_destroy(&c35_rhs186);
  sf_mex_destroy(&c35_lhs186);
}

static real_T c35_mpower(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_a)
{
  real_T c35_b_a;
  real_T c35_c_a;
  real_T c35_ak;
  real_T c35_d_a;
  real_T c35_e_a;
  real_T c35_b;
  c35_b_a = c35_a;
  c35_c_a = c35_b_a;
  c35_eml_scalar_eg(chartInstance);
  c35_ak = c35_c_a;
  c35_d_a = c35_ak;
  c35_eml_scalar_eg(chartInstance);
  c35_e_a = c35_d_a;
  c35_b = c35_d_a;
  return c35_e_a * c35_b;
}

static void c35_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_eml_error(SFc35_simulationInstanceStruct *chartInstance)
{
  int32_T c35_i118;
  static char_T c35_cv1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c35_u[30];
  const mxArray *c35_y = NULL;
  int32_T c35_i119;
  static char_T c35_cv2[4] = { 's', 'q', 'r', 't' };

  char_T c35_b_u[4];
  const mxArray *c35_b_y = NULL;
  for (c35_i118 = 0; c35_i118 < 30; c35_i118++) {
    c35_u[c35_i118] = c35_cv1[c35_i118];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 10, 0U, 1U, 0U, 2, 1, 30),
                FALSE);
  for (c35_i119 = 0; c35_i119 < 4; c35_i119++) {
    c35_b_u[c35_i119] = c35_cv2[c35_i119];
  }

  c35_b_y = NULL;
  sf_mex_assign(&c35_b_y, sf_mex_create("y", c35_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c35_y, 14, c35_b_y));
}

static void c35_b_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static boolean_T c35_eml_use_refblas(SFc35_simulationInstanceStruct
  *chartInstance)
{
  return FALSE;
}

static void c35_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance, int32_T
  c35_k, real_T c35_A_data[17], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[9])
{
  int32_T c35_b_k;
  int32_T c35_b_ldb;
  int32_T c35_i120;
  int32_T c35_c_k;
  real_T c35_alpha1;
  int32_T c35_c_ldb;
  real_T c35_beta1;
  char_T c35_TRANSB;
  char_T c35_TRANSA;
  ptrdiff_t c35_m_t;
  ptrdiff_t c35_n_t;
  int32_T c35_var;
  ptrdiff_t c35_k_t;
  ptrdiff_t c35_lda_t;
  int32_T c35_b_var;
  ptrdiff_t c35_ldb_t;
  ptrdiff_t c35_ldc_t;
  double * c35_alpha1_t;
  double * c35_Aia0_t;
  double * c35_Bib0_t;
  double * c35_beta1_t;
  double * c35_Cic0_t;
  c35_b_k = c35_k;
  c35_b_ldb = c35_ldb;
  for (c35_i120 = 0; c35_i120 < 9; c35_i120++) {
    c35_C[c35_i120] = 0.0;
  }

  c35_c_k = c35_b_k;
  c35_alpha1 = 1.0;
  c35_c_ldb = c35_b_ldb;
  c35_beta1 = 0.0;
  c35_TRANSB = 'N';
  c35_TRANSA = 'N';
  c35_m_t = (ptrdiff_t)(1);
  c35_n_t = (ptrdiff_t)(9);
  c35_var = c35_c_k;
  c35_k_t = (ptrdiff_t)(c35_var);
  c35_lda_t = (ptrdiff_t)(1);
  c35_b_var = c35_c_ldb;
  c35_ldb_t = (ptrdiff_t)(c35_b_var);
  c35_ldc_t = (ptrdiff_t)(1);
  c35_alpha1_t = (double *)(&c35_alpha1);
  c35_Aia0_t = (double *)(&c35_A_data[0]);
  c35_Bib0_t = (double *)(&c35_B[0]);
  c35_beta1_t = (double *)(&c35_beta1);
  c35_Cic0_t = (double *)(&c35_C[0]);
  dgemm(&c35_TRANSA, &c35_TRANSB, &c35_m_t, &c35_n_t, &c35_k_t, c35_alpha1_t,
        c35_Aia0_t, &c35_lda_t, c35_Bib0_t, &c35_ldb_t, c35_beta1_t, c35_Cic0_t,
        &c35_ldc_t);
}

static void c35_check_forloop_overflow_error(SFc35_simulationInstanceStruct
  *chartInstance, boolean_T c35_overflow)
{
  int32_T c35_i121;
  static char_T c35_cv3[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c35_u[34];
  const mxArray *c35_y = NULL;
  int32_T c35_i122;
  static char_T c35_cv4[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c35_b_u[23];
  const mxArray *c35_b_y = NULL;
  for (c35_i121 = 0; c35_i121 < 34; c35_i121++) {
    c35_u[c35_i121] = c35_cv3[c35_i121];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 10, 0U, 1U, 0U, 2, 1, 34),
                FALSE);
  for (c35_i122 = 0; c35_i122 < 23; c35_i122++) {
    c35_b_u[c35_i122] = c35_cv4[c35_i122];
  }

  c35_b_y = NULL;
  sf_mex_assign(&c35_b_y, sf_mex_create("y", c35_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c35_y, 14, c35_b_y));
}

static void c35_c_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_mldivide(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_A[81], real_T c35_B[81], real_T c35_Y[81])
{
  int32_T c35_i123;
  real_T c35_b_A[81];
  int32_T c35_info;
  int32_T c35_ipiv[9];
  int32_T c35_b_info;
  int32_T c35_c_info;
  int32_T c35_d_info;
  int32_T c35_i124;
  int32_T c35_b_i;
  int32_T c35_c_i;
  int32_T c35_ip;
  int32_T c35_j;
  int32_T c35_b_j;
  real_T c35_temp;
  int32_T c35_i125;
  real_T c35_c_A[81];
  int32_T c35_i126;
  real_T c35_d_A[81];
  for (c35_i123 = 0; c35_i123 < 81; c35_i123++) {
    c35_b_A[c35_i123] = c35_A[c35_i123];
  }

  c35_b_eml_matlab_zgetrf(chartInstance, c35_b_A, c35_ipiv, &c35_info);
  c35_b_info = c35_info;
  c35_c_info = c35_b_info;
  c35_d_info = c35_c_info;
  if (c35_d_info > 0) {
    c35_eml_warning(chartInstance);
  }

  for (c35_i124 = 0; c35_i124 < 81; c35_i124++) {
    c35_Y[c35_i124] = c35_B[c35_i124];
  }

  for (c35_b_i = 1; c35_b_i < 10; c35_b_i++) {
    c35_c_i = c35_b_i - 1;
    if (c35_ipiv[c35_c_i] != c35_c_i + 1) {
      c35_ip = c35_ipiv[c35_c_i];
      for (c35_j = 1; c35_j < 10; c35_j++) {
        c35_b_j = c35_j - 1;
        c35_temp = c35_Y[c35_c_i + 9 * c35_b_j];
        c35_Y[c35_c_i + 9 * c35_b_j] = c35_Y[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
          c35_ip, 1, 9, 1, 0) + 9 * c35_b_j) - 1];
        c35_Y[(c35_ip + 9 * c35_b_j) - 1] = c35_temp;
      }
    }
  }

  for (c35_i125 = 0; c35_i125 < 81; c35_i125++) {
    c35_c_A[c35_i125] = c35_b_A[c35_i125];
  }

  c35_c_eml_xtrsm(chartInstance, c35_c_A, c35_Y);
  for (c35_i126 = 0; c35_i126 < 81; c35_i126++) {
    c35_d_A[c35_i126] = c35_b_A[c35_i126];
  }

  c35_d_eml_xtrsm(chartInstance, c35_d_A, c35_Y);
}

static void c35_realmin(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_eps(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_eml_matlab_zgetrf(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_b_A[81], int32_T c35_ipiv[9], int32_T *c35_info)
{
  int32_T c35_i127;
  for (c35_i127 = 0; c35_i127 < 81; c35_i127++) {
    c35_b_A[c35_i127] = c35_A[c35_i127];
  }

  c35_b_eml_matlab_zgetrf(chartInstance, c35_b_A, c35_ipiv, c35_info);
}

static int32_T c35_eml_ixamax(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_n, real_T c35_x[81], int32_T c35_ix0)
{
  int32_T c35_idxmax;
  int32_T c35_b_n;
  int32_T c35_b_ix0;
  int32_T c35_c_n;
  int32_T c35_c_ix0;
  int32_T c35_ix;
  real_T c35_b_x;
  real_T c35_c_x;
  real_T c35_d_x;
  real_T c35_y;
  real_T c35_e_x;
  real_T c35_f_x;
  real_T c35_b_y;
  real_T c35_smax;
  int32_T c35_d_n;
  int32_T c35_b;
  int32_T c35_b_b;
  boolean_T c35_overflow;
  int32_T c35_k;
  int32_T c35_b_k;
  int32_T c35_a;
  real_T c35_g_x;
  real_T c35_h_x;
  real_T c35_i_x;
  real_T c35_c_y;
  real_T c35_j_x;
  real_T c35_k_x;
  real_T c35_d_y;
  real_T c35_s;
  c35_b_n = c35_n;
  c35_b_ix0 = c35_ix0;
  c35_c_n = c35_b_n;
  c35_c_ix0 = c35_b_ix0;
  if (c35_c_n < 1) {
    c35_idxmax = 0;
  } else {
    c35_idxmax = 1;
    if (c35_c_n > 1) {
      c35_ix = c35_c_ix0;
      c35_b_x = c35_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_ix, 1, 81, 1, 0) - 1];
      c35_c_x = c35_b_x;
      c35_d_x = c35_c_x;
      c35_y = muDoubleScalarAbs(c35_d_x);
      c35_e_x = 0.0;
      c35_f_x = c35_e_x;
      c35_b_y = muDoubleScalarAbs(c35_f_x);
      c35_smax = c35_y + c35_b_y;
      c35_d_n = c35_c_n;
      c35_b = c35_d_n;
      c35_b_b = c35_b;
      c35_overflow = (c35_b_b > 2147483646);
      if (c35_overflow) {
        c35_check_forloop_overflow_error(chartInstance, TRUE);
      }

      for (c35_k = 2; c35_k <= c35_d_n; c35_k++) {
        c35_b_k = c35_k;
        c35_a = c35_ix + 1;
        c35_ix = c35_a;
        c35_g_x = c35_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_ix, 1, 81, 1, 0) - 1];
        c35_h_x = c35_g_x;
        c35_i_x = c35_h_x;
        c35_c_y = muDoubleScalarAbs(c35_i_x);
        c35_j_x = 0.0;
        c35_k_x = c35_j_x;
        c35_d_y = muDoubleScalarAbs(c35_k_x);
        c35_s = c35_c_y + c35_d_y;
        if (c35_s > c35_smax) {
          c35_idxmax = c35_b_k;
          c35_smax = c35_s;
        }
      }
    }
  }

  return c35_idxmax;
}

static void c35_eml_xger(SFc35_simulationInstanceStruct *chartInstance, int32_T
  c35_m, int32_T c35_n, real_T c35_alpha1, int32_T c35_ix0, int32_T c35_iy0,
  real_T c35_A[81], int32_T c35_ia0, real_T c35_b_A[81])
{
  int32_T c35_i128;
  for (c35_i128 = 0; c35_i128 < 81; c35_i128++) {
    c35_b_A[c35_i128] = c35_A[c35_i128];
  }

  c35_b_eml_xger(chartInstance, c35_m, c35_n, c35_alpha1, c35_ix0, c35_iy0,
                 c35_b_A, c35_ia0);
}

static void c35_eml_warning(SFc35_simulationInstanceStruct *chartInstance)
{
  int32_T c35_i129;
  static char_T c35_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c35_u[27];
  const mxArray *c35_y = NULL;
  for (c35_i129 = 0; c35_i129 < 27; c35_i129++) {
    c35_u[c35_i129] = c35_varargin_1[c35_i129];
  }

  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", c35_u, 10, 0U, 1U, 0U, 2, 1, 27),
                FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c35_y));
}

static void c35_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_A[81], real_T c35_B[81], real_T c35_b_B[81])
{
  int32_T c35_i130;
  int32_T c35_i131;
  real_T c35_b_A[81];
  for (c35_i130 = 0; c35_i130 < 81; c35_i130++) {
    c35_b_B[c35_i130] = c35_B[c35_i130];
  }

  for (c35_i131 = 0; c35_i131 < 81; c35_i131++) {
    c35_b_A[c35_i131] = c35_A[c35_i131];
  }

  c35_c_eml_xtrsm(chartInstance, c35_b_A, c35_b_B);
}

static void c35_below_threshold(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_b_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81], real_T c35_b_B[81])
{
  int32_T c35_i132;
  int32_T c35_i133;
  real_T c35_b_A[81];
  for (c35_i132 = 0; c35_i132 < 81; c35_i132++) {
    c35_b_B[c35_i132] = c35_B[c35_i132];
  }

  for (c35_i133 = 0; c35_i133 < 81; c35_i133++) {
    c35_b_A[c35_i133] = c35_A[c35_i133];
  }

  c35_d_eml_xtrsm(chartInstance, c35_b_A, c35_b_B);
}

static void c35_d_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_b_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B_data[17], int32_T c35_B_sizes[1], real_T c35_C
  [9], real_T c35_b_C[9])
{
  int32_T c35_i134;
  int32_T c35_i135;
  real_T c35_b_A[81];
  int32_T c35_b_B_sizes;
  int32_T c35_loop_ub;
  int32_T c35_i136;
  real_T c35_b_B_data[17];
  for (c35_i134 = 0; c35_i134 < 9; c35_i134++) {
    c35_b_C[c35_i134] = c35_C[c35_i134];
  }

  for (c35_i135 = 0; c35_i135 < 81; c35_i135++) {
    c35_b_A[c35_i135] = c35_A[c35_i135];
  }

  c35_b_B_sizes = c35_B_sizes[0];
  c35_loop_ub = c35_B_sizes[0] - 1;
  for (c35_i136 = 0; c35_i136 <= c35_loop_ub; c35_i136++) {
    c35_b_B_data[c35_i136] = c35_B_data[c35_i136];
  }

  c35_d_eml_xgemm(chartInstance, c35_b_A, c35_b_B_data, *(int32_T (*)[1])&
                  c35_b_B_sizes, c35_b_C);
}

static void c35_b_below_threshold(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_e_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_tanh(SFc35_simulationInstanceStruct *chartInstance, real_T
                     c35_x[9], real_T c35_b_x[9])
{
  int32_T c35_i137;
  for (c35_i137 = 0; c35_i137 < 9; c35_i137++) {
    c35_b_x[c35_i137] = c35_x[c35_i137];
  }

  c35_b_tanh(chartInstance, c35_b_x);
}

static void c35_f_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_g_eml_scalar_eg(SFc35_simulationInstanceStruct *chartInstance)
{
}

static void c35_c_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_k, real_T c35_A_data[153], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[81], real_T c35_b_C[81])
{
  int32_T c35_i138;
  int32_T c35_b_A_sizes[2];
  int32_T c35_A;
  int32_T c35_b_A;
  int32_T c35_loop_ub;
  int32_T c35_i139;
  real_T c35_b_A_data[153];
  int32_T c35_i140;
  real_T c35_b_B[81];
  for (c35_i138 = 0; c35_i138 < 81; c35_i138++) {
    c35_b_C[c35_i138] = c35_C[c35_i138];
  }

  c35_b_A_sizes[0] = 9;
  c35_b_A_sizes[1] = c35_A_sizes[1];
  c35_A = c35_b_A_sizes[0];
  c35_b_A = c35_b_A_sizes[1];
  c35_loop_ub = c35_A_sizes[0] * c35_A_sizes[1] - 1;
  for (c35_i139 = 0; c35_i139 <= c35_loop_ub; c35_i139++) {
    c35_b_A_data[c35_i139] = c35_A_data[c35_i139];
  }

  for (c35_i140 = 0; c35_i140 < 81; c35_i140++) {
    c35_b_B[c35_i140] = c35_B[c35_i140];
  }

  c35_e_eml_xgemm(chartInstance, c35_k, c35_b_A_data, c35_b_A_sizes, c35_b_B,
                  c35_ldb, c35_b_C);
}

static const mxArray *c35_h_sf_marshallOut(void *chartInstanceVoid, void
  *c35_inData)
{
  const mxArray *c35_mxArrayOutData = NULL;
  int32_T c35_u;
  const mxArray *c35_y = NULL;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_mxArrayOutData = NULL;
  c35_u = *(int32_T *)c35_inData;
  c35_y = NULL;
  sf_mex_assign(&c35_y, sf_mex_create("y", &c35_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c35_mxArrayOutData, c35_y, FALSE);
  return c35_mxArrayOutData;
}

static int32_T c35_l_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId)
{
  int32_T c35_y;
  int32_T c35_i141;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), &c35_i141, 1, 6, 0U, 0, 0U, 0);
  c35_y = c35_i141;
  sf_mex_destroy(&c35_u);
  return c35_y;
}

static void c35_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c35_mxArrayInData, const char_T *c35_varName, void *c35_outData)
{
  const mxArray *c35_b_sfEvent;
  const char_T *c35_identifier;
  emlrtMsgIdentifier c35_thisId;
  int32_T c35_y;
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)chartInstanceVoid;
  c35_b_sfEvent = sf_mex_dup(c35_mxArrayInData);
  c35_identifier = c35_varName;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_y = c35_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c35_b_sfEvent),
    &c35_thisId);
  sf_mex_destroy(&c35_b_sfEvent);
  *(int32_T *)c35_outData = c35_y;
  sf_mex_destroy(&c35_mxArrayInData);
}

static uint8_T c35_m_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_b_is_active_c35_simulation, const char_T
  *c35_identifier)
{
  uint8_T c35_y;
  emlrtMsgIdentifier c35_thisId;
  c35_thisId.fIdentifier = c35_identifier;
  c35_thisId.fParent = NULL;
  c35_y = c35_n_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c35_b_is_active_c35_simulation), &c35_thisId);
  sf_mex_destroy(&c35_b_is_active_c35_simulation);
  return c35_y;
}

static uint8_T c35_n_emlrt_marshallIn(SFc35_simulationInstanceStruct
  *chartInstance, const mxArray *c35_u, const emlrtMsgIdentifier *c35_parentId)
{
  uint8_T c35_y;
  uint8_T c35_u0;
  sf_mex_import(c35_parentId, sf_mex_dup(c35_u), &c35_u0, 1, 3, 0U, 0, 0U, 0);
  c35_y = c35_u0;
  sf_mex_destroy(&c35_u);
  return c35_y;
}

static void c35_b_eml_matlab_zgetrf(SFc35_simulationInstanceStruct
  *chartInstance, real_T c35_A[81], int32_T c35_ipiv[9], int32_T *c35_info)
{
  int32_T c35_i142;
  int32_T c35_j;
  int32_T c35_b_j;
  int32_T c35_a;
  int32_T c35_jm1;
  int32_T c35_b;
  int32_T c35_mmj;
  int32_T c35_b_a;
  int32_T c35_c;
  int32_T c35_b_b;
  int32_T c35_jj;
  int32_T c35_c_a;
  int32_T c35_jp1j;
  int32_T c35_d_a;
  int32_T c35_b_c;
  int32_T c35_i143;
  int32_T c35_i144;
  int32_T c35_i145;
  real_T c35_b_A[81];
  int32_T c35_e_a;
  int32_T c35_jpiv_offset;
  int32_T c35_f_a;
  int32_T c35_c_b;
  int32_T c35_jpiv;
  int32_T c35_g_a;
  int32_T c35_d_b;
  int32_T c35_c_c;
  int32_T c35_e_b;
  int32_T c35_jrow;
  int32_T c35_h_a;
  int32_T c35_f_b;
  int32_T c35_jprow;
  int32_T c35_ix0;
  int32_T c35_iy0;
  int32_T c35_b_ix0;
  int32_T c35_b_iy0;
  int32_T c35_c_ix0;
  int32_T c35_c_iy0;
  int32_T c35_ix;
  int32_T c35_iy;
  int32_T c35_k;
  real_T c35_temp;
  int32_T c35_i_a;
  int32_T c35_j_a;
  int32_T c35_b_jp1j;
  int32_T c35_k_a;
  int32_T c35_d_c;
  int32_T c35_l_a;
  int32_T c35_g_b;
  int32_T c35_i146;
  int32_T c35_m_a;
  int32_T c35_h_b;
  int32_T c35_n_a;
  int32_T c35_i_b;
  boolean_T c35_overflow;
  int32_T c35_b_i;
  int32_T c35_c_i;
  real_T c35_x;
  real_T c35_y;
  real_T c35_z;
  int32_T c35_j_b;
  int32_T c35_e_c;
  int32_T c35_o_a;
  int32_T c35_f_c;
  int32_T c35_p_a;
  int32_T c35_g_c;
  int32_T c35_m;
  int32_T c35_n;
  int32_T c35_d_ix0;
  int32_T c35_d_iy0;
  int32_T c35_ia0;
  real_T c35_d7;
  c35_realmin(chartInstance);
  c35_eps(chartInstance);
  for (c35_i142 = 0; c35_i142 < 9; c35_i142++) {
    c35_ipiv[c35_i142] = 1 + c35_i142;
  }

  *c35_info = 0;
  for (c35_j = 1; c35_j < 9; c35_j++) {
    c35_b_j = c35_j;
    c35_a = c35_b_j - 1;
    c35_jm1 = c35_a;
    c35_b = c35_b_j;
    c35_mmj = 9 - c35_b;
    c35_b_a = c35_jm1;
    c35_c = c35_b_a * 10;
    c35_b_b = c35_c + 1;
    c35_jj = c35_b_b;
    c35_c_a = c35_jj + 1;
    c35_jp1j = c35_c_a;
    c35_d_a = c35_mmj;
    c35_b_c = c35_d_a;
    c35_i143 = 0;
    for (c35_i144 = 0; c35_i144 < 9; c35_i144++) {
      for (c35_i145 = 0; c35_i145 < 9; c35_i145++) {
        c35_b_A[c35_i145 + c35_i143] = c35_A[c35_i145 + c35_i143];
      }

      c35_i143 += 9;
    }

    c35_e_a = c35_eml_ixamax(chartInstance, c35_b_c + 1, c35_b_A, c35_jj) - 1;
    c35_jpiv_offset = c35_e_a;
    c35_f_a = c35_jj;
    c35_c_b = c35_jpiv_offset;
    c35_jpiv = c35_f_a + c35_c_b;
    if (c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_jpiv, 1, 81, 1, 0) - 1] != 0.0)
    {
      if (c35_jpiv_offset != 0) {
        c35_g_a = c35_b_j;
        c35_d_b = c35_jpiv_offset;
        c35_c_c = c35_g_a + c35_d_b;
        c35_ipiv[c35_b_j - 1] = c35_c_c;
        c35_e_b = c35_jm1 + 1;
        c35_jrow = c35_e_b;
        c35_h_a = c35_jrow;
        c35_f_b = c35_jpiv_offset;
        c35_jprow = c35_h_a + c35_f_b;
        c35_ix0 = c35_jrow;
        c35_iy0 = c35_jprow;
        c35_b_ix0 = c35_ix0;
        c35_b_iy0 = c35_iy0;
        c35_c_ix0 = c35_b_ix0;
        c35_c_iy0 = c35_b_iy0;
        c35_ix = c35_c_ix0;
        c35_iy = c35_c_iy0;
        for (c35_k = 1; c35_k < 10; c35_k++) {
          c35_temp = c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_ix, 1, 81, 1, 0)
            - 1];
          c35_A[c35_ix - 1] = c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_iy, 1,
            81, 1, 0) - 1];
          c35_A[c35_iy - 1] = c35_temp;
          c35_i_a = c35_ix + 9;
          c35_ix = c35_i_a;
          c35_j_a = c35_iy + 9;
          c35_iy = c35_j_a;
        }
      }

      c35_b_jp1j = c35_jp1j;
      c35_k_a = c35_mmj;
      c35_d_c = c35_k_a;
      c35_l_a = c35_jp1j;
      c35_g_b = c35_d_c - 1;
      c35_i146 = c35_l_a + c35_g_b;
      c35_m_a = c35_b_jp1j;
      c35_h_b = c35_i146;
      c35_n_a = c35_m_a;
      c35_i_b = c35_h_b;
      if (c35_n_a > c35_i_b) {
        c35_overflow = FALSE;
      } else {
        c35_overflow = (c35_i_b > 2147483646);
      }

      if (c35_overflow) {
        c35_check_forloop_overflow_error(chartInstance, TRUE);
      }

      for (c35_b_i = c35_b_jp1j; c35_b_i <= c35_i146; c35_b_i++) {
        c35_c_i = c35_b_i;
        c35_x = c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_c_i, 1, 81, 1, 0) - 1];
        c35_y = c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_jj, 1, 81, 1, 0) - 1];
        c35_z = c35_x / c35_y;
        c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_c_i, 1, 81, 1, 0) - 1] = c35_z;
      }
    } else {
      *c35_info = c35_b_j;
    }

    c35_j_b = c35_b_j;
    c35_e_c = 9 - c35_j_b;
    c35_o_a = c35_jj;
    c35_f_c = c35_o_a;
    c35_p_a = c35_jj;
    c35_g_c = c35_p_a;
    c35_m = c35_mmj;
    c35_n = c35_e_c;
    c35_d_ix0 = c35_jp1j;
    c35_d_iy0 = c35_f_c + 9;
    c35_ia0 = c35_g_c + 10;
    c35_d7 = -1.0;
    c35_b_eml_xger(chartInstance, c35_m, c35_n, c35_d7, c35_d_ix0, c35_d_iy0,
                   c35_A, c35_ia0);
  }

  if (*c35_info == 0) {
    if (!(c35_A[80] != 0.0)) {
      *c35_info = 9;
    }
  }
}

static void c35_b_eml_xger(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_m, int32_T c35_n, real_T c35_alpha1, int32_T c35_ix0, int32_T
  c35_iy0, real_T c35_A[81], int32_T c35_ia0)
{
  int32_T c35_b_m;
  int32_T c35_b_n;
  int32_T c35_b_ix0;
  int32_T c35_b_iy0;
  int32_T c35_b_ia0;
  int32_T c35_c_m;
  int32_T c35_c_n;
  int32_T c35_c_ix0;
  int32_T c35_c_iy0;
  int32_T c35_c_ia0;
  int32_T c35_d_m;
  int32_T c35_d_n;
  int32_T c35_d_ix0;
  int32_T c35_d_iy0;
  int32_T c35_d_ia0;
  int32_T c35_ixstart;
  int32_T c35_a;
  int32_T c35_jA;
  int32_T c35_jy;
  int32_T c35_e_n;
  int32_T c35_j;
  real_T c35_yjy;
  real_T c35_temp;
  int32_T c35_ix;
  int32_T c35_b;
  int32_T c35_i147;
  int32_T c35_b_a;
  int32_T c35_b_b;
  int32_T c35_i148;
  int32_T c35_c_a;
  int32_T c35_c_b;
  int32_T c35_d_a;
  int32_T c35_d_b;
  boolean_T c35_overflow;
  int32_T c35_ijA;
  int32_T c35_b_ijA;
  int32_T c35_e_a;
  int32_T c35_f_a;
  int32_T c35_g_a;
  c35_b_m = c35_m;
  c35_b_n = c35_n;
  c35_b_ix0 = c35_ix0;
  c35_b_iy0 = c35_iy0;
  c35_b_ia0 = c35_ia0;
  c35_c_m = c35_b_m;
  c35_c_n = c35_b_n;
  c35_c_ix0 = c35_b_ix0;
  c35_c_iy0 = c35_b_iy0;
  c35_c_ia0 = c35_b_ia0;
  c35_d_m = c35_c_m;
  c35_d_n = c35_c_n;
  c35_d_ix0 = c35_c_ix0;
  c35_d_iy0 = c35_c_iy0;
  c35_d_ia0 = c35_c_ia0;
  c35_ixstart = c35_d_ix0;
  c35_a = c35_d_ia0 - 1;
  c35_jA = c35_a;
  c35_jy = c35_d_iy0;
  c35_e_n = c35_d_n;
  for (c35_j = 1; c35_j <= c35_e_n; c35_j++) {
    c35_yjy = c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_jy, 1, 81, 1, 0) - 1];
    if (c35_yjy != 0.0) {
      c35_temp = -c35_yjy;
      c35_ix = c35_ixstart;
      c35_b = c35_jA + 1;
      c35_i147 = c35_b;
      c35_b_a = c35_d_m;
      c35_b_b = c35_jA;
      c35_i148 = c35_b_a + c35_b_b;
      c35_c_a = c35_i147;
      c35_c_b = c35_i148;
      c35_d_a = c35_c_a;
      c35_d_b = c35_c_b;
      if (c35_d_a > c35_d_b) {
        c35_overflow = FALSE;
      } else {
        c35_overflow = (c35_d_b > 2147483646);
      }

      if (c35_overflow) {
        c35_check_forloop_overflow_error(chartInstance, TRUE);
      }

      for (c35_ijA = c35_i147; c35_ijA <= c35_i148; c35_ijA++) {
        c35_b_ijA = c35_ijA;
        c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_b_ijA, 1, 81, 1, 0) - 1] =
          c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_b_ijA, 1, 81, 1, 0) - 1] +
          c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_ix, 1, 81, 1, 0) - 1] *
          c35_temp;
        c35_e_a = c35_ix + 1;
        c35_ix = c35_e_a;
      }
    }

    c35_f_a = c35_jy + 9;
    c35_jy = c35_f_a;
    c35_g_a = c35_jA + 9;
    c35_jA = c35_g_a;
  }
}

static void c35_c_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81])
{
  real_T c35_alpha1;
  char_T c35_DIAGA;
  char_T c35_TRANSA;
  char_T c35_UPLO;
  char_T c35_SIDE;
  ptrdiff_t c35_m_t;
  ptrdiff_t c35_n_t;
  ptrdiff_t c35_lda_t;
  ptrdiff_t c35_ldb_t;
  double * c35_Aia0_t;
  double * c35_Bib0_t;
  double * c35_alpha1_t;
  c35_below_threshold(chartInstance);
  c35_alpha1 = 1.0;
  c35_DIAGA = 'U';
  c35_TRANSA = 'N';
  c35_UPLO = 'L';
  c35_SIDE = 'L';
  c35_m_t = (ptrdiff_t)(9);
  c35_n_t = (ptrdiff_t)(9);
  c35_lda_t = (ptrdiff_t)(9);
  c35_ldb_t = (ptrdiff_t)(9);
  c35_Aia0_t = (double *)(&c35_A[0]);
  c35_Bib0_t = (double *)(&c35_B[0]);
  c35_alpha1_t = (double *)(&c35_alpha1);
  dtrsm(&c35_SIDE, &c35_UPLO, &c35_TRANSA, &c35_DIAGA, &c35_m_t, &c35_n_t,
        c35_alpha1_t, c35_Aia0_t, &c35_lda_t, c35_Bib0_t, &c35_ldb_t);
}

static void c35_d_eml_xtrsm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B[81])
{
  real_T c35_alpha1;
  char_T c35_DIAGA;
  char_T c35_TRANSA;
  char_T c35_UPLO;
  char_T c35_SIDE;
  ptrdiff_t c35_m_t;
  ptrdiff_t c35_n_t;
  ptrdiff_t c35_lda_t;
  ptrdiff_t c35_ldb_t;
  double * c35_Aia0_t;
  double * c35_Bib0_t;
  double * c35_alpha1_t;
  c35_below_threshold(chartInstance);
  c35_alpha1 = 1.0;
  c35_DIAGA = 'N';
  c35_TRANSA = 'N';
  c35_UPLO = 'U';
  c35_SIDE = 'L';
  c35_m_t = (ptrdiff_t)(9);
  c35_n_t = (ptrdiff_t)(9);
  c35_lda_t = (ptrdiff_t)(9);
  c35_ldb_t = (ptrdiff_t)(9);
  c35_Aia0_t = (double *)(&c35_A[0]);
  c35_Bib0_t = (double *)(&c35_B[0]);
  c35_alpha1_t = (double *)(&c35_alpha1);
  dtrsm(&c35_SIDE, &c35_UPLO, &c35_TRANSA, &c35_DIAGA, &c35_m_t, &c35_n_t,
        c35_alpha1_t, c35_Aia0_t, &c35_lda_t, c35_Bib0_t, &c35_ldb_t);
}

static void c35_d_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  real_T c35_A[81], real_T c35_B_data[17], int32_T c35_B_sizes[1], real_T c35_C
  [9])
{
  real_T c35_b;
  real_T c35_b_b;
  real_T c35_flt;
  boolean_T c35_p;
  real_T c35_c_b;
  real_T c35_d_b;
  real_T c35_b_flt;
  boolean_T c35_b_p;
  int32_T c35_i149;
  int32_T c35_i150;
  int32_T c35_i151;
  int32_T c35_ic;
  int32_T c35_b_ic;
  int32_T c35_ar;
  int32_T c35_ib;
  int32_T c35_b_ib;
  real_T c35_temp;
  int32_T c35_ia;
  int32_T c35_c_ic;
  int32_T c35_a;
  int32_T c35_b_a;
  boolean_T guard1 = FALSE;
  if (c35_eml_use_refblas(chartInstance)) {
  } else {
    c35_b_below_threshold(chartInstance);
  }

  c35_e_eml_scalar_eg(chartInstance);
  c35_b = (real_T)c35_B_sizes[0];
  c35_b_b = c35_b;
  c35_flt = c35_b_b;
  c35_p = (9.0 == c35_flt);
  guard1 = FALSE;
  if (c35_p) {
    c35_c_b = (real_T)c35_B_sizes[0];
    c35_d_b = c35_c_b;
    c35_b_flt = c35_d_b;
    c35_b_p = (9.0 == c35_b_flt);
    if (c35_b_p) {
      for (c35_i149 = 0; c35_i149 < 9; c35_i149++) {
        c35_C[c35_i149] = 0.0;
        c35_i150 = 0;
        for (c35_i151 = 0; c35_i151 < 9; c35_i151++) {
          c35_C[c35_i149] += c35_A[c35_i150 + c35_i149] * c35_B_data[c35_i151];
          c35_i150 += 9;
        }
      }
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    for (c35_ic = 1; c35_ic < 10; c35_ic++) {
      c35_b_ic = c35_ic - 1;
      c35_C[c35_b_ic] = 0.0;
    }

    c35_ar = 0;
    for (c35_ib = 1; c35_ib < 10; c35_ib++) {
      c35_b_ib = c35_ib;
      if (c35_B_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_b_ib, 1, c35_B_sizes[0],
           1, 0) - 1] != 0.0) {
        c35_temp = c35_B_data[_SFD_EML_ARRAY_BOUNDS_CHECK("", c35_b_ib, 1,
          c35_B_sizes[0], 1, 0) - 1];
        c35_ia = c35_ar;
        for (c35_c_ic = 1; c35_c_ic < 10; c35_c_ic++) {
          c35_b_ic = c35_c_ic - 1;
          c35_a = c35_ia + 1;
          c35_ia = c35_a;
          c35_C[c35_b_ic] += c35_temp * c35_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
            c35_ia, 1, 81, 1, 0) - 1];
        }
      }

      c35_b_a = c35_ar + 9;
      c35_ar = c35_b_a;
    }
  }
}

static void c35_b_tanh(SFc35_simulationInstanceStruct *chartInstance, real_T
  c35_x[9])
{
  int32_T c35_k;
  real_T c35_b_k;
  real_T c35_b_x;
  real_T c35_c_x;
  for (c35_k = 0; c35_k < 9; c35_k++) {
    c35_b_k = 1.0 + (real_T)c35_k;
    c35_b_x = c35_x[(int32_T)c35_b_k - 1];
    c35_c_x = c35_b_x;
    c35_c_x = muDoubleScalarTanh(c35_c_x);
    c35_x[(int32_T)c35_b_k - 1] = c35_c_x;
  }
}

static void c35_e_eml_xgemm(SFc35_simulationInstanceStruct *chartInstance,
  int32_T c35_k, real_T c35_A_data[153], int32_T c35_A_sizes[2], real_T c35_B[81],
  int32_T c35_ldb, real_T c35_C[81])
{
  int32_T c35_b_k;
  int32_T c35_b_ldb;
  int32_T c35_c_k;
  real_T c35_alpha1;
  int32_T c35_c_ldb;
  real_T c35_beta1;
  char_T c35_TRANSB;
  char_T c35_TRANSA;
  ptrdiff_t c35_m_t;
  ptrdiff_t c35_n_t;
  int32_T c35_var;
  ptrdiff_t c35_k_t;
  ptrdiff_t c35_lda_t;
  int32_T c35_b_var;
  ptrdiff_t c35_ldb_t;
  ptrdiff_t c35_ldc_t;
  double * c35_alpha1_t;
  double * c35_Aia0_t;
  double * c35_Bib0_t;
  double * c35_beta1_t;
  double * c35_Cic0_t;
  c35_b_k = c35_k;
  c35_b_ldb = c35_ldb;
  c35_c_k = c35_b_k;
  c35_alpha1 = 1.0;
  c35_c_ldb = c35_b_ldb;
  c35_beta1 = 0.0;
  c35_TRANSB = 'N';
  c35_TRANSA = 'N';
  c35_m_t = (ptrdiff_t)(9);
  c35_n_t = (ptrdiff_t)(9);
  c35_var = c35_c_k;
  c35_k_t = (ptrdiff_t)(c35_var);
  c35_lda_t = (ptrdiff_t)(9);
  c35_b_var = c35_c_ldb;
  c35_ldb_t = (ptrdiff_t)(c35_b_var);
  c35_ldc_t = (ptrdiff_t)(9);
  c35_alpha1_t = (double *)(&c35_alpha1);
  c35_Aia0_t = (double *)(&c35_A_data[0]);
  c35_Bib0_t = (double *)(&c35_B[0]);
  c35_beta1_t = (double *)(&c35_beta1);
  c35_Cic0_t = (double *)(&c35_C[0]);
  dgemm(&c35_TRANSA, &c35_TRANSB, &c35_m_t, &c35_n_t, &c35_k_t, c35_alpha1_t,
        c35_Aia0_t, &c35_lda_t, c35_Bib0_t, &c35_ldb_t, c35_beta1_t, c35_Cic0_t,
        &c35_ldc_t);
}

static void init_dsm_address_info(SFc35_simulationInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c35_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1531140892U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2237912724U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1058641105U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(599793920U);
}

mxArray *sf_c35_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("pBMiY0EJcj4eioiEh38iGC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(9);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(9);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(9);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(9);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c35_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c35_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c35_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[10],T\"Sigma\",},{M[1],M[33],T\"cov_i\",},{M[1],M[5],T\"state_est\",},{M[1],M[32],T\"state_i\",},{M[8],M[0],T\"is_active_c35_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c35_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc35_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc35_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           35,
           1,
           1,
           13,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_simulationMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simulationMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simulationMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"Sigma");
          _SFD_SET_DATA_PROPS(1,1,1,0,"state_predict");
          _SFD_SET_DATA_PROPS(2,1,1,0,"Sigma_predict");
          _SFD_SET_DATA_PROPS(3,1,1,0,"z");
          _SFD_SET_DATA_PROPS(4,2,0,1,"state_est");
          _SFD_SET_DATA_PROPS(5,1,1,0,"covl");
          _SFD_SET_DATA_PROPS(6,1,1,0,"statel");
          _SFD_SET_DATA_PROPS(7,10,0,0,"nsonar");
          _SFD_SET_DATA_PROPS(8,1,1,0,"j");
          _SFD_SET_DATA_PROPS(9,10,0,0,"i");
          _SFD_SET_DATA_PROPS(10,2,0,1,"state_i");
          _SFD_SET_DATA_PROPS(11,2,0,1,"cov_i");
          _SFD_SET_DATA_PROPS(12,10,0,0,"kl");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,4,0,0,0,0,0,4,2);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1536);
        _SFD_CV_INIT_EML_IF(0,1,0,234,243,255,297);
        _SFD_CV_INIT_EML_IF(0,1,1,255,265,279,297);
        _SFD_CV_INIT_EML_IF(0,1,2,313,326,-1,341);
        _SFD_CV_INIT_EML_IF(0,1,3,421,443,-1,1467);

        {
          static int condStart[] = { 316, 323 };

          static int condEnd[] = { 319, 326 };

          static int pfixExpr[] = { 0, 1, -2 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,316,326,2,0,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 424, 432 };

          static int condEnd[] = { 428, 443 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,1,424,443,2,2,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_d_sf_marshallOut,(MexInFcnForType)
            c35_d_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c35_e_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_c_sf_marshallOut,(MexInFcnForType)
            c35_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c35_e_sf_marshallOut,(MexInFcnForType)
          c35_e_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c35_e_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c35_e_sf_marshallOut,(MexInFcnForType)
          c35_e_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_b_sf_marshallOut,(MexInFcnForType)
            c35_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_sf_marshallOut,(MexInFcnForType)
            c35_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 9;
          dimVector[1]= 9;
          _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c35_d_sf_marshallOut,(MexInFcnForType)
            c35_d_sf_marshallIn);
        }

        {
          real_T *c35_z;
          real_T *c35_j;
          real_T (*c35_Sigma)[81];
          real_T (*c35_state_predict)[9];
          real_T (*c35_Sigma_predict)[81];
          real_T (*c35_state_est)[9];
          real_T (*c35_covl)[4];
          real_T (*c35_statel)[2];
          real_T (*c35_state_i)[2];
          real_T (*c35_cov_i)[4];
          c35_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
          c35_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
          c35_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c35_statel = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
          c35_covl = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 3);
          c35_state_est = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S,
            2);
          c35_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c35_Sigma_predict = (real_T (*)[81])ssGetInputPortSignal
            (chartInstance->S, 1);
          c35_state_predict = (real_T (*)[9])ssGetInputPortSignal
            (chartInstance->S, 0);
          c35_Sigma = (real_T (*)[81])ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, *c35_Sigma);
          _SFD_SET_DATA_VALUE_PTR(1U, *c35_state_predict);
          _SFD_SET_DATA_VALUE_PTR(2U, *c35_Sigma_predict);
          _SFD_SET_DATA_VALUE_PTR(3U, c35_z);
          _SFD_SET_DATA_VALUE_PTR(4U, *c35_state_est);
          _SFD_SET_DATA_VALUE_PTR(5U, *c35_covl);
          _SFD_SET_DATA_VALUE_PTR(6U, *c35_statel);
          _SFD_SET_DATA_VALUE_PTR(7U, &chartInstance->c35_nsonar);
          _SFD_SET_DATA_VALUE_PTR(8U, c35_j);
          _SFD_SET_DATA_VALUE_PTR(9U, &chartInstance->c35_i);
          _SFD_SET_DATA_VALUE_PTR(10U, *c35_state_i);
          _SFD_SET_DATA_VALUE_PTR(11U, *c35_cov_i);
          _SFD_SET_DATA_VALUE_PTR(12U, chartInstance->c35_kl);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simulationMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "enD2oeLM1OTy9hAeoXHiqB";
}

static void sf_opaque_initialize_c35_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc35_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c35_simulation((SFc35_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c35_simulation((SFc35_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c35_simulation(void *chartInstanceVar)
{
  enable_c35_simulation((SFc35_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c35_simulation(void *chartInstanceVar)
{
  disable_c35_simulation((SFc35_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c35_simulation(void *chartInstanceVar)
{
  sf_c35_simulation((SFc35_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c35_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c35_simulation
    ((SFc35_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c35_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c35_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c35_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c35_simulation((SFc35_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c35_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c35_simulation(S);
}

static void sf_opaque_set_sim_state_c35_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c35_simulation(S, st);
}

static void sf_opaque_terminate_c35_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc35_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c35_simulation((SFc35_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc35_simulation((SFc35_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c35_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c35_simulation((SFc35_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c35_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     i kl nsonar
   */
  const char_T *rtParamNames[] = { "i", "kl", "nsonar" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for i*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for kl*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for nsonar*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      35);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,35,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,35,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,35);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,35,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,35,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 6; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,35);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2498740637U));
  ssSetChecksum1(S,(4276164908U));
  ssSetChecksum2(S,(1786057152U));
  ssSetChecksum3(S,(769854035U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c35_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c35_simulation(SimStruct *S)
{
  SFc35_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc35_simulationInstanceStruct *)utMalloc(sizeof
    (SFc35_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc35_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c35_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c35_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c35_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c35_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c35_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c35_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c35_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c35_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c35_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c35_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c35_simulation;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c35_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c35_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c35_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c35_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c35_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
