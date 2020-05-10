#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <limits>

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/OptimizationStopCondition.h>
#include <Bpp/Numeric/AutoParameter.h>

#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>


#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace bpp;

// #define p_vdel 0.65


// constants used when parsing files
#define MAXLLINE 100000
#define MAXLNAME 1000

// constants used for numerical calculations
#define QuasiZero 0.0001
#define OPTIM_TOLERANCE 0.00001
#define OPTIM_MAX_EVAL 1000
// #define OPTIM_MIN_EVAL 10
#define OPTIM_MIN_EVAL 2
#define LARGE_S_VALUE 500
#define VERY_LARGE_S_VALUE 10000
#define INTEGRAL_NB_RECTANGLES_EXTREME 1000
#define INTEGRAL_NB_RECTANGLES_MEDIAN 1000
#define GAMMA_INTEGRAL_BREAK 0.5
#define ARBITRARY_THETA 1.
#define PI 3.14159265358979323846


// constants used to parametrize the DFEM
#define SCALED_BETA_BOUNDARY 25
#define MAX_ORIENTATION_ERROR 0.1


// constants used to calculate likelihod surface and CI
#define LNL_TOL 2   // tolerable number of log-likelihood units around the maximum (used to calculate CI around parameters)
#define GSH_SURF_RANGE 0.5  // range for Gamma shape in Gamma likelihood surface calculation
#define GSH_SURF_NB 50    // number of values for Gamma shape in Gamma likelihood surface calculation
#define logGM_SURF_RANGE 0.5  // range for log Gamma mean in Gamma likelihood surface calculation
#define logGM_SURF_NB 50  // number of values for log Gamma mean in Gamma likelihood surface calculation

#define ARBITRARY_LENGTH 100000
#define ARBITRARY_LAMBDA 1.

int debug = 0;

// shared shape default initial value
#define DEFAULT_INITIAL_SHARED_SHAPE 0.2

struct MK_data
{
  char name[MAXLNAME + 1];
  int nb_gene;      // number of chromosomes sampled
  double Ls_poly, Ln_poly;  // number of S and NS sites for polymorphism data
  double* specS, * specN;    // specS[j] <- number of synonymous SNP at observed frequency j in the sample (1<=j<=nb_SFScat)
  int nb_SFScat;      // nb_gene/2 if folded; nb_gene-1 if unfolded
  double Ls_div, Ln_div;  // number of S and NS sites for divergence data
  double fixS, fixN;    // number of fixed S and N differences
  bool div_data_available;  // divergence data available
  bool folded;      // folded/unfolded SFS data
  vector< vector< vector<double> > > exp_count;
};


struct model
{
  string name;
  vector<string> param_name;
  map<string, bool> optimized;
  map<string, double> fixed_value;
  vector< vector<double> > constraint;
  vector<bool> shared;
  vector<string> reparametrization;
};

struct model current_model;


struct SFS_div
{
  vector<double> specS, specN;
  double divS, divN;
};


struct parameter_point
{
  // current model
  string model_name;
  // DFEM
  std::map<string, double> DFEM_param; // relevant DFEM parameters + orientation error
  // demography, data biases
  double theta;       // 4Ne.mu
  vector<double> ri;      // Eyre-Walker et al's (2006) factors
  // DFEM discretization
  vector<double> discr_bound;   // discretization boundaries
  vector<double> discr_mass;    // probability mass of intervals definedby discr_bound
  // data-dependent variables
  double lnL;       // log likelihood
  // struct SFS_div expected;              //expected SFS and divergence
  double alpha;       // expected proportion of substitutions whose 4Ne.s > neut_thresh
  double alpha_down, alpha_up;    // CI boundaries for alpha
  double omegaA;      // alpha dN/dS
  double omegaA_down, omegaA_up;  // CI boundaries for omegaA
  double omegaNA;     // (1-alpha) dN/dS
  double omegaNA_down, omegaNA_up;  // CI boundaries for omegaNA
  vector<double> discr_DFE_fixed; // discrete distribution of fitness effects of non-syno subst
};


struct fixed_parameters
{
  int nb;
  vector<string> model_name;
  vector<string> param_name;
  vector<double> value;
};


struct shared_parameter_point
{
  double negGshape;
  double posGshape;
  double pos_prop;
};


struct options
{
  vector<bool> do_model;            // which model to do (same index as in implemented_model_names)
  string model_name;
  bool use_divergence_data;         // true: use divergence data if any; false: only use SFS data
  bool use_divergence_parameter;    // true: include an extra parameter controlling dN (classical EW-like approach)
                                    // false: predict dN as well as SFS from parametrized DFEM
  bool use_syno_orientation_error;  // true: allow some percentage of misorientation specific to synonymous SNPs (in addition to the ri's)
                                    // false: assume equally probable misorientation at syno and non-syno SNPs (captured by the ri's)
  bool use_poisson;                 // use a Poisson or a multinomial likelihood model
  int nb_random_start;              // number of (additional) random starting avalues in likelihood optimization procedure
  double nearly_neutral_threshold;  // minimum Ne.s value defined as truly adaptive
  double FWW_threshold;             // threshold for low allele frequency removal in Fay et al (2001) alpha estimation
  struct fixed_parameters fixed_p;  // list of paraemeters that should be kept constant (not optimized) and their value
  bool do_separate, do_shared;
  vector<double> fixed_shared_negGshape;
  bool fold;
  double p_vdel;
};


/* global variales used throughut */

vector<string> implemented_model_names;
struct options opt;
string current_DFEM_fit;
map<string, vector<double> > constraints;
map<string, string> reparametrization;

int global_current_species;
vector<double> global_integral_S_point;
vector<double> global_integral_S_width;
map<double, double> S_width;

double lsqrtpi;
double* log_facto;
double** positive_dN_corrective_term;
double** negative_dN_corrective_term;

int nb_random_start_values;

double ancient_to_recent_Ne_ratio = 1.;


/* fin */
void fin(const char* s){printf("%s\n", s); exit(0);}

/* check_alloc */
void* check_alloc(int nb_blocks, size_t block_size)
{
  void* tab;
  tab = calloc(nb_blocks, block_size);
  if (tab == NULL)
    fin("Not enough memory");
  return tab;
}

/* printtime */
void printtime(double rtime)
{
  int d, h, m;
  int s;

  d = (int)rtime / 86400;
  rtime -= d * 86400;
  if (d)
    printf("%d day", d);
  if (d > 1)
    printf("s");
  if (d)
    printf(" ");

  h = (int)rtime / 3600;
  rtime -= h * 3600;
  if (h)
    printf("%d hour", h);
  if (h > 1)
    printf("s");
  if (h)
    printf(" ");

  m = (int)rtime / 60;
  rtime -= m * 60;
  if (m)
    printf("%d minute", m);
  if (m > 1)
    printf("s");
  if (m)
    printf(" ");

  s = (int)rtime;
  if (s)
    printf("%d second", s);
  if (s > 1)
    printf("s");

  if (d == 0 && h == 0 && m == 0 && s == 0)
    printf("less than 1 second");

  printf("\n");
}


/* precalculate_log_facto */
void precalculate_log_facto(int maxn)
{
  log_facto[0] = log_facto[1] = 0;
  for (int i = 2; i < maxn; i++)
  {
    log_facto[i] = log_facto[i - 1] + log(i);
  }
}


/* set_implemented_models*/

void set_implemented_models()
{
//  implemented_model_names.push_back("ExpoZero");
//  implemented_model_names.push_back("ExpoExpo");
//  implemented_model_names.push_back("ExpoGamma");
  implemented_model_names.push_back("GammaZero");
//  implemented_model_names.push_back("GammaExpo");
//  implemented_model_names.push_back("GammaGamma");
//  implemented_model_names.push_back("DisplGamma");
//  implemented_model_names.push_back("ScaledBeta");
//  implemented_model_names.push_back("FGMBesselK");
  implemented_model_names.push_back("ReflectedGamma");
}


/* set_reparametrization */

void set_reparametrization()
{
  string none = "none";
  string log = "log";

  reparametrization["negGmean"] = log;
  reparametrization["posGmean"] = log;
  reparametrization["BKscale"] = log;
  reparametrization["posGshape"] = none;
  reparametrization["negGshape"] = none;
  reparametrization["BKsigma"] = none;
  reparametrization["pos_prop"] = none;
  reparametrization["med_prop"] = none;
  reparametrization["s0"] = none;
  reparametrization["Ba"] = none;
  reparametrization["Bb"] = none;
  reparametrization["BKm"] = none;
  reparametrization["BKnsz2"] = none;
  reparametrization["syno_orientation_error"] = none;
}


/* set_constraints */

void set_constraints()
{
  // Gamma or Expo negative mean
  vector<double> constr_mean; constr_mean.push_back(0.1); constr_mean.push_back(pow(10., 20.));
  constraints["negGmean"] = constr_mean;

  // Gamma or Expo positive mean
  vector<double> constr_posmean; constr_posmean.push_back(QuasiZero); constr_posmean.push_back(pow(10., 4.));
  constraints["posGmean"] = constr_posmean;
  constraints["BKscale"] = constr_posmean;

  // Gamma shape
  vector<double> constr_shape; constr_shape.push_back(0.05); constr_shape.push_back(100.);
  constraints["posGshape"] = constr_shape;
  constraints["negGshape"] = constr_shape;
  constraints["BKsigma"] = constr_shape;

  // unconstrained proportion
  vector<double> constr_prop; constr_prop.push_back(0.); constr_prop.push_back(1.);
  constraints["med_prop"] = constr_prop;

  // displacement
  vector<double> constr_displ; constr_displ.push_back(-100.); constr_displ.push_back(0.);
  constraints["s0"] = constr_displ;

  // Beta parameters
  vector<double> constr_Beta; constr_Beta.push_back(QuasiZero); constr_Beta.push_back(10000.);
  constraints["Ba"] = constr_Beta;
  constraints["Bb"] = constr_Beta;

  // BesselK specifics
  vector<double> constr_BKm; constr_BKm.push_back(1.); constr_BKm.push_back(10.);
  vector<double> constr_BKnsz2; constr_BKnsz2.push_back(1.); constr_BKnsz2.push_back(10000.);
  constraints["BKm"] = constr_BKm;
  constraints["BKnsz2"] = constr_BKnsz2;

  // constrained proportion
  vector<double> constr_prop2; constr_prop2.push_back(0.); constr_prop2.push_back(MAX_ORIENTATION_ERROR);
  constraints["syno_orientation_error"] = constr_prop2;
  constraints["pos_prop"] = constr_prop2;
}


/***************  FROM NUMERICAL RECIPES:   ****************/
/***************   NUMERICAL INTEGRATION    ****************/

#define FUNC(x) ((*func)(x))
#define EPS 1.0e-6
// #define JMAX 20
#define JMAX 100
#define JMAXP (JMAX + 1)
#define K 5
double c[K + 1], d[K + 1]; // global variables used by polint

double trapzd(double (* func)(double), double a, double b, int n)
{
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if (n == 1)
  {
    return s = 0.5 * (b - a) * (FUNC(a) + FUNC(b));
  }
  else
  {
    for (it = 1, j = 1; j < n - 1; j++)
    {
      it <<= 1;
    }
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for (sum = 0.0, j = 1; j <= it; j++, x += del)
    {
      sum += FUNC(x);
    }
    s = 0.5 * (s + (b - a) * sum / tnm);
    return s;
  }
}


void my_polint(double xa[], double ya[], int n, double x, double* y, double* dy)
{
  int i, m, ns = 1;
  double den, dif, dift, ho, hp, w;


  dif = fabs(x - xa[1]);

  for (i = 1; i <= n; i++)
  {
    if ( (dift = fabs(x - xa[i])) < dif)
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
//	*y=ya[ns--];

  *y = ya[ns];
  ns--;

  for (m = 1; m < n; m++)
  {
    for (i = 1; i <= n - m; i++)
    {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if ( (den = ho - hp) == 0.0)
        fin("Error in routine my_polint");
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
//		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));

    if (2 * ns < n - m)
      *dy = c[ns + 1];
    else
    {
      *dy = d[ns]; ns--;
    }
    *y += (*dy);
  }
}

/* Returns the one-dimension integral of f  from a to b */

double qromb(double (* func)(double), double a, double b)
{
  double ss, dss;
  double s[JMAXP + 1], h[JMAXP + 1];
  int j;

  h[1] = 1.0;
  for (j = 1; j <= JMAX; j++)
  {
    double prov = trapzd(func, a, b, j);
    s[j] = prov;
    if (j >= K)
    {
      my_polint(&h[j - K], &s[j - K], K, 0.0, &ss, &dss);
      if (fabs(dss) < EPS * fabs(ss))
        return ss;
    }
    s[j + 1] = s[j];
    h[j + 1] = 0.25 * h[j];
  }
  fin("Too many steps in routine qromb");
  return 0.0;
}


/***************   FROM JEAN-PIERRE MOREAU WEB SITE:   ****************/
/***************   MODIFIED BESSEL FUNCTION 2nd KIND   ****************/


double BESSJ0 (double X)
{
/***********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
 ************************************************************************/
  const double
    P1 = 1.0, P2 = -0.1098628627E-2, P3 = 0.2734510407E-4,
    P4 = -0.2073370639E-5, P5 = 0.2093887211E-6,
    Q1 = -0.1562499995E-1, Q2 = 0.1430488765E-3, Q3 = -0.6911147651E-5,
    Q4 = 0.7621095161E-6, Q5 = -0.9349451520E-7,
    R1 = 57568490574.0, R2 = -13362590354.0, R3 = 651619640.7,
    R4 = -11214424.18, R5 = 77392.33017, R6 = -184.9052456,
    S1 = 57568490411.0, S2 = 1029532985.0, S3 = 9494680.718,
    S4 = 59272.64853, S5 = 267.8532712, S6 = 1.0;
  double
    AX, FR, FS, Z, FP, FQ, XX, Y, TMP;

  if (X == 0.0)
    return 1.0;
  AX = fabs(X);
  if (AX < 8.0)
  {
    Y = X * X;
    FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
    FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * S6))));
    TMP = FR / FS;
  }
  else
  {
    Z = 8. / AX;
    Y = Z * Z;
    XX = AX - 0.785398164;
    FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
    FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
    TMP = sqrt(0.636619772 / AX) * (FP * cos(XX) - Z * FQ * sin(XX));
  }
  return TMP;
}

double Sign(double X, double Y)
{
  if (Y < 0.0)
    return -fabs(X);
  else
    return fabs(X);
}

double BESSJ1 (double X)
{
/**********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 1, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
 ***********************************************************************/
  const double
    P1 = 1.0, P2 = 0.183105E-2, P3 = -0.3516396496E-4, P4 = 0.2457520174E-5,
    P5 = -0.240337019E-6,  P6 = 0.636619772,
    Q1 = 0.04687499995, Q2 = -0.2002690873E-3, Q3 = 0.8449199096E-5,
    Q4 = -0.88228987E-6, Q5 = 0.105787412E-6,
    R1 = 72362614232.0, R2 = -7895059235.0, R3 = 242396853.1,
    R4 = -2972611.439,   R5 = 15704.48260,  R6 = -30.16036606,
    S1 = 144725228442.0, S2 = 2300535178.0, S3 = 18583304.74,
    S4 = 99447.43394,    S5 = 376.9991397,  S6 = 1.0;

  double AX, FR, FS, Y, Z, FP, FQ, XX, TMP;

  AX = fabs(X);
  if (AX < 8.0)
  {
    Y = X * X;
    FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
    FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * S6))));
    TMP = X * (FR / FS);
  }
  else
  {
    Z = 8.0 / AX;
    Y = Z * Z;
    XX = AX - 2.35619491;
    FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
    FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
    TMP = sqrt(P6 / AX) * (cos(XX) * FP - Z * sin(XX) * FQ) * Sign(S6, X);
  }
  return TMP;
}


double BESSY0 (double X)
{
/* --------------------------------------------------------------------
      This subroutine calculates the Second Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
   --------------------------------------------------------------------- */
  const double
    P1 = 1.0, P2 = -0.1098628627E-2, P3 = 0.2734510407E-4,
    P4 = -0.2073370639E-5, P5 = 0.2093887211E-6,
    Q1 = -0.1562499995E-1, Q2 = 0.1430488765E-3, Q3 = -0.6911147651E-5,
    Q4 = 0.7621095161E-6, Q5 = -0.9349451520E-7,
    R1 = -2957821389.0, R2 = 7062834065.0, R3 = -512359803.6,
    R4 = 10879881.29,  R5 = -86327.92757, R6 = 228.4622733,
    S1 = 40076544269.0, S2 = 745249964.8, S3 = 7189466.438,
    S4 = 47447.26470,   S5 = 226.1030244, S6 = 1.0;
  double FS, FR, Z, FP, FQ, XX, Y;
  if (X == 0.0)
    return -1e30;
  if (X < 8.0)
  {
    Y = X * X;
    FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
    FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * S6))));
    return FR / FS + 0.636619772 * BESSJ0(X) * log(X);
  }
  else
  {
    Z = 8.0 / X;
    Y = Z * Z;
    XX = X - 0.785398164;
    FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
    FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
    return sqrt(0.636619772 / X) * (FP * sin(XX) + Z * FQ * cos(XX));
  }
}

double BESSY1 (double X)
{
/* ---------------------------------------------------------------------
      This subroutine calculates the Second Kind Bessel Function of
      order 1, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
   ---------------------------------------------------------------------- */
  const double
    P1 = 1.0, P2 = 0.183105E-2, P3 = -0.3516396496E-4,
    P4 = 0.2457520174E-5, P5 = -0.240337019E-6,
    Q1 = 0.04687499995, Q2 = -0.2002690873E-3, Q3 = 0.8449199096E-5,
    Q4 = -0.88228987E-6, Q5 = 0.105787412E-6,
    R1 = -0.4900604943E13, R2 = 0.1275274390E13, R3 = -0.5153438139E11,
    R4 = 0.7349264551E9,  R5 = -0.4237922726E7,  R6 = 0.8511937935E4,
    S1 = 0.2499580570E14, S2 = 0.4244419664E12, S3 = 0.3733650367E10,
    S4 = 0.2245904002E8,  S5 = 0.1020426050E6,  S6 = 0.3549632885E3, S7 = 1.0;
  double FR, FS, Z, FP, FQ, XX, Y;
  if (X == 0.0)
    return -1e30;
  if (X < 8.0)
  {
    Y = X * X;
    FR = R1 + Y * (R2 + Y * (R3 + Y * (R4 + Y * (R5 + Y * R6))));
    FS = S1 + Y * (S2 + Y * (S3 + Y * (S4 + Y * (S5 + Y * (S6 + Y * S7)))));
    return X * (FR / FS) + 0.636619772 * (BESSJ1(X) * log(X) - 1.0 / X);
  }
  else
  {
    Z = 8. / X;
    Y = Z * Z;
    XX = X - 2.356194491;
    FP = P1 + Y * (P2 + Y * (P3 + Y * (P4 + Y * P5)));
    FQ = Q1 + Y * (Q2 + Y * (Q3 + Y * (Q4 + Y * Q5)));
    return sqrt(0.636619772 / X) * (sin(XX) * FP + Z * cos(XX) * FQ);
  }
}

double BESSY (int N, double X)
{
/* -----------------------------------------------------------------
      This subroutine calculates the second kind Bessel Function of
      integer order N, for any real X. We use here the classical
      recursive formula.
   ------------------------------------------------------------------ */
  double TOX, BY, BYM, BYP; int J;
  if (N == 0)
    return BESSY0(X);
  if (N == 1)
    return BESSY1(X);
  if (X == 0.0)
    return -1e30;
  TOX = 2.0 / X;
  BY  = BESSY1(X);
  BYM = BESSY0(X);
  for (J = 1; J < N; J++)
  {
    BYP = J * TOX * BY - BYM;
    BYM = BY;
    BY  = BYP;
  }
  return BY;
}


/*******************************************/
/***********   EXPECTATIONS     ************/
/*******************************************/

/* WARNING: global variables used by the series of functions below: */

bool global_folded;
int global_j, global_n, global_nbcat;
double global_L, global_theta, global_lambda, global_Gmean, global_Gshape, global_posGmean, global_posGshape, global_pos_prop, global_Bb, global_Ba, global_med_prop, global_s0, global_BKm, global_BKsigma, global_BKscale, global_BKnsz2, global_S;

vector< vector< vector<double> > > global_exp_count;


/*  onearg_corrective_fct */
/* to be integrated to calculate TB's corrective term for dN */
double onearg_corrective_fct(double z)
{
  double S = global_S;
  int n = global_n;

  if (fabs(S) <= QuasiZero)
    return 1. / global_n;
  return (1 - exp(-S * (1 - z))) / (z * (1 - z) * (1 - exp(-S))) * pow(z, n);
}

/* precalculate_dN_corrective_term */
/* Pre-calculates (for integer values only: approx) the factor necessary to correct */
/* for polymorphic sites at which we have sampled a single type of alleles */
/* (thus seen as divergence sites) */

void precalculate_dN_corrective_term()
{
  for (double S = -LARGE_S_VALUE - QuasiZero; S < 0.; S += 1)
  {
    global_S = S;
    negative_dN_corrective_term[global_current_species][abs((int)round(S))] = qromb(&onearg_corrective_fct, QuasiZero, 1 - QuasiZero);
  }


  for (double S = 1. - QuasiZero; S < VERY_LARGE_S_VALUE; S += 1)
  {
    global_S = S;
    positive_dN_corrective_term[global_current_species][abs((int)round(S))] = qromb(&onearg_corrective_fct, QuasiZero, 1 - QuasiZero);
  }
}


/* onearg_sojourn_time_density */
/* Density of sojourn time spent at frequency x for a mutation under population selection coefficient S, with argument S passed as a global variable */

double onearg_sojourn_time_density(double x)
{
  double S = global_S;

  return 2 * (1 - exp(-S * (1 - x))) / (x * (1 - x) * (1 - exp(-S)));
}


/* onearg_sampling_j_out_of_n */
/* Binomial probability B(n, j, x), with arguments j and n passed as global variables */

double onearg_sampling_j_out_of_n(double x)
{
  int j = global_j;
  int n = global_n;
  double logx = log(x);
  double log1x = log(1 - x);

  double l1 = log_facto[n] - log_facto[j] - log_facto[n - j];

  return exp(l1 + j * logx + (n - j) * log1x);
}


/* sojourn_sampling_product */
/* Correspond to the product of the above 2 functions. */
/* (equivalent to the H*Q product in Eyre-Walker et al 2006) */
/* NB: makes use of global variables global_n, global_j, global_S */

double sojourn_sampling_product(double x)
{
  return onearg_sojourn_time_density(x) * onearg_sampling_j_out_of_n(x);
}


/* onearg_S_unfolded_expected_allele_count */

double onearg_S_unfolded_expected_allele_count(double S)
{
  // neutral case
  if (fabs(S) < QuasiZero)
    return global_L * global_theta / global_j;

  // selected case
  global_S = S;
  double inte = qromb(&sojourn_sampling_product, QuasiZero, 1. - QuasiZero);

  return global_L * global_theta * inte / 2.;
}


/* onearg_S_folded_expected_allele_count */

double onearg_S_folded_expected_allele_count(double S)
{
  if (global_n % 2 == 0 && global_j == global_n / 2)
    return onearg_S_unfolded_expected_allele_count(S);

  double fold1 = onearg_S_unfolded_expected_allele_count(S);
  global_j = global_n - global_j;
  double fold2 = onearg_S_unfolded_expected_allele_count(S);
  global_j = global_n - global_j;

  return fold1 + fold2;
}


/* onearg_negative_cst_density(S) */
/*
   double onearg_negative_cst_density(double mx){

   double S=-global_Gmean;
   if(mx>=S && mx-S<=S_width[mx]/2.) return 1./S_width[mx];
   if(mx<S && S-mx<S_width[mx]/2.) return 1./S_width[mx];

   return 0;

   }
 */


/* onearg_negative_Expo_density */
/* Density of a negative Exponential distribution, with Gmean passed as global variables */
/* argument mx must be negative or zero is returned */

double onearg_negative_Expo_density(double mx)
{
  if (mx > 0.)
    return 0.;
  double x = -mx;
  return exp(-x / global_Gmean) / global_Gmean;
}


/* onearg_negative_Expo_positive_Expo_density */
/* Density of a negative-Gammapositive-exponential distribution, with Gmean, Gshape, posGmean and pos_prop passed as global variables */

double onearg_negative_Expo_positive_Expo_density(double mx)
{
  if (mx > 0.)
    return global_pos_prop * exp(-mx / global_posGmean) / global_posGmean;
  else
  {
    double x = -mx;
    return (1 - global_pos_prop) * exp(-x / global_Gmean) / global_Gmean;
  }
}


/* onearg_negative_Expo_positive_Gamma_density */
/* Density of a negative-Gammapositive-Gamma distribution, with negGmean, negGshape, posGmean, posGshape, and pos_prop passed as global variables */

double onearg_negative_Expo_positive_Gamma_density(double mx)
{
  double lp;

  if (mx > 0.)
  {
    lp = (global_posGshape - 1) * log(mx) - (mx * global_posGshape / global_posGmean) + global_posGshape * log(global_posGshape / global_posGmean);
    lp -= lgamma(global_posGshape);
    return global_pos_prop * exp(lp);
  }
  else
  {
    double x = -mx;
    return (1 - global_pos_prop) * exp(-x / global_Gmean) / global_Gmean;
  }
}


/* onearg_negative_Gamma_density */
/* Density of a negative Gamma distribution, with Gmean and Gshape passed as global variables */
/* argument mx must be negative or zero is returned */

double onearg_negative_Gamma_density(double mx)
{
  double lp;

  if (mx > 0.)
    return 0.;
  double x = -mx;

  lp = (global_Gshape - 1) * log(x) - (x * global_Gshape / global_Gmean) + global_Gshape * log(global_Gshape / global_Gmean);
  lp -= lgamma(global_Gshape);
  return exp(lp);
}


/* onearg_reflected_Gamma_density */
/* Density of a negative Gamma distribution, with Gmean and Gshape passed as global variables */
/* argument mx can be negative or positive */

double onearg_reflected_Gamma_density(double mx)
{
  double lp;

  double x = fabs(mx);

  lp = (global_Gshape - 1) * log(x) - (x * global_Gshape / global_Gmean) + global_Gshape * log(global_Gshape / global_Gmean);
  lp -= lgamma(global_Gshape);
  return exp(lp) / (1 + exp(mx));
}


/* onearg_negative_Gamma_positive_expo_density */
/* Density of a negative-Gammapositive-exponential distribution, with Gmean, Gshape, posGmean and pos_prop passed as global variables */

double onearg_negative_Gamma_positive_Expo_density(double mx)
{
  double lp;

  if (mx > 0.)
    return global_pos_prop * exp(-mx / global_posGmean) / global_posGmean;
  else
  {
    double x = -mx;
    lp = (global_Gshape - 1) * log(x) - (x * global_Gshape / global_Gmean) + global_Gshape * log(global_Gshape / global_Gmean);
    lp -= lgamma(global_Gshape);
    return (1. - global_pos_prop) * exp(lp);
  }
}


/* onearg_negative_Gamma_positive_Gamma_density */
/* Density of a negative-Gammapositive-Gamma distribution, with negGmean, negGshape, posGmean, posGshape, and pos_prop passed as global variables */

double onearg_negative_Gamma_positive_Gamma_density(double mx)
{
  double lp;

  if (mx > 0.)
  {
    lp = (global_posGshape - 1) * log(mx) - (mx * global_posGshape / global_posGmean) + global_posGshape * log(global_posGshape / global_posGmean);
    lp -= lgamma(global_posGshape);
    return global_pos_prop * exp(lp);
  }
  else
  {
    double x = -mx;
    lp = (global_Gshape - 1) * log(x) - (x * global_Gshape / global_Gmean) + global_Gshape * log(global_Gshape / global_Gmean);
    lp -= lgamma(global_Gshape);
    return (1. - global_pos_prop) * exp(lp);
  }
}


/* onearg_negative_displaced_Gamma_density */
/* Density of a negative displaced Gamma distribution, with Gmean, Gshape and s0 (displacement) passed as global variables */
/* argument mx must be below s0 or zero is returned */

double onearg_negative_displaced_Gamma_density(double mx)
{
  double lp;

  if (mx > global_s0)
    return 0.;
  double x = -mx;

  lp = (global_Gshape - 1) * log(x - global_s0) - ((x - global_s0) * global_Gshape / global_Gmean) + global_Gshape * log(global_Gshape / global_Gmean);
  lp -= lgamma(global_Gshape);
  return exp(lp);
}


/* onearg_scaled_Beta_density */
/* Density of a centered, scaled Gamma distribution, with Bb, Ba, and med_prop passed as global variables */
/* argument mx must be between -SCALED_BETA_BOUNDARY and +SCALED_BETA_BOUNDARY or zero is returned */

double onearg_scaled_Beta_density(double x)
{
  if (x < -SCALED_BETA_BOUNDARY)
    return 0.;
  if (x > SCALED_BETA_BOUNDARY)
    return 0.;

  return global_med_prop * gsl_ran_beta_pdf((x / (2 * SCALED_BETA_BOUNDARY)) + 0.5, global_Ba, global_Bb) / (2 * SCALED_BETA_BOUNDARY);
}


/* onearg_BesselK_density */
/* Density of the FGM-associated DFEM derived by Lourenco et al 2011 Evolution */

double onearg_BesselK_density(double x)
{
  double lp;

  if (x == 0.)
    return (onearg_BesselK_density(-QuasiZero) + onearg_BesselK_density(QuasiZero)) / 2.;

  x = x / global_BKscale;

  double m = global_BKm;
  double nsz2 = global_BKnsz2;
  double s = global_BKsigma;
  // double z=global_BKz;

  double s2 = s * s;
  // double z2=z*z;
  double x2 = x * x;

  lp = -log(2.) * (m + 1.) / 2.;
  lp -= x * nsz2 / 4.;
  lp += log(nsz2) / 2.;
  double A = (1. + 4 / (nsz2 * s2)) / x2;
  lp += log(A) * (1. - m) / 4.;
  lp -= m * log(s);
  lp -= lsqrtpi;
  lp -= log(gsl_sf_gamma(m / 2.));
  double B = 1. / (nsz2 * x2 * (nsz2 + 4. / s2));
  double C = 1. / (4 * sqrt(B));
  if (C > 50.)
    return 0.;
  lp += log(gsl_sf_bessel_Knu((m - 1.) / 2., C));

  return exp(lp) / global_BKscale;
}


/* precalculate_integral_S_width */

vector<double> precalculate_integral_S_width()
{
  double width;
  vector<double> v_width;
  int k = global_integral_S_point.size();

  for (int i = 0; i < k; i++)
  {
    if (i == 0)
      width = global_integral_S_point[i + 1] - global_integral_S_point[i];
    else if (i == k - 1)
      width = global_integral_S_point[i] - global_integral_S_point[i - 1];
    else
      width = (global_integral_S_point[i + 1] - global_integral_S_point[i - 1]) / 2.;
    v_width.push_back(width);
    S_width[global_integral_S_point[i]] = width;
  }

  return v_width;
}


/* precalculate_integral_S_points */
/* Returns  a vector of values of S for which expected counts will be pre-calculated */
/* and used for numerical integration of the likelihood */

vector<double> precalculate_integral_S_points()
{
  vector<double> v;

  double inte_break_S = GAMMA_INTEGRAL_BREAK;
  double large = LARGE_S_VALUE;
  double very_large = VERY_LARGE_S_VALUE;
  int k1 = INTEGRAL_NB_RECTANGLES_EXTREME;
  int k2 = INTEGRAL_NB_RECTANGLES_MEDIAN;

  // large negative: from -LARGE_S_VALUE to -GAMMA_INTEGRAL_BREAK
  double a1 = (large - inte_break_S) / k1;
  for (int i = 0; i < k1; i++)
  {
    v.push_back(-large + (i + 0.5) * a1);
  }

  // medium values: from -GAMMA_INTEGRAL_BREAK to GAMMA_INTEGRAL_BREAK
  double a2 = 2. * inte_break_S / k2;
  for (int i = 0; i < k2; i++)
  {
    v.push_back(-inte_break_S + (i + 0.5) * a2);
  }

  // large positive: from GAMMA_INTEGRAL_BREAK to LARGE_S_VALUE
  for (int i = 0; i < k1; i++)
  {
    v.push_back(inte_break_S + (i + 0.5) * a1);
  }

  // very large positive: from LARGE_S_VALUE to VERY_LARGE_S_VALUE
  double a3 = (very_large - large) / k1;
  for (int i = 0; i < k1; i++)
  {
    v.push_back(large + (i + 0.5) * a3);
  }

  return v;
}


/* precalculate_expected_allele_counts */

vector<double> precalculate_expected_allele_counts(int j)
{
  vector<double> v;

  global_j = j;

  for (size_t i = 0; i < global_integral_S_point.size(); i++)
  {
    if (global_folded)
      v.push_back(onearg_S_folded_expected_allele_count(global_integral_S_point[i]));
    else
      v.push_back(onearg_S_unfolded_expected_allele_count(global_integral_S_point[i]));
  }
  return v;
}


/* negative_density */
double negative_density(string model_name, double S)
{
  if (model_name == "ExpoZero")
    return onearg_negative_Expo_density(S);
  else if (model_name == "ExpoExpo")
    return onearg_negative_Expo_positive_Expo_density(S);
  else if (model_name == "ExpoGamma")
    return onearg_negative_Expo_positive_Gamma_density(S);
  else if (model_name == "GammaZero")
    return onearg_negative_Gamma_density(S);
  else if (model_name == "GammaExpo")
    return onearg_negative_Gamma_positive_Expo_density(S);
  else if (model_name == "GammaGamma")
    return onearg_negative_Gamma_positive_Gamma_density(S);
  else if (model_name == "DisplGamma")
    return onearg_negative_displaced_Gamma_density(S);
  else if (model_name == "ScaledBeta")
    return onearg_scaled_Beta_density(S);
  else if (model_name == "FGMBesselK")
    return onearg_BesselK_density(S);
  else if (model_name == "ReflectedGamma")
    return onearg_reflected_Gamma_density(S);
  else
    return -1.;
}


/* expected_allele_count_product */
/* Product of DFEM density by S-conditional folded expected allele count */
/* (equivalent to D*H*Q in Eyre-Walker et al 2006)*/
/* NB: makes use of global variables global_n, global_j, global_mean, global_shape, global_posGmean, global_pos_prop, global_s0 */

double expected_allele_count_product(double S, int S_rank, int j)
{
  global_S = S;

  double a = negative_density(current_DFEM_fit, S);

  double b = global_exp_count[global_current_species][j][S_rank];

  return a * b;
}


/* expected_allele_count */
/* Expected number of SNPs of the j vs n-j sort in a sample */
/* of size n taken from a sequence of length L in a WF population */
/* of population mutation rate theta. The DFEM is specificied by global */
/* variable current_DFEM_fit (string) and defined by (some of) parameters */
/* Gmean, Gshape, posGmean, pos_prop and s0*/
/* Integration is entirely numerical */

double expected_allele_count(double L, double theta, map<string, double> param, int n, int j)
{
  global_j = j;
  global_n = n;
  global_L = L;
  global_theta = theta;

  std::map<string, double>::iterator it;
  for (it = param.begin(); it != param.end(); it++)
  {
    if (it->first == "negGmean")
      global_Gmean = it->second;
    else if (it->first == "negGshape")
      global_Gshape = it->second;
    else if (it->first == "posGmean")
      global_posGmean = it->second;
    else if (it->first == "posGshape")
      global_posGshape = it->second;
    else if (it->first == "pos_prop")
      global_pos_prop = it->second;
    else if (it->first == "s0")
      global_s0 = it->second;
    else if (it->first == "Ba")
      global_Ba = it->second;
    else if (it->first == "Bb")
      global_Bb = it->second;
    else if (it->first == "med_prop")
      global_med_prop = it->second;
    else if (it->first == "BKnsz2")
      global_BKnsz2 = it->second;
    else if (it->first == "BKm")
      global_BKm = it->second;
    else if (it->first == "BKsigma")
      global_BKsigma = it->second;
    else if (it->first == "BKscale")
      global_BKscale = it->second;
  }

  int k = global_integral_S_point.size();
  double inte_exp_count = 0.;

  for (int i = 0; i < k; i++)
  {
    inte_exp_count += expected_allele_count_product(global_integral_S_point[i], i, j) * global_integral_S_width[i];
  }

  return inte_exp_count;
}


/* S_expected_divergence */

double S_expected_divergence(double L, double lambda, double theta, int n, double S)
{
  // neutral case
  if (fabs(S) < QuasiZero)
    return L * (lambda + theta / n);

  // selected case
  double dN_correction;
  if (S > 0)
    dN_correction = positive_dN_corrective_term[global_current_species][abs((int)round(S))];
  else
    dN_correction = negative_dN_corrective_term[global_current_species][abs((int)round(S))];

  return L * (lambda * S / (1 - exp(-S)) + theta * dN_correction);
}


/* onearg_S_expected_divergence */

double onearg_S_expected_divergence(double S)
{
  // neutral case
  if (fabs(S) < QuasiZero)
    return global_L * (global_lambda + global_theta / global_n);

  // selected case
  double dN_correction;
  if (S > 0)
    dN_correction = positive_dN_corrective_term[global_current_species][abs((int)round(S))];
  else
    dN_correction = negative_dN_corrective_term[global_current_species][abs((int)round(S))];

  return global_L * (global_lambda * S / (1 - exp(-S)) + global_theta * dN_correction);
}


/* expected_div_product */

double expected_div_product(double S)
{
  double a = negative_density(current_DFEM_fit, S);
  double b = onearg_S_expected_divergence(S);

  return a * b;
}


/* expected_divergence */
/* Expected number of fixed differences; DFEM is specified by the global variable current_DFEM_fit */
/* L: sequence length; lambda:divergence time*mutation rate; theta: 2*pop size*mutation rate; */
/* Gmean, Gshape, posGmean, pos_prop, s0: parameters of the DFEM; n:sample size */

double expected_divergence(double L, double lambda, double theta, map<string, double> param, int n)
{
  global_L = L;
  global_lambda = lambda;
  global_theta = theta * ancient_to_recent_Ne_ratio;
  global_n = n;

  std::map<string, double>::iterator it;

  for (it = param.begin(); it != param.end(); it++)
  {
    if (it->first == "negGmean")
      global_Gmean = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "negGshape")
      global_Gshape = it->second;
    else if (it->first == "posGmean")
      global_posGmean = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "posGshape")
      global_posGshape = it->second;
    else if (it->first == "pos_prop")
      global_pos_prop = it->second;
    else if (it->first == "s0")
      global_s0 = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "Ba")
      global_Ba = it->second;
    else if (it->first == "Bb")
      global_Bb = it->second;
    else if (it->first == "med_prop")
      global_med_prop = it->second / ancient_to_recent_Ne_ratio;
    else if (it->first == "BKnsz2")
      global_BKnsz2 = it->second;
    else if (it->first == "BKm")
      global_BKm = it->second;
    else if (it->first == "BKsigma")
      global_BKsigma = it->second;
    else if (it->first == "BKscale")
      global_BKscale = it->second *= ancient_to_recent_Ne_ratio;
  }

  int k = global_integral_S_point.size();
  double exp_div = 0.;

  for (int i = 0; i < k; i++)
  {
    exp_div += expected_div_product(global_integral_S_point[i]) * global_integral_S_width[i];
  }

  return exp_div;
}


/* expected_divergence_truncated */
/* Expected number of fixed differences; DFEM is specified by the global variable current_DFEM_fit */
/* L: sequence length; lambda:divergence time*mutation rate; theta: 2*pop size*mutation rate; */
/* Gmean, Gshape, posGmean, pos_prop, s0: parameters of the DFEM; n:sample size */
/* only values of S below minS are considered */

double expected_divergence_truncated(double L, double lambda, double theta, map<string, double> param, int n, double minS)
{
  global_L = L;
  global_lambda = lambda;
  global_theta = theta * ancient_to_recent_Ne_ratio;
  global_n = n;

  std::map<string, double>::iterator it;
  for (it = param.begin(); it != param.end(); it++)
  {
    if (it->first == "negGmean")
      global_Gmean = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "negGshape")
      global_Gshape = it->second;
    else if (it->first == "posGmean")
      global_posGmean = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "posGshape")
      global_posGshape = it->second;
    else if (it->first == "pos_prop")
      global_pos_prop = it->second;
    else if (it->first == "s0")
      global_s0 = it->second * ancient_to_recent_Ne_ratio;
    else if (it->first == "Ba")
      global_Ba = it->second;
    else if (it->first == "Bb")
      global_Bb = it->second;
    else if (it->first == "med_prop")
      global_med_prop = it->second / ancient_to_recent_Ne_ratio;
    else if (it->first == "BKnsz2")
      global_BKnsz2 = it->second;
    else if (it->first == "BKm")
      global_BKm = it->second;
    else if (it->first == "BKsigma")
      global_BKsigma = it->second;
    else if (it->first == "BKscale")
      global_BKscale = it->second * ancient_to_recent_Ne_ratio;
  }

  int k = global_integral_S_point.size();
  double exp_div = 0.;

  for (int i = 0; i < k; i++)
  {
    if (global_integral_S_point[i] > minS)
      continue;
    exp_div += expected_div_product(global_integral_S_point[i]) * global_integral_S_width[i];
  }

  return exp_div;
}


/* Gamma_Welch_expected_div */
/* Expected dN based on Welch 2008 JME equation 23 */

double Gamma_Welch_expected_div(double L, double lambda, double Gmean, double Gshape)
{
  return L * lambda * pow(Gshape, Gshape + 1.) * (1 + (pow(0.666666, Gshape) / Gshape)) * pow(Gmean, -Gshape);
}


/*****************************************************/
/*****  Likelihood calculation and optimisation  *****/
/*****************************************************/

// global variables used to avoid alloc/free at each likelihood calculation
vector<double> prov_exp_specS, exp_specS, exp_specN;


/* create_model */

struct model create_model(string modname)
{
  struct model m;

  m.name = modname;
  if (modname == "ExpoZero")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
  }
  else if (modname == "ExpoExpo")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("posGmean"); m.optimized["posGmean"] = true;
    m.param_name.push_back("pos_prop"); m.optimized["pos_prop"] = true;
  }
  else if (modname == "ExpoGamma")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("posGmean"); m.optimized["posGmean"] = true;
    m.param_name.push_back("posGshape"); m.optimized["posGshape"] = true;
    m.param_name.push_back("pos_prop"); m.optimized["pos_prop"] = true;
  }
  else if (modname == "GammaZero")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("negGshape"); m.optimized["negGshape"] = true;
  }
  else if (modname == "GammaExpo")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("negGshape"); m.optimized["negGshape"] = true; /* m.optimized["negGshape"]=false; m.fixed_value["negGshape"]=5.;*/
    m.param_name.push_back("posGmean"); m.optimized["posGmean"] = true;
    m.param_name.push_back("pos_prop"); m.optimized["pos_prop"] = true;
  }
  else if (modname == "GammaGamma")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("negGshape"); m.optimized["negGshape"] = true;
    m.param_name.push_back("posGmean"); m.optimized["posGmean"] = true;
    m.param_name.push_back("posGshape"); m.optimized["posGshape"] = true;
    m.param_name.push_back("pos_prop"); m.optimized["pos_prop"] = true;
  }
  else if (modname == "DisplGamma")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("negGshape"); m.optimized["negGshape"] = true;
    m.param_name.push_back("s0"); m.optimized["s0"] = true;
  }
  else if (modname == "ScaledBeta")
  {
    m.param_name.push_back("Ba"); m.optimized["Ba"] = true;
    m.param_name.push_back("Bb"); m.optimized["Bb"] = true;
    m.param_name.push_back("med_prop"); m.optimized["med_prop"] = true;
  }
  else if (modname == "FGMBesselK")
  {
    m.param_name.push_back("BKm"); m.optimized["BKm"] = true;
    m.param_name.push_back("BKnsz2"); m.optimized["BKnsz2"] = true;
    m.param_name.push_back("BKsigma"); m.optimized["BKsigma"] = true;
    m.param_name.push_back("BKscale"); m.optimized["BKscale"] = true;
  }
  else if (modname == "ReflectedGamma")
  {
    m.param_name.push_back("negGmean"); m.optimized["negGmean"] = true;
    m.param_name.push_back("negGshape"); m.optimized["negGshape"] = true;
  }

  if (!global_folded && opt.use_syno_orientation_error)
  {
    m.param_name.push_back("syno_orientation_error"); m.optimized["syno_orientation_error"] = true;
  }

  for (size_t i = 0; i < m.param_name.size(); i++)
  {
    m.reparametrization.push_back(reparametrization[m.param_name[i]]);
    m.constraint.push_back(constraints[m.param_name[i]]);
  }

  return m;
}

/* manage_fixed_params */
void manage_fixed_params(vector<struct model>* model_list, struct fixed_parameters fixed_p)
{
  if (fixed_p.nb == 0)
    return;

  for (size_t i = 0; i < (*model_list).size(); i++)
  {
    for (int j = 0; j < fixed_p.nb; j++)
    {
      if ((*model_list)[i].name == fixed_p.model_name[j])
      {
        (*model_list)[i].optimized[fixed_p.param_name[j]] = false;
        if (reparametrization[fixed_p.param_name[j]] == "log")
          (*model_list)[i].fixed_value[fixed_p.param_name[j]] = log(fixed_p.value[j]);
        else
          (*model_list)[i].fixed_value[fixed_p.param_name[j]] = fixed_p.value[j];
        // cout <<"Model " <<fixed_p.model_name[j]  <<": parameter " <<fixed_p.param_name[j] <<" fixed to " <<fixed_p.value[j] <<endl;
      }
    }
  }
}


/* random_starting_point */

map<string, double> random_starting_point(struct model* m)
{
  map<string, double> pp;
  vector<string> param = m->param_name;

  for (size_t i = 0; i < param.size(); i++)
  {
    double c_down = constraints[param[i]][0];
    double c_up = constraints[param[i]][1];
    double ran = -1.;

    double dran = drand48();

    if (reparametrization[param[i]] == "none")
      ran = c_down + dran * (c_up - c_down);
    else if (reparametrization[param[i]] == "log")
    {
      if (c_down == 0.)
        c_down = QuasiZero;
      ran = exp(log(c_down) + dran * (log(c_up) - log(c_down)));
    }

    if (ran >= c_up)
      ran = c_up - QuasiZero;
    if (ran <= c_down)
      ran = c_down + QuasiZero;

    pp[param[i]] = ran;

// cout <<"starting " <<param[i] <<": " <<ran <<endl;
  }

  return pp;
}


/* create_parameter_point */

struct parameter_point create_parameter_point(string model)
{
  struct parameter_point pp;

  pp.discr_bound.push_back(-100.); pp.discr_bound.push_back(-25); pp.discr_bound.push_back(-10.); pp.discr_bound.push_back(-2.); pp.discr_bound.push_back(0.);
  pp.discr_bound.push_back(2.); pp.discr_bound.push_back(10.); pp.discr_bound.push_back(25); pp.discr_bound.push_back(100.);

  pp.model_name = model;

  if (model == "ExpoZero")
  {
    pp.DFEM_param["negGmean"] = 1000.;
  }
  else if (model == "ExpoExpo")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["posGmean"] = 1000.;
    pp.DFEM_param["pos_prop"] = 0.5;
  }
  else if (model == "ExpoGamma")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["posGmean"] = 1000.;
    pp.DFEM_param["posGshape"] = 1.;
    pp.DFEM_param["pos_prop"] = 0.5;
  }
  else if (model == "GammaZero")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["negGshape"] = 1.;
  }
  else if (model == "GammaExpo")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["negGshape"] = 1.;
    pp.DFEM_param["posGmean"] = 1000.;
    pp.DFEM_param["pos_prop"] = 0.5;
  }
  else if (model == "GammaGamma")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["negGshape"] = 1.;
    pp.DFEM_param["posGmean"] = 1000.;
    pp.DFEM_param["posGshape"] = 1.;
    pp.DFEM_param["pos_prop"] = 0.5;
  }
  else if (model == "DisplGamma")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["negGshape"] = 1.;
    pp.DFEM_param["s0"] = -1.;
  }
  else if (model == "ScaledBeta")
  {
    pp.DFEM_param["Ba"] = 1.;
    pp.DFEM_param["Bb"] = 1.;
    pp.DFEM_param["med_prop"] = 0.5;
  }
  else if (model == "FGMBesselK")
  {
    pp.DFEM_param["BKm"] = 2.;
    pp.DFEM_param["BKnsz2"] = 100.;
    pp.DFEM_param["BKsigma"] = 0.5;
    pp.DFEM_param["BKscale"] = 1000.;
  }
  else if (model == "ReflectedGamma")
  {
    pp.DFEM_param["negGmean"] = 1000.;
    pp.DFEM_param["negGshape"] = 1.;
  }

  if (!global_folded && opt.use_syno_orientation_error) pp.DFEM_param["syno_orientation_error"] = 0.01;

  return pp;
}


/* minus_SFS_loglikelihood */
/* SFS likelihood function: */
/* data must contain the observed spectra */
/* exp_pS must contain the expected synonymous spectrum */
/* exp_pN must contain the expected non-synonymous spectrum */
/* exp_dS must contain the expected synonymous divergence */
/* exp_dN must contain the expected non-synonymous divergence */

double minus_SFS_loglikelihood(struct MK_data data, vector<double> exp_pS, vector<double> exp_pN, double exp_dS, double exp_dN)
{
  int nc = data.nb_SFScat;
  double l, lnL;
  int k;

  if (opt.use_poisson) // calculate likelihood (poisson)
  {
    lnL = 0.;
    for (int j = 1; j <= nc; j++)
    {
      l = exp_pS[j];
      if (l <= 0)
        continue;
      k = (int)round(data.specS[j]);
      lnL += k * log(l) - log_facto[k] - l;
      l = exp_pN[j];
      if (l <= 0)
        continue;
      k = (int)round(data.specN[j]);
      lnL += k * log(l) - log_facto[k] - l;
    }
    if (data.div_data_available)
    {
      l = exp_dS;
      k = (int)round(data.fixS);
      if (l != 0)
        lnL += k * log(l) - log_facto[k] - l;
      l = exp_dN;
      k = (int)round(data.fixN);
      if (l != 0)
        lnL += k * log(l) - log_facto[k] - l;
    }
  }
  else    // calculate likelihood (multinomial)
  {
    int totobs = 0.;
    double totexp = 0.;
    for (int j = 1; j <= nc; j++)
    {
      totobs += (int)round(data.specS[j]) + (int)round(data.specN[j]);
      totexp += exp_pS[j] + exp_pN[j];
    }
    if (data.div_data_available)
      totobs += (int)round(data.fixS) + (int)round(data.fixN);
    totexp += exp_dS + exp_dN;
    lnL = log_facto[totobs];
    for (int j = 1; j <= nc; j++)
    {
      l = exp_pS[j] / totexp;
      if (l <= 0)
        continue;
      k = (int)round(data.specS[j]);
      lnL += k * log(l) - log_facto[k];
      l = exp_pN[j] / totexp;
      if (l <= 0)
        continue;
      k = (int)round(data.specN[j]);
      lnL += k * log(l) - log_facto[k];
    }
    if (data.div_data_available)
    {
      l = exp_dS / totexp;
      k = (int)round(data.fixS);
      if (l != 0)
        lnL += k * log(l) - log_facto[k];
      l = exp_dN / totexp;
      k = (int)round(data.fixN);
      if (l != 0)
        lnL += k * log(l) - log_facto[k];
    }
  }

  return -lnL;
}


/* optimize_saturated */
double optimize_saturated(struct MK_data* data)
{
  int nc = data->nb_SFScat;
  for (int j = 1; j <= nc; j++)
  {
    exp_specS[j] = data->specS[j];
    exp_specN[j] = data->specN[j];
  }

  return minus_SFS_loglikelihood(*data, exp_specS, exp_specN, data->fixS, data->fixN);
}


/* optimize_neutral */

double optimize_neutral(struct MK_data* data, double* MLE_oer, double* MLE_neut, struct SFS_div* expectations)
{
  int nc = data->nb_SFScat;
  int ng = data->nb_gene;
  double tot_syn, tot_nsyn, tot_syn_sites, tot_nsyn_sites, f;

  // f0
  tot_syn = tot_nsyn = tot_syn_sites = tot_nsyn_sites = 0;
  for (int j = 1; j <= nc; j++)
  {
    tot_syn += data->specS[j];
    tot_nsyn += data->specN[j];
    tot_syn_sites += data->Ls_poly;
    tot_nsyn_sites += data->Ln_poly;
    if (!opt.use_divergence_parameter)
    {
      tot_syn += data->fixS;
      tot_nsyn += data->fixN;
      tot_syn_sites += data->Ls_div;
      tot_nsyn_sites += data->Ln_div;
    }
  }
  f = tot_nsyn / (tot_nsyn + tot_syn);
  if (MLE_neut)
    *MLE_neut = f * tot_syn_sites / tot_nsyn_sites;

  // expected SFS
  for (int j = 1; j <= nc; j++)
  {
    exp_specS[j] = (1. - f) * (data->specS[j] + data->specN[j]);
    exp_specN[j] = f * (data->specS[j] + data->specN[j]);
  }

  // expected divergence
  double exp_dS, exp_dN;

  if (opt.use_divergence_parameter)
  {
    exp_dN = data->fixN;
    exp_dS = data->fixS;
  }
  else
  {
    exp_dN = data->fixS * tot_nsyn / tot_syn;
    exp_dS = data->fixS;
    double coeff = (data->fixN + data->fixS) / (exp_dN + exp_dS);
    exp_dN *= coeff;
    exp_dS *= coeff;
  }


  // calculate likelihood, optimizing orientation locally if needed
  if (data->folded == true || !opt.use_syno_orientation_error)
  {
    if (expectations)
    {
      expectations->divS = exp_dS; expectations->divN = exp_dN;
      expectations->specS = exp_specS; expectations->specN = exp_specN;
    }
    return minus_SFS_loglikelihood(*data, exp_specS, exp_specN, exp_dS, exp_dN);
  }
  else// my ugly optimization of orientation error
  {
    double maxL = minus_SFS_loglikelihood(*data, exp_specS, exp_specN, exp_dS, exp_dN);
    double best_oer = 0.;
    vector<double> or_exp_specS(nc + 2);
    int nbval = 100;
    for (int i = 1; i <= nbval; i++)
    {
      double oer = i * (MAX_ORIENTATION_ERROR / (double)nbval);
      for (int j = 1; j <= nc; j++)
      {
        or_exp_specS[j] = exp_specS[j] * (1. - oer) + exp_specS[ng - j] * oer;
      }
      double L = minus_SFS_loglikelihood(*data, or_exp_specS, exp_specN, exp_dS, exp_dN);
      if (L < maxL)
      {
        maxL = L; best_oer = oer;
      }
    }
    if (MLE_oer)
      *MLE_oer = best_oer;
    if (expectations)
    {
      expectations->divS = exp_dS; expectations->divN = exp_dN;
      expectations->specS = exp_specS; expectations->specN = exp_specN;
    }
    return maxL;
  }
}


/* expected_SFS_div */
/* Calculate expected SFS and divergence, for a given data set and set of parameters */
/* Current model and options are specified by global variables */
/* Optimal theta and rj are returned if the last 2 args are not passed as NULL */

struct SFS_div expected_SFS_div(struct MK_data data, map<string, double> param, double* the, vector<double>* rj)
{
  int nc = data.nb_SFScat;
  int ng = data.nb_gene;


  /* 1. SFS */

  // calculate neutral SFS expectations assuming theta=1 and rj's=1
  for (int j = 1; j <= nc; j++)
  {
    prov_exp_specS[j] = data.Ls_poly / j;
    if (data.folded && (ng % 2 != 0 || j != ng / 2)) prov_exp_specS[j] += data.Ls_poly / (ng - j);
  }

  // apply orientation error rate if required
  double oer;
  if (data.folded == false && opt.use_syno_orientation_error)
  {
    oer = param["syno_orientation_error"];
    for (int j = 1; j <= nc; j++)
    {
      exp_specS[j] = prov_exp_specS[j] * (1. - oer) + prov_exp_specS[ng - j] * oer;
    }
  }
  else
  {
    for (int j = 1; j <= nc; j++)
    {
      exp_specS[j] = prov_exp_specS[j];
    }
  }

  // calculate selected SFS expectations assuming theta=1 and rj's=1
  for (int j = 1; j <= nc; j++)
  {
    double provd = expected_allele_count(data.Ln_poly, 1., param, ng, j);
    exp_specN[j] = provd * (1. - opt.p_vdel);
  }

  // estimate ML theta and rj's
  double theta = (data.specS[1] + data.specN[1]) / (exp_specS[1] + exp_specN[1]);
  vector<double> r(nc + 1);
  r[1] = 1.;
  for (int j = 2; j <= nc; j++)
  {
    r[j] = (data.specS[j] + data.specN[j]) / (exp_specS[j] + exp_specN[j]);
    r[j] /= theta;
  }
  if (the) *the = theta;
  if (rj) *rj = r;

  // modify SFS expectations
  exp_specN[1] *= theta;
  exp_specS[1] *= theta;
  for (int j = 2; j <= nc; j++)
  {
    exp_specN[j] *= theta * r[j];
    exp_specS[j] *= theta * r[j];
  }

  // account for fixed proportion of very deleterious mutations
  // for(int j=2;j<=nc;j++)
  //  exp_specN[j]*=(1.-opt.p_vdel);

  /* 2. divergence */

  double exp_divS = 0.; double exp_divN = 0.;
  if (data.div_data_available)
  {
    if (opt.use_divergence_parameter)
    {
      exp_divS = data.fixS;
      exp_divN = data.fixN;
    }
    else
    {
      double prov_lambda = ((double)data.fixS / (double)data.Ls_div) - theta / ng;
      exp_divS = prov_lambda * data.Ls_div;
      exp_divN = (1 - opt.p_vdel) * expected_divergence(data.Ln_div, prov_lambda, theta, param, ng);

      double coeff = (data.fixS + data.fixN) / (exp_divS + exp_divN);
      exp_divS *= coeff;
      exp_divN *= coeff;
    }
  }

  /* 3. fill SFS_div object and return */

  struct SFS_div expected;

  expected.specS.push_back(-1.);
  for (int j = 1; j <= nc; j++)
  {
    expected.specS.push_back(exp_specS[j]);
  }

  expected.specN.push_back(-1.);
  for (int j = 1; j <= nc; j++)
  {
    expected.specN.push_back(exp_specN[j]);
  }

  if (data.div_data_available)
  {
    expected.divS = exp_divS;
    expected.divN = exp_divN;
  }

// if(debug){ printf("exp_specS: "); for(int j=1;j<=nc;j++) printf("%f ", exp_specS[j]); printf("\n");}
// if(debug){ printf("exp_specN: "); for(int j=1;j<=nc;j++) printf("%f ", exp_specN[j]); printf("\n");}

  return expected;
}


/* minus_loglikelihood */
/* SFS likelihood function, with theta and rj's automatically set at their ML value given */
/* Gmean, Gshape, posGmean, pos_prop and/or s0. */
/* Optimal theta and rj are returned if the last 2 args are not passed as NULL */

double minus_loglikelihood(struct MK_data data, map<string, double> param, double* the, vector<double>* rj)
{
  struct SFS_div expected = expected_SFS_div(data, param, the, rj);
  double mln = minus_SFS_loglikelihood(data, expected.specS, expected.specN, expected.divS, expected.divN);
  return mln;
}


/* bpp_SFS_minus_lnl */
/* bio++ likelihood function */

class bpp_SFS_minus_lnl : public virtual Function,
  public AbstractParametrizable
{
private:
  struct MK_data* data;
  double lnL;
  static constexpr double prec = 0.1;

public:
  // constructor: store data d in private object and create parameters as required by model m
  bpp_SFS_minus_lnl(struct MK_data* d, struct model* m) : AbstractParametrizable(""), data(d)
  {
    for (size_t i = 0; i < m->param_name.size(); i++)
    {
      if (m->optimized[m->param_name[i]] == false) continue;
      double val, cd, cu;
      shared_ptr<IntervalConstraint> ic;
      if (m->reparametrization[i] == "log")
      {
        cd = log(m->constraint[i][0]);
        cu = log(m->constraint[i][1]);
      }
      else
      {
        cd = m->constraint[i][0];
        cu = m->constraint[i][1];
      }
      val = (cd + cu) / 2.;
      ic.reset(new IntervalConstraint(cd, cu, true, true, prec));
      addParameter_(new Parameter(m->param_name[i], val, ic));
    }
  }

  bpp_SFS_minus_lnl* clone() const { return new bpp_SFS_minus_lnl(*this); }

public:
  // modify parameters; function matchParametersValue does call  fireParameterChanged
  void setParameters(const ParameterList& pl)
  throw (ParameterNotFoundException, ConstraintException, Exception)
  {
    matchParametersValues(pl);
  }

  // returns current likelihood
  double getValue() const throw (Exception) { return lnL; }

  // fireParameterChanged: called whenever a param is modified; calculates likelihood and store in lnL
  void fireParameterChanged(const ParameterList& pl)
  {
    map<string, double> p;
    vector<string> plnames = current_model.param_name;
    for (size_t i = 0; i < plnames.size(); i++)
    {
      if (current_model.optimized[plnames[i]] == false)
        p[plnames[i]] = current_model.fixed_value[plnames[i]];
      else
        p[plnames[i]] = getParameterValue(plnames[i]);
      if (current_model.reparametrization[i] == "log") p[plnames[i]] = exp(p[plnames[i]]);
    }

    lnL = minus_loglikelihood(*data, p, NULL, NULL);
  }
};


/* SFS_max_likelihood */
/* Likelihood optimizing function: creates bio++ likelihood function, bio++ reparametrizer and bio++ optimizer, */
/* initializes parameters, launches optimization, returns optimal values (through pl) */

double SFS_max_likelihood(bpp_SFS_minus_lnl lnL_fct, map<string, double> p_init, ParameterList* pl)
{
  ParameterList init = lnL_fct.getParameters();

  vector<string> plnames = current_model.param_name;
  for (size_t i = 0; i < plnames.size(); i++)
  {
    if (current_model.optimized[plnames[i]] == false)
      continue;
    double in;
    if (current_model.reparametrization[i] == "log")
      in = log(p_init[plnames[i]]);
    else
      in = p_init[plnames[i]];
    init.setParameterValue(plnames[i], in);
  }

  lnL_fct.setParameters(init);
  ReparametrizationFunctionWrapper rpf(&lnL_fct, false);
  ThreePointsNumericalDerivative tpnd(&rpf);
  tpnd.setParametersToDerivate(rpf.getParameters().getParameterNames());
  AbstractOptimizer* optim = 0;
  // optim = new BfgsMultiDimensions(&tpnd);
  optim = new PseudoNewtonOptimizer(&tpnd);
  optim->setVerbose(0);
  optim->setProfiler(0);
  optim->setMessageHandler(0);
  optim->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);

  FunctionStopCondition mystop(optim);
  mystop.setTolerance(OPTIM_TOLERANCE);
  optim->setStopCondition(mystop);

  optim->setMaximumNumberOfEvaluations(OPTIM_MAX_EVAL);

  double ml = 0;

  ParameterList init_transformed = rpf.getParameters();
  try
  {
    optim->init(init_transformed);
    ml = optim->optimize();
  }
  catch (Exception& e)
  {
    cout << "ERROR! " << e.what() << endl; return 1.;
  }

  *pl = lnL_fct.getParameters();
  if (optim->getNumberOfEvaluations() < OPTIM_MIN_EVAL)
    ml = -1.;

  delete optim;

  return -ml;
}


/* set_initial_values */

vector< map<string, double> > set_initial_values(struct model* m)
{
  vector< map<string, double> > v_init_p;

  // basic values
  map<string, double> basic;
  basic["negGshape"] = 0.3; basic["negGmean"] = 10000.;
  basic["posGshape"] = 1.; basic["posGmean"] = 100.;  basic["pos_prop"] = 0.025;
  basic["s0"] = -1.;
  basic["Bb"] = 3;  basic["Ba"] = 0.5;  basic["med_prop"] = 0.2;
  basic["BKm"] = 2; basic["BKnsz2"] = 100; basic["BKsigma"] = 0.5; basic["BKscale"] = 1000; basic["syno_orientation_error"] = 0.01;

  // chosen starting values
  if (current_DFEM_fit == "ExpoZero")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean"] = 100.;
    v_init_p.push_back(basic); v_init_p[1]["negGmean"] = 1000.;
    v_init_p.push_back(basic); v_init_p[2]["negGmean"] = 10000.;
    v_init_p.push_back(basic); v_init_p[3]["negGmean"] = 100000.;
  }
  else if (current_DFEM_fit == "ExpoExpo")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean"] = 100.; v_init_p[0]["posGmean"] = 1.;
    v_init_p.push_back(basic); v_init_p[1]["negGmean"] = 1000.; v_init_p[1]["posGmean"] = 10.;
    v_init_p.push_back(basic); v_init_p[2]["negGmean"] = 10000.; v_init_p[2]["posGmean"] = 100.;
    v_init_p.push_back(basic); v_init_p[3]["negGmean"] = 100000.; v_init_p[3]["posGmean"] = 1000.;
  }
  else if (current_DFEM_fit == "ExpoGamma")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean"] = 100.; v_init_p[0]["posGmean"] = 1.; v_init_p[0]["posGshape"] = 5.;
    v_init_p.push_back(basic); v_init_p[1]["negGmean"] = 1000.; v_init_p[1]["posGmean"] = 10.; v_init_p[1]["posGshape"] = 1.;
    v_init_p.push_back(basic); v_init_p[2]["negGmean"] = 10000.; v_init_p[2]["posGmean"] = 100.; v_init_p[2]["posGshape"] = 0.5;
    v_init_p.push_back(basic); v_init_p[3]["negGmean"] = 100000.; v_init_p[3]["posGmean"] = 1000.; v_init_p[3]["posGshape"] = 0.2;
  }
  else if (current_DFEM_fit == "GammaZero")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean"] = 100000000.; v_init_p[0]["negGshape"] = 0.1;
    v_init_p.push_back(basic); v_init_p[1]["negGmean"] = 100000.; v_init_p[1]["negGshape"] = 0.5;
    v_init_p.push_back(basic); v_init_p[2]["negGmean"] = 10000000.; v_init_p[2]["negGshape"] = 0.2;
  }
  else if (current_DFEM_fit == "GammaExpo")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean=1000."]; v_init_p[0]["negGshape"] = 1.; v_init_p[0]["posGmean"] = 1.; v_init_p[0]["pos_prop"] = 0.001;
    v_init_p.push_back(basic); v_init_p[1]["negGmean=100000."]; v_init_p[1]["negGshape"] = 0.5; v_init_p[1]["posGmean"] = 10.; v_init_p[1]["pos_prop"] = 0.001;
    v_init_p.push_back(basic); v_init_p[2]["negGmean=10000000."]; v_init_p[2]["negGshape"] = 0.2; v_init_p[2]["posGmean"] = 100.; v_init_p[2]["pos_prop"] = 0.001;
    v_init_p.push_back(basic); v_init_p[3]["negGmean=10000000."]; v_init_p[3]["negGshape"] = 0.2; v_init_p[3]["posGmean"] = 1000.; v_init_p[3]["pos_prop"] = 0.001;
  }
  else if (current_DFEM_fit == "GammaGamma")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean=1000."]; v_init_p[0]["negGshape"] = 1.; v_init_p[0]["posGmean"] = 1.; v_init_p[0]["posGshape"] = 5.;
    v_init_p.push_back(basic); v_init_p[1]["negGmean=100000."]; v_init_p[1]["negGshape"] = 0.5; v_init_p[1]["posGmean"] = 10.; v_init_p[1]["posGshape"] = 1.;
    v_init_p.push_back(basic); v_init_p[2]["negGmean=10000000."]; v_init_p[2]["negGshape"] = 0.2; v_init_p[2]["posGmean"] = 100.; v_init_p[2]["posGshape"] = 0.5;
    v_init_p.push_back(basic); v_init_p[3]["negGmean=10000000."]; v_init_p[3]["negGshape"] = 0.2; v_init_p[3]["posGmean"] = 1000.; v_init_p[3]["posGshape"] = 0.2;
  }
  else if (current_DFEM_fit == "DisplGamma")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean=1000."]; v_init_p[0]["negGshape"] = 1.; v_init_p[0]["s0"] = -2;
    v_init_p.push_back(basic); v_init_p[1]["negGmean=1000."]; v_init_p[1]["negGshape"] = 1.; v_init_p[1]["s0"] = -0.2;
    v_init_p.push_back(basic); v_init_p[2]["negGmean=100000."]; v_init_p[2]["negGshape"] = 0.5; v_init_p[2]["s0"] = -0.5;
    v_init_p.push_back(basic); v_init_p[3]["negGmean=10000000."]; v_init_p[3]["negGshape"] = 0.2; v_init_p[3]["s0"] = -0.1;
  }
  else if (current_DFEM_fit == "ScaledBeta")
  {
    v_init_p.push_back(basic); v_init_p[0]["Ba"] = 0.01; v_init_p[0]["Bb"] = 3.;
    v_init_p.push_back(basic); v_init_p[1]["Ba"] = 0.1; v_init_p[1]["Bb"] = 1.;
    v_init_p.push_back(basic); v_init_p[2]["Ba"] = 0.1; v_init_p[2]["Bb"] = 5.;
    v_init_p.push_back(basic); v_init_p[3]["Ba"] = 1.; v_init_p[3]["Bb"] = 10.;
  }
  else if (current_DFEM_fit == "FGMBesselK")
  {
    v_init_p.push_back(basic); v_init_p[0]["BKm"] = 1; v_init_p[0]["BKnsz2"] = 500.;
    v_init_p.push_back(basic); v_init_p[1]["BKm"] = 2; v_init_p[1]["BKnsz2"] = 500.;
    v_init_p.push_back(basic); v_init_p[2]["BKm"] = 1; v_init_p[2]["BKnsz2"] = 5000.;
    v_init_p.push_back(basic); v_init_p[3]["BKm"] = 2; v_init_p[3]["BKnsz2"] = 5000.;
  }
  else if (current_DFEM_fit == "ReflectedGamma")
  {
    v_init_p.push_back(basic); v_init_p[0]["negGmean"] = 100000000.; v_init_p[0]["negGshape"] = 0.1;
    v_init_p.push_back(basic); v_init_p[1]["negGmean"] = 100000.; v_init_p[1]["negGshape"] = 0.5;
    v_init_p.push_back(basic); v_init_p[2]["negGmean"] = 10000000.; v_init_p[2]["negGshape"] = 0.2;
  }

  // special case: fixed shape
  if (!m->optimized["negGshape"])
  {
    double adjusted_negGmean;
    if (m->fixed_value["negGshape"] < 0.12)
      adjusted_negGmean = 1000000.;
    else if (m->fixed_value["negGshape"] < 0.2)
      adjusted_negGmean = 100000.;
    else if (m->fixed_value["negGshape"] < 0.35)
      adjusted_negGmean = 10000.;
    else if (m->fixed_value["negGshape"] < 0.5)
      adjusted_negGmean = 1000.;
    else
      adjusted_negGmean = 100.;
    v_init_p.clear();
    v_init_p.push_back(basic);
    v_init_p[0]["negGmean"] = adjusted_negGmean;
  }


  // random starting values
  for (int i = 0; i < nb_random_start_values; i++)
  {
    v_init_p.push_back(random_starting_point(m));
  }

  return v_init_p;
}


/* optimize_DFEM */
/* Performs optimization: generate (sets of) starting values, launches optimizer, recovers optimum */

void optimize_DFEM(struct MK_data* data, struct model* m, struct parameter_point* optimum)
{
// cout <<"aa" <<endl;

  bpp_SFS_minus_lnl lnL_fct(data, m);

// cout <<"bb" <<endl;

  // starting points

  vector< map<string, double> > v_init_p = set_initial_values(m);
// cout <<v_init_p.size() <<" starting points" <<endl;

// cout <<"cc" <<endl;

  // optimization

  ParameterList pl, opt_pl;
  double lnL, maxlnL = -std::numeric_limits<double>::max();

  for (size_t i = 0; i < v_init_p.size(); i++)
  {
// cout <<"Starting values " <<i+1 <<" out of " <<v_init_p.size() <<endl;
    lnL = SFS_max_likelihood(lnL_fct, v_init_p[i], &pl);
// cout <<endl <<"lnL:" <<lnL <<endl;
    if (lnL > maxlnL && lnL < 0.)
    {
      maxlnL = lnL; opt_pl = pl;
    }
  }

// cout <<"dd" <<endl;

  if (maxlnL == -std::numeric_limits<double>::max())
    fin("Sorry: optimization failed; try more starting values\n");


  // recover optimal values

  vector<string> param_names = m->param_name;

  for (size_t i = 0; i < param_names.size(); i++)
  {
    if (m->optimized[param_names[i]] == false)
      optimum->DFEM_param[param_names[i]] = m->fixed_value[param_names[i]];
    else
      optimum->DFEM_param[param_names[i]] = opt_pl.getParameterValue(param_names[i]);
    if (m->reparametrization[i] == "log")
      optimum->DFEM_param[param_names[i]] = exp(optimum->DFEM_param[param_names[i]]);
  }

  // theta and ri's

  double opt_theta;
  vector<double> opt_r;
  minus_loglikelihood(*data, optimum->DFEM_param, &opt_theta, &opt_r);

  optimum->theta = opt_theta; optimum->ri = opt_r; optimum->lnL = maxlnL;
}


/* DFEM_mass */
/* returns the distribution mass between b1 and b2 for the current DFEM, parameters from param */
/* if b1>b2 then: if b2 is negative, return the mass between -inf and b2 */
/*                if b2 is positive, return the mass between b1 and +inf */

double DFEM_mass(double b1, double b2, struct parameter_point param)
{
// cout <<"DFEM_mass:" <<b1 <<" " <<b2 <<endl;


  if (b1 < 0. && b2 > 0.)
    return DFEM_mass(b1, 0., param) + DFEM_mass(0., b2, param);

  if (param.model_name == "ExpoZero")
  {
    double sh = 1.;
    double m = param.DFEM_param["negGmean"];
    double gsg = gsl_sf_gamma(sh);
    if (b1 > b2 && b2 < 0.)
      return gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg;
    else if (b1 > b2 && b2 >= 0.)
      return 0.;
    else if (b1 >= 0)
      return 0.;
    else
      return (gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg;
  }

  if (param.model_name == "ExpoExpo")
  {
    double sh = 1.;
    double m = param.DFEM_param["negGmean"];
    double psh = 1.;
    double pm = param.DFEM_param["posGmean"];
    double pp = param.DFEM_param["pos_prop"];
    double gsg_neg = gsl_sf_gamma(sh);
    double gsg_pos = gsl_sf_gamma(psh);
    if (b1 > b2 && b2 < 0.)
      return (1 - pp) * gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg_neg;
    else if (b1 > b2 && b2 >= 0.)
      return pp * gsl_sf_gamma_inc(psh, b1 * psh / pm) / gsg_pos;
    else if (b1 >= 0)
      return pp * (gsl_sf_gamma_inc(psh, b1 * psh / pm) - gsl_sf_gamma_inc(psh, b2 * psh / pm)) / gsg_pos;
    else
      return (1 - pp) * ((gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg_neg);
  }

  if (param.model_name == "ExpoGamma")
  {
    double sh = 1.;
    double m = param.DFEM_param["negGmean"];
    double psh = param.DFEM_param["posGshape"];
    double pm = param.DFEM_param["posGmean"];
    double pp = param.DFEM_param["pos_prop"];
    double gsg_neg = gsl_sf_gamma(sh);
    double gsg_pos = gsl_sf_gamma(psh);
    if (b1 > b2 && b2 < 0.)
      return (1 - pp) * gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg_neg;
    else if (b1 > b2 && b2 >= 0.)
      return pp * gsl_sf_gamma_inc(psh, b1 * psh / pm) / gsg_pos;
    else if (b1 >= 0)
      return pp * (gsl_sf_gamma_inc(psh, b1 * psh / pm) - gsl_sf_gamma_inc(psh, b2 * psh / pm)) / gsg_pos;
    else
      return (1 - pp) * ((gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg_neg);
  }

  if (param.model_name == "GammaZero")
  {
    double sh = param.DFEM_param["negGshape"];
    double m = param.DFEM_param["negGmean"];
    double gsg = gsl_sf_gamma(sh);
    if (b1 > b2 && b2 < 0.)
      return gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg;
    else if (b1 > b2 && b2 >= 0.)
      return 0.;
    else if (b1 >= 0)
      return 0.;
    else
      return (gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg;
  }

  if (param.model_name == "GammaExpo")
  {
    double sh = param.DFEM_param["negGshape"];
    double m = param.DFEM_param["negGmean"];
    double psh = 1.;
    double pm = param.DFEM_param["posGmean"];
    double pp = param.DFEM_param["pos_prop"];
    double gsg_neg = gsl_sf_gamma(sh);
    double gsg_pos = gsl_sf_gamma(psh);
    if (b1 > b2 && b2 < 0.)
      return (1 - pp) * gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg_neg;
    else if (b1 > b2 && b2 >= 0.)
      return pp * gsl_sf_gamma_inc(psh, b1 * psh / pm) / gsg_pos;
    else if (b1 >= 0)
      return pp * (gsl_sf_gamma_inc(psh, b1 * psh / pm) - gsl_sf_gamma_inc(psh, b2 * psh / pm)) / gsg_pos;
    else
      return (1 - pp) * ((gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg_neg);
  }

  if (param.model_name == "GammaGamma")
  {
    double sh = param.DFEM_param["negGshape"];
    double m = param.DFEM_param["negGmean"];
    double psh = param.DFEM_param["posGshape"];
    double pm = param.DFEM_param["posGmean"];
    double pp = param.DFEM_param["pos_prop"];
    double gsg_neg = gsl_sf_gamma(sh);
    double gsg_pos = gsl_sf_gamma(psh);
    if (b1 > b2 && b2 < 0.)
      return (1 - pp) * gsl_sf_gamma_inc(sh, -b2 * sh / m) / gsg_neg;
    else if (b1 > b2 && b2 >= 0.)
      return pp * gsl_sf_gamma_inc(psh, b1 * psh / pm) / gsg_pos;
    else if (b1 >= 0)
      return pp * (gsl_sf_gamma_inc(psh, b1 * psh / pm) - gsl_sf_gamma_inc(psh, b2 * psh / pm)) / gsg_pos;
    else
      return (1 - pp) * ((gsl_sf_gamma_inc(sh, -b2 * sh / m) - gsl_sf_gamma_inc(sh, -b1 * sh / m)) / gsg_neg);
  }

  if (param.model_name == "DisplGamma")
  {
    double s0 = param.DFEM_param["s0"];
    param.model_name = "GammaZero";
    return DFEM_mass(b1 + s0, b2 + s0, param);
  }

  if (param.model_name == "ScaledBeta")
  {
    double mp = param.DFEM_param["med_prop"];
    double Ba = param.DFEM_param["Ba"];
    double Bb = param.DFEM_param["Bb"];
    if (b1 > b2 && b1 <= SCALED_BETA_BOUNDARY)
      return 1 - mp;
    if (b1 > b2 && b2 >= SCALED_BETA_BOUNDARY)
      return 0.;
    if (b1 > b2 && b2 < 0.)
      return DFEM_mass(-SCALED_BETA_BOUNDARY, b2, param);
    if (b1 > b2 && b2 >= 0.)
      return DFEM_mass(b1, SCALED_BETA_BOUNDARY, param);
    if (b1 < -SCALED_BETA_BOUNDARY && b2 > SCALED_BETA_BOUNDARY)
      return 1.;
    if (b1 < -SCALED_BETA_BOUNDARY)
      return DFEM_mass(-SCALED_BETA_BOUNDARY, b2, param);
    if (b2 > SCALED_BETA_BOUNDARY)
      return DFEM_mass(b1, SCALED_BETA_BOUNDARY, param);
    ;

    double scb1 = (b1 / (2 * SCALED_BETA_BOUNDARY)) + 0.5;
    double scb2 = (b2 / (2 * SCALED_BETA_BOUNDARY)) + 0.5;
    if (b1 > b2 && b2 < 0.)
      return mp * gsl_cdf_beta_P(scb2, Ba, Bb);
    else if (b1 > b2 && b2 >= 0.)
      return mp * gsl_cdf_beta_Q(scb1, Ba, Bb);
    else
      return mp * (gsl_cdf_beta_P(scb2, Ba, Bb) - gsl_cdf_beta_P(scb1, Ba, Bb));
  }

  if (param.model_name == "FGMBesselK") // numerical discretization for this one
  {
    global_BKm = param.DFEM_param["BKm"];
    global_BKnsz2 = param.DFEM_param["BKnsz2"];
    global_BKsigma = param.DFEM_param["BKsigma"];
    global_BKscale = param.DFEM_param["BKscale"];

    if (b2 > b1)
    {
      double inte = 0.;
      int kk = global_integral_S_point.size();

      for (int i = 0; i < kk; i++)
      {
        if (global_integral_S_point[i] < b1)
          continue;
        if (global_integral_S_point[i] > b2)
          continue;

        inte += onearg_BesselK_density(global_integral_S_point[i]) * global_integral_S_width[i];
      }
      return inte;
    }
  }

  if (param.model_name == "ReflectedGamma")// numerical discretization for this one too
  {
    global_Gshape = param.DFEM_param["negGshape"];
    global_Gmean = param.DFEM_param["negGmean"];
    if (b2 > b1)
    {
      double inte = 0.;
      int kk = global_integral_S_point.size();

      for (int i = 0; i < kk; i++)
      {
        if (global_integral_S_point[i] < b1)
          continue;
        if (global_integral_S_point[i] > b2)
          continue;

        inte += onearg_reflected_Gamma_density(global_integral_S_point[i]) * global_integral_S_width[i];
      }
      return inte;
    }
  }

  return 0.;
}


/* calculate_alpha_basic */

vector<double> calculate_alpha_basic(struct MK_data* data)
{
  vector<double> v;

  if (data->div_data_available == false)
  {
    v.push_back(-1.); v.push_back(-1.); v.push_back(-1.); return v;
  }

  double tot_poly_N = 0., tot_poly_S = 0.;
  for (int j = 1; j <= global_nbcat; j++)
  {
    tot_poly_N += data->specN[j];
    tot_poly_S += data->specS[j];
  }

  double alpha = 1. - (double)(tot_poly_N * data->fixS) / (double)(tot_poly_S * data->fixN);
  double omegaA = (alpha * data->fixN * data->Ls_div) / (data->fixS * data->Ln_div);
  double omegaNA = ((1 - alpha) * data->fixN * data->Ls_div) / (data->fixS * data->Ln_div);

  v.push_back(alpha); v.push_back(omegaA); v.push_back(omegaNA);

  return v;
}


/* calculate_alpha_fww */

vector<double> calculate_alpha_fww(struct MK_data* data, double alfreq_thresh)
{
  vector<double> v;

  if (data->div_data_available == false)
  {
    v.push_back(-1.); v.push_back(-1.); v.push_back(-1.); return v;
  }

  double tot_poly_N = 0., tot_poly_S = 0.;
  for (int j = 1; j <= global_nbcat; j++)
  {
    if ((double)j / (double)global_n < alfreq_thresh)
      continue;
    tot_poly_N += data->specN[j];
    tot_poly_S += data->specS[j];
  }
  double alpha = 1. - (double)(tot_poly_N * data->fixS) / (double)(tot_poly_S * data->fixN);
  double omegaA = (alpha * data->fixN * data->Ls_div) / (data->fixS * data->Ln_div);
  double omegaNA = ((1 - alpha) * data->fixN * data->Ls_div) / (data->fixS * data->Ln_div);

  v.push_back(alpha); v.push_back(omegaA);  v.push_back(omegaNA);

  return v;
}


/* calculate_alpha_CI */

void calculate_alpha_CI(struct MK_data* data, struct parameter_point* optimal, double threshold)
{
  optimal->alpha_down = optimal->alpha_up = 0.;

  double increment = 0.001;

  if (opt.use_divergence_parameter)
  {
    double opt_lnL = optimal->lnL;
    double opt_lambda = ((double)data->fixS / (double)data->Ls_div) - optimal->theta / data->nb_gene;
    struct SFS_div opt_exp = expected_SFS_div(*data, optimal->DFEM_param, NULL, NULL);
    double opt_divN = opt_exp.divN;
    double exp_div_N3 = expected_divergence_truncated(data->Ln_div, opt_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene, opt.nearly_neutral_threshold);

    // down boundary
    double prev_diff = 0.;
    double cur_divN = opt_divN;
    bool problem = false;
    while (1)
    {
      cur_divN *= (1. - increment);
      double cur_lnL = -minus_SFS_loglikelihood(*data, opt_exp.specS, opt_exp.specN, opt_exp.divS, cur_divN);
      double cur_diff = opt_lnL - cur_lnL;
      if (cur_diff >= threshold && prev_diff <= threshold)
        break;
      if (cur_divN < opt_divN / 5.)
      {
        problem = true; break;
      }
      prev_diff = cur_diff;
    }
    if (!problem)
      optimal->alpha_down = (cur_divN - exp_div_N3) / data->fixN;

    // up boundary
    prev_diff = 0.;
    cur_divN = opt_divN;
    while (1)
    {
      cur_divN *= (1. + increment);
      double cur_lnL = -minus_SFS_loglikelihood(*data, opt_exp.specS, opt_exp.specN, opt_exp.divS, cur_divN);
      double cur_diff = opt_lnL - cur_lnL;
      if (cur_diff >= threshold && prev_diff <= threshold)
        break;
      if (cur_divN > opt_divN * 3.)
      {
        problem = true; break;
      }
      prev_diff = cur_diff;
    }
    if (!problem)
      optimal->alpha_up = (cur_divN - exp_div_N3) / data->fixN;
  }

  optimal->omegaA_down = optimal->omegaA * (optimal->alpha_down / optimal->alpha);
  optimal->omegaA_up = optimal->omegaA * (optimal->alpha_up / optimal->alpha);

  optimal->omegaNA_up = (optimal->omegaA + optimal->omegaNA) * (1 - optimal->alpha_down);
  optimal->omegaNA_down = (optimal->omegaA + optimal->omegaNA) * (1 - optimal->alpha_up);
}


/* calculate_alpha_DFEM */
/* Calculates alpha, omegaA and discretized distributino of fitness effect of non-syno substitutions */
/* for a given model and a given data set using MLE of parameters as passed in param optimal */

double calculate_alpha_DFEM(struct MK_data* data, struct parameter_point* optimal)
{
  double totN;
  double opt_lambda;

  if (data->div_data_available == false)
  {
    data->Ln_div = ARBITRARY_LENGTH;
    opt_lambda = ARBITRARY_LAMBDA;
    totN = expected_divergence(data->Ln_div, opt_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene);
  }
  else if (opt.use_divergence_parameter)
  {
    totN = data->fixN;
    opt_lambda = ((double)data->fixS / (double)data->Ls_div) - optimal->theta / data->nb_gene;
  }
  else
  {
    double prov_lambda = ((double)data->fixS / (double)data->Ls_div) - optimal->theta / data->nb_gene;
    double prov_totN = expected_divergence(data->Ln_div, prov_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene);
    double coeff = (data->fixS + data->fixN) / (prov_lambda * data->Ls_div + prov_totN);
    opt_lambda = prov_lambda * coeff;
    totN = prov_totN * coeff;
  }


  double exp_div_N1 = expected_divergence_truncated(data->Ln_div, opt_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene, -opt.nearly_neutral_threshold);
  double exp_div_N2 = expected_divergence_truncated(data->Ln_div, opt_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene, 0);
  double exp_div_N3 = expected_divergence_truncated(data->Ln_div, opt_lambda, optimal->theta, optimal->DFEM_param, data->nb_gene, opt.nearly_neutral_threshold);


  double prop_fix_sdel = exp_div_N1 / totN;
  double prop_fix_wdel = (exp_div_N2 - exp_div_N1) / totN;
  double prop_fix_wadv = (exp_div_N3 - exp_div_N2) / totN;
  double prop_fix_sadv = (totN - exp_div_N3) / totN;

  optimal->alpha = prop_fix_sadv;
  optimal->omegaA = (prop_fix_sadv * totN * data->Ls_div) / (data->fixS * data->Ln_div);

  optimal->omegaNA = ((1. - prop_fix_sadv) * totN * data->Ls_div) / (data->fixS * data->Ln_div);
  optimal->discr_DFE_fixed.push_back(prop_fix_sdel);
  optimal->discr_DFE_fixed.push_back(prop_fix_wdel);
  optimal->discr_DFE_fixed.push_back(prop_fix_wadv);
  optimal->discr_DFE_fixed.push_back(prop_fix_sadv);

  return optimal->alpha;
}


/* console_output */

void console_output(vector<struct parameter_point> v_optimum, vector<struct model> v_model, double sat_lnL, double neut_lnL, double neut_oer, double neut_neut, vector<double> alpha_basic, vector<double> alpha_FWW, int nbcat)
{
  // likelihood
  printf("\nmax likelihood:\n");
// cout <<"nbcat:" <<nbcat <<endl;
  int nbparam_neut = nbcat + 2;
  if (opt.use_divergence_parameter)
    nbparam_neut++;
  if (opt.use_syno_orientation_error)
    nbparam_neut++;
  printf("  Neutral:\t%.2f\t[%d parameters]\n", neut_lnL, nbparam_neut);
  for (size_t i = 0; i < v_optimum.size(); i++)
  {
    int nbparam = v_model[i].param_name.size() + nbcat + 1;
    if (opt.use_divergence_parameter)
      nbparam++;
    printf("  %s:\t%.2f\t[%d parameters]\n", v_optimum[i].model_name.c_str(), v_optimum[i].lnL, nbparam);
  }
  int nbparam_sat = 2 * nbcat + 2;
  printf("  Saturated:\t%.2f\t[%d parameters]\n", sat_lnL, nbparam_sat);

  // DFEM
  cout << endl << "Estimated DFEM:" << endl;
  cout << "  " << "Neutral:" << endl;
  if (opt.use_syno_orientation_error)
    cout << "    syno_orientation_error:" << neut_oer << endl;
  cout << "    [-inf]:" << 1. - neut_neut << " " << "0:" << neut_neut << endl;
  for (size_t i = 0; i < v_optimum.size(); i++)
  {
    cout << "  " << v_optimum[i].model_name << ":" << endl;
    std::map<string, double>::iterator it;
    cout << "    ";
    for (it = v_optimum[i].DFEM_param.begin(); it != v_optimum[i].DFEM_param.end(); it++)
    {
      cout << it->first << ":" << it->second << " ";
    }
    cout << endl;

    int nb_bound = v_optimum[i].discr_bound.size();
    cout << "    [-inf," << v_optimum[i].discr_bound[0] << "]:" << v_optimum[i].discr_mass[0] << "  ";
    for (int j = 1; j < nb_bound; j++)
    {
      cout << "[" << v_optimum[i].discr_bound[j - 1] << "," << v_optimum[i].discr_bound[j] << "]:" << v_optimum[i].discr_mass[j] << "  ";
      if (v_optimum[i].discr_bound[j] == 0. || (j > 1 && v_optimum[i].discr_bound[j - 1] > 0 && v_optimum[i].discr_bound[j - 2] < 0))
        cout << endl << "    ";
    }
    cout << "[" << v_optimum[i].discr_bound[nb_bound - 1] << ",+inf]:" << v_optimum[i].discr_mass[nb_bound] << endl;
  }

  // alpha
  cout << endl << "Estimated alpha:" << endl;
  cout << "  basic:\t" << alpha_basic[0] << endl;
  cout << "  FWW[0.15]:\t" << alpha_FWW[0] << endl;
  for (size_t i = 0; i < v_optimum.size(); i++)
  {
    cout << "  " << v_optimum[i].model_name << ":\t" << v_optimum[i].alpha /* <<" [" <<v_optimum[i].alpha_down <<"," <<v_optimum[i].alpha_up <<"]"*/ << endl;
  }

  // omegaA
  cout << endl << "Estimated omegaA:" << endl;
  cout << "  basic:\t" << alpha_basic[1] << endl;
  cout << "  FWW[0.15]:\t" << alpha_FWW[1] << endl;
  for (size_t i = 0; i < v_optimum.size(); i++)
  {
    cout << "  " << v_optimum[i].model_name << ":\t" << v_optimum[i].omegaA /* <<" [" <<v_optimum[i].omegaA_down <<"," <<v_optimum[i].omegaA_up <<"]"*/ << endl;
  }

  // omegaA
  cout << endl << "Estimated omegaNA:" << endl;
  cout << "  basic:\t" << alpha_basic[2] << endl;
  cout << "  FWW[0.15]:\t" << alpha_FWW[2] << endl;
  for (size_t i = 0; i < v_optimum.size(); i++)
  {
    cout << "  " << v_optimum[i].model_name << ":\t" << v_optimum[i].omegaNA /* <<" [" <<v_optimum[i].omegaNA_down <<"," <<v_optimum[i].omegaNA_up <<"]" */ << endl;
  }
}


/* file_output */

void file_output(struct MK_data* data, vector<struct parameter_point> v_optimum, vector<struct model> v_model, double sat_lnL, double neut_lnL, double neut_oer, double neut_neut, vector<double> alpha_basic, vector<double> alpha_FWW, FILE* out)
{
  vector<struct model> mod;
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    mod.push_back(create_model(implemented_model_names[i]));
  }

  int nbcat = data->nb_SFScat;

  struct SFS_div observed;
  for (int i = 0; i <= nbcat; i++)
  {
    observed.specS.push_back(data->specS[i]);
  }
  for (int i = 0; i <= nbcat; i++)
  {
    observed.specN.push_back(data->specN[i]);
  }
  observed.divS = data->fixS;
  observed.divN = data->fixN;

  double totpS, totpN;
  totpS = totpN = 0;
  for (int i = 1; i <= nbcat; i++)
  {
    totpS += observed.specS[i];
    totpN += observed.specN[i];
  }

  struct SFS_div neutral_expected, expected;


  // header
  fprintf(out, "dataset,pn,Lpn,ps,Lps,dn,Ldn,ds,Lds,model,#param,lnL,alpha,alpha_down,alpha_up,omegaA,omegaA_down,omegaA_up,");
  fprintf(out, "prop_subst_sdel,prop_subst_wdel,prop_subst_wadv,prop_subst_sadv,Neutral:f0,");
  if (opt.use_syno_orientation_error)
    fprintf(out, "Neutral:syno_orientation_error,");
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    string model_name = implemented_model_names[i];
    for (size_t j = 0; j < mod[i].param_name.size(); j++)
    {
      fprintf(out, "%s:%s,", model_name.c_str(), mod[i].param_name[j].c_str());
    }
  }
  fprintf(out, "theta,");
  for (int j = 2; j <= nbcat; j++)
  {
    fprintf(out, "r%d,", j);
  }
  for (int j = 1; j <= nbcat; j++)
  {
    fprintf(out, "obs_specS%d,", j);
  }
  for (int j = 1; j <= nbcat; j++)
  {
    fprintf(out, "exp_specS%d,", j);
  }
  for (int j = 1; j <= nbcat; j++)
  {
    fprintf(out, "obs_specN%d,", j);
  }
  for (int j = 1; j <= nbcat; j++)
  {
    fprintf(out, "exp_specN%d,", j);
  }
  fprintf(out, "obs_dS,exp_dS,obs_dN,exp_dN\n");


  // neutral
  int nbparam_neut = nbcat + 2;
  if (opt.use_divergence_parameter)
    nbparam_neut++;
  if (opt.use_syno_orientation_error)
    nbparam_neut++;
  optimize_neutral(data, NULL, NULL, &neutral_expected);
  fprintf(out, "%s,", data->name);
  fprintf(out, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,", totpN, data->Ln_poly, totpS, data->Ls_poly, data->fixN, data->Ln_div, data->fixS, data->Ls_div);
  fprintf(out, "Neutral,%d,%.4f,%.4f,NA, NA,%.4f,NA,NA,NA,NA,NA,NA,%f,", nbparam_neut, neut_lnL, alpha_basic[0], alpha_basic[1], neut_neut);
  if (opt.use_syno_orientation_error)
    fprintf(out, "%.4f,", neut_oer);
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    string model_name = implemented_model_names[i];
    for (size_t j = 0; j < mod[i].param_name.size(); j++)
    {
      fprintf(out, "NA,");
    }
  }
  fprintf(out, "NA,");
  for (int j = 2; j <= nbcat; j++)
  {
    fprintf(out, "NA,");
  }

  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specS[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", neutral_expected.specS[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specN[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", neutral_expected.specN[j]);
  }
  fprintf(out, "%.4f,%.4f,%.4f,%.4f\n", observed.divS, neutral_expected.divS, observed.divN, neutral_expected.divN);


  // other
  for (size_t k = 0; k < v_optimum.size(); k++)
  {
    int nbparam = v_model[k].param_name.size() + nbcat + 1;
    if (opt.use_divergence_parameter)
      nbparam++;
    current_DFEM_fit = v_model[k].name;
    expected = expected_SFS_div(*data, v_optimum[k].DFEM_param, NULL, NULL);
    fprintf(out, "%s,", data->name);
    fprintf(out, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,", totpN, data->Ln_poly, totpS, data->Ls_poly, data->fixN, data->Ln_div, data->fixS, data->Ls_div);
    fprintf(out, "%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,", v_optimum[k].model_name.c_str(), nbparam, v_optimum[k].lnL, v_optimum[k].alpha, v_optimum[k].alpha_down, v_optimum[k].alpha_up, v_optimum[k].omegaA, v_optimum[k].omegaA_down, v_optimum[k].omegaA_up, v_optimum[k].discr_DFE_fixed[0], v_optimum[k].discr_DFE_fixed[1], v_optimum[k].discr_DFE_fixed[2], v_optimum[k].discr_DFE_fixed[3]);
    fprintf(out, "NA,");
    if (opt.use_syno_orientation_error)
      fprintf(out, "NA,");
    for (size_t i = 0; i < implemented_model_names.size(); i++)
    {
      string model_name = implemented_model_names[i];
      for (size_t j = 0; j < mod[i].param_name.size(); j++)
      {
        if (model_name == v_model[k].name)
          fprintf(out, "%f,", v_optimum[k].DFEM_param[mod[i].param_name[j]]);
        else
          fprintf(out, "NA,");
      }
    }
    fprintf(out, "%f,", v_optimum[k].theta);
    for (int j = 2; j <= nbcat; j++)
    {
      fprintf(out, "%f,", v_optimum[k].ri[j]);
    }

    for (int j = 1; j < nbcat; j++)
    {
      fprintf(out, "%.4f,", observed.specS[j]);
    }
    for (int j = 1; j < nbcat; j++)
    {
      fprintf(out, "%.4f,", expected.specS[j]);
    }
    for (int j = 1; j < nbcat; j++)
    {
      fprintf(out, "%.4f,", observed.specN[j]);
    }
    for (int j = 1; j < nbcat; j++)
    {
      fprintf(out, "%.4f,", expected.specN[j]);
    }
    fprintf(out, "%.4f,%.4f,%.4f,%.4f\n", observed.divS, expected.divS, observed.divN, expected.divN);
  }

  // saturated
  int nbparam_sat = 2 * nbcat + 2;
  fprintf(out, "%s,", data->name);
  fprintf(out, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,", totpN, data->Ln_poly, totpS, data->Ls_poly, data->fixN, data->Ln_div, data->fixS, data->Ls_div);
  fprintf(out, "Saturated,%d,%.4f,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,", nbparam_sat, sat_lnL);
  fprintf(out, "NA,");
  if (opt.use_syno_orientation_error)
    fprintf(out, "NA,");
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    string model_name = implemented_model_names[i];
    for (size_t j = 0; j < mod[i].param_name.size(); j++)
    {
      fprintf(out, "NA,");
    }
  }
  fprintf(out, "NA,");
  for (int j = 2; j <= nbcat; j++)
  {
    fprintf(out, "NA,");
  }

  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specS[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specS[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specN[j]);
  }
  for (int j = 1; j < nbcat; j++)
  {
    fprintf(out, "%.4f,", observed.specN[j]);
  }
  fprintf(out, "%.4f,%.4f,%.4f,%.4f\n", observed.divS, observed.divS, observed.divN, observed.divN);
}


/* fold */
/* Transform an unfolded sfs into a folded one (double*) */
/* Unfolded sfs has size n (i.e., elements from 1 to n-1) */
/* Folded version returned in the same array */
/* Returns the size of folded SFS */

int fold(double* sfs, int n)
{
  if (n % 2 == 1)
  {
    for (int i = 1; i <= n / 2; i++)
    {
      sfs[i] += sfs[n - i];
    }
    for (int i = n / 2 + 1; i < n; i++)
    {
      sfs[i] = -1.;
    }
  }
  else
  {
    for (int i = 1; i < n / 2; i++)
    {
      sfs[i] += sfs[n - i];
    }
    for (int i = n / 2 + 1; i < n; i++)
    {
      sfs[i] = -1.;
    }
  }

  return n / 2;
}


/* read_dofe */
/* Reads dofe input file and return a data set structure */
/* Two formats covered: 2006 (poly only) and 2009 (poly+ div) */

struct MK_data* read_dofe(FILE* in, int* nb_dataset, bool use_div_data)
{
  char line[MAXLLINE + 1], * tok;
  int ret, nbclass;
  bool folded;


// pass 1 : folded/unfolded, count loci
  folded = true;
  char* cret = fgets(line, MAXLLINE, in);
  if (cret == NULL) return NULL;
  int nbl = 0;
  while (fgets(line, MAXLLINE, in))
  {
    if (line[0] == '\0') break;
    if (strncmp(line, "#unfolded", 7) == 0) {folded = false; continue;}
    if (line[0] != '\n') nbl++;
  }


  MK_data* dataset = (MK_data*)calloc(nbl, sizeof(struct MK_data));
  if (dataset == NULL) {printf("Not enough memory\n"); exit(0);}
  rewind(in);


// pass 2 : read data
  cret = fgets(line, MAXLLINE, in);
  if (cret == NULL) return NULL;
  nbl = 0;
  while (fgets(line, MAXLLINE, in))
  {
    if (line[0] == '\n') continue;
    if (line[0] == '#') continue;
    // read locus name
    sprintf(dataset[nbl].name, "%s", strtok(line, "\t\n"));
    dataset[nbl].folded = folded;
    // read nb genes, deduce nb categories
    ret = sscanf(strtok(NULL, "\t\n"), "%d", &(dataset[nbl].nb_gene));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    if (folded) dataset[nbl].nb_SFScat = dataset[nbl].nb_gene / 2;
    else dataset[nbl].nb_SFScat = dataset[nbl].nb_gene - 1;
    // calculate nb categories and allocate spectra
    nbclass = dataset[nbl].nb_SFScat;
    dataset[nbl].specS = (double*)calloc(nbclass + 1, sizeof(double));
    dataset[nbl].specN = (double*)calloc(nbclass + 1, sizeof(double));
    // read Ln_poly
    ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].Ln_poly));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    // read N spectrum
    for (int i = 1; i <= nbclass; i++)
    {
      ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].specN[i]));
      if (ret != 1) {printf("Bad input file\n"); exit(0);}
    }
    // read Ls_poly
    ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].Ls_poly));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    // read S spectrum
    for (int i = 1; i <= nbclass; i++)
    {
      ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].specS[i]));
      if (ret != 1) {printf("Bad input file\n"); exit(0);}
    }
    // check presence div data
    tok = strtok(NULL, "\t\n");
    if (tok == NULL || !use_div_data)
    {
      dataset[nbl].Ln_div = dataset[nbl].Ls_div = dataset[nbl].fixN = dataset[nbl].fixS = -1.;
      dataset[nbl].div_data_available = false;
      nbl++;
      continue;
    }
    dataset[nbl].div_data_available = true;
    // read Ln_div
    ret = sscanf(tok, "%le", &(dataset[nbl].Ln_div));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    // read fixN
    ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].fixN));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    // read Ls_div
    ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].Ls_div));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    // read fixS
    ret = sscanf(strtok(NULL, "\t\n"), "%le", &(dataset[nbl].fixS));
    if (ret != 1) {printf("Bad input file\n"); exit(0);}
    nbl++;
  }

  *nb_dataset = nbl;


  return dataset;
}


/* usage */
void usage()
{
  cout << endl << "usage: grapes -in input_file -out output_file -model model_name  [options] " << endl << endl;

  cout << "options:" << endl;
  // cout <<"[-nearly_neutral float]           maximal \"nearly-neutral\" Ns            default=5." <<endl;
  // cout <<"[-FWW_threshold float]            min allele freq in FWW alpha           default=0.15" <<endl;
  // cout <<"[-multinomial]               multinomial likelihood calculation     default=Poisson" <<endl;
  // cout <<"[-no_div_data]                    don't use divergence data at all       default=false" <<endl;
  // cout <<"[-no_div_param]                   use div. data only to estimate DFE     default=false" <<endl;
  cout << "[-fold]                           fold SFS                               default=false" << endl;
  cout << "[-no_syn_orient_error]            equal syn/non-syn misorientation       default=false" << endl;
  // cout <<"[-anc_to_rec_Ne_ratio float]      divergence/polymorphism Ne ratio       default=1." <<endl;
  cout << "[-nb_rand_start int]              nb optimization random start values    default=0" << endl;
  cout << "[-no_separate]                    do not perform separate optimization   default=false" << endl;
  cout << "[-no_shared]                      do not perform shared optimization     default=false" << endl;
  cout << "[-shared_shape float,float...]    use fixed shape(s)                     default=false" << endl;
  cout << "[-p_lethal float]                 prop. of lethal non-syn mutations      default=0." << endl;
  // cout <<"[-fixed_param string]             fixed (non-optimized) param file name  default=none" <<endl <<endl;

  cout << "implemented models: GammaZero ReflectedGamma" << endl;

  exit(0);
}


/* read_fixed_param_file */

void read_fixed_param_file(char* fname, struct options* opt)
{
  FILE* f = fopen(fname, "r");
  if (f == NULL)
  {
    printf("Can't read file: %s\n", fname); exit(0);
  }
  char line[MAXLLINE + 1];

  while (fgets(line, MAXLLINE, f))
  {
    char* s1, * s2, * s3;
    if (line[0] == '\n')
      continue;
    if (line[0] == '#')
      continue;
    s1 = strtok(line, ",");
    s2 = strtok(NULL, ",");
    s3 = strtok(NULL, ",");
    if (s1 == NULL || s2 == NULL || s3 == NULL)
      continue;
    double val;
    if (sscanf(s3, "%le", &val) != 1)
      continue;
    opt->fixed_p.value.push_back(val);
    string st1(s1);
    opt->fixed_p.model_name.push_back(st1);
    string st2(s2);
    opt->fixed_p.param_name.push_back(st2);
    opt->fixed_p.nb++;
  }
}


/* read_fixed_shapes */

vector<double> read_fixed_shapes(char* s)
{
  char* prov;
  double fi;
  vector<double> vfi;
  int ret;

  prov = strtok(s, ",");
  if (!prov)
    fin("error -shared_shape option: comma-separated floats expected");

  while (prov)
  {
    ret = sscanf(prov, "%le", &fi);
    if (ret != 1)
      fin("error -shared_shape option: comma-separated floats expected");
    vfi.push_back(fi);
    prov = strtok(NULL, ",");
  }

  return vfi;
}


/* set_default_options */

void set_default_options(struct options* opt)
{
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    opt->do_model.push_back(false);
  }

  opt->use_divergence_data = false;
  opt->use_divergence_parameter = true;
  opt->use_syno_orientation_error = true;
  opt->use_poisson = true;
  opt->nb_random_start = 0;
  opt->nearly_neutral_threshold = 5.;
  opt->FWW_threshold = 0.15;
  opt->fixed_p.nb = 0;
  opt->do_separate = true;
  opt->do_shared = true;
  opt->fold = false;
  opt->p_vdel = 0.;
  // opt->fixed_shared_negGshape=-1.;
}


/* get_options */

void  get_options(int argc, char** argv, struct options* opt, string* infile_name, string* outfile_name)
{
  string argument, argument1;

  for (int i = 1; i < argc; i++)
  {
    argument = argv[i];

    if (argument == "-in" && i != argc - 1)
      (*infile_name) = argv[i + 1];
    if (argument == "-out" && i != argc - 1)
      (*outfile_name) = argv[i + 1];

    if (argument == "-model" && i != argc - 1)
    {
      argument1 = argv[i + 1];
      for (size_t j = 0; j < implemented_model_names.size(); j++)
      {
        if (argument1 == implemented_model_names[j] || argument1 == "all")
          opt->do_model[j] = true;
      }
    }

    if (argument == "-no_div_data")
      opt->use_divergence_data = false;
    if (argument == "-no_div_param")
      opt->use_divergence_parameter = false;
    if (argument == "-no_syn_orient_error")
      opt->use_syno_orientation_error = false;
    if (argument == "-multinomial")
      opt->use_poisson = false;
    if (argument == "-nearly_neutral_threshold" && i != argc - 1)
      sscanf(argv[i + 1], "%le", &(opt->nearly_neutral_threshold));
    if (argument == "-FWW_threshold" && i != argc - 1)
      sscanf(argv[i + 1], "%le", &(opt->FWW_threshold));
    if (argument == "-nb_rand_start" && i != argc - 1)
      sscanf(argv[i + 1], "%d", &(opt->nb_random_start));
    if (argument == "-anc_to_rec_Ne_ratio" && i != argc - 1)
      sscanf(argv[i + 1], "%le", &ancient_to_recent_Ne_ratio);
    if (argument == "-no_separate")
      opt->do_separate = false;
    if (argument == "-no_shared")
      opt->do_shared = false;
    if (argument == "-shared_shape" && i != argc - 1)
      opt->fixed_shared_negGshape = read_fixed_shapes(argv[i + 1]);
    if (argument == "-fold")
      opt->fold = true;
    if (argument == "-p_lethal" && i != argc - 1)
      sscanf(argv[i + 1], "%le", &opt->p_vdel);

//    if(argument=="-fixed_param" && i!=argc-1){
//      read_fixed_param_file(argv[i+1], opt);
//    }
// no fixed param in multi_grapes
  }
}


/* chek_options */

void check_options(struct options* opt, string infile_name, string outfile_name)
{
  string err_message = "";
  string warn_message = "";

  if (infile_name.size() == 0)
    err_message += "Error: missing input file name\n";


  if (outfile_name.size() == 0)
    warn_message += "Warning: missing output file name - no file will be written\n";


  bool do_some_model = false;
  int nb_model = 0;

  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    if (opt->do_model[i])
    {
      do_some_model = true; nb_model++; opt->model_name = implemented_model_names[i];
    }
  }

  if (do_some_model == false)
    err_message += "Error: no selected model\n";

  if (nb_model > 1)
    err_message += "Error: a single model should be passed\n";

  if (opt->FWW_threshold < 0.)
    err_message += "Error: negative FWW_threshold\n";

  if (ancient_to_recent_Ne_ratio <= 0.)
    err_message += "Error: non-positive anc_to_rec_Ne_ratio\n";

  if (opt->use_divergence_data == false)
    opt->use_divergence_parameter = false;

  if (!(opt->do_separate || opt->do_shared))
    err_message += "Error: neither separate nor shared analysis requested\n";

  if (err_message.size() > 0)
  {
    cout << err_message;
    usage();
  }

  if (warn_message.size() > 0)
  {
    cout << warn_message;
  }

  return;
}


/* one_species_max_likelihood */
/* Calculates the likelihood of data set dataset for model model */
/* Global variable opt needs to be set; it defines many aspects of */
/* the calculation: model, fixed params, etc... */

struct parameter_point one_species_max_likelihood(MK_data* dataset)
{
  global_n = dataset->nb_gene;
  global_L = dataset->Ln_poly;
  global_nbcat = dataset->nb_SFScat;

  // create model list
  vector<struct model> model_list;
  for (size_t i = 0; i < implemented_model_names.size(); i++)
  {
    if (opt.do_model[i])
      model_list.push_back(create_model(implemented_model_names[i]));
  }

  // manage fixed parameters
  manage_fixed_params(&model_list, opt.fixed_p);

  // optimize likelihood
  vector<struct parameter_point> v_optimal;
  for (size_t i = 0; i < model_list.size(); i++)
  {
    current_model = model_list[i];
    current_DFEM_fit = model_list[i].name;
    struct parameter_point opt_prov = create_parameter_point(model_list[i].name);
    optimize_DFEM(dataset, &model_list[i], &opt_prov);
    v_optimal.push_back(opt_prov);
  }


  return v_optimal[0];
}


/* median_negGshape */

double median_negGshape(vector<struct parameter_point> v_opt)
{
  double min;
  unsigned int rg_min;
  vector<double> basic, sorted;

  for (size_t i = 0; i < v_opt.size(); i++)
  {
    basic.push_back(v_opt[i].DFEM_param["negGshape"]);
  }

  for (size_t i = 0; i < basic.size(); i++)
  {
    min = 99.;
    rg_min = -1;
    for (size_t j = 0; j < basic.size(); j++)
    {
      if (basic[j] < min)
      {
        min = basic[j];
        rg_min = i;
      }
    }
    sorted.push_back(min);
    basic[rg_min] = 100.;
  }

  int ns2 = basic.size() / 2;
  if (basic.size() % 2 == 1)
    return sorted[ns2];
  else
    return (sorted[ns2 - 1] + sorted[ns2]) / 2.;
}


/* reset_opt_fixed_p */
void reset_opt_fixed_p()
{
  opt.fixed_p.nb = 0;
  opt.fixed_p.model_name.clear();
  opt.fixed_p.param_name.clear();
  opt.fixed_p.value.clear();
}


/* multi_max_likelihood */

vector<struct parameter_point> multi_max_likelihood(struct MK_data* dataset, int nb_dataset, struct shared_parameter_point shared_pp, double* multi_max_lnL)
{
  // set shared param as fixed
  reset_opt_fixed_p();
  if (shared_pp.negGshape > 0.)
  {
    opt.fixed_p.model_name.push_back(opt.model_name);
    opt.fixed_p.param_name.push_back("negGshape");
    opt.fixed_p.value.push_back(shared_pp.negGshape);
    opt.fixed_p.nb++;
  }
  if (shared_pp.posGshape > 0.)
  {
    opt.fixed_p.model_name.push_back(opt.model_name);
    opt.fixed_p.param_name.push_back("posGshape");
    opt.fixed_p.value.push_back(shared_pp.posGshape);
    opt.fixed_p.nb++;
  }
  if (shared_pp.pos_prop > 0.)
  {
    opt.fixed_p.model_name.push_back(opt.model_name);
    opt.fixed_p.param_name.push_back("pos_prop");
    opt.fixed_p.value.push_back(shared_pp.pos_prop);
    opt.fixed_p.nb++;
  }

  // call one_species likelihood
  vector<struct parameter_point> v_opt;
  for (int k = 0; k < nb_dataset; k++)
  {
    global_current_species = k;
    v_opt.push_back(one_species_max_likelihood(dataset + k));
  }

  // calculate joint likelihood
  double mml = 0.;
  for (int k = 0; k < nb_dataset; k++)
  {
    mml += v_opt[k].lnL;
  }
  *multi_max_lnL = mml;

// cout <<shared_pp.negGshape <<":" <<mml <<endl;

  return v_opt;
}


/* bpp_SFS_minus_lnl_shared */
/* bio++ likelihood function, multiple species, shared_parameters */

class bpp_SFS_minus_lnl_shared : public virtual Function,
  public AbstractParametrizable
{
private:
  struct MK_data* data;
  int nb_dataset;
  double lnL;
  static constexpr double prec = 0.1;

public:
  bpp_SFS_minus_lnl_shared(struct MK_data* d, int nbd, struct shared_parameter_point shared_pp) : AbstractParametrizable(""), data(d), nb_dataset(nbd)
  {
    if (shared_pp.negGshape <= 0. && shared_pp.posGshape <= 0.  && shared_pp.pos_prop < 0.)
    {
      cout << "no shared param to optimize" << endl;
      exit(0);
    }

    if (shared_pp.negGshape > 0.)
    {
      double val, cd, cu;
      shared_ptr<IntervalConstraint> ic;
      cd = constraints["negGshape"][0];
      cu = constraints["negGshape"][1];
      val = (cd + cu) / 2.;
      ic.reset(new IntervalConstraint(cd, cu, true, true, prec));
      addParameter_(new Parameter("negGshape", val, ic));
    }

    if (shared_pp.posGshape > 0.)
    {
      double val, cd, cu;
      shared_ptr<IntervalConstraint> ic;
      cd = constraints["posGshape"][0];
      cu = constraints["posGshape"][1];
      val = (cd + cu) / 2.;
      ic.reset(new IntervalConstraint(cd, cu, true, true, prec));
      addParameter_(new Parameter("posGshape", val, ic));
    }

    if (shared_pp.pos_prop > 0.)
    {
      double val, cd, cu;
      shared_ptr<IntervalConstraint> ic;
      cd = constraints["pos_prop"][0];
      cu = constraints["pos_prop"][1];
      val = (cd + cu) / 2.;
      ic.reset(new IntervalConstraint(cd, cu, true, true, prec));
      addParameter_(new Parameter("pos_prop", val, ic));
    }
  }

  bpp_SFS_minus_lnl_shared* clone() const { return new bpp_SFS_minus_lnl_shared(*this); }

public:
  void setParameters(const ParameterList& pl)
  throw (ParameterNotFoundException, ConstraintException, Exception)
  {
    matchParametersValues(pl);
  }

  double getValue() const throw (Exception) { return lnL; }

  void fireParameterChanged(const ParameterList& pl)
  {
    struct shared_parameter_point my_shared_pp;

    if (hasParameter("negGshape")) my_shared_pp.negGshape = getParameterValue("negGshape"); else my_shared_pp.negGshape = -1.;
    if (hasParameter("posGshape")) my_shared_pp.posGshape = getParameterValue("posGshape"); else my_shared_pp.posGshape = -1.;
    if (hasParameter("pos_prop")) my_shared_pp.pos_prop = getParameterValue("pos_prop"); else my_shared_pp.pos_prop = -1.;

    double minus_lnL;
    multi_max_likelihood(data, nb_dataset, my_shared_pp, &minus_lnL);
    lnL = -minus_lnL;
  }
};


/* optimize_DFEM_shared */

vector< struct parameter_point > optimize_DFEM_shared(struct MK_data* data, int nb_dataset, struct shared_parameter_point* shared_pp, double* multi_ml)
{
  vector< struct parameter_point > v_opt;
  bpp_SFS_minus_lnl_shared lnL_fct(data, nb_dataset, *shared_pp);

  // starting points

  ParameterList init = lnL_fct.getParameters();
  if (init.hasParameter("negGshape"))
    init.setParameterValue("negGshape", shared_pp->negGshape);
  if (init.hasParameter("posGshape"))
    init.setParameterValue("posGshape", shared_pp->posGshape);
  if (init.hasParameter("pos_prop"))
    init.setParameterValue("pos_prop", shared_pp->pos_prop);

  lnL_fct.setParameters(init);

  // create relevant objects

  ReparametrizationFunctionWrapper rpf(&lnL_fct, false);

  ThreePointsNumericalDerivative tpnd(&rpf);
  tpnd.setParametersToDerivate(rpf.getParameters().getParameterNames());

  AbstractOptimizer* optim = 0;
  optim = new PseudoNewtonOptimizer(&tpnd);

  optim->setVerbose(1);
  optim->setProfiler(0);
  optim->setMessageHandler(0);
  optim->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);

  FunctionStopCondition mystop(optim);
  mystop.setTolerance(OPTIM_TOLERANCE);
  optim->setStopCondition(mystop);
  optim->setMaximumNumberOfEvaluations(OPTIM_MAX_EVAL);

  // optimize

  ParameterList init_transformed = rpf.getParameters();
  optim->init(init_transformed);
  optim->optimize();

  // collect results
  if (optim->getNumberOfEvaluations() < OPTIM_MIN_EVAL)
  {
    *multi_ml = -1.; return v_opt;
  }
  ParameterList finit = lnL_fct.getParameters();
  if (finit.hasParameter("negGshape"))
    shared_pp->negGshape = finit.getParameterValue("negGshape");
  if (finit.hasParameter("posGshape"))
    shared_pp->posGshape = finit.getParameterValue("posGshape");
  if (finit.hasParameter("pos_prop"))
    shared_pp->pos_prop = finit.getParameterValue("pos_prop");

  v_opt = multi_max_likelihood(data, nb_dataset, *shared_pp, multi_ml);

  delete optim;

  return v_opt;
}


/*************************************/
/************    MAIN    *************/
/*************************************/


int main(int argc, char** argv)
{
  /* 1. initializations */

  time_t t1, t2;
  t1 = time(&t1);
  srand48(1234);
  lsqrtpi = log(sqrt(PI));
  set_implemented_models();
  set_constraints();
  set_reparametrization();


  /* 2. manage options */

  string infile_name, outfile_name;

  set_default_options(&opt);
  get_options(argc, argv, &opt, &infile_name, &outfile_name);
  check_options(&opt, infile_name, outfile_name);
  nb_random_start_values = opt.nb_random_start;

  /* 3. read data file, pre-process data */

  FILE* input = fopen(infile_name.c_str(), "r");
  if (input == NULL)
  {
    printf("Can't read file: %s\n", infile_name.c_str()); exit(0);
  }
  int nb_dataset = 0;
  struct MK_data* dataset = read_dofe(input, &nb_dataset, opt.use_divergence_data);
  if (dataset == NULL)
  {
    printf("Bad input file: %s\n", infile_name.c_str()); exit(0);
  }

  global_folded = dataset[0].folded; // NB: same folded status for all datasets

  if (!global_folded & opt.fold)
  {
    for (int k = 0; k < nb_dataset; k++)
    {
      dataset[k].nb_SFScat = fold(dataset[k].specS, dataset[k].nb_gene);
      fold(dataset[k].specN, dataset[k].nb_gene);
      dataset[k].folded = true;
    }
    global_folded = true;
  }
  if (global_folded)
    opt.use_syno_orientation_error = false;

  cout << infile_name << ":" << endl;
  printf("%d data sets found\n", nb_dataset);
  printf("\n");

  /* 4. precalculations */

  cout << "Precalculations..." << endl;

  // factorial logarithms
  int maxtot = 0;
  for (int k = 0; k < nb_dataset; k++)
  {
    int tot = 0;
    for (int j = 1; j <= dataset[k].nb_SFScat; j++)
    {
      tot += (int)round(dataset[k].specN[j]) + (int)round(dataset[k].specS[j]);
    }
    if (dataset[k].div_data_available)
      tot += (int)round(dataset[k].fixN) + (int)round(dataset[k].fixS);
    if (tot > maxtot)
      maxtot = tot;
  }
  log_facto = (double*)calloc(maxtot + 3, sizeof(double));
  precalculate_log_facto(maxtot + 2);

  t2 = time(&t2);
  cout << endl << "log-facto done: ";
  printtime(difftime(t2, t1));

  // corrective term for dN (observing a substitution when the site is actually not fixed)
  positive_dN_corrective_term = (double**)check_alloc(nb_dataset, sizeof(double*));
  negative_dN_corrective_term = (double**)check_alloc(nb_dataset, sizeof(double*));
  for (int k = 0; k < nb_dataset; k++)
  {
    positive_dN_corrective_term[k] = (double*)check_alloc((int)VERY_LARGE_S_VALUE + 1, sizeof(double));
    negative_dN_corrective_term[k] = (double*)check_alloc((int)LARGE_S_VALUE + 1, sizeof(double));
  }
  for (int k = 0; k < nb_dataset; k++)
  {
    global_current_species = k;
    global_n = dataset[k].nb_gene;
    global_L = dataset[k].Ln_poly;
    precalculate_dN_corrective_term();
  }

  t2 = time(&t2);
  cout << endl << "corrective dN done: ";
  printtime(difftime(t2, t1));


  // expected counts for specific values of S (to be reused many times in numerical integration)
  // beware!! theta=ARBITRARY_THETA=1. assumed in precalculations!
  // Need to be consistent when calling function Gamma_folded_expected_allele_count (argument theta must be 1.)


  global_integral_S_point = precalculate_integral_S_points();
  global_integral_S_width = precalculate_integral_S_width();

  t2 = time(&t2);
  cout << endl << "integral done: ";
  printtime(difftime(t2, t1));

  vector< vector<double> > sp_exp_count;

  for (int k = 0; k < nb_dataset; k++)
  {
    global_theta = ARBITRARY_THETA;
    vector<double> zerov;
    sp_exp_count.clear();
    sp_exp_count.push_back(zerov);

    global_current_species = k;
    global_n = dataset[k].nb_gene;
    global_L = dataset[k].Ln_poly;
    for (int j = 1; j <= dataset[k].nb_SFScat; j++)
    {
      sp_exp_count.push_back(precalculate_expected_allele_counts(j));
    }

    global_exp_count.push_back(sp_exp_count);
  }

  t2 = time(&t2);
  cout << endl << "all done: ";
  printtime(difftime(t2, t1));


  // resize global expected count vectors

  int max_nbcat = dataset[0].nb_SFScat;
  for (int k = 1; k < nb_dataset; k++)
  {
    if (dataset[k].nb_SFScat > max_nbcat)
      max_nbcat = dataset[k].nb_SFScat;
  }

  prov_exp_specS.resize(max_nbcat + 1, 0);
  exp_specS.resize(max_nbcat + 1, 0);
  exp_specN.resize(max_nbcat + 1, 0);


  /* 5. likelihood calculation/maximisation */

  // step 0: saturated model likelihood

  vector<double> saturated_lnL;
  for (int k = 0; k < nb_dataset; k++)
  {
    saturated_lnL.push_back(-optimize_saturated(dataset + k));
  }

  // step 1: species after species

  vector< struct parameter_point > v_opt_sep;
  double multi_maxlnL_sep = 0.;

  if (opt.do_separate)
  {
    cout << "Likelihood maximization, separate..." << endl;

    for (int k = 0; k < nb_dataset; k++)
    {
      global_current_species = k;
      v_opt_sep.push_back(one_species_max_likelihood(&dataset[k]));
      cout << "  " << string(dataset[k].name) << "  ";
      cout << v_opt_sep[k].DFEM_param["negGshape"] << "," << v_opt_sep[k].DFEM_param["negGmean"] << endl;
    }
    for (int k = 0; k < nb_dataset; k++)
    {
      multi_maxlnL_sep += v_opt_sep[k].lnL;
    }
  }

  // step 2: shared parameters

  struct shared_parameter_point shared_pp;
  vector<struct shared_parameter_point> v_shared_pp;
  vector< vector< struct parameter_point > >  v_opt_shared;
  double multi_maxlnL_shared = 0.;
  vector<double> v_multi_maxlnL_shared;

  if (opt.do_shared)
  {
    if (opt.fixed_shared_negGshape.size() > 0)   // fixed shared shape(s)
    {
      shared_pp.posGshape = -1.;
      shared_pp.pos_prop = -1.;
      for (size_t i = 0; i < opt.fixed_shared_negGshape.size(); i++)
      {
        shared_pp.negGshape = opt.fixed_shared_negGshape[i];
        cout << "Likelihood calculation, shared shape=" << opt.fixed_shared_negGshape[i]  << "..." << endl;
        v_opt_shared.push_back(multi_max_likelihood(dataset, nb_dataset, shared_pp, &multi_maxlnL_shared));
        v_multi_maxlnL_shared.push_back(multi_maxlnL_shared);
        v_shared_pp.push_back(shared_pp);
      }
    }
    else// optimize shared shape
    {
      if (opt.do_separate)
        shared_pp.negGshape = median_negGshape(v_opt_sep);
      else
        shared_pp.negGshape = DEFAULT_INITIAL_SHARED_SHAPE;
      shared_pp.posGshape = -1.;
      shared_pp.pos_prop = -1.;
      cout << "Likelihood maximization, shared..." << endl;
      v_opt_shared.push_back(optimize_DFEM_shared(dataset, nb_dataset, &shared_pp, &multi_maxlnL_shared));
      v_multi_maxlnL_shared.push_back(multi_maxlnL_shared);
      v_shared_pp.push_back(shared_pp);
    }
  }


  /* 6. file output */

  vector<double> totSNP_S, totSNP_N;
  for (int k = 0; k < nb_dataset; k++)
  {
    double nS = 0.; double nN = 0.;
    for (int j = 1; j <= dataset[k].nb_SFScat; j++)
    {
      nN += dataset[k].specN[j];
      nS += dataset[k].specS[j];
    }
    totSNP_S.push_back(nS);
    totSNP_N.push_back(nN);
  }

  if (outfile_name.size())
  {
    FILE* output = fopen(outfile_name.c_str(), "w");
    // header line
    fprintf(output, "species,nb_gene,nb_SNP_S,nb_SNP_N,saturated_lnL");
    for (int j = 0; j < max_nbcat; j++)
    {
      fprintf(output, ",obs_S%d", j + 1);
    }
    for (int j = 0; j < max_nbcat; j++)
    {
      fprintf(output, ",obs_N%d", j + 1);
    }
    if (opt.do_separate)
    {
      fprintf(output, ",sep_lnL,sep_shape,sep_mean,sep_theta");
      for (int j = 0; j < max_nbcat; j++)
      {
        fprintf(output, ",exp_sep_S%d", j + 1);
      }
      for (int j = 0; j < max_nbcat; j++)
      {
        fprintf(output, ",exp_sep_N%d", j + 1);
      }
      for (int j = 0; j < max_nbcat; j++)
      {
        fprintf(output, ",sep_r%d", j + 1);
      }
    }
    if (opt.do_shared)
    {
      for (size_t i = 0; i < v_opt_shared.size(); i++)
      {
        fprintf(output, ",shared_lnL,shared_shape,shared_mean,shared_theta");
        for (int j = 0; j < max_nbcat; j++)
        {
          fprintf(output, ",exp_sha_S%d", j + 1);
        }
        for (int j = 0; j < max_nbcat; j++)
        {
          fprintf(output, ",exp_sha_N%d", j + 1);
        }
      }
    }
    fprintf(output, "\n");
    // result lines (one per species)
    struct SFS_div expected;
    for (int k = 0; k < nb_dataset; k++)
    {
      fprintf(output, "%s,%d,%d,%d,%f", dataset[k].name, dataset[k].nb_gene, (int)totSNP_S[k], (int)totSNP_N[k], saturated_lnL[k]);
      for (int i = 0; i < dataset[k].nb_SFScat; i++)
      {
        fprintf(output, ",%f", dataset[k].specS[i + 1]);
      }
      for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
      {
        fprintf(output, ",NA");
      }
      for (int i = 0; i < dataset[k].nb_SFScat; i++)
      {
        fprintf(output, ",%f", dataset[k].specN[i + 1]);
      }
      for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
      {
        fprintf(output, ",NA");
      }
      if (opt.do_separate)
      {
        fprintf(output, ",%f,%f,%f,%f", v_opt_sep[k].lnL, v_opt_sep[k].DFEM_param["negGshape"], v_opt_sep[k].DFEM_param["negGmean"], v_opt_sep[k].theta);
        global_current_species = k;
        expected = expected_SFS_div(dataset[k], v_opt_sep[k].DFEM_param, NULL, NULL);
        for (int i = 0; i < dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",%f", expected.specS[i + 1]);
        }
        for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",NA");
        }
        for (int i = 0; i < dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",%f", expected.specN[i + 1]);
        }
        for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",NA");
        }
        for (int i = 0; i < dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",%f", v_opt_sep[k].ri[i + 1]);
        }
        for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
        {
          fprintf(output, ",NA");
        }
      }

      if (opt.do_shared)
      {
        for (size_t i = 0; i < v_opt_shared.size(); i++)
        {
          fprintf(output, ",%f,%f,%f,%f", v_opt_shared[i][k].lnL, v_opt_shared[i][k].DFEM_param["negGshape"], v_opt_shared[i][k].DFEM_param["negGmean"], v_opt_shared[i][k].theta);
          global_current_species = k;
          expected = expected_SFS_div(dataset[k], v_opt_shared[i][k].DFEM_param, NULL, NULL);
          for (int i = 0; i < dataset[k].nb_SFScat; i++)
          {
            fprintf(output, ",%f", expected.specS[i + 1]);
          }
          for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
          {
            fprintf(output, ",NA");
          }
          for (int i = 0; i < dataset[k].nb_SFScat; i++)
          {
            fprintf(output, ",%f", expected.specN[i + 1]);
          }
          for (int i = 0; i < max_nbcat - dataset[k].nb_SFScat; i++)
          {
            fprintf(output, ",NA");
          }
        }
      }
      fprintf(output, "\n");
    }
  }

  /* 7. console output */

  cout << endl << endl;

  if (opt.do_separate)
  {
    cout << "Separate:" << endl;
    cout << "Maximum likelihood: " << multi_maxlnL_sep << endl;
    cout << "DFE shape:";
    for (int k = 0; k < nb_dataset; k++)
    {
      cout << " " << v_opt_sep[k].DFEM_param["negGshape"];
    }
    cout << endl;
    cout << "mean S:";
    for (int k = 0; k < nb_dataset; k++)
    {
      cout << " " << v_opt_sep[k].DFEM_param["negGmean"];
    }
    cout << endl;
    cout << "theta:";
    for (int k = 0; k < nb_dataset; k++)
    {
      cout << " " << v_opt_sep[k].theta;
    }
    cout << endl;
    cout << endl;
  }

  if (opt.do_shared)
    cout << "Shared:" << endl;
  for (size_t i = 0; i < v_opt_shared.size(); i++)
  {
    cout << "Maximum likelihood: " << v_multi_maxlnL_shared[i] << endl;
    cout << "DFE shape: " << v_shared_pp[i].negGshape << endl;
    cout << "mean S:";
    for (int k = 0; k < nb_dataset; k++)
    {
      cout << " " << v_opt_shared[i][k].DFEM_param["negGmean"];
    }
    cout << endl;
    cout << "theta:";
    for (int k = 0; k < nb_dataset; k++)
    {
      cout << " " << v_opt_shared[i][k].theta;
    }
    cout << endl << endl;
  }


  t2 = time(&t2);
  cout << endl << "Running time: ";
  printtime(difftime(t2, t1));
}
