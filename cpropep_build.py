# file "cpropep_build.py"
import os
from cffi import FFI
from glob import glob

ffibuilder = FFI()

exclude_c_files = ['test.c']

cpropep_libs = ['libnum', 'libthermo', 'libcpropep', 'libcompat']
inc_dir = [('pypropep/cpropep/' + d + '/include/') for d in cpropep_libs]

MAX_PRODUCT  = 400
MAX_ELEMENT  = 15
MAX_COMP     = 20

src_files = []
for l in cpropep_libs:
    prefix = 'pypropep/cpropep/' + l + '/src/'
    src_files += glob(prefix + '*.c')
    for ef in exclude_c_files:
        if (prefix + ef) in src_files:
            print('Excluding {} from source list'.format(prefix+ef))
            src_files.remove(prefix + ef)

inc_files = ''
for i in inc_dir:
    files = glob(i + '/*.h')
    for f in files:
        inc_files += '#include "%s"\n' % (os.path.basename(f))

ffibuilder.set_source("pypropep.cpropep._cpropep",
    inc_files,
    sources=src_files,
    include_dirs=inc_dir)

# TODO:Find a way to scrape #defines from headers rather than hard coding const
ffibuilder.cdef("""
//**** libcpropep/type.h ****//
typedef enum
{
  GAS,
  CONDENSED,
  STATE_LAST,
  ...
} state_t;

typedef enum
{
  TP,          /* assign temperature and pressure */
  HP,          /* assign enthalpy and pressure */
  SP,           /* assign entropy and pressure */
  ...
} problem_t;

typedef enum
{
  SUBSONIC_AREA_RATIO,
  SUPERSONIC_AREA_RATIO,
  PRESSURE,
  ...
} exit_condition_t;

typedef struct _performance_prop
{
  double ae_at;   /* Exit aera / Throat aera              */
  double a_dotm;  /* Exit aera / mass flow rate (m/s/atm) */
  double cstar;   /* Characteristic velocity              */
  double cf;      /* Coefficient of thrust                */
  double Ivac;    /* Specific impulse (vacuum)            */
  double Isp;     /* Specific impulse                     */
  ...;
} performance_prop_t;

typedef struct _composition
{
  short  ncomp;              /* Number of different component */
  short  molecule[20];      /* Molecule code                 */
  double coef[20];          /* Moles of molecule             */
  double density;            /* Density of propellant         */
  ...;
} composition_t;

typedef struct _product
{
  int   element_listed;                 /* true if element have been listed */
  int   product_listed;                 /* true if product have been listed */
  int   isequil;                        /* true if equilibrium is ok        */

  /* coefficient matrix for the gases */
  unsigned short A[15][400];

  short  n_element;                        /* n. of different element        */
  short  element[15];                     /* element list                   */
  short  n[STATE_LAST];                    /* n. of species for each state   */
  short  n_condensed;                      /* n. of total possible condensed */
  short  species[STATE_LAST][400];         /* possible species in each state */
  double coef[STATE_LAST][400];            /* coef. of each molecule         */
  ...;
} product_t;

typedef struct _iteration_var
{
  double n;                        /* mol/g of the mixture                  */
  double ln_n;                     /* ln(n)                                 */
  double sumn;                     /* sum of all the nj                     */
  double delta_ln_n;               /* delta ln(n) in the iteration process  */
  double delta_ln_T;               /* delta ln(T) in the iteration process  */
  double delta_ln_nj[400];         /* delta ln(nj) in the iteration process */
  double ln_nj[400];               /* ln(nj) nj are the individual mol/g    */
  ...;
} iteration_var_t;

typedef struct _equilib_prop
{
  double P;    /* Pressure (atm)              */
  double T;    /* Temperature (K)             */
  double H;    /* Enthalpy (kJ/kg)            */
  double U;    /* Internal energy (kJ/kg)     */
  double G;    /* Gibbs free energy (kJ/kg)   */
  double S;    /* Entropy (kJ/(kg)(K))        */
  double M;    /* Molar mass (g/mol)          */
  double dV_P; /* (d ln(V) / d ln(P))t        */
  double dV_T; /* (d ln(V) / d ln(T))p        */
  double Cp;   /* Specific heat (kJ/(kg)(K))  */
  double Cv;   /* Specific heat (kJ/(kg)(K))  */
  double Isex; /* Isentropic exponent (gamma) */
  double Vson; /* Sound speed (m/s)           */
  ...;
} equilib_prop_t;


typedef struct _new_equilibrium
{
  int equilibrium_ok;  /* true if the equilibrium have been compute */
  int properties_ok;   /* true if the properties have been compute  */
  int performance_ok;  /* true if the performance have been compute */

  //temporarily
  double entropy;

  iteration_var_t    itn;
  composition_t      propellant;
  product_t          product;
  equilib_prop_t     properties;
  performance_prop_t performance;
  ...;
} equilibrium_t;

//**** libthermo/load.h ****//
int load_thermo(char *filename);
int load_propellant(char *filename);

//**** libthermo/thermo.h *****//
typedef struct _thermo
{
  char    name[19];
  char    comments[57];
  int     nint;         /* number of different temperature interval */
  char    id[7];        /* identification code */
  int     elem[5];
  int     coef[5];
  state_t state;
  double  weight;       /* molecular weight */
  float   heat;         /* heat of formation at 298.15 K  (J/mol)  */
  double  dho;          /* HO(298.15) - HO(0) */
  float   range[4][2];  /* temperature range */
  int     ncoef[4];     /* number of coefficient for Cp0/R   */
  int     ex[4][8];     /* exponent in empirical equation */

  double param[4][9];

  /* for species with data at only one temperature */
  /* especially condensed                          */
  float temp;
  float enth;

} thermo_t;

typedef struct _propellant{
  char  name[120]; /* name of the propellant */
  int   elem[6];   /* element in the molecule (atomic number) max 6 */
  int   coef[6];   /* stochiometric coefficient of this element
		                  (0 for none) */
  float heat;      /* heat of formation in Joule/gram */
  float density;   /* density in g/cubic cm */
} propellant_t;

extern propellant_t	*propellant_list;
extern thermo_t	    *thermo_list;

extern const float molar_mass[];
extern const char symb[][3];

extern unsigned long num_thermo;
extern unsigned long num_propellant;

int thermo_search(char *str);
int propellant_search(char *str);
int atomic_number(char *symbole);
int propellant_search_by_formula(char *str);
double enthalpy_0(int sp, float T);
double entropy_0(int sp, float T);
double entropy(int sp, state_t st, double ln_nj_n, float T, float P);
double specific_heat_0(int sp, float T);
double mixture_specific_heat_0(equilibrium_t *e, double temp);
int temperature_check(int sp, float T);
double transition_temperature(int sp, float T);
double propellant_enthalpy(equilibrium_t *e);
double product_enthalpy(equilibrium_t *e);
double product_entropy(equilibrium_t *e);
double propellant_mass(equilibrium_t *e);
int compute_density(composition_t *c);
double gibbs_0(int sp, float T);
double gibbs(int sp, state_t st, double nj_n_n, float T, float P);
double heat_of_formation(int molecule);
double propellant_molar_mass(int molecule);

//**** libcpropep/equillibrium.h ****//
int reset_element_list(equilibrium_t *e);
int initialize_equilibrium(equilibrium_t *e);
int reset_equilibrium(equilibrium_t *e);
int copy_equilibrium(equilibrium_t *dest, equilibrium_t *src);
int compute_thermo_properties(equilibrium_t *e);
int set_state(equilibrium_t *e, double T, double P);
int add_in_propellant(equilibrium_t *e, int sp, double mol);
int equilibrium(equilibrium_t *equil, problem_t P);
double product_molar_mass(equilibrium_t *e);

//**** libcpropep/performance.h ****//
int frozen_performance(equilibrium_t *e, exit_condition_t exit_type,
                       double value);
int shifting_performance(equilibrium_t *e, exit_condition_t exit_type,
                         double value);
    """)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
