#ifndef conversion_h
#define conversion_h


/* pressure units */
enum
{
  ATM,
  PSI,
  BAR,
  KPA
};

/* Transform calories to joules */
#define CAL_TO_JOULE       4.1868

/* Transform pound/(cubic inch) to gram/(cubic centimeter) */
#define LBS_IN3_TO_G_CM3  27.679905

/* Transform different pressure units */

#define ATM_TO_PA         101325.0
#define ATM_TO_PSI        14.695949
#define ATM_TO_BAR         1.01325

#define BAR_TO_PSI        14.503774

#define BAR_TO_ATM         0.98692327
#define PSI_TO_ATM         0.068045964
#define KPA_TO_ATM         0.0098692327

/* Length */

#define M_TO_CM             100.0
#define M_TO_IN              39.370079

#define IN_TO_M               0.0254

/* Surface */

#define M2_TO_CM2         10000.0
#define M2_TO_IN2          1550.0031

/* Volume */

#define M3_TO_CM3       1000000.0
#define M3_TO_IN3         61023.744

/* Mass flow */

#define KG_S_TO_LB_S     2.2046226


/* force */

/* newton to pound-force */
#define N_TO_LBF         0.22480894
#define LBF_TO_N         4.4482216

#endif
