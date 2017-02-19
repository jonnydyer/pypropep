
/* Python header files */
#include "Python.h"

#if PY_MAJOR_VERSION == 2
#define PyLong_AsLong PyInt_AsLong
#define PyLong_FromLong PyInt_FromLong
#define PyUnicode_Check3 PyString_Check
#define PyUnicode_FromFormat PyString_FromFormat
#define PyUnicode_FromString PyString_FromString
#define PyUnicode_FromStringAndSize PyString_FromStringAndSize
#define PyUnicode_Type PyString_Type
#define PyUnicode_AsUTF8 PyString_AsString
#define PyVarObject_HEAD_INIT(p, b) PyObject_HEAD_INIT(p) 0,
#define OB_REFCNT ob_refcnt
#define OB_TYPE ob_type
#else
#define PyUnicode_Check3 PyUnicode_Check
#define OB_REFCNT ob_base.ob_refcnt
#define OB_TYPE ob_base.ob_type
#if PY_MINOR_VERSION == 0 || PY_MINOR_VERSION == 1 || PY_MINOR_VERSION == 2
#define PyUnicode_AsUTF8 _PyUnicode_AsString
#endif
#endif

/* libthermo header files */

#include "load.h"
#include "thermo.h"

