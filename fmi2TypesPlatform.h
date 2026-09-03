#ifndef FMI2_TYPES_PLATFORM_H
#define FMI2_TYPES_PLATFORM_H

#define fmi2TypesPlatform "default"
#define fmi2Version "2.0"

typedef void* fmi2Component;
typedef void* fmi2ComponentEnvironment;
typedef void* fmi2FMUstate;
typedef unsigned int fmi2ValueReference;
typedef double fmi2Real;
typedef int fmi2Integer;
typedef int fmi2Boolean;
typedef const char* fmi2String;
typedef char fmi2Byte;

#ifndef fmi2True
#define fmi2True 1
#endif
#ifndef fmi2False
#define fmi2False 0
#endif

#endif
