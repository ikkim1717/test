#ifndef FMI2_FUNCTION_TYPES_H
#define FMI2_FUNCTION_TYPES_H

#include "fmi2TypesPlatform.h"

typedef enum {
    fmi2OK,
    fmi2Warning,
    fmi2Discard,
    fmi2Error,
    fmi2Fatal,
    fmi2Pending
} fmi2Status;

typedef enum {
    fmi2ModelExchange,
    fmi2CoSimulation
} fmi2Type;

typedef enum {
    fmi2DoStepStatus,
    fmi2PendingStatus,
    fmi2LastSuccessfulTime,
    fmi2Terminated
} fmi2StatusKind;

typedef void (*fmi2CallbackLogger)(
    fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message,
    ...);

typedef void* (*fmi2CallbackAllocateMemory)(size_t nobj, size_t size);
typedef void (*fmi2CallbackFreeMemory)(void* obj);
typedef void (*fmi2StepFinished)(fmi2ComponentEnvironment componentEnvironment, fmi2Status status);

typedef struct {
    fmi2CallbackLogger logger;
    fmi2CallbackAllocateMemory allocateMemory;
    fmi2CallbackFreeMemory freeMemory;
    fmi2StepFinished stepFinished;
    fmi2ComponentEnvironment componentEnvironment;
} fmi2CallbackFunctions;

#endif
