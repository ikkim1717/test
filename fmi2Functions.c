#include "fmi2Functions.h"
#include "fsm_suspension_core_from_m_original_strict.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FSM_FMU_GUID "8b861c20-2a5f-4b9d-9ee7-f5085a002001"
#define FSM_MODEL_IDENTIFIER "fsm_suspension_fmu"

#define VR_IN_BASE 0u
#define VR_OUT_BASE 5u

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    char *instanceName;
    fmi2CallbackFunctions cb;
    fmi2Boolean loggingOn;
    fmi2Boolean initialized;
    fmi2Real time;
    fmi2Real inputs[FSM_N_INPUTS];
    fmi2Real outputs[FSM_N_OUTPUTS];
    FsmState state;
} FmuInstance;

static void* fmu_alloc(const fmi2CallbackFunctions *cb, size_t nobj, size_t size)
{
    if (cb && cb->allocateMemory)
        return cb->allocateMemory(nobj, size);
    return calloc(nobj, size);
}

static void fmu_free(const fmi2CallbackFunctions *cb, void *p)
{
    if (!p)
        return;
    if (cb && cb->freeMemory)
        cb->freeMemory(p);
    else
        free(p);
}

static char* fmu_strdup(const fmi2CallbackFunctions *cb, const char *s)
{
    size_t n;
    char *p;
    if (!s) s = "";
    n = strlen(s) + 1u;
    p = (char*)fmu_alloc(cb, n, sizeof(char));
    if (p)
        memcpy(p, s, n);
    return p;
}

static FmuInstance* asInstance(fmi2Component c)
{
    return (FmuInstance*)c;
}

static void logMessage(FmuInstance *inst, fmi2Status status, const char *category, const char *message)
{
    if (inst && inst->loggingOn && inst->cb.logger) {
        inst->cb.logger(inst->cb.componentEnvironment,
            inst->instanceName ? inst->instanceName : FSM_MODEL_IDENTIFIER,
            status,
            category ? category : "log",
            message ? message : "");
    }
}

static fmi2Status unsupported(FmuInstance *inst, const char *what)
{
    logMessage(inst, fmi2Error, "unsupported", what);
    return fmi2Error;
}

const char* fmi2GetTypesPlatform(void)
{
    return fmi2TypesPlatform;
}

const char* fmi2GetVersion(void)
{
    return fmi2Version;
}

fmi2Status fmi2SetDebugLogging(fmi2Component c, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[])
{
    FmuInstance *inst = asInstance(c);
    (void)nCategories;
    (void)categories;
    if (!inst)
        return fmi2Error;
    inst->loggingOn = loggingOn;
    return fmi2OK;
}

fmi2Component fmi2Instantiate(fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGUID,
    fmi2String fmuResourceLocation, const fmi2CallbackFunctions* functions, fmi2Boolean visible, fmi2Boolean loggingOn)
{
    FmuInstance *inst;
    size_t k;

    (void)fmuResourceLocation;
    (void)visible;

    if (fmuType != fmi2CoSimulation)
        return NULL;
    if (fmuGUID && strcmp(fmuGUID, FSM_FMU_GUID) != 0)
        return NULL;

    inst = (FmuInstance*)fmu_alloc(functions, 1u, sizeof(FmuInstance));
    if (!inst)
        return NULL;
    memset(inst, 0, sizeof(*inst));

    if (functions)
        inst->cb = *functions;
    inst->loggingOn = loggingOn;

    inst->instanceName = fmu_strdup(functions, instanceName ? instanceName : FSM_MODEL_IDENTIFIER);
    if (!inst->instanceName) {
        fmu_free(functions, inst);
        return NULL;
    }

    for (k = 0; k < FSM_N_INPUTS; k++)
        inst->inputs[k] = 0.0;
    for (k = 0; k < FSM_N_OUTPUTS; k++)
        inst->outputs[k] = 0.0;

    fsm_state_reset(&inst->state);
    inst->initialized = fmi2False;
    inst->time = 0.0;

    return inst;
}

void fmi2FreeInstance(fmi2Component c)
{
    FmuInstance *inst = asInstance(c);
    fmi2CallbackFunctions cb;
    if (!inst)
        return;
    cb = inst->cb;
    fmu_free(&cb, inst->instanceName);
    fmu_free(&cb, inst);
}

fmi2Status fmi2SetupExperiment(fmi2Component c, fmi2Boolean toleranceDefined, fmi2Real tolerance,
    fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime)
{
    FmuInstance *inst = asInstance(c);
    (void)toleranceDefined;
    (void)tolerance;
    (void)stopTimeDefined;
    (void)stopTime;
    if (!inst)
        return fmi2Error;
    inst->time = startTime;
    return fmi2OK;
}

fmi2Status fmi2EnterInitializationMode(fmi2Component c)
{
    FmuInstance *inst = asInstance(c);
    if (!inst)
        return fmi2Error;
    fsm_state_reset(&inst->state);
    memset(inst->outputs, 0, sizeof(inst->outputs));
    return fmi2OK;
}

fmi2Status fmi2ExitInitializationMode(fmi2Component c)
{
    FmuInstance *inst = asInstance(c);
    if (!inst)
        return fmi2Error;
    inst->initialized = fmi2True;
    return fmi2OK;
}

fmi2Status fmi2Terminate(fmi2Component c)
{
    FmuInstance *inst = asInstance(c);
    if (!inst)
        return fmi2Error;
    inst->initialized = fmi2False;
    return fmi2OK;
}

fmi2Status fmi2Reset(fmi2Component c)
{
    FmuInstance *inst = asInstance(c);
    if (!inst)
        return fmi2Error;
    fsm_state_reset(&inst->state);
    memset(inst->outputs, 0, sizeof(inst->outputs));
    inst->time = 0.0;
    inst->initialized = fmi2False;
    return fmi2OK;
}

fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[])
{
    FmuInstance *inst = asInstance(c);
    size_t i;
    if (!inst || !vr || !value)
        return fmi2Error;
    for (i = 0; i < nvr; i++) {
        if (vr[i] < FSM_N_INPUTS) {
            value[i] = inst->inputs[vr[i]];
        } else if (vr[i] >= VR_OUT_BASE && vr[i] < VR_OUT_BASE + FSM_N_OUTPUTS) {
            value[i] = inst->outputs[vr[i] - VR_OUT_BASE];
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[])
{
    FmuInstance *inst = asInstance(c);
    size_t i;
    if (!inst || !vr || !value)
        return fmi2Error;
    for (i = 0; i < nvr; i++) {
        if (vr[i] < FSM_N_INPUTS) {
            if (!isfinite(value[i]))
                return fmi2Error;
            inst->inputs[vr[i]] = value[i];
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

fmi2Status fmi2DoStep(fmi2Component c, fmi2Real currentCommunicationPoint,
    fmi2Real communicationStepSize, fmi2Boolean noSetFMUStatePriorToCurrentPoint)
{
    FmuInstance *inst = asInstance(c);
    FsmStatus status;
    int nSteps;
    int k;
    double target;

    (void)noSetFMUStatePriorToCurrentPoint;

    if (!inst)
        return fmi2Error;

    if (communicationStepSize <= 0.0 || !isfinite(communicationStepSize)) {
        logMessage(inst, fmi2Error, "error", "communicationStepSize must be positive and finite.");
        return fmi2Error;
    }

    nSteps = (int)floor(communicationStepSize / FSM_DT + 0.5);
    if (nSteps < 1)
        nSteps = 1;

    target = (double)nSteps * FSM_DT;
    if (fabs(target - communicationStepSize) > 1e-9) {
        logMessage(inst, fmi2Error, "error", "communicationStepSize must be an integer multiple of 0.001 seconds.");
        return fmi2Error;
    }

    fsm_status_ok(&status);

    for (k = 0; k < nSteps; k++) {
        if (!fsm_step(&inst->state,
                inst->inputs[0],
                inst->inputs[1],
                inst->inputs[2],
                inst->inputs[3],
                inst->inputs[4],
                inst->outputs,
                &status)) {
            logMessage(inst, fmi2Error, "error", status.message[0] ? status.message : "fsm_step failed.");
            return fmi2Error;
        }
    }

    inst->time = currentCommunicationPoint + communicationStepSize;
    return fmi2OK;
}

fmi2Status fmi2CancelStep(fmi2Component c)
{
    return unsupported(asInstance(c), "fmi2CancelStep is not supported.");
}

fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[], const fmi2Real value[])
{
    (void)vr; (void)nvr; (void)order; (void)value;
    return unsupported(asInstance(c), "Input derivatives are not supported.");
}

fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr,
    const fmi2Integer order[], fmi2Real value[])
{
    (void)vr; (void)nvr; (void)order; (void)value;
    return unsupported(asInstance(c), "Output derivatives are not supported.");
}

fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind s, fmi2Status* value)
{
    (void)s;
    if (!c || !value)
        return fmi2Error;
    *value = fmi2OK;
    return fmi2OK;
}

fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind s, fmi2Real* value)
{
    FmuInstance *inst = asInstance(c);
    if (!inst || !value)
        return fmi2Error;
    if (s == fmi2LastSuccessfulTime) {
        *value = inst->time;
        return fmi2OK;
    }
    return fmi2Discard;
}

fmi2Status fmi2GetIntegerStatus(fmi2Component c, const fmi2StatusKind s, fmi2Integer* value)
{
    (void)c; (void)s; (void)value;
    return fmi2Discard;
}

fmi2Status fmi2GetBooleanStatus(fmi2Component c, const fmi2StatusKind s, fmi2Boolean* value)
{
    (void)c; (void)s; (void)value;
    return fmi2Discard;
}

fmi2Status fmi2GetStringStatus(fmi2Component c, const fmi2StatusKind s, fmi2String* value)
{
    (void)c; (void)s; (void)value;
    return fmi2Discard;
}

fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "Integer variables are not used.");
}

fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "Integer variables are not used.");
}

fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "Boolean variables are not used.");
}

fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "Boolean variables are not used.");
}

fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "String variables are not used.");
}

fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[])
{
    (void)vr; (void)nvr; (void)value;
    return unsupported(asInstance(c), "String variables are not used.");
}

fmi2Status fmi2GetFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
    (void)FMUstate;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2SetFMUstate(fmi2Component c, fmi2FMUstate FMUstate)
{
    (void)FMUstate;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate)
{
    (void)FMUstate;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t* size)
{
    (void)FMUstate; (void)size;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2SerializeFMUstate(fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size)
{
    (void)FMUstate; (void)serializedState; (void)size;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2DeSerializeFMUstate(fmi2Component c, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate)
{
    (void)serializedState; (void)size; (void)FMUstate;
    return unsupported(asInstance(c), "FMU state serialization is not supported.");
}

fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown_ref[], size_t nUnknown,
    const fmi2ValueReference vKnown_ref[], size_t nKnown, const fmi2Real dvKnown[], fmi2Real dvUnknown[])
{
    (void)vUnknown_ref; (void)nUnknown; (void)vKnown_ref; (void)nKnown; (void)dvKnown; (void)dvUnknown;
    return unsupported(asInstance(c), "Directional derivatives are not supported.");
}
