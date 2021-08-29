#ifndef __CRAY_PM_HANDLER_H
#define __CRAY_PM_HANDLER_H

#ifdef __cplusplus
extern "C"
{
#endif

#define MAX_NUM_COUNTERS 12
#define CRAY_FREQ_INDEX 7
#define NUM_COUNTERS 12

struct system_info_t;

typedef enum { ENERGY, POWER, CPU_ENERGY, CPU_POWER, MEMORY_ENERGY, MEMORY_POWER, POWER_CAP, RAW_SCAN_HZ, FRESHNESS, GENERATION, VERSION, STARTUP} cray_counter_type;

struct pm_counter {
    cray_counter_type type;
    char *pm_filename;
    int pm_file;
    double measurement;
};

struct cray_measurement
{
    double node_energy;
    double node_power;
    double node_measured_power;
    double cpu_energy;
    double cpu_power;
    double cpu_measured_power;
    double memory_energy;
    double memory_power;
    double memory_measured_power;
};

struct system_cray_info {
    struct pm_counter *counters;
    int num_counters;
};

int init_cray_pm_counters (struct system_info_t * system_info);
int finalize_cray_pm_counters (struct system_info_t * system_info);
int cray_open_pm_file (char *counter_name);
double cray_read_pm_counter(int pm_file);
int get_cray_measurement (struct cray_measurement *cm, struct system_info_t * system_info);
int compute_cray_total_measurements (struct cray_measurement *cm, struct cray_measurement *end, struct cray_measurement *start, double total_time);

#ifdef __cplusplus
}
#endif

#endif