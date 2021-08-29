/* Original code from rapl-read.c by Vince Weaver, http://web.eece.maine.edu/~vweaver/projects/rapl/
* Modified by Ivana Marincic -- imarincic @ uchicago.edu -- July 2017 */
#ifndef __MSR_HANDLER_H
#define __MSR_HANDLER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>

#include <sys/syscall.h>
#include <linux/perf_event.h>

#define MSR_RAPL_POWER_UNIT     0x606

/*
 * Platform specific RAPL Domains.
 * Note that PP1 RAPL Domain is supported on 062A only
 * And DRAM RAPL Domain is supported on 062D only
 */
/* Package RAPL Domain */
#define MSR_PKG_POWER_LIMIT    0x610
#define MSR_PKG_ENERGY_STATUS       0x611
#define MSR_PKG_PERF_STATUS     0x613
#define MSR_PKG_POWER_INFO      0x614

/* PP0 RAPL Domain */
#define MSR_PP0_POWER_LIMIT     0x638
#define MSR_PP0_ENERGY_STATUS       0x639
#define MSR_PP0_POLICY          0x63A
#define MSR_PP0_PERF_STATUS     0x63B

/* PP1 RAPL Domain, may reflect to uncore devices */
#define MSR_PP1_POWER_LIMIT     0x640
#define MSR_PP1_ENERGY_STATUS       0x641
#define MSR_PP1_POLICY          0x642

/* DRAM RAPL Domain */
#define MSR_DRAM_POWER_LIMIT        0x618
#define MSR_DRAM_ENERGY_STATUS      0x619
#define MSR_DRAM_PERF_STATUS        0x61B
#define MSR_DRAM_POWER_INFO     0x61C

/* PSYS RAPL Domain */
//#define MSR_PLATFORM_ENERGY_STATUS  0x64d
#define MSR_PLATFORM_ENERGY_COUNTER  0x64d
#define MSR_PLATFORM_POWER_LIMIT 0x65C

#define IA32_THERM_STATUS 0x19C
#define IA32_MPERF 0xE7
#define IA32_APERF 0xE8

/* RAPL UNIT BITMASK */
#define POWER_UNIT_OFFSET   0
#define POWER_UNIT_MASK     0x0F

#define ENERGY_UNIT_OFFSET  0x08
#define ENERGY_UNIT_MASK    0x1F00

#define TIME_UNIT_OFFSET    0x10
#define TIME_UNIT_MASK      0xF000

#define CPU_SANDYBRIDGE     42
#define CPU_SANDYBRIDGE_EP  45
#define CPU_IVYBRIDGE       58
#define CPU_IVYBRIDGE_EP    62
#define CPU_HASWELL     60  // 69,70 too?
#define CPU_HASWELL_EP      63
#define CPU_BROADWELL       61  // 71 too?
#define CPU_BROADWELL_EP    79
#define CPU_BROADWELL_DE    86
#define CPU_SKYLAKE     78
#define CPU_SKYLAKE_HS      94
#define CPU_KNIGHTS_LANDING 87
#define CPU_KABYLAKE        142
#define CPU_KABYLAKE_2      158

#define MAX_CPUS    1024
#define MAX_PACKAGES    16
#define MAX_MSRS 25

#define DEFAULT_PKG_POW 215.0
#define DEFAULT_CORE_POW 0.0 //when disabled
#define DEFAULT_SHORT 258.0
#define DEFAULT_SECONDS_LONG 1.0
#define DEFAULT_SECONDS_SHORT 0.009765625000
#define DEFAULT_CORE_SECONDS 0.000976562500 //when disabled
#define MAX_WATTS 300.0
#define MIN_WATTS 50.0
#define NUM_ZONES 5 //PACKAGE, CORE, UNCORE, PLATFORM, DRAM (system/architecture dependent)
#define ZONE_NAME_LEN 10

#define PACKAGE_INDEX 0 //defining this because package is most commonly used zone

#define WATTS_LONG_START_BITS 0
#define WATTS_LONG_END_BITS 14
#define WATTS_SHORT_START_BITS 32
#define WATTS_SHORT_END_BITS 46
#define SECONDS_LONG_START_BITS 17
#define SECONDS_LONG_END_BITS 23
#define SECONDS_SHORT_START_BITS 49
#define SECONDS_SHORT_END_BITS 55
//these include the clamping bits too
#define ENABLED_LONG_START_BITS 15
#define ENABLED_LONG_END_BITS 16
#define ENABLED_SHORT_START_BITS 47
#define ENABLED_SHORT_END_BITS 48

#define NUM_RAPL_DOMAINS    5

#define BUFSIZE 500

struct system_info_t;

struct msr_info {
    int msr;
    int package_id;
    int cpu_id;
    double thermal_spec_power;
    double minimum_power;
    double maximum_power;
    double maximum_time_window;
};

struct msr_energy {
    int msr;
    int package_id;
    int cpu_id;
    double last_energy;
    double total_energy;
    int num_overflows;
    double cpu_energy_units;
    double dram_energy_units;
};

typedef enum zone_labels { PACKAGE, CORE, UNCORE, PLATFORM, DRAM} zone_label_t;

struct msr_pcap {
    int msr;
    int package_id;
    int cpu_id;
    zone_label_t zone_label;
    int enabled_long;
    int clamped_long;
    int enabled_short;
    int clamped_short;
    double watts_long;
    double watts_short;
    double seconds_long;
    double seconds_short;
};

struct msr_perf {
    int msr;
    int package_id;
    int cpu_id;
    double throttled_time;
};

struct msr_policy {
    int msr;
    int package_id;
    int cpu_id;
    int policy;
};

struct rapl_energy {
    double package;
    double pp0;
    double pp1;
    double dram;
    double platform;
};

struct system_msr_info {
    int error_state;
    /* general info */
    int cpu_model;
    int total_cores;
    /* package info */
    int total_packages;
    int package_map[MAX_PACKAGES];
    int package_fd[MAX_PACKAGES];

    /* msr units */
    double power_units;
    double time_units;
    double cpu_energy_units[MAX_PACKAGES];
    double dram_energy_units[MAX_PACKAGES];

    /* msr info */
    int msr_nums[5];
    int msrs[5][MAX_MSRS];

    /* list of specific msr groups */
    struct msr_info *info_msrs;
    struct msr_energy *energy_msrs;
    struct msr_pcap *pcap_msrs;
    struct msr_perf *perf_msrs;
    struct msr_policy *policy_msrs;

    int num_zones;
};

void init_msrs (struct system_info_t *system_info);
int finalize_msrs (struct system_info_t * system_info);
int rapl_set_power_cap (char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short, struct system_info_t * system_info, int enable);

int rapl_read_energy (struct rapl_energy * re, struct system_info_t * system_info);
int rapl_compute_total_power (struct rapl_energy *re, struct rapl_energy *energy, double time);
int rapl_compute_total_energy (struct rapl_energy *re, struct rapl_energy *end, struct rapl_energy *start);

int rapl_get_power_cap (struct msr_pcap *pcap, char *zone_name, struct system_info_t * system_info);
int rapl_get_power_cap_info(char *zone_name, double *min, double *max,
    double *thermal_spec, double *max_time_window, struct system_info_t * system_info);

#ifdef __cplusplus
}
#endif

#endif
