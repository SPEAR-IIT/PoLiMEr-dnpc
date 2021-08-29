/* Original code for accessing MSRs from rapl-read.c by Vince Weaver, http://web.eece.maine.edu/~vweaver/projects/rapl/
* Original code for writing to MSRs from raplcap-msr.c by Connor Imes, https://github.com/connorimes/raplcap/blob/master/msr/raplcap-msr.c
* Heavily modified by Ivana Marincic -- imarincic @ uchicago.edu -- July 2017 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include <sys/syscall.h>
#include <linux/perf_event.h>

#include "PoLiMEr.h"
#include "PoLiLog.h"
#include "msr-handler.h"

static int short_term_supported (int msr);
static int verify_power_limits(double watts, int enable);
static int get_msr_for_zone_name(char *zone_name, int get_pcap);

static int detect_cpu(void);
static int verify_model(int model);
static int detect_packages (struct system_info_t *system_info);

static void get_msr_units(struct system_info_t *system_info, int package);

static int open_msr(int core);
static long long read_msr(int fd, int msr_address);

static int read_msr_info (struct msr_info *msr_info, struct system_info_t *system_info, int package_id);
static int read_msr_pcap (struct msr_pcap *msr_pcap, struct system_info_t *system_info, int package_id);
static int read_msr_perf (struct msr_perf *msr_perf, struct system_info_t *system_info, int package_id);
static int read_msr_policy (struct msr_policy *msr_policy, struct system_info_t *system_info, int package_id);
static int read_msr_energy (struct msr_energy *msr_energy, struct system_info_t * system_info, int package_id);

static int set_msr_pcap(struct msr_pcap *pcap, struct system_info_t * system_info, int package_id);
static uint64_t to_msr_power(double watts, double power_units);
static int write_msr(int fd, int msr_address, uint64_t data);
static uint64_t replace_bits(uint64_t msrval, uint64_t data, uint8_t first, uint8_t last);
static uint64_t get_bits(uint64_t msrval, uint8_t first, uint8_t last);
static uint64_t to_msr_time(double seconds, double time_units);
static double from_msr_time(uint64_t y, uint64_t f, double time_units);

static int short_term_supported (int msr);

static uint64_t log2_u64(uint64_t y);
static uint64_t pow2_u64(uint64_t y);

int sandybridge_energy_msrs[3] = {MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int sandybridge_pcap_msrs[3] = {MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int sandybridge_perf_msrs[1] = {-1};
int sandybridge_info_msrs[1] = {MSR_PKG_POWER_INFO};
int sandybridge_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int sandybridge_msr_nums[5] = {3,3,1,0,2};

int sandybridge_ep_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int sandybridge_ep_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int sandybridge_ep_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int sandybridge_ep_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int sandybridge_ep_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int sandybridge_ep_msr_nums[5] = {4,4,2,2,2};

int ivybridge_energy_msrs[3] = {MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int ivybridge_pcap_msrs[3] = {MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int ivybridge_perf_msrs[1] = {-1};
int ivybridge_info_msrs[1] = {MSR_PKG_POWER_INFO};
int ivybridge_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int ivybridge_msr_nums[5] = {3,3,1,0,2};

int ivybridge_ep_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int ivybridge_ep_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int ivybridge_ep_perf_msrs[1] = {MSR_PKG_PERF_STATUS};
int ivybridge_ep_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int ivybridge_ep_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int ivybridge_ep_msr_nums[5] = {4,4,2,1,2};

int haswell_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int haswell_pcap_msrs[3] = {MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int haswell_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int haswell_info_msrs[1] = {MSR_PKG_POWER_INFO};
int haswell_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int haswell_msr_nums[5] = {4,3,1,2,2};

int haswell_ep_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int haswell_ep_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int haswell_ep_perf_msrs[2] = {MSR_DRAM_PERF_STATUS};
int haswell_ep_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int haswell_ep_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int haswell_ep_msr_nums[5] = {4,4,2,1,2};

int broadwell_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int broadwell_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int broadwell_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int broadwell_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int broadwell_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int broadwell_msr_nums[5] = {4,4,2,2,2};

int broadwell_ep_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int broadwell_ep_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int broadwell_ep_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int broadwell_ep_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int broadwell_ep_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int broadwell_ep_msr_nums[5] = {4,4,2,2,2};

int broadwell_de_energy_msrs[4] = {MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int broadwell_de_pcap_msrs[4] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int broadwell_de_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int broadwell_de_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int broadwell_de_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int broadwell_de_msr_nums[5] = {4,4,2,2,2};

int skylake_energy_msrs[5] = {MSR_PLATFORM_ENERGY_COUNTER, MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int skylake_pcap_msrs[5] = {MSR_PLATFORM_POWER_LIMIT, MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int skylake_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int skylake_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int skylake_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int skylake_msr_nums[5] = {5,5,2,2,2};

int skylake_hs_energy_msrs[5] = {MSR_PLATFORM_ENERGY_COUNTER, MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int skylake_hs_pcap_msrs[5] = {MSR_PLATFORM_POWER_LIMIT, MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int skylake_hs_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int skylake_hs_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int skylake_hs_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int skylake_hs_msr_nums[5] = {5,5,2,2,2};

int kabylake_energy_msrs[5] = {MSR_PLATFORM_ENERGY_COUNTER, MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int kabylake_pcap_msrs[5] = {MSR_PLATFORM_POWER_LIMIT, MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int kabylake_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int kabylake_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int kabylake_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int kabylake_msr_nums[5] = {5,5,2,2,2};

int kabylake_2_energy_msrs[5] = {MSR_PLATFORM_ENERGY_COUNTER, MSR_DRAM_ENERGY_STATUS, MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_PP1_ENERGY_STATUS};
int kabylake_2_pcap_msrs[5] = {MSR_PLATFORM_POWER_LIMIT, MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT, MSR_PP1_POWER_LIMIT};
int kabylake_2_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int kabylake_2_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int kabylake_2_policy_msrs[2] = {MSR_PP0_POLICY, MSR_PP1_POLICY};
int kabylake_2_msr_nums[5] = {5,5,2,2,2};


int knights_landing_energy_msrs[3] = {MSR_PKG_ENERGY_STATUS, MSR_PP0_ENERGY_STATUS, MSR_DRAM_ENERGY_STATUS};
int knights_landing_pcap_msrs[3] = {MSR_DRAM_POWER_LIMIT, MSR_PKG_POWER_LIMIT, MSR_PP0_POWER_LIMIT};
int knights_landing_perf_msrs[2] = {MSR_PKG_PERF_STATUS, MSR_DRAM_PERF_STATUS};
int knights_landing_info_msrs[2] = {MSR_PKG_POWER_INFO, MSR_DRAM_POWER_INFO};
int knights_landing_policy_msrs[1] = {-1};
int knights_landing_msr_nums[5] = {3,3,2,2,0};

void init_msrs (struct system_info_t * system_info)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    system_info->sysmsr = malloc(sizeof(struct system_msr_info));

    system_info->sysmsr->error_state = 1;

    system_info->sysmsr->info_msrs = 0;
    system_info->sysmsr->energy_msrs = 0;
    system_info->sysmsr->pcap_msrs = 0;
    system_info->sysmsr->perf_msrs = 0;
    system_info->sysmsr->policy_msrs = 0;
    system_info->sysmsr->num_zones = 0;

    system_info->sysmsr->cpu_model = detect_cpu();

    if (system_info->sysmsr->cpu_model < 0)
    {
        poli_log(ERROR, NULL, "Something went wrong with determining CPU. Will not use RAPL interface.");
        return;
    }
    else
        system_info->sysmsr->error_state = 0;

    detect_packages(system_info);

    int i, j;
    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < MAX_MSRS; j++)
            system_info->sysmsr->msrs[i][j] = 0;
    }



    //TODO
    system_info->sysmsr->num_zones = 3;

    switch (system_info->sysmsr->cpu_model)
    {
        case CPU_SANDYBRIDGE:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = sandybridge_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = sandybridge_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = sandybridge_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = sandybridge_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = sandybridge_policy_msrs[j];

            system_info->sysmsr->msrs[3][0] = sandybridge_perf_msrs[0];
            break;
        case CPU_SANDYBRIDGE_EP:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = sandybridge_ep_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = sandybridge_ep_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = sandybridge_ep_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = sandybridge_ep_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = sandybridge_ep_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = sandybridge_ep_policy_msrs[j];
            break;
        case CPU_IVYBRIDGE:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = ivybridge_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = ivybridge_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = ivybridge_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = ivybridge_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = ivybridge_policy_msrs[j];

            system_info->sysmsr->msrs[3][0] = ivybridge_perf_msrs[0];
            break;
        case CPU_IVYBRIDGE_EP:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = ivybridge_ep_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = ivybridge_ep_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = ivybridge_ep_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = ivybridge_ep_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = ivybridge_ep_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = ivybridge_ep_policy_msrs[j];
            break;
        case CPU_HASWELL:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = haswell_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = haswell_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = haswell_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = haswell_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = haswell_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = haswell_policy_msrs[j];
            break;
        case CPU_HASWELL_EP:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = haswell_ep_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = haswell_ep_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = haswell_ep_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = haswell_ep_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = haswell_ep_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = haswell_ep_policy_msrs[j];
            break;
        case CPU_BROADWELL:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = broadwell_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = broadwell_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = broadwell_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = broadwell_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = broadwell_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = broadwell_policy_msrs[j];
            break;
        case CPU_BROADWELL_EP:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = broadwell_ep_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = broadwell_ep_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = broadwell_ep_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = broadwell_ep_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = broadwell_ep_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = broadwell_ep_policy_msrs[j];
            break;
        case CPU_BROADWELL_DE:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = broadwell_de_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = broadwell_de_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = broadwell_de_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = broadwell_de_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = broadwell_de_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = broadwell_de_policy_msrs[j];
            break;
        case CPU_SKYLAKE:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = skylake_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = skylake_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = skylake_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = skylake_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = skylake_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = skylake_policy_msrs[j];
            break;
        case CPU_SKYLAKE_HS:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = skylake_hs_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = skylake_hs_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = skylake_hs_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = skylake_hs_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = skylake_hs_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = skylake_hs_policy_msrs[j];
            break;
        case CPU_KABYLAKE:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = kabylake_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = kabylake_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = kabylake_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = kabylake_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = kabylake_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = kabylake_policy_msrs[j];
            break;
        case CPU_KABYLAKE_2:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = kabylake_2_msr_nums[i];
            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = kabylake_2_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = kabylake_2_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = kabylake_2_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = kabylake_2_perf_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[4]; j++)
                system_info->sysmsr->msrs[4][j] = kabylake_2_policy_msrs[j];
            break;
        case CPU_KNIGHTS_LANDING:
            for (i = 0; i < 5; i++)
                system_info->sysmsr->msr_nums[i] = knights_landing_msr_nums[i];

            for (j = 0; j < system_info->sysmsr->msr_nums[0]; j++)
                system_info->sysmsr->msrs[0][j] = knights_landing_energy_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[1]; j++)
                system_info->sysmsr->msrs[1][j] = knights_landing_pcap_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[2]; j++)
                system_info->sysmsr->msrs[2][j] = knights_landing_info_msrs[j];
            for (j = 0; j < system_info->sysmsr->msr_nums[3]; j++)
                system_info->sysmsr->msrs[3][j] = knights_landing_perf_msrs[j];

            system_info->sysmsr->msrs[4][0] = knights_landing_policy_msrs[0];
            break;
        default:
            break;
    }

    system_info->sysmsr->energy_msrs = calloc(system_info->sysmsr->msr_nums[0] * system_info->sysmsr->total_packages, sizeof(struct msr_energy));

    system_info->sysmsr->pcap_msrs = calloc(system_info->sysmsr->msr_nums[1] * system_info->sysmsr->total_packages, sizeof(struct msr_pcap));

    system_info->sysmsr->info_msrs = calloc(system_info->sysmsr->msr_nums[2] * system_info->sysmsr->total_packages, sizeof(struct msr_info));

    system_info->sysmsr->perf_msrs = calloc(system_info->sysmsr->msr_nums[3] * system_info->sysmsr->total_packages, sizeof(struct msr_perf));

    system_info->sysmsr->policy_msrs = calloc(system_info->sysmsr->msr_nums[4] * system_info->sysmsr->total_packages, sizeof(struct msr_policy));

    int package;

    for (package = 0; package < system_info->sysmsr->total_packages; package++)
    {
        int cpu_id = system_info->sysmsr->package_map[package];
        int fd = open_msr(cpu_id);
        if (fd < 0)
        {
            poli_log(ERROR, NULL, "Failed to open any MSR file. There won't be any measurements using RAPL.");
            system_info->sysmsr->error_state = 1;
            return;
        }
        else
            system_info->sysmsr->error_state = 0;

        system_info->sysmsr->package_fd[package] = fd;
        get_msr_units(system_info, package);
        int msr;

        for (msr = 0; msr < system_info->sysmsr->msr_nums[0]; msr++)
        {
            struct msr_energy *emsr = &system_info->sysmsr->energy_msrs[msr + (package * system_info->sysmsr->msr_nums[0])];
            emsr->msr = system_info->sysmsr->msrs[0][msr];
            emsr->package_id = package;
            emsr->cpu_id = cpu_id;
            emsr->last_energy = 0.0;
            emsr->num_overflows = 0;
            emsr->cpu_energy_units = system_info->sysmsr->cpu_energy_units[package];
            emsr->dram_energy_units = system_info->sysmsr->dram_energy_units[package];
        }

        for (msr = 0; msr < system_info->sysmsr->msr_nums[1]; msr++)
        {
            struct msr_pcap *pcapmsr = &system_info->sysmsr->pcap_msrs[msr + (package * system_info->sysmsr->msr_nums[1])];
            pcapmsr->msr = system_info->sysmsr->msrs[1][msr];
            pcapmsr->package_id = package;
            pcapmsr->cpu_id = cpu_id;
            pcapmsr->enabled_long = 0;
            pcapmsr->clamped_long = 0;
            pcapmsr->enabled_short = 0;
            pcapmsr->enabled_long = 0;
            pcapmsr->watts_long = 0.0;
            pcapmsr->watts_short = 0.0;
            pcapmsr->seconds_long = 0.0;
            pcapmsr->seconds_short = 0.0;
        }

        for (msr = 0; msr < system_info->sysmsr->msr_nums[2]; msr++)
        {
            struct msr_info *imsr = &system_info->sysmsr->info_msrs[msr + (package * system_info->sysmsr->msr_nums[2])];
            imsr->msr = system_info->sysmsr->msrs[2][msr];
            imsr->package_id = package;
            imsr->cpu_id = cpu_id;
            imsr->thermal_spec_power = 0.0;
            imsr->minimum_power = -1.0;
            imsr->maximum_power = 0.0;
            imsr->maximum_time_window = 0.0;
        }

        for (msr = 0; msr < system_info->sysmsr->msr_nums[3]; msr++)
        {
            struct msr_perf *perfmsr = &system_info->sysmsr->perf_msrs[msr + package * system_info->sysmsr->msr_nums[3]];
            perfmsr->msr = system_info->sysmsr->msrs[3][msr];
            perfmsr->package_id = package;
            perfmsr->cpu_id = cpu_id;
            perfmsr->throttled_time = 0.0;
        }

        for (msr = 0; msr < system_info->sysmsr->msr_nums[4]; msr++)
        {
            struct msr_policy *polmsr = &system_info->sysmsr->policy_msrs[msr + package * system_info->sysmsr->msr_nums[4]];
            polmsr->msr = system_info->sysmsr->msrs[4][msr];
            polmsr->package_id = package;
            polmsr->cpu_id = cpu_id;
            polmsr->policy = -1;
        }
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);
}

int finalize_msrs (struct system_info_t * system_info)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    if (!system_info->sysmsr->error_state)
    {
        int package;
        for (package = 0; package < system_info->sysmsr->total_packages; package++)
            if (system_info->sysmsr->package_fd[package])
                close(system_info->sysmsr->package_fd[package]);
    }

    if(system_info->sysmsr)
        free(system_info->sysmsr);

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

int rapl_read_energy(struct rapl_energy *re, struct system_info_t * system_info)
{
    re->package = -1.0;
    re->pp0 = -1.0;
    re->pp1 = -1.0;
    re->dram = -1.0;
    re->platform = -1.0;

    if (system_info->sysmsr->error_state)
    {
        poli_log(WARNING, NULL, "RAPL Interface couldn't be set up. Energy readings are not possible.");
        return 0;
    }

    //TODO
    int i, package_id = 0; //fixing package to be 0, this needs to be changed for other platforms
    int num_energy_msrs = system_info->sysmsr->msr_nums[0];

    for (i = 0; i < num_energy_msrs; i++)
    {
        struct msr_energy *emsr = &system_info->sysmsr->energy_msrs[package_id * num_energy_msrs + i];
        read_msr_energy(emsr, system_info, 0);
        if (emsr->msr == MSR_PKG_ENERGY_STATUS)
            re->package = emsr->total_energy;
        else if (emsr->msr == MSR_PP0_ENERGY_STATUS)
            re->pp0 = emsr->total_energy;
        else if (emsr->msr == MSR_PP1_ENERGY_STATUS)
            re->pp1 = emsr->total_energy;
        else if (emsr->msr == MSR_DRAM_ENERGY_STATUS)
            re->dram = emsr->total_energy;
        else if (emsr->msr == MSR_PLATFORM_ENERGY_COUNTER)
            re->platform = emsr->total_energy;
        else
            poli_log(ERROR, NULL, "%s: Unrecognized energy MSR %#010X", __FUNCTION__, emsr->msr);
    }

    if (re->package == -1.0 && re->pp0 == -1.0 && re->pp1 == -1.0 && re->dram == -1.0 && re->platform == -1.0)
        poli_log(ERROR, NULL, "%s: wasn't able to get any energy measurments!", __FUNCTION__);

    return 0;
}

static int verify_power_limits(double watts, int enable)
{
    double minwatts = MIN_WATTS;
    if (!enable)
        minwatts = 0.0;
    if (watts > MAX_WATTS || watts < minwatts )
    {
        poli_log(ERROR, NULL, "Invalid power cap value %lf. Cannot be less than %lf and larger than %lf",
            watts, MIN_WATTS, MAX_WATTS);
        return 0;
    }

    return 1;
}

int rapl_compute_total_energy (struct rapl_energy *re, struct rapl_energy *end, struct rapl_energy *start)
{
    re->package = -1;
    re->pp0 = -1;
    re->pp1 = -1;
    re->dram = -1;
    re->platform = -1;

    double endpackage, endpp0, endpp1, enddram, endplatform = 0.0;
    if (end)
    {
        endpackage = end->package;
        endpp0 = end->pp0;
        endpp1 = end->pp1;
        enddram = end->dram;
        endplatform = end->platform;
    }

    if (endpackage > -1 && start->package > -1 && (endpackage >= start->package))
        re->package = endpackage - start->package;
    if (endpp0 > -1 && start->pp0 > -1 && (endpp0 >= start->pp0))
        re->pp0 = endpp0 - start->pp0;
    if (endpp1 > -1 && start->pp1 > -1 && (endpp1 >= start->pp1))
        re->pp1 = endpp1 - start->pp1;
    if (enddram > -1 && start->dram > -1 && (enddram >= start->dram))
        re->dram = enddram - start->dram;
    if (endplatform > -1 && start->platform > -1 && (endplatform >= start->platform))
        re->platform = endplatform - start->platform;
    return 0;
}

int rapl_compute_total_power (struct rapl_energy *re, struct rapl_energy *energy, double time)
{
    re->package = -1;
    re->pp0 = -1;
    re->pp1 = -1;
    re->dram = -1;
    re->platform = -1;
    if (energy->package > -1)
        re->package = energy->package / time;
    if (energy->pp0 > -1)
        re->pp0 = energy->pp0 / time;
    if (energy->pp1 > -1)
        re->pp1 = energy->pp1 / time;
    if (energy->dram > -1)
        re->dram = energy->dram / time;
    if (energy->platform > -1)
        re->platform = energy->platform / time;
    return 0;
}

static int get_msr_for_zone_name(char *zone_name, int get_pcap)
{
    int msr_address = -1;

    if (!strcmp(zone_name, "PACKAGE"))
    {
        if (get_pcap)
            msr_address = MSR_PKG_POWER_LIMIT;
        else
            msr_address = MSR_PKG_POWER_INFO;
    }
    else if (!strcmp(zone_name, "CORE"))
        msr_address = MSR_PP0_POWER_LIMIT;
    else if (!strcmp(zone_name, "UNCORE"))
        msr_address = MSR_PP1_POWER_LIMIT;
    else if (!strcmp(zone_name, "DRAM"))
    {
        if (get_pcap)
            msr_address = MSR_DRAM_POWER_LIMIT;
        else
            msr_address = MSR_DRAM_POWER_INFO;
    }
    else if (!strcmp(zone_name, "PLATFORM"))
        msr_address = MSR_PLATFORM_POWER_LIMIT;
    else
        poli_log(ERROR, NULL, "%s: Unsupported zone for power capping: %s", __FUNCTION__, zone_name);

    return msr_address;
}

static int rapl_init_power_cap(struct msr_pcap *pcap, char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short, int enable)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    int msr_address = get_msr_for_zone_name(zone_name, 1);

    if (msr_address == -1)
        pcap = 0;
    else
    {
        pcap->msr = msr_address;
        if (pcap->msr != MSR_PKG_POWER_LIMIT && pcap->msr != MSR_PP0_POWER_LIMIT &&
        pcap->msr != MSR_PP1_POWER_LIMIT && pcap->msr != MSR_DRAM_POWER_LIMIT &&
        pcap->msr != MSR_PLATFORM_POWER_LIMIT)
        {
            poli_log(ERROR, NULL, "%s: Invalid msr %#010X", __FUNCTION__, pcap->msr);
            return 1;
        }
        pcap->watts_long = watts_long;
        pcap->seconds_long = seconds_long;
        pcap->enabled_long = enable;
        pcap->clamped_long = enable;

        if (!short_term_supported(pcap->msr))
        {
            pcap->enabled_short = 0;
            pcap->clamped_short = 0;
            pcap->watts_short = 0.0;
            pcap->seconds_short = 0.0;
        }
        else
        {
            pcap->enabled_short = enable;
            pcap->clamped_short = enable;
            pcap->watts_short = watts_short;
            pcap->seconds_short = seconds_short;
        }
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

int rapl_set_power_cap(char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short, struct system_info_t * system_info, int enable)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    int ret = 1;

    if (system_info->sysmsr->error_state)
    {
        poli_log(WARNING, NULL, "RAPL Interface couldn't be set up. Setting power cap is not possible.");
        return ret;
    }

    if (!verify_power_limits(watts_long, enable) || !verify_power_limits(watts_short, enable))
    {
        poli_log(WARNING, NULL, "%s: The requested power cap is invalid. Will reset system to default values...", __FUNCTION__);
        if (strcmp(zone_name, "PACKAGE") == 0)
        {
            watts_long = DEFAULT_PKG_POW;
            watts_short = DEFAULT_SHORT;
            seconds_long = DEFAULT_SECONDS_LONG;
            seconds_short = DEFAULT_SECONDS_SHORT;
            enable = 1;
        }
        else if (strcmp(zone_name, "CORE") == 0)
        {
            watts_long = DEFAULT_CORE_POW;
            seconds_long = DEFAULT_SECONDS_LONG;
            enable = 0;
        }
    }

    if (enable > 1)
        enable = 1;
    else if (enable < 0)
        enable = 0;

    struct msr_pcap pcap;

    if (rapl_init_power_cap(&pcap, zone_name, watts_long, watts_short, seconds_long, seconds_short, enable) == 0)
    {
        int package = 0; //TODO figure out what to do with multiple packages
        pcap.package_id = package;
        ret = set_msr_pcap(&pcap, system_info, package);
        if (ret != 0)
            poli_log(ERROR, NULL, "Something went wrong with setting a power cap!");
    }
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    return ret;
}

static int short_term_supported (int msr)
{
    int supported = 0;
    switch (msr)
    {
        case MSR_PLATFORM_POWER_LIMIT:
            supported = 1;
            break;
        case MSR_PKG_POWER_LIMIT:
            supported = 1;
            break;
        default:
            break;
    }
    return supported;
}


int rapl_get_power_cap_info(char *zone_name, double *min, double *max,
    double *thermal_spec, double *max_time_window, struct system_info_t * system_info)
{
    if (system_info->sysmsr->error_state)
    {
        poli_log(WARNING, NULL, "RAPL Interface couldn't be set up. Getting power cap info is not possible.");
        return 1;
    }

    int msr_address = get_msr_for_zone_name(zone_name, 0);
    int ret = 0;
    if (msr_address != MSR_PKG_POWER_INFO && msr_address != MSR_DRAM_POWER_INFO)
    {
        poli_log(ERROR, NULL, "Power capping info is only allowed for PACKAGE and DRAM zones!");
        ret = 1;
    }
    else if (msr_address != -1)
    {
        int package = 0; //TODO

        uint64_t data;
        if (pread(system_info->sysmsr->package_fd[package], &data, sizeof(uint64_t), msr_address) != sizeof(uint64_t))
        {
            poli_log(ERROR, NULL, "%s: Something went wrong with getting power cap info for msr %#010X : %s", __FUNCTION__, msr_address, strerror(errno));
            ret = 1;
        }
        else
        {
            *thermal_spec = system_info->sysmsr->power_units * (double) (data & 0x7fff);
            *min = system_info->sysmsr->power_units * (double) ((data >> 16) & 0x7fff);
            *max = system_info->sysmsr->power_units * (double) ((data >> 32) & 0x7fff);
            *max_time_window = system_info->sysmsr->time_units * (double) ((data >> 48) & 0x7fff);
        }
    }

    if (ret)
    {
        *thermal_spec = -1;
        *min = -1;
        *max = -1;
        *max_time_window = -1;
    }

    return ret;
}

int rapl_get_power_cap(struct msr_pcap *pcap, char *zone_name, struct system_info_t * system_info)
{
    if (system_info->sysmsr->error_state)
    {
        poli_log(WARNING, NULL, "RAPL Interface couldn't be set up. Setting power cap is not possible.");
        return 1;
    }

    int msr_address = get_msr_for_zone_name(zone_name, 1);
    int ret = 0;
    if (msr_address != -1)
    {
        pcap->msr = msr_address;
        int package = 0; //TODO

        uint64_t data;
        if (pread(system_info->sysmsr->package_fd[package], &data, sizeof(uint64_t), pcap->msr) != sizeof(uint64_t))
        {
            poli_log(ERROR, NULL, "%s: Something went wrong with getting power cap for msr %#010X : %s", __FUNCTION__, pcap->msr, strerror(errno));
            ret = 1;
        }
        else
        {
            if (msr_address == MSR_PKG_POWER_LIMIT)
                pcap->zone_label = PACKAGE;
            else if (msr_address == MSR_PP0_POWER_LIMIT)
                pcap->zone_label = CORE;
            else if (msr_address == MSR_PP1_POWER_LIMIT)
                pcap->zone_label = UNCORE;
            else if (msr_address == MSR_PLATFORM_POWER_LIMIT)
                pcap->zone_label = PLATFORM;
            else if (msr_address == MSR_DRAM_POWER_LIMIT)
                pcap->zone_label = DRAM;

            // for now clamping will be forced to be the same as enabled bit
            uint64_t enabled_long = get_bits(data, ENABLED_LONG_START_BITS, ENABLED_LONG_END_BITS);
            pcap->enabled_long = (enabled_long == 0x3) ? 1 : 0;
            pcap->clamped_long = pcap->enabled_long;
            pcap->watts_long = get_bits(data, WATTS_LONG_START_BITS, WATTS_LONG_END_BITS) * system_info->sysmsr->power_units;
            pcap->seconds_long = from_msr_time(get_bits(data, SECONDS_LONG_START_BITS, SECONDS_LONG_END_BITS - 2),
                get_bits(data, SECONDS_LONG_END_BITS - 1, SECONDS_LONG_END_BITS), system_info->sysmsr->time_units);

            if (!short_term_supported(pcap->msr))
            {
                pcap->enabled_short = 0;
                pcap->clamped_short = 0;
                pcap->watts_short = 0;
                pcap->seconds_short = 0;
            }
            else
            {
                uint64_t enabled_short = get_bits(data, ENABLED_SHORT_START_BITS, ENABLED_SHORT_END_BITS);
                pcap->enabled_short = (enabled_short == 0x3) ? 1 : 0;
                pcap->clamped_short = pcap->enabled_short;
                pcap->watts_short = get_bits(data, WATTS_SHORT_START_BITS, WATTS_SHORT_END_BITS) * system_info->sysmsr->power_units;
                pcap->seconds_short = from_msr_time(get_bits(data, SECONDS_SHORT_START_BITS, SECONDS_SHORT_END_BITS - 2),
                    get_bits(data, SECONDS_SHORT_END_BITS - 1, SECONDS_SHORT_END_BITS), system_info->sysmsr->time_units);
            }
        }
    }

    return ret;
}

//from raplcap-msr.c
static uint64_t get_bits(uint64_t msrval, uint8_t first, uint8_t last)
{
  assert(first <= last);
  assert(last < 64);
  return (msrval >> first) & (((uint64_t) 1 << (last - first + 1)) - 1);
}

static int open_msr(int core)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    char msr_filename[BUFSIZE];
    int fd;
    int ret = 0;
    sprintf(msr_filename, "/dev/cpu/%d/msr_safe", core);
    fd = open(msr_filename, O_RDWR);
    if ( fd < 0 )
    {
        if ( errno == ENXIO )
        {
            poli_log(ERROR, NULL, "rdmsr: No CPU %d", core);
            ret = 2;
        }
        else if ( errno == EIO )
        {
            poli_log(ERROR, NULL, "rdmsr: CPU %d doesn't support MSRs", core);
            ret = 3;
        }
        else
        {
            poli_log(WARNING, NULL, "Couldn't open the msr_safe file. Trying regular msr...");
            sprintf(msr_filename, "/dev/cpu/%d/msr", core);
            fd = open(msr_filename, O_RDWR);
            if ( fd < 0)
            {
                if ( errno == ENXIO )
                {
                    poli_log(ERROR, NULL, "rdmsr: No CPU %d", core);
                    ret = 2;
                }
                else if ( errno == EIO )
                {
                    poli_log(ERROR, NULL, "rdmsr: CPU %d doesn't support MSRs", core);
                    ret = 3;
                }
                else
                {
                    poli_log(ERROR, NULL, "Failed to open %s: %s",msr_filename, strerror(errno));
                    ret = 127;
                }
            }
        }
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);
    if (!ret)
        return fd;
    return ret;
}

static void get_msr_units(struct system_info_t *system_info, int package)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);
    /* Calculate the units used */
    long long result = read_msr(system_info->sysmsr->package_fd[package], MSR_RAPL_POWER_UNIT);

    system_info->sysmsr->power_units = 1.0 / pow2_u64(result & 0xf);
    system_info->sysmsr->time_units = 1.0 / pow2_u64((result >> 16) & 0xf);

    system_info->sysmsr->cpu_energy_units[package] = 1.0 / pow2_u64((result >> 8) & 0x1f);

    /* On Haswell EP and Knights Landing */
    /* The DRAM units differ from the CPU ones */
    if ((system_info->sysmsr->cpu_model == CPU_HASWELL_EP) || (system_info->sysmsr->cpu_model == CPU_KNIGHTS_LANDING))
        system_info->sysmsr->dram_energy_units[package] = 1.0 / pow2_u64(16);
    else
        system_info->sysmsr->dram_energy_units[package] = system_info->sysmsr->cpu_energy_units[package];

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);
}

static int write_msr (int fd, int msr_address, uint64_t data)
{
    off_t msr = (off_t) msr_address;
    assert(msr >= 0);
    if (pwrite(fd, &data, sizeof(uint64_t), msr) == sizeof(uint64_t))
        return 0;
    poli_log(ERROR, NULL, "Something went wrong with writing to msr %#010X : %s", msr_address, strerror(errno));
    return 1;
}

//based on raplcap-msr.c
static int set_msr_pcap(struct msr_pcap *pcap, struct system_info_t * system_info, int package_id)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    uint64_t msrval_enable, msrval;

    if (pread(system_info->sysmsr->package_fd[package_id], &msrval, sizeof(uint64_t), pcap->msr) != sizeof(uint64_t) )
    {
        poli_log(ERROR, NULL, "%s: Couldn't read MSR at address %#010X", __FUNCTION__, pcap->msr);
        return 0;
    }

    const uint64_t enabled_long_bits = (pcap->enabled_long) ? 0x3 : 0x0;
    const uint64_t enabled_short_bits = (pcap->enabled_short) ? 0x3 : 0x0;

    msrval = replace_bits(msrval, enabled_long_bits, ENABLED_LONG_START_BITS, ENABLED_LONG_END_BITS);

    if (pcap->enabled_short)
        msrval = replace_bits(msrval, enabled_short_bits, ENABLED_SHORT_START_BITS, ENABLED_SHORT_END_BITS);

    if (write_msr(system_info->sysmsr->package_fd[package_id], pcap->msr, msrval))
        poli_log(ERROR, NULL, "%s: Something went wrong with enabling the MSR", __FUNCTION__);

    msrval = replace_bits(msrval, to_msr_power(pcap->watts_long, system_info->sysmsr->power_units), WATTS_LONG_START_BITS, WATTS_LONG_END_BITS);

    if (pcap->seconds_long > 0)
        msrval = replace_bits(msrval, to_msr_time(pcap->seconds_long, system_info->sysmsr->time_units), SECONDS_LONG_START_BITS, SECONDS_LONG_END_BITS);

    if (pcap->enabled_short && pcap->clamped_short)
    {
        msrval = replace_bits(msrval, to_msr_power(pcap->watts_short, system_info->sysmsr->power_units), WATTS_SHORT_START_BITS, WATTS_SHORT_END_BITS);
        if (pcap->seconds_short > 0)
            msrval = replace_bits(msrval, to_msr_time(pcap->seconds_short, system_info->sysmsr->time_units), SECONDS_SHORT_START_BITS, SECONDS_SHORT_END_BITS);
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return write_msr(system_info->sysmsr->package_fd[package_id], pcap->msr, msrval);
}

//from raplcap-msr.c
static uint64_t to_msr_power(double watts, double power_units)
{
    assert(watts >= 0);
    assert(power_units > 0);
    // Lower bound is 0, but upper bound is limited by what fits in 15 bits
    static const uint64_t MSR_POWER_MAX = 0x7FFF;
    uint64_t p = (uint64_t) (watts / power_units);
    if (p > MSR_POWER_MAX)
        p = MSR_POWER_MAX;
    return p;
}

//from raplcap-msr.c
static uint64_t replace_bits(uint64_t msrval, uint64_t data, uint8_t first, uint8_t last)
{
    assert(first <= last);
    assert(last < 64);
    const uint64_t mask = (((uint64_t) 1 << (last - first + 1)) - 1) << first;
    return (msrval & ~mask) | ((data << first) & mask);
}

static uint64_t log2_u64(uint64_t y)
{
  // log2(y); returns 0 for y = 0
  uint8_t ret = 0;
  while (y >>= 1) {
    ret++;
  }
  return ret;
}

//from raplcap-msr.c
static uint64_t pow2_u64(uint64_t y)
{
  // 2^y
  return ((uint64_t) 1) << y;
}

//from raplcap-msr.c
static double from_msr_time(uint64_t y, uint64_t f, double time_units)
{
  return pow2_u64(y) * ((4 + f) / 4.0) * time_units;
}

//from raplcap-msr.c
static uint64_t to_msr_time(double seconds, double time_units)
{
    assert(seconds > 0);
    assert(time_units > 0);
    // Seconds cannot be shorter than the smallest time unit - log2 would get a negative value and overflow "y".
    // They also cannot be larger than 2^2^5-1 so that log2 doesn't produce a value that uses more than 5 bits for "y".
    // Clamping prevents values outside the allowable range, but precision can still be lost in the conversion.
    static const double MSR_TIME_MIN = 1.0;
    static const double MSR_TIME_MAX = (double) 0xFFFFFFFF;
    double t = seconds / time_units;
    if (t < MSR_TIME_MIN)
        t = MSR_TIME_MIN;
    else if (t > MSR_TIME_MAX)
        t = MSR_TIME_MAX;

    // y = log2((4*t)/(4+f)); we can ignore "f" since t >= 1 and 0 <= f <= 3; we can also drop the real part of "t"
    const uint64_t y = log2_u64((uint64_t) t);
    // f = (4*t)/(2^y)-4; the real part of "t" only matters for t < 4, otherwise it's insignificant in computing "f"
    const uint64_t f = (((uint64_t) (4 * t)) / pow2_u64(y)) - 4;
    return ((y & 0x1F) | ((f & 0x3) << 5));
}


static long long read_msr (int fd, int msr_address)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    uint64_t data;
    if ( pread(fd, &data, sizeof(uint64_t), msr_address) != sizeof(uint64_t) )
    {
        poli_log(ERROR, NULL, "Couldn't read MSR at address %#010X: %s", msr_address, strerror(errno));
        return -1;
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return (long long)data;
}


static int read_msr_info (struct msr_info *msr_info, struct system_info_t *system_info, int package_id)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    if ((msr_info->msr != MSR_PKG_POWER_INFO) && (msr_info->msr != MSR_DRAM_POWER_INFO))
    {
        poli_log(ERROR, NULL, "%s: The requested msr at address %#010X is not valid!", __FUNCTION__, msr_info->msr);
        return 1;
    }

    long long result = read_msr(system_info->sysmsr->package_fd[package_id], msr_info->msr);
    if (!result)
    {
        poli_log(ERROR, NULL, "Something went wrong with reading the info msr %#010X", msr_info->msr);
        return 1;
    }
    msr_info->thermal_spec_power = system_info->sysmsr->power_units * ((double) (result & 0x7fff));
    msr_info->minimum_power = system_info->sysmsr->power_units * ((double) ((result >> 16) & 0x7fff));
    msr_info->maximum_power = system_info->sysmsr->power_units * ((double) ((result >> 32) & 0x7fff));
    msr_info->maximum_time_window = system_info->sysmsr->time_units * ((double) ((result >> 48) & 0x7fff));

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);
    return 0;
}

static int read_msr_pcap (struct msr_pcap *msr_pcap, struct system_info_t *system_info, int package_id)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    if ((msr_pcap->msr != MSR_PKG_POWER_LIMIT) && (msr_pcap->msr != MSR_PP0_POWER_LIMIT)
        && (msr_pcap->msr != MSR_PP1_POWER_LIMIT) && (msr_pcap->msr != MSR_DRAM_POWER_LIMIT) &&
        (msr_pcap->msr != MSR_PLATFORM_POWER_LIMIT))
    {
        poli_log(ERROR, NULL, "%s: The requested msr at address %#010X is not valid!", __FUNCTION__, msr_pcap->msr);
        return 1;
    }

    long long result = read_msr(system_info->sysmsr->package_fd[package_id], msr_pcap->msr);
    if (result < 0)
    {
        poli_log(ERROR, NULL, "%s: Something went wrong with reading the power limit msr %#010X", __FUNCTION__, msr_pcap->msr);
        return 1;
    }

    msr_pcap->watts_long = system_info->sysmsr->power_units * (double) ((result >> WATTS_LONG_START_BITS) & 0x7FFF);
    msr_pcap->seconds_long = system_info->sysmsr->time_units * (double) ((result >> SECONDS_LONG_START_BITS) & 0x007F);
    msr_pcap->enabled_long = (int) (result & (1LL << ENABLED_LONG_START_BITS));
    msr_pcap->clamped_long = (int) (result & (1LL << ENABLED_LONG_END_BITS));

    if (short_term_supported(msr_pcap->msr))
    {
        msr_pcap->watts_short = system_info->sysmsr->power_units * (double) ((result >> WATTS_SHORT_START_BITS) & 0x7FFF);
        msr_pcap->seconds_short = system_info->sysmsr->time_units * (double) ((result >> SECONDS_SHORT_START_BITS) & 0x007F);
        msr_pcap->enabled_short = (int) (result & (1LL << ENABLED_SHORT_START_BITS));
        msr_pcap->clamped_short = (int) (result & (1LL << ENABLED_SHORT_END_BITS));
    }
    else
    {
        msr_pcap->watts_short = 0.0;
        msr_pcap->seconds_short = 0.0;
        msr_pcap->enabled_short = 0.0;
        msr_pcap->clamped_short = 0.0;
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);
    return 0;
}

static int read_msr_perf (struct msr_perf *msr_perf, struct system_info_t *system_info, int package_id)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    if ((msr_perf->msr != MSR_PKG_PERF_STATUS) && (msr_perf->msr != MSR_DRAM_PERF_STATUS) &&
        (msr_perf->msr != MSR_PP0_PERF_STATUS))
    {
        poli_log(ERROR, NULL, "%s: The requested msr at address %#010X is not valid!", __FUNCTION__, msr_perf->msr);
        return 1;
    }

    long long result = read_msr(system_info->sysmsr->package_fd[package_id], msr_perf->msr);
    if (!result)
    {
        poli_log(ERROR, NULL, "Something went wrong with reading the perf msr %#010X", msr_perf->msr);
        return 1;
    }
    msr_perf->throttled_time = ((double) result) * system_info->sysmsr->time_units;

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

static int read_msr_policy (struct msr_policy *msr_policy, struct system_info_t *system_info, int package_id)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    if ((msr_policy->msr != MSR_PP0_POLICY) && (msr_policy->msr != MSR_PP1_POLICY))
    {
        poli_log(ERROR, NULL, "%s: The requested msr at address %d is not valid!", __FUNCTION__, msr_policy->msr);
        return 1;
    }
    long long result = read_msr(system_info->sysmsr->package_fd[package_id], msr_policy->msr);
    if (!result)
    {
        poli_log(ERROR, NULL, "Something went wrong with reading the power limit msr");
        return 1;
    }
    msr_policy->policy = (int) result & 0x001f;

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

static int read_msr_energy (struct msr_energy *msr_energy, struct system_info_t * system_info, int package_id)
{
    if ((msr_energy->msr != MSR_PKG_ENERGY_STATUS) && (msr_energy->msr != MSR_PP0_ENERGY_STATUS) &&
        (msr_energy->msr != MSR_PP1_ENERGY_STATUS) && (msr_energy->msr != MSR_DRAM_ENERGY_STATUS) &&
        (msr_energy->msr != MSR_PLATFORM_ENERGY_COUNTER))
    {
        poli_log(ERROR, NULL, "%s: The requested msr at address %#010X is not valid!", __FUNCTION__, msr_energy->msr);
        return 1;
    }
    int fd = system_info->sysmsr->package_fd[package_id];
    uint64_t data;
    if (pread(fd, &data, sizeof(uint64_t), msr_energy->msr) != sizeof(uint64_t))
    {
        poli_log(ERROR, NULL, "%s: Couldn't read MSR at address %#010X", __FUNCTION__, msr_energy->msr);
        return 0;
    }

    data &= 0xFFFFFFFF;

    if (data < msr_energy->last_energy)
        msr_energy->num_overflows++;

    msr_energy->last_energy = (double) data;
    if (msr_energy->msr == MSR_DRAM_ENERGY_STATUS)
        msr_energy->total_energy = (double) ((data + msr_energy->num_overflows * (uint64_t) UINT32_MAX) * msr_energy->dram_energy_units);
    else
        msr_energy->total_energy = (double) ((data + msr_energy->num_overflows * (uint64_t) UINT32_MAX) * msr_energy->cpu_energy_units);

    return 0;
}

static int detect_cpu(void)
{

    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    int family, model = -1;
    char buffer[BUFSIZE], vendor[BUFSIZE];

    FILE *cpuinfo;
    cpuinfo = fopen("/proc/cpuinfo", "r");
    if (cpuinfo==NULL)
    {
        poli_log(ERROR, NULL, "Couldn't access cpuinfo! %s", strerror(errno));
        return -1;
    }

    int done = 0;
    int verified_vendor, verified_cpufam, model_found = 0;

    while ((fgets(buffer,BUFSIZE,cpuinfo) != NULL) && !done)
    {
        char *token, *value;
        token = strtok(buffer, "\t:");
        value = strtok(NULL, "\t:");

        if (!verified_vendor)
        {
            if (!strncmp(token,"vendor_id", 8))
            {
                char *s = strtok(value, " ");
                if (strncmp(s,"GenuineIntel", 12))
                {
                    poli_log(ERROR, NULL, "%s not an Intel chip",value);
                    done = 1;
                }
                else
                    verified_vendor = 1;
            }
        }
        if (!verified_cpufam)
        {
            if (!strncmp(token,"cpu family",10))
            {
                sscanf(value, "%d",&family);
                if (family != 6)
                {
                    poli_log(ERROR, NULL, "Wrong CPU family %d",family);
                    done = 1;
                }
                else
                    verified_cpufam = 1;
            }
        }
        if (!model_found)
        {
            if (!strncmp(token,"model",5))
            {
                sscanf(value,"%d",&model);
                model_found = 1;
            }
        }
        if (verified_vendor && verified_cpufam && model_found && !done)
            done = 1;
        else if (done && !verified_vendor && !verified_cpufam)
        {
            poli_log(ERROR, NULL, "Wrong CPU family or not an Intel chip. RAPL Interface won't be used.");
            return -1;
        }
    }

    fclose(cpuinfo);

    if (!verify_model(model))
    {
        poli_log(ERROR, NULL, "Your CPU is not currently supported.RAPL Interface won't be used.");
        return -1;
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return model;
}

static int verify_model(int model)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    switch(model)
    {
        case CPU_SANDYBRIDGE:
            break;
        case CPU_SANDYBRIDGE_EP:
            break;
        case CPU_IVYBRIDGE:
            break;
        case CPU_IVYBRIDGE_EP:
            break;
        case CPU_HASWELL:
            break;
        case CPU_HASWELL_EP:
            break;
        case CPU_BROADWELL:
            break;
        case CPU_BROADWELL_EP:
            break;
        case CPU_SKYLAKE:
            break;
        case CPU_SKYLAKE_HS:
            break;
        case CPU_KABYLAKE:
            break;
        case CPU_KABYLAKE_2:
            break;
        case CPU_KNIGHTS_LANDING:
            break;
        default:
            poli_log(ERROR, NULL, "Unsupported model %d",model);
            model = 0;
            break;
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return model;
}

static int detect_packages (struct system_info_t *system_info)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    char filename[BUFSIZE];
    FILE *package_id_file;
    int package;
    int i;

    /* Initalize the array */
    for (i = 0; i < MAX_PACKAGES; i++)
        system_info->sysmsr->package_map[i] = -1;

    system_info->sysmsr->total_packages = 0;
    system_info->sysmsr->total_cores = 0;

    for (i = 0; i < MAX_CPUS; i++)
    {
        sprintf(filename, "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", i);
        package_id_file = fopen(filename, "r");
        if (package_id_file == NULL) break;
        if (fscanf(package_id_file, "%d", &package) < 1)
        {
            poli_log(ERROR, NULL, "Was not able to read package ID: %s! Setting package ID to 0.", strerror(errno));
            package = 0;
        }

        fclose(package_id_file);

        if (system_info->sysmsr->package_map[package] == -1)
        {
            system_info->sysmsr->total_packages++;
            system_info->sysmsr->package_map[package] = i;
        }
    }
    system_info->sysmsr->total_cores = i;

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}
