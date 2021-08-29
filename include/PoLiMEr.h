#ifndef __POLIMER_H
#define __POLIMER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <sys/types.h>
#include <sys/time.h>
#include <stdint.h>
#include <signal.h>
#include <time.h>

#ifndef _NOMPI
#include "mpi.h"
#endif

#include "msr-handler.h"

#ifdef _CRAY
#include "cray_pm-handler.h"
#endif

#ifdef _BGQ
#include "bgq-handler.h"
#endif


// Maximum number of user-specified tags
#define MAX_TAGS     10000
// Maximum number of polling records
#define MAX_POLL_SAMPLES 500000
#define POLL_INTERVAL 0.5
#define INITIAL_TIMER_DELAY 100000

struct monitor_t {
    int imonitor;
    int world_rank;
    int world_size;
    int node_rank;
    int node_size;
    int color;
    int num_threads;
    int thread_id;
    char *jobid;
#ifndef _NOMPI
    MPI_Comm mynode_comm;
    char my_host[MPI_MAX_PROCESSOR_NAME];
#else
    char *my_host;
#endif
};

struct poller_t {
    int time_counter;
#ifdef _BENCH
    int time_counter_em;
#endif
#ifndef _TIMER_OFF
    struct sigaction sa;
    struct itimerval timer;
    volatile int timer_on;
#endif
};


struct energy_reading {
  struct rapl_energy rapl_energy;
#ifdef _CRAY
  struct cray_measurement cray_meas;
#elif _BGQ
  struct bgq_measurement bgq_meas;
#endif
};

struct poli_tag {
    int id;
    char *tag_name;
    int monitor_id;
    int monitor_rank;

    struct energy_reading start_energy;
    struct energy_reading end_energy;
    struct energy_reading total_energy;
    struct energy_reading total_power;

    double start_time;
    double end_time;
    struct timeval start_timestamp;
    struct timeval end_timestamp;
    int start_timer_count;
    int end_timer_count;
    int closed;
};

typedef enum pcap_flags { DEFAULT, USER_SET, SYSTEM_RESET, INTERNAL, INITIAL } pcap_flag_t;

struct pcap_tag {
    int id;
    int monitor_id;
    int monitor_rank;
    char zone[ZONE_NAME_LEN];
    double watts_long;
    double watts_short;
    double seconds_long;
    double seconds_short;
    int enabled;
    double wtime;
    struct timeval timestamp;
    pcap_flag_t pcap_flag; //to have some idea if system reset, user set or controlled by library
    struct poli_tag active_poli_tags[MAX_TAGS]; //stores all active tags
    int num_active_poli_tags;
    int start_timer_count;
};

struct pcap_info {
    int monitor_id;
    int monitor_rank;
    zone_label_t zone_label;
    char zone[ZONE_NAME_LEN];
    int enabled_long;
    int enabled_short;
    int clamped_long;
    int clamped_short;
    double watts_long;
    double watts_short;
    double seconds_long;
    double seconds_short;
    double thermal_spec;
    double min;
    double max;
    double max_time_window;
};


struct frequency {
    double freq;
#ifdef _CRAY
    double cray_freq;
#endif
};

struct system_poll_info {
    int counter;
    double pkg_pcap;
    struct pcap_info pcap_info_list[NUM_ZONES];
    struct energy_reading last_energy;
    struct energy_reading current_energy;
    struct energy_reading computed_power;
    struct frequency freq;

    double wtime;
    double poll_iter_time;
};

struct system_info_t {

    /* record initial start times */
    struct timeval initial_start_time;
    double initial_mpi_wtime;

    /*record initial energy*/
    struct energy_reading initial_energy;
    struct energy_reading final_energy;

    struct poli_tag *poli_tag_list;
    int poli_opentag_tracker;
    int poli_closetag_tracker;
    struct pcap_tag *pcap_tag_list;
    struct pcap_info *current_pcap_list; //stores PACKAGE, CORE, DRAM in that order

#ifndef _TIMER_OFF
    struct system_poll_info *system_poll_list;
#endif
#ifdef _BENCH
    struct system_poll_info *system_poll_list_em;
#endif

    int num_poli_tags;
    int num_open_tags;
    int num_closed_tags;
    int num_pcap_tags;

    int cur_freq_file;

    /* add all system-dependent structs here*/
    struct system_msr_info *sysmsr;
#ifdef _CRAY
    struct system_cray_info *syscray;
#endif
};

/******************************************************************************/
/*                      MODIFIED                                              */
/******************************************************************************/

// energy reading memory operations
struct energy_reading *poli_malloc_energy_reading(void);
void poli_free_energy_reading(struct energy_reading *er);

/* poli_get_current_energy - gets the current energy reading 
   returns: 0 if no errors, 1 otherwise*/
int poli_get_current_energy(struct energy_reading *current_energy);

/* poli_get_current_power - gets the current energy power
   returns: 0 if no errors, 1 otherwise*/
int poli_get_current_power(struct energy_reading *current_power);

/******************************************************************************/
/*                      INITIALIZATION                                        */
/******************************************************************************/
/* emom_init - initalizes the necessary MPI and Energymon environments. must be called after MPI_Init() in application
   returns: 0 if everything went well*/
int poli_init(void);

/*                      END OF INITIALIZATION                                 */

/******************************************************************************/
/*                      EMON TAGS                                             */
/******************************************************************************/

/* start_poli_tag - starts an poli tag
   input: tag name
   returns: 0 if no errors, 1 otherwise*/
int poli_start_tag(char *tag_name);

/* end_poli_tag - ends an poli tag
   input: name of tag to end
   returns: 0 if no errors, 1 otherwises*/
int poli_end_tag(char *tag_name);


/*                      END OF EMON TAGS                                      */


/******************************************************************************/
/*                     SETTING POWER CAPS                                     */
/******************************************************************************/

/* poli_set_power_cap - sets a general power cap (uses PACKAGE zone as default)
   input: desired power in watts
   returns: 0 if no errors, 1 otherwise*/
int poli_set_power_cap (double watts);

int poli_set_power_cap_with_params(char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short);

/* poli_reset_system - resets the system power caps to default values, which also creates a pcap tag
   returns: 0 if no errors, 1 otherwise*/
int poli_reset_system(void);

/*                    END OF SETTING POWER CAPS                               */

/******************************************************************************/
/*                    GETTING POWER CAPS                                      */
/******************************************************************************/

/* poli_get_power_cap - returns the package power cap as default
   input: pointer to double holding the resulting watts
   returns: 0 if successful, -1 otherwise*/
int poli_get_power_cap (double *watts);

/* poli_get_power_cap - returns power cap (watts_long) last recorded by a process (not executing system call for better performance)
   input: the name of the zone for which power cap is requested, pointer to double holding the resulting watts
   returns: 0 if successful, -1 otherwise*/
int poli_get_power_cap_for_param (char *zone_name, char *param, double *result);

int poli_get_power_cap_limits (char *zone_name, double *min, double *max);

void poli_print_power_cap_info (void);
/* print_power_cap_info - prints all information related to currently recorded power caps*/
void poli_print_power_cap_info_verbose (void);


/*               END OF GETTING POWER CAPS                                    */

/******************************************************************************/
/*              FREQUENCY                                                     */
/******************************************************************************/

int poli_get_current_frequency (double *freq);
int poli_print_frequency_info (void);

/*               END OF FREQUENNCY                                            */

/******************************************************************************/
/*              HELPERS                                                       */
/******************************************************************************/

#ifndef _NOMPI
int poli_am_monitor (void);
int poli_get_monitor_id (int *id);
int poli_get_monitor_node (int *node_id);
int poli_get_subcommunicator (MPI_Comm *comm);
/* get_num_monitors - returns the total number of monitors */
int poli_get_num_monitors (void);
#endif

int coordsToInt (int *coords, int dim);


/*                  END OF HELPERS                                            */

/******************************************************************************/
/*                  FINALIZING/CLEANUP                                        */
/******************************************************************************/

/* poli_finalize - manages output of results, and performs general cleanup
   returns: 0*/
int poli_finalize(void);
/*                   END OF FINALIZING/CLEANUP                                */

#ifdef __cplusplus
}
#endif

#endif
