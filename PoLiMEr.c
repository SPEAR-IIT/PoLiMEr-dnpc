#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <math.h>

#ifndef _NOMPI
#include <mpi.h>
#endif
#ifndef _NOOMP
#include <omp.h>
#endif

#ifdef _PMI
#include <pmi.h>
#endif

#include "PoLiMEr.h"
#include "PoLiLog.h"

#include "msr-handler.h"

#ifdef _CRAY
#include "cray_pm-handler.h"
#endif

#ifdef _BGQ
#include "bgq-handler.h"
#endif

struct monitor_t *monitor = 0;
struct poller_t *poller = 0;

//the main struct holding the entire system states used througout the program
struct system_info_t *system_info = 0;

/* THIS IS ONLY FOR KNL! */
//we need to know what the power capping zones are
char *zone_names[3] = {"PACKAGE", "CORE", "DRAM"};
//for easier string manipulation
int zone_names_len[3] = {7,4,4};

static void init_system_info (void);

#ifndef _NOMPI
static int get_comm_split_color (struct monitor_t * monitor);
#ifdef _PMI
static int get_comm_split_color_pmi (struct monitor_t * monitor);
#endif
static int get_comm_split_color_hostname (struct monitor_t * monitor);
#endif

static void init_power_interfaces (struct system_info_t * system_info);
static void finalize_power_interfaces (struct system_info_t * system_info);

static struct poli_tag *get_poli_tag_for_end_time_counter(int counter);
static int start_poli_tag_no_sync (char *tag_name);
static int end_poli_tag_no_sync (char *tag_name);
static struct poli_tag *find_poli_tag_for_name (char *tag_name);
static struct poli_tag *get_poli_tag_for_start_time_counter(int counter);
static int end_existing_poli_tag (struct poli_tag *this_poli_tag);

static int init_pcap_tag (char *zone, double watts_long, double watts_short, double seconds_long, double seconds_short, pcap_flag_t pcap_flag);
static struct pcap_tag *get_pcap_for_time_counter(int counter);

/* get_system_power_cap_for_zone - returns power cap (watts_long) after executing a system call for the specified zone
   input: name of zone requested
   returns: the power in watts*/
static int get_system_power_cap_for_zone (int zone_index);
static int get_system_power_caps (void);

#ifndef _TIMER_OFF
static int setup_timer (void);
static int stop_timer (void);
static void timer_handler (int signum);
#endif

static int compute_current_power(struct system_poll_info * info, double time, struct system_info_t * system_info);
static int get_current_frequency (struct system_poll_info * info);
static int read_cpufreq (double *freq);

static void poli_sync (void);
static void poli_sync_node (void);
static double get_time (void);
static struct energy_reading read_current_energy (struct system_info_t * system_info);
static void get_timestamp(double time_from_start, char *time_str_buffer, size_t buff_len);
static FILE * open_file (char *filename);
int coordsToInt (int *coords, int dim);

static void init_energy_reading (struct energy_reading *reading);

/* get_zone_index - returns the index in zone_names[] of a zone
   input: the zone name
   returns: the index, or -1 if error*/
static int get_zone_index (char *zone_name);

static int file_handler (void);
static int poli_tags_to_file (void);
static int pcap_tags_to_file (void);
static int polling_info_to_file (void);

/******************************************************************************/
/*                      MODIFIED                                              */
/******************************************************************************/
static void copy_energy_reading(struct energy_reading *to, struct energy_reading *from)
{
  to->rapl_energy.package = from->rapl_energy.package;
  to->rapl_energy.dram = from->rapl_energy.dram;
  to->rapl_energy.pp0 = from->rapl_energy.pp0;
  to->rapl_energy.pp1 = from->rapl_energy.pp1;
  to->rapl_energy.platform = from->rapl_energy.platform;
}

int poli_get_current_energy(struct energy_reading *current_energy)
{
  struct energy_reading energy;
  energy = read_current_energy(system_info);
  copy_energy_reading(current_energy, &(energy));

  return 0;
}

#ifndef _TIMER_OFF
int poli_get_current_power(struct energy_reading *current_power)
{
  if (poller->timer_on && monitor->imonitor) {
    if (poller->time_counter < MAX_POLL_SAMPLES) {
      struct system_poll_info *info = &system_info->system_poll_list[poller->time_counter];
      struct energy_reading last_energy;
      struct energy_reading energy;
      struct energy_reading power;
      struct rapl_energy diff;
      double time;

      energy = read_current_energy(system_info); 

      if (poller->time_counter > 0) {
        last_energy = system_info->system_poll_list[poller->time_counter-1].current_energy;
        //last_energy = system_info->system_poll_list[poller->time_counter].current_energy;
        // may not correspond to a polling time
        //time = ((double) POLL_INTERVAL);
        // this should be more accurate
        time = get_time() - system_info->system_poll_list[poller->time_counter-1].wtime;
      } else {
        init_energy_reading(&last_energy);
        time = system_info->initial_mpi_wtime;
      }

      rapl_compute_total_energy(&(diff), &(energy.rapl_energy), &(last_energy.rapl_energy));
      rapl_compute_total_power(&(power.rapl_energy), &(diff), time);
      copy_energy_reading(current_power, &(power));

      return 0;
    }
  }

  return 1;
}
#endif

/******************************************************************************/
/*                      INITIALIZATION                                        */
/******************************************************************************/

int poli_init (void)
{
    monitor = malloc(sizeof(struct monitor_t));
#ifndef _NOMPI
    // get current MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &monitor->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &monitor->world_rank);

    // determine host and color based on host
    memset(monitor->my_host, '\0', sizeof(monitor->my_host));
    get_comm_split_color(monitor);

    // set up new MPI environment
    MPI_Comm_split(MPI_COMM_WORLD, monitor->color, monitor->world_rank, &monitor->mynode_comm);
    MPI_Comm_size(monitor->mynode_comm, &monitor->node_size);
    MPI_Comm_rank(monitor->mynode_comm, &monitor->node_rank);
#else

    monitor->world_size = 1;
    monitor->node_size = 0;
    monitor->color = 0;
    monitor->world_rank = 0;
    monitor->node_rank = 0;
/*
#ifdef _COBALT
    monitor->my_host = getenv("COBALT_PARTNAME");
#else
    monitor->my_host = getenv("HOSTNAME");
#endif
    if (monitor->my_host == NULL)
        strcpy(monitor->my_host, "host");
*/
    monitor->my_host = getenv("HOSTNAME");
    if (monitor->my_host == NULL)
        strcpy(monitor->my_host, "host");

#endif // end of _NOMPI

#ifndef _NOOMP
    monitor->num_threads = omp_get_num_threads();
    monitor->thread_id = omp_get_thread_num();
#else
    monitor->num_threads = 0;
    monitor->thread_id = 0;
#endif

    if (monitor->node_rank == 0 && monitor->thread_id == 0)
        monitor->imonitor = 1;
/*
#ifdef _COBALT
    monitor->jobid = getenv("COBALT_JOBID");
#else
    monitor->jobid = getenv("PoLi_JOBNAME");
    if (monitor->jobid == NULL)
        strcpy(monitor->jobid, "myjob");
#endif
*/
    monitor->jobid = getenv("PoLi_JOBNAME");
    if (monitor->jobid == NULL)
        strcpy(monitor->jobid, "myjob");

    if (monitor->imonitor)
    {
        poller = malloc(sizeof(struct poller_t));
        poller->time_counter = 0;

        //initialize the main struct
        init_system_info();
        //initialize all power monitoring and control interfaces
        init_power_interfaces(system_info);

        // record start time
        system_info->initial_mpi_wtime = get_time();
        gettimeofday(&system_info->initial_start_time, NULL);

        // may not be relevant but here it is assumed that the system is at default settings
        if (get_system_power_caps() != 0)
            poli_log(ERROR, monitor, "Couldn't get power caps on init!");

        start_poli_tag_no_sync("application_summary");
        // record energy
        system_info->initial_energy = read_current_energy(system_info);
    }

    poli_sync();

#ifndef _TIMER_OFF
    //setup and start timer
    if (monitor->imonitor)
        setup_timer();
#endif

    poli_log(TRACE, monitor,   "Finishing %s\n", __FUNCTION__);

    return 0;
}

static void init_system_info (void)
{
    //allocate memory and set dummy values
    system_info = malloc(sizeof(struct system_info_t));

    system_info->poli_tag_list = 0;
    system_info->pcap_tag_list = 0;
    system_info->current_pcap_list = 0;

#ifndef _TIMER_OFF
    system_info->system_poll_list = 0;
#endif

    system_info->num_poli_tags = 0;
    system_info->num_open_tags = 0;
    system_info->num_closed_tags = 0;

    system_info->poli_opentag_tracker = -1;
    system_info->poli_closetag_tracker = -1;

    system_info->num_pcap_tags = 0;

    // allocate list of poli tags (power measurements)
    system_info->poli_tag_list = calloc(MAX_TAGS, sizeof(struct poli_tag));

    // allocate list of power cap tags (create a tag each time we specifically set a power cap)
    system_info->pcap_tag_list = calloc(MAX_TAGS, sizeof(struct pcap_tag));

    //
    system_info->current_pcap_list = calloc(NUM_ZONES, sizeof(struct pcap_info));

#ifndef _TIMER_OFF
    //allocate list keeping the poll info
    system_info->system_poll_list = calloc(MAX_POLL_SAMPLES, sizeof(struct system_poll_info));
#endif

#ifdef _BENCH
    system_info->system_poll_list_em = calloc(MAX_POLL_SAMPLES, sizeof(struct system_poll_info));
#endif

    char *freq_path = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq";
    system_info->cur_freq_file = open(freq_path, O_RDONLY);

    if (system_info->cur_freq_file < 0)
    {
        poli_log(ERROR, monitor,   "Failed to open file to read frequency at %s!\n Error code: %s\n Trying cpuinfo_cur_frequency...", freq_path, strerror(errno));
        system_info->cur_freq_file = open("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_cur_freq", O_RDONLY);
        if (system_info->cur_freq_file < 0)
            poli_log(ERROR, monitor,   "Unable to access frequency on your system. Make sure you have read permission on cpuinfo_cur_freq.\n Error code: %s", strerror(errno));
    }
}


static void init_power_interfaces (struct system_info_t * system_info)
{
    //initialize the msr environment to read from/write to msrs
    init_msrs(system_info);
#ifdef _CRAY
    //initialize the cray environment to read the power monitoring counters
    init_cray_pm_counters (system_info);
#endif
}

#ifndef _NOMPI
static int get_comm_split_color (struct monitor_t * monitor)
{
#ifdef _PMI
    return get_comm_split_color_pmi (monitor);
#elif _BGQ
    return get_comm_split_color_bgq (monitor);
#else
    return get_comm_split_color_hostname (monitor);
#endif
    return 0;
}

#ifdef _PMI
/* This does not work when on the debug queue*/
static int get_comm_split_color_pmi (struct monitor_t * monitor)
{
    // identify which rank is on what node
    int coord[4];
    pmi_mesh_coord_t xyz;
    int nid;
    PMI_Get_nid(monitor->world_rank, &nid);
    PMI_Get_meshcoord((pmi_nid_t) nid, &xyz);
    coord[0] = xyz.mesh_x;
    coord[1] = xyz.mesh_y;
    coord[2] = xyz.mesh_z;
    // the nid helps to identify one of the four nodes hosted on a Aries router.
    coord[3] = nid;

    // set sub communicator color to nid
    monitor->color = coord[3];

    snprintf(monitor->my_host, sizeof(monitor->my_host), "%d", monitor->color);
    return 0;
}
#else
#ifndef _BGQ
static int get_comm_split_color_hostname (struct monitor_t * monitor)
{
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int resultlen;
    MPI_Get_processor_name(hostname, &resultlen);

    //the following is a hack for JLSE
    //hostnames should be of the form knlxx.ftm.alcf.anl.gov

    char *token;
    token = strtok(hostname, " \t.\n");
    strncpy(monitor->my_host, token, 5);
    char *hostnum = strtok(monitor->my_host, "knl");
    char clean_hostnum[strlen(hostnum)];
    memset(clean_hostnum, '\0', sizeof(clean_hostnum));
    int i;
    int prev_zero = 0;
    int j = 0;
    for (i = 0; i < strlen(hostnum); i++)
    {
        if (i == 0)
        {
            if (hostnum[i] == '0')
                prev_zero = 1;
        }

        if (hostnum[i] != '0')
        {
            prev_zero = 0;
            clean_hostnum[j] = hostnum[i];
            j++;
        }
        else
        {
            if (!prev_zero)
            {
                clean_hostnum[j] = hostnum[i];
                j++;
            }
        }
    }

    monitor->color = atoi(clean_hostnum);

    return 0;
}
#endif
#endif
#endif

/*                          END OF INITIALIZATION                             */

/******************************************************************************/
/*                      EMON TAGS                                             */
/******************************************************************************/

static struct poli_tag *find_poli_tag_for_name (char *tag_name)
{
    int tag_num;
    for (tag_num = 0; tag_num < system_info->num_poli_tags; tag_num++)
    {
        struct poli_tag *this_poli_tag = &system_info->poli_tag_list[tag_num];
        if (strcmp(this_poli_tag->tag_name, tag_name) == 0 && !this_poli_tag->closed) {
            return this_poli_tag;
        }
    }
    return 0;
}

static struct poli_tag *get_poli_tag_for_start_time_counter(int counter)
{
    if (monitor->imonitor)
    {
        int i;
        for (i = 0; i < system_info->num_poli_tags; i++)
        {
            struct poli_tag *this_poli_tag = &system_info->poli_tag_list[i];
            if (this_poli_tag->start_timer_count == counter)
                return this_poli_tag;
        }
    }
    return 0;
}

static struct poli_tag *get_poli_tag_for_end_time_counter(int counter)
{
    if (monitor->imonitor)
    {
        int i;
        for (i = 0; i < system_info->num_poli_tags; i++)
        {
            struct poli_tag *this_poli_tag = &system_info->poli_tag_list[i];
            if (this_poli_tag->end_timer_count == counter)
                return this_poli_tag;
        }
    }
    return 0;
}

int poli_start_tag (char *tag_name)
{
    poli_sync_node();
    poli_log(TRACE, monitor,   "Entering %s %s\n", __FUNCTION__, tag_name);
    int ret = start_poli_tag_no_sync(tag_name);
    poli_log(TRACE, monitor,   "Finishing %s %s\n", __FUNCTION__, tag_name);
    return ret;
}

static int start_poli_tag_no_sync (char *tag_name)
{
    if (monitor->imonitor)
    {
        int num_tags = system_info->num_poli_tags;
        struct poli_tag *new_poli_tag = &system_info->poli_tag_list[num_tags];
        new_poli_tag->id = num_tags;
        system_info->poli_opentag_tracker = num_tags;
        new_poli_tag->tag_name = tag_name;
        new_poli_tag->monitor_id = monitor->color;
        new_poli_tag->monitor_rank = monitor->world_rank;

        new_poli_tag->start_energy = read_current_energy(system_info);

        new_poli_tag->start_time = get_time();
        gettimeofday(&(new_poli_tag->start_timestamp), NULL);
        new_poli_tag->start_timer_count = poller->time_counter;

        system_info->num_poli_tags++;
        system_info->num_open_tags++;
    }
    return 0;
}

int poli_end_tag (char *tag_name)
{
    poli_sync_node();
    poli_log(TRACE, monitor,   "Entering %s %s\n", __FUNCTION__, tag_name);
    int ret = end_poli_tag_no_sync(tag_name);
    poli_log(TRACE, monitor,   "Finishing %s %s\n", __FUNCTION__, tag_name);
    return ret;
}

static int end_poli_tag_no_sync (char *tag_name)
{
    int ret = 0;
    if (monitor->imonitor)
    {
        struct poli_tag *this_poli_tag = &system_info->poli_tag_list[system_info->poli_opentag_tracker];
        if (system_info->poli_opentag_tracker < 1) //1 because we don't want application_summary to be prematurely closed
        {
            if (strcmp(this_poli_tag->tag_name, tag_name) != 0)
            {
                poli_log(WARNING, monitor, "You attempted to close tag %s, but no tags were opened! This tag will be omitted", tag_name);
                return -1;
            }
        }

        //if (strcmp(this_poli_tag->tag_name, tag_name) != 0) //this enables interleaving tags
        //    this_poli_tag = find_poli_tag_for_name(tag_name);
        if (this_poli_tag == 0)
        {
            poli_log(ERROR, monitor,   "%s: Failed to find poli tag to end: %s\n", __FUNCTION__, tag_name);
            return 1;
        }
        ret = end_existing_poli_tag(this_poli_tag);
    }
    return ret;
}

/*note: having a barrier here can cause problems if we have lots of ranks per node..*/
static int end_existing_poli_tag (struct poli_tag *this_poli_tag)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor,   "Entering %s", __FUNCTION__);

        this_poli_tag->end_energy = read_current_energy(system_info);

        this_poli_tag->end_time = get_time();
        gettimeofday(&(this_poli_tag->end_timestamp), NULL);
        this_poli_tag->end_timer_count = poller->time_counter;
        this_poli_tag->closed = 1;
        system_info->poli_closetag_tracker = system_info->poli_opentag_tracker;
        system_info->num_closed_tags--; //yes, decrement
        system_info->poli_closetag_tracker = system_info->poli_opentag_tracker;
        system_info->poli_opentag_tracker--;

        poli_log(TRACE, monitor, "Finishing %s", __FUNCTION__);
    }
    return 0;
}

/*                      END OF EMON TAGS                                      */

/******************************************************************************/
/*                     POWER CAP TAGS                                         */
/******************************************************************************/

static int init_pcap_tag (char *zone, double watts_long, double watts_short, double seconds_long, double seconds_short, pcap_flag_t pcap_flag)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor, "Entering %s", __FUNCTION__);

        /* check if requested zone is an allowed parameter */
        int i = get_zone_index(zone);
        if ( i < 0)
        {
            poli_log(ERROR, monitor, "Invalid zone name: %s !", zone);
            return 1;
        }

        struct pcap_tag *new_pcap_tag = &system_info->pcap_tag_list[system_info->num_pcap_tags];

        new_pcap_tag->id = system_info->num_pcap_tags;
        new_pcap_tag->monitor_id = monitor->color;
        new_pcap_tag->monitor_rank = monitor->world_rank;

        memset(new_pcap_tag->zone, '\0', ZONE_NAME_LEN);
        strncpy(new_pcap_tag->zone, zone_names[i], zone_names_len[i]);

        new_pcap_tag->watts_long = watts_long;
        new_pcap_tag->watts_short = watts_short;
        new_pcap_tag->seconds_long = seconds_long;
        new_pcap_tag->seconds_short = seconds_short;

        new_pcap_tag->wtime = get_time();
        gettimeofday(&(new_pcap_tag->timestamp), NULL);

        new_pcap_tag->start_timer_count = poller->time_counter;
        new_pcap_tag->pcap_flag = pcap_flag;

        int found_num = 0;
        int tag_num;

        for (tag_num = 0; tag_num < system_info->num_poli_tags; tag_num++)
        {
            struct poli_tag current_poli_tag = system_info->poli_tag_list[tag_num];
            if (current_poli_tag.closed == 0)
            {
                new_pcap_tag->active_poli_tags[found_num] = current_poli_tag;
                found_num++;
            }
        }

        new_pcap_tag->num_active_poli_tags = found_num;

        system_info->num_pcap_tags++;

        poli_log(TRACE, monitor,   "Finishing %s", __FUNCTION__);

    }

    return 0;
}

static struct pcap_tag *get_pcap_for_time_counter(int counter)
{
    poli_log(TRACE, monitor,   "Entering %s", __FUNCTION__);
    int tag_num;

    for (tag_num = 0; tag_num < system_info->num_pcap_tags; tag_num++)
    {
        struct pcap_tag *tag = &system_info->pcap_tag_list[tag_num];
        if (tag->start_timer_count == counter)
            return tag;
    }
    return 0;
}

/*                     END OF POWER CAP TAGS                                  */

/******************************************************************************/
/*                     SETTING POWER CAPS                                     */
/******************************************************************************/


//todo make this more template like
int poli_set_power_cap (double watts)
{
    return poli_set_power_cap_with_params(zone_names[PACKAGE_INDEX], watts, watts, DEFAULT_SECONDS_LONG, DEFAULT_SECONDS_SHORT);
}

int poli_set_power_cap_with_params(char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor, "Entering %s", __FUNCTION__);

        int i = get_zone_index(zone_name);
        if (i < 0)
        {
            poli_log(ERROR, monitor, "Something went wrong with getting index for zone %s. Are you sure the zone name is valid?", zone_name);
            return 1;
        }

        if (rapl_set_power_cap(zone_name, watts_long, watts_short, seconds_long, seconds_short, system_info, 1) != 0)
        {
            poli_log(ERROR, monitor,   "%s: Something went wrong with setting rapl power cap!", __FUNCTION__);
            return 1;
        }

        if (init_pcap_tag(zone_name, watts_long, watts_short, seconds_long, seconds_short, USER_SET) != 0)
        {
            poli_log(ERROR, monitor,   "%s: Something went wrong with initializing a new power cap tag!", __FUNCTION__);
            return 1;
        }

        struct pcap_info *info = &system_info->current_pcap_list[i];

        info->monitor_id = monitor->color;
        info->monitor_rank = monitor->world_rank;
        memset(info->zone, '\0', ZONE_NAME_LEN);
        strncpy(info->zone, zone_name, zone_names_len[i]);
        if (watts_long > 0)
        {
            info->enabled_long = 1;
            info->clamped_long = 1;
        }
        else
        {
            info->enabled_long = 0;
            info->clamped_long = 0;
        }
        if (watts_short > 0)
        {
            info->enabled_short = 1;
            info->clamped_short = 1;
        }
        else
        {
            info->enabled_short = 0;
            info->clamped_short = 0;
        }
        info->watts_long = watts_long;
        info->watts_short = watts_short;
        info->seconds_long = seconds_long;
        info->seconds_short = seconds_short;

        poli_log(TRACE, monitor, "Finishing %s", __FUNCTION__);
    }
    return 0;
}

int poli_reset_system (void)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor, "Entering %s", __FUNCTION__);

        if (rapl_set_power_cap("PACKAGE", (double) DEFAULT_PKG_POW, (double) DEFAULT_SHORT, (double) DEFAULT_SECONDS_LONG, (double) DEFAULT_SECONDS_SHORT, system_info, 1) ||
            rapl_set_power_cap("CORE", (double) DEFAULT_CORE_POW, 0, (double) DEFAULT_CORE_SECONDS, 0, system_info, 0))
        {
            poli_log(ERROR, monitor,   "%s: Something went wrong with setting power caps. Returning...\n", __FUNCTION__);
            return 1;
        }

        /* Set up new pcap tags to indicate change in power caps */
        if (init_pcap_tag("PACKAGE", (double) DEFAULT_PKG_POW, (double) DEFAULT_SHORT, (double) DEFAULT_SECONDS_LONG, (double) DEFAULT_SECONDS_SHORT, SYSTEM_RESET) != 0 ||
            init_pcap_tag("CORE", (double) DEFAULT_CORE_POW, 0, (double) DEFAULT_SECONDS_LONG, 0, SYSTEM_RESET) != 0)
        {
            poli_log(ERROR, monitor,   "%s: Something went wrong with initializing a new power cap tag!\n", __FUNCTION__);
            return 1;
        }

        get_system_power_caps(); //to reset the system_info->current_pcap_list

        poli_log(TRACE, monitor, "Finishing %s", __FUNCTION__);
    }

    return 0;
}

/*                    END OF SETTING POWER CAPS                               */

/******************************************************************************/
/*                    GETTING POWER CAPS                                      */
/******************************************************************************/

int poli_get_power_cap_for_param (char *zone_name, char *param, double *result)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor,   "Entering %s", __FUNCTION__);

        int found = 0;
        int pkg_exists = 0;
        /* check if the package power cap is there, in case zone_name is invalid*/
        struct pcap_info *pkg_info = &system_info->current_pcap_list[PACKAGE_INDEX];
        if (pkg_info != 0)
            pkg_exists = 1;
        /* look for power cap of zone requested */
        int i = get_zone_index(zone_name);
        if (i > -1)
            found = 1;

        if ( !found && pkg_exists)
        {
            i = PACKAGE_INDEX;
            poli_log(WARNING, monitor, "Power cap for requested zone %s could not be found. Defaulting to zone PACKAGE", zone_name);
        }
        else if ( !found && !pkg_exists)
        {
            poli_log(ERROR, monitor, "%s: Could not find power cap parameter info for zone %s", zone_name, __FUNCTION__);
            *result = 0.0;
        }
        else
        {
            if (result != NULL)
            {
                struct pcap_info *info = &system_info->current_pcap_list[i];
                if (strcmp(param, "watts_long") == 0)
                    *result = info->watts_long;
                else if (strcmp(param, "watts_short") == 0)
                    *result = info->watts_short;
                else if (strcmp(param, "seconds_long") == 0)
                    *result = info->seconds_long;
                else if (strcmp(param, "seconds_short") == 0)
                    *result = info->seconds_short;
                else if (strcmp(param, "enabled_long") == 0)
                    *result = (double) info->enabled_long;
                else if (strcmp(param, "enabled_short") == 0)
                    *result = (double) info->enabled_short;
                else if (strcmp(param, "clamped_long") == 0)
                    *result = (double) info->clamped_long;
                else if (strcmp(param, "clamped_short") == 0)
                    *result = (double) info->clamped_short;
                else
                {
                    poli_log(WARNING, monitor, "The parameter: %s was not recognized. Returning watts_long for zone %s as default.", param, zone_name);
                    *result = info->watts_long;
                }
            }
        }

        poli_log(TRACE, monitor, "Finishing %s", __FUNCTION__);
    }

#ifndef _NOMPI
    MPI_Bcast(result, 1, MPI_DOUBLE, 0, monitor->mynode_comm); //this is necessary for user to retrieve value
#endif

    return 0;
}

int poli_get_power_cap (double *watts)
{
    return poli_get_power_cap_for_param(zone_names[PACKAGE_INDEX], "watts_long", watts);
}

//TODO currently supporting only RAPL
static int get_system_power_cap_for_zone (int zone_index)
{
    if (monitor->imonitor)
    {
        struct msr_pcap pcap;
        int pret = rapl_get_power_cap(&pcap, zone_names[zone_index], system_info);
        if (pret != 0)
        {
            poli_log(ERROR, monitor, "%s: Something went wrong with getting RAPL power cap", __FUNCTION__);
            return -1;
        }

        struct pcap_info *info = &system_info->current_pcap_list[zone_index];

        info->monitor_id = monitor->color;
        info->monitor_rank = monitor->world_rank;
        memset(info->zone, '\0', ZONE_NAME_LEN);
        strncpy(info->zone, zone_names[zone_index], zone_names_len[zone_index]);
        info->zone_label = pcap.zone_label;
        info->enabled_long = pcap.enabled_long;
        info->watts_long = pcap.watts_long;
        info->watts_short = pcap.watts_short;
        info->seconds_long = pcap.seconds_long;
        info->seconds_short = pcap.seconds_short;
        info->enabled_short = pcap.enabled_short;
        info->clamped_long = pcap.clamped_long;
        info->clamped_short = pcap.clamped_short;
    }

    return 0;
}

//TODO this is also specific to RAPL
static int get_system_power_caps (void)
{
    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor,   "Entering %s", __FUNCTION__);

        int i;
        for (i = 0; i < system_info->sysmsr->num_zones; i++)
            get_system_power_cap_for_zone(i);

        poli_log(TRACE, monitor,   "Finishing %s", __FUNCTION__);
    }

    return 0;
}

int poli_get_power_cap_limits (char* zone_name, double *min, double *max)
{
    if (monitor->imonitor)
    {
        int index = get_zone_index(zone_name);
        struct pcap_info *current_pcap = &system_info->current_pcap_list[index];

        if (!current_pcap->min || !current_pcap->max)
            rapl_get_power_cap_info(zone_name, &current_pcap->min, &current_pcap->max,
                &current_pcap->thermal_spec, &current_pcap->max_time_window, system_info);

        *min = current_pcap->min;
        *max = current_pcap->max;
    }

#ifndef _NOMPI
    MPI_Bcast(min, 1, MPI_DOUBLE, 0, monitor->mynode_comm);
    MPI_Bcast(max, 1, MPI_DOUBLE, 0, monitor->mynode_comm);
#endif

    return 0;
}

void poli_print_power_cap_info (void)
{
    if (monitor->imonitor)
    {
        printf("************************************************************\n");
        printf("                       POWER CAP INFO                       \n");
        printf("                       RANK: %d NODE: %s                    \n", monitor->world_rank, monitor->my_host);
        printf("------------------------------------------------------------\n");

        struct pcap_info *current_pcap = &system_info->current_pcap_list[PACKAGE_INDEX];
        if (current_pcap != 0)
        {
            printf("\tzone: %s\n", current_pcap->zone);
            printf("\twatts_long: %lf\n", current_pcap->watts_long);
            printf("\tseconds_long: %lf\n", current_pcap->seconds_long);
            printf("\tenabled_long: %d\n", current_pcap->enabled_long);
            printf("\tclamped_long: %d\n", current_pcap->clamped_long);
            printf("\twatts_short: %lf\n", current_pcap->watts_short);
            printf("\tseconds_short: %lf\n", current_pcap->seconds_short);
            printf("\tenabled_short: %d\n", current_pcap->enabled_short);
            printf("\tclamped_short: %d\n", current_pcap->clamped_short);
            printf("\tthermal specification: %lf\n", current_pcap->thermal_spec);
            printf("\tmax power cap %lf\n", current_pcap->max);
            printf("\tmin power cap %lf\n", current_pcap->min);
            printf("\tmax time window %lf\n", current_pcap->max_time_window);
        }
        else
            printf("Power cap for zone PACKAGE on rank %d has not yet been recorded.\n", monitor->world_rank);

        printf("************************************************************\n");
    }
}

void poli_print_power_cap_info_verbose (void)
{
    if (monitor->imonitor)
    {
        int i;
        printf("************************************************************\n");
        printf("                POWER CAP INFO VERBOSE                      \n");
        printf("                RANK: %d NODE: %s                           \n", monitor->world_rank, monitor->my_host);
        printf("------------------------------------------------------------\n");
        for (i = 0; i < system_info->sysmsr->num_zones; i++)
        {
            struct pcap_info *current_pcap = &system_info->current_pcap_list[i];
            if (current_pcap != 0)
            {
                printf("\tzone: %s\n", current_pcap->zone);
                printf("\twatts_long: %lf\n", current_pcap->watts_long);
                printf("\tseconds_long: %lf\n", current_pcap->seconds_long);
                printf("\tenabled_long: %d\n", current_pcap->enabled_long);
                printf("\tclamped_long: %d\n", current_pcap->clamped_long);
                int short_supported = 0;
                if (current_pcap->zone_label == PACKAGE || current_pcap->zone_label == PLATFORM)
                    short_supported = 1;
                if (short_supported)
                {
                    printf("\twatts_short: %lf\n", current_pcap->watts_short);
                    printf("\tseconds_short: %lf\n", current_pcap->seconds_short);
                    printf("\tenabled_short: %d\n", current_pcap->enabled_short);
                    printf("\tclamped_short: %d\n", current_pcap->clamped_short);
                }
                printf("\tthermal specification: %lf\n", current_pcap->thermal_spec);
                printf("\tmax power cap %lf\n", current_pcap->max);
                printf("\tmin power cap %lf\n", current_pcap->min);
                printf("\tmax time window %lf\n", current_pcap->max_time_window);
            }
            else
                printf("Power cap for zone %s on rank %d has not yet been recorded.\n", zone_names[i], monitor->world_rank);
        }
        printf("************************************************************\n");
    }
}


/*               END OF GETTING POWER CAPS                                    */


/******************************************************************************/
/*              FREQUENCY                                                     */
/******************************************************************************/

int poli_get_current_frequency (double *freq)
{
    int ret = 0;
    if (monitor->imonitor)
    {
        struct system_poll_info info;
        get_current_frequency(&info);
        double cpufreq = info.freq.freq;
        if (cpufreq == 0.0)
        {
            poli_log(ERROR, monitor,   "Unable to get frequency from cpufreq");
            ret = 1;
        }
        else
        {
            poli_log(TRACE, monitor,   "Frequency obtained from cpufreq: %lf KHz\n", cpufreq);
            (*freq) = cpufreq;
            return 0;
        }

#ifdef _CRAY
        double crayfreq = info.freq.cray_freq;
        if (crayfreq == 0.0)
        {
            poli_log(ERROR, monitor,   "Unable to get frequency from Cray stack.");
            ret = 1;
        }
        else
        {
            poli_log(TRACE, monitor,   "Frequency obtained from Cray stack: %lf KHz\n", crayfreq);
            (*freq) = crayfreq;
            return 0;
        }
#endif
    }
    return ret;
}

int poli_print_frequency_info (void)
{
    struct system_poll_info info;
    get_current_frequency(&info);
    double cpufreq = info.freq.freq;
    printf("************************************************************\n");
    printf("                     FREQUENCY INFO                         \n");
    printf("                     RANK: %d NODE: %s                      \n", monitor->world_rank, monitor->my_host);
    printf("------------------------------------------------------------\n");
    printf("Frequency obtained from cpufreq: %lf KHz\n", cpufreq);
#ifdef _CRAY
    double crayfreq = info.freq.cray_freq;
    printf("Frequency obtained from Cray stack: %lf KHz\n", crayfreq);
#endif
    printf("************************************************************\n");
    return 0;
}

static int get_current_frequency (struct system_poll_info * info)
{
    int fret = read_cpufreq(&(info->freq.freq));
    if (fret)
    {
        poli_log(ERROR, monitor,   "%s: Something went wrong with getting frequency form cpufreq", __FUNCTION__);
        info->freq.freq = 0.0;
    }
#ifdef _CRAY
    info->freq.cray_freq = cray_read_pm_counter(system_info->syscray->counters[CRAY_FREQ_INDEX].pm_file);
#endif
    return 0;
}

static int read_cpufreq (double *freq)
{
    if (monitor->imonitor)
    {
        char buff[200];
        memset(buff, '\0', sizeof(buff));

        int hz = 0;
        if (system_info->cur_freq_file > 0)
        {
            if (poller->time_counter > 0)
            {
                if(lseek(system_info->cur_freq_file, 0, SEEK_SET) < 0)
                {
                    poli_log(ERROR, monitor,   "Attempt to set file descriptor to beginning of frequency file failed! %s", strerror(errno));
                    return 1;
                }
            }
            #ifdef _TIMER_OFF
            if (lseek(system_info->cur_freq_file, 0, SEEK_SET) < 0) {
                poli_log(ERROR, monitor,   "Attempt to set file descriptor to beginning of frequency file failed! %s", strerror(errno));
                return 1;
            }
            #endif
            if (read(system_info->cur_freq_file, buff, sizeof(buff)))
            {
                char *token;
                token = strtok(buff, " \t\n");
                hz = atoi(token);
            }

            (*freq) = (double) (hz / 1000.0);
        }
        else
            (*freq) = 0.0;
    }
    return 0;
}

/*               END OF FREQUENNCY                                            */

/******************************************************************************/
/*              TIMER                                                         */
/******************************************************************************/

#ifndef _TIMER_OFF
static int setup_timer (void)
{
    if (monitor->imonitor)
    {
        poller->timer_on = 1;

        memset(&poller->sa, 0, sizeof(poller->sa));
        poller->sa.sa_handler = &timer_handler;

        int status = sigaction(SIGALRM, &poller->sa, NULL);
        if (0 != status)
        {
            poli_log(ERROR, monitor,   "Failed to set SIGACTION: %s", strerror(errno));
            return 1;
        }

        poller->timer.it_value.tv_sec = 0;
        poller->timer.it_value.tv_usec = INITIAL_TIMER_DELAY;

        poller->timer.it_interval.tv_sec = 0;
        poller->timer.it_interval.tv_usec = POLL_INTERVAL * 1000000;

        status = setitimer(ITIMER_REAL, &poller->timer, NULL);
        if (0 != status)
        {
            poli_log(ERROR, monitor,   "Failed to set timer: %s", strerror(errno));
            return 1;
        }
    }
    return 0;
}

static int stop_timer (void)
{
    poller->timer_on = 0;
    sigaction(SIGALRM, &poller->sa, NULL);

    poller->timer.it_value.tv_sec = 0;
    poller->timer.it_value.tv_usec = 0;

    poller->timer.it_interval.tv_sec = 0;
    poller->timer.it_interval.tv_usec = 0;

    setitimer(ITIMER_REAL, &poller->timer, NULL);
    return 0;
}

static void init_energy_reading (struct energy_reading *reading)
{
    struct rapl_energy rapl_energy = {0};
    reading->rapl_energy = rapl_energy;
#ifdef _CRAY
    struct cray_measurement cray_meas = {0};
    reading->cray_meas = cray_meas;
#else
#ifdef _BGQ
    struct bgq_measurement bgq_meas = {0};
    reading->bgq_meas = bgq_meas;
#endif
#endif
    return;
}

static void timer_handler (int signum)
{
    if (poller->timer_on && monitor->imonitor)
    {
        if (poller->time_counter < MAX_POLL_SAMPLES)
        {
            double start_iter_time = get_time();
            double pcap;

            struct system_poll_info *info = &system_info->system_poll_list[poller->time_counter];

            struct energy_reading last_energy;

            if (poller->time_counter > 0)
                last_energy = system_info->system_poll_list[poller->time_counter-1].current_energy;
            else
                init_energy_reading(&last_energy);

            int zone;
            for (zone = 0; zone < system_info->sysmsr->num_zones; zone++)
            {
                info->pcap_info_list[zone] = system_info->current_pcap_list[zone];
            }

            info->counter = poller->time_counter;

            info->pkg_pcap = system_info->current_pcap_list[PACKAGE_INDEX].watts_long;
            get_current_frequency(info);
            info->current_energy = read_current_energy(system_info);
            info->last_energy = last_energy;

            if (poller->time_counter == 0)
                compute_current_power(info, system_info->initial_mpi_wtime, system_info);
            else
                compute_current_power(info, ((double) POLL_INTERVAL), system_info);

            info->wtime = get_time();
            info->poll_iter_time = info->wtime - start_iter_time;

            poller->time_counter++;
        }
    }
    return;
}
#endif


#ifdef _BENCH
void timer_handler_em (void)
{
    if (monitor->imonitor)
    {
        if (poller->time_counter_em < MAX_POLL_SAMPLES)
        {
            double start_iter_time = get_time();
            double pcap;

            struct system_poll_info *info = &system_info->system_poll_list_em[time_counter_em];

            struct energy_reading last_energy;

            if (poller->time_counter > 0)
                last_energy = system_info->system_poll_list_em[time_counter_em-1].current_energy;
            else
            {
                struct rapl_energy rapl_energy = {0};
                last_energy.rapl_energy = rapl_energy;
            #ifdef _CRAY
                struct cray_measurement cray_meas = {0};
                last_energy.cray_meas = cray_meas;
            #endif
            }

            int zone;
            for (zone = 0; zone < system_info->sysmsr->num_zones; zone++)
            {
                info->pcap_info_list[zone] = system_info->current_pcap_list[zone];
            }

            info->counter = time_counter_em;

            info->pkg_pcap = system_info->current_pcap_list[PACKAGE_INDEX].watts_long;
            get_current_frequency(info);
            info->current_energy = read_current_energy(system_info);
            info->last_energy = last_energy;

            if (time_counter == 0)
                compute_current_power(info, system_info->initial_mpi_wtime, system_info);
            else
                compute_current_power(info, ((double) POLL_INTERVAL), system_info);

            info->current_energy.rapl_energy.package, info->current_energy.rapl_energy.dram, info->computed_power.rapl_energy.package, info->computed_power.rapl_energy.dram,
            info->current_energy.cray_meas.node_energy, info->current_energy.cray_meas.cpu_energy, info->current_energy.cray_meas.memory_energy,
            info->computed_power.cray_meas.node_power, info->computed_power.cray_meas.cpu_power, info->computed_power.cray_meas.memory_power);

            info->wtime = get_time();
            info->poll_iter_time = info->wtime - start_iter_time;

            time_counter_em++;
        }
    }

    return;
}
#endif

/*               END OF TIMER                                                 */

/******************************************************************************/
/*              HELPERS                                                       */
/******************************************************************************/

static double get_time (void)
{
#ifndef _NOMPI
    return MPI_Wtime();
#else
    return (double) time(NULL);
#endif
}

static int compute_current_power (struct system_poll_info * info, double time, struct system_info_t * system_info)
{
    if (time == system_info->initial_mpi_wtime)
        rapl_compute_total_power(&(info->computed_power.rapl_energy), &(info->current_energy.rapl_energy), time);
    else
    {
        struct rapl_energy diff;
        rapl_compute_total_energy(&(diff), &(info->current_energy.rapl_energy), &(info->last_energy.rapl_energy));
        rapl_compute_total_power(&(info->computed_power.rapl_energy), &(diff), time);
    }
#ifdef _CRAY
    compute_cray_total_measurements(&(info->computed_power.cray_meas), &(info->current_energy.cray_meas), &(info->last_energy.cray_meas), time);
#elif _BGQ
    compute_bgq_total_measurements(&(info->computed_power.bgq_meas), &(info->current_energy.bgq_meas), &(info->last_energy.bgq_meas), time);
#endif
    return 0;
}

static struct energy_reading read_current_energy (struct system_info_t * system_info)
{
    struct energy_reading current_energy;
    rapl_read_energy(&(current_energy.rapl_energy), system_info);
#ifdef _CRAY
    get_cray_measurement(&(current_energy.cray_meas), system_info);
#elif _BGQ
    init_bgq_measurement(&(current_energy.bgq_meas));
    get_bgq_measurement(&(current_energy.bgq_meas), system_info);
#endif
    return current_energy;
}

static void poli_sync (void)
{
#ifndef _NOMPI
    int finalized;
    MPI_Finalized(&finalized);
    if (monitor->world_size > 1 && !finalized)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    return;
}

static void poli_sync_node (void)
{
#ifndef _NOMPI
    int finalized;
    MPI_Finalized(&finalized);
    if (monitor->node_size > 1 && !finalized)
        MPI_Barrier(monitor->mynode_comm);
#endif
    return;
}


static void get_timestamp(double time_from_start, char *time_str_buffer, size_t buff_len)
{
    double frac, intpart;
    int64_t fracpart;
    struct timeval cur_time;

    char time_str[8192];
    memset(time_str, '\0', 8192);
    char *time_str_ptr;

    frac = modf(time_from_start, &intpart);
    cur_time.tv_sec = system_info->initial_start_time.tv_sec + (int64_t) intpart;
    fracpart = (int32_t) (frac * 1000000.0);

    if ((system_info->initial_start_time.tv_usec + fracpart) >= 1e6) {
        cur_time.tv_usec = (system_info->initial_start_time.tv_usec + fracpart - 1e6);
        cur_time.tv_sec++;
    }
    else {
        cur_time.tv_usec = system_info->initial_start_time.tv_usec + fracpart;
    }

    time_t rawtime = (time_t)cur_time.tv_sec;
    struct tm *timeinfo = localtime(&rawtime);

    strftime (time_str_buffer, buff_len, "%Y-%m-%d %H:%M:%S", timeinfo);

    return;

}

static FILE * open_file (char *filename)
{
    FILE * fp;
    char *prefix;
    prefix = getenv("PoLi_PREFIX");
    char file[1000];
    if (prefix != NULL)
        sprintf(file, "%s%s_%s_%s.txt", prefix, filename, monitor->my_host, monitor->jobid);
    else
        sprintf(file, "%s_%s_%s.txt", filename, monitor->my_host, monitor->jobid);

    fp = fopen(file, "w");
    if (!fp)
    {
        poli_log(ERROR, monitor,   "Failed to open file %s: %s", file, strerror(errno));
        fclose(fp);
        return NULL;
    }
    return fp;
}

int get_zone_index (char *zone_name)
{
    int i;
    for (i = 0; i < system_info->sysmsr->num_zones; i++)
    {
        if (strncmp(zone_name, zone_names[i], zone_names_len[i]) == 0)
            return i;
    }
    return -1;
}

#ifndef _NOMPI
int poli_am_monitor (void)
{
    if (monitor->imonitor)
        return 1;
    else
        return 0;
}

int poli_get_monitor_id (int *id)
{
    if (monitor->imonitor)
        (*id) = monitor->world_rank;
    else
        (*id) = -1;
    return 0;
}

int poli_get_monitor_node (int *node_id)
{
    if (monitor->imonitor)
        (*node_id) = atoi(monitor->my_host);
    else
        (*node_id) = -1;
    return 0;
}

int poli_get_subcommunicator (MPI_Comm *comm)
{
    *comm = monitor->mynode_comm;
    return 0;
}

int poli_get_num_monitors (void)
{
    int num_monitors;
    MPI_Allreduce(&num_monitors, &monitor->imonitor, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return num_monitors;
}
#endif

/*
coordsToInt - composes an integer out of the coordinates on an Aries router
input: the coordinates to convert and the number of coordinates
returns: the resulting integer
*/
// (x, y, z) == (2, 3, 4) becomes res = 234
int coordsToInt (int *coords, int dim)
{
    int i, res = 0;

    for ( i = 0; i < dim; i++ )
        res += coords[i] * pow (10.0, (double)i);

    return res;
}

static int compute_power_from_tag(struct poli_tag *tag, double time)
{
    rapl_compute_total_energy(&(tag->total_energy.rapl_energy), &(tag->end_energy.rapl_energy), &(tag->start_energy.rapl_energy));
    rapl_compute_total_power(&(tag->total_power.rapl_energy), &(tag->total_energy.rapl_energy), time);
#ifdef _CRAY
    compute_cray_total_measurements(&(tag->total_energy.cray_meas), &(tag->end_energy.cray_meas), &(tag->start_energy.cray_meas), time);
#elif _BGQ
    compute_bgq_total_measurements(&(tag->total_energy.bgq_meas), &(tag->end_energy.bgq_meas), &(tag->start_energy.bgq_meas), time);
#endif

    return 0;
}

/*                  END OF HELPERS                                            */

/******************************************************************************/
/*                  FINALIZING/CLEANUP                                        */
/******************************************************************************/

static int file_handler (void)
{
    int ret = 0;
    if (monitor->imonitor)
    {
        if (polling_info_to_file() != 0)
        {
            ret = 1;
            poli_log(ERROR, monitor,   "Something went wrong with writing polling output to file");
        }
        if (system_info->num_poli_tags > 0)
        {
            if (poli_tags_to_file() != 0)
            {
                ret = (ret || 1);
                poli_log(ERROR, monitor,   "Something went wrong with writing energy tags to file\n");
            }
        }
        if (system_info->num_pcap_tags > 0)
        {
            if (pcap_tags_to_file() != 0)
            {
                ret = (ret || 1);
                poli_log(ERROR, monitor,   "Something went wrong with \n");
            }
        }
    }
    return ret;
}

int poli_tags_to_file (void)
{
    if (monitor->imonitor)
    {
        FILE *fp = open_file("PoLiMEr_energy-tags");
        if (fp == NULL)
            return 1;

#ifndef _HEADER_OFF
        fprintf(fp, "Tag Name\tTimestamp\tStart Time (s)\tEnd Time (s)\tTotal Time (s)\t");
        if (!system_info->sysmsr->error_state)
        {
            fprintf(fp, "Total RAPL pkg E (J)\tTotal RAPL PP0 E (J)\tTotal RAPL PP1 E (J)\tTotal RAPL platform E (J)\tTotal RAPL dram E (J)\t");
            fprintf(fp, "Total RAPL pkg P (W)\tTotal RAPL PP0 P (W)\tTotal RAPL PP1 P (W)\tTotal RAPL platform P (W)\tTotal RAPL dram P (W)");
        }
#ifdef _CRAY
        fprintf(fp, "\tTotal Cray node E (J)\tTotal Cray cpu E (J)\tTotal Cray memory E (J)\t");
        fprintf(fp, "Total Cray node P (W)\tTotal Cray cpu P (W)\tTotal Cray memory P (W)\t");
        fprintf(fp, "Total Cray node calc P (W)\tTotal Cray cpu calc P (W)\tTotal Cray memory calc P (W)");
#endif
#ifdef _BGQ
        write_bgq_header(&fp);
#endif
        fprintf(fp, "\n");
#endif
        int tag_num;
        for (tag_num = 0; tag_num < system_info->num_poli_tags; tag_num++)
        {
            struct poli_tag *tag = &system_info->poli_tag_list[tag_num];

            double total_time = tag->end_time - tag->start_time;
            double start_offset, end_offset = 0.0;
            if (strcmp(tag->tag_name, "application_summary") == 0 && total_time < 0)
            {
                total_time = get_time() - system_info->initial_mpi_wtime;
                end_offset = total_time;
            }

            compute_power_from_tag(tag, total_time);

            if (strcmp(tag->tag_name, "application_summary") == 0)
            {
                start_offset = 0.0;
                if (end_offset == 0.0)
                    end_offset = tag->end_time - tag->start_time;
            }
            else
            {
                start_offset = tag->start_time - system_info->initial_mpi_wtime;
                end_offset = tag->end_time - system_info->initial_mpi_wtime;
            }

            char time_str_buffer[20];
            get_timestamp(start_offset, time_str_buffer, sizeof(time_str_buffer));

            fprintf(fp, "%s\t%s\t%lf\t%lf\t%lf\t", tag->tag_name, time_str_buffer, start_offset, end_offset, total_time);

            if (!system_info->sysmsr->error_state)
            {
                struct rapl_energy total_energy = tag->total_energy.rapl_energy;
                struct rapl_energy total_power = tag->total_power.rapl_energy;
                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t", total_energy.package, total_energy.pp0, total_energy.pp1, total_energy.platform, total_energy.dram);
                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf", total_power.package, total_power.pp0, total_power.pp1, total_power.platform, total_power.dram);
            }
#ifdef _CRAY
            struct cray_measurement total_measurements = tag->total_energy.cray_meas;
            fprintf(fp, "\t%lf\t%lf\t%lf\t", total_measurements.node_energy, total_measurements.cpu_energy, total_measurements.memory_energy);
            fprintf(fp, "%lf\t%lf\t%lf\t", total_measurements.node_power, total_measurements.cpu_power, total_measurements.memory_power);
            fprintf(fp, "%lf\t%lf\t%lf\n", total_measurements.node_measured_power, total_measurements.cpu_measured_power, total_measurements.memory_measured_power);
#else
#ifdef _BGQ
            struct bgq_measurement bgq_meas = tag->total_energy.bgq_meas;
            printf("\t%lf\n", bgq_meas.card_power);
            printf("%lf\n", bgq_meas.cpu);
            printf("%lf\n", bgq_meas.dram);
            printf("%lf\n", bgq_meas.optics);
            printf("%lf\n", bgq_meas.pci);
            printf("%lf\n", bgq_meas.network);
            printf("%lf\n", bgq_meas.link_chip);
            printf("%lf\n", bgq_meas.sram);
            write_bgq_output(&fp, &bgq_meas);
#endif
            fprintf(fp, "\n");
#endif
        }
        fclose(fp);
    }

    return 0;
}

int pcap_tags_to_file (void)
{
    if (monitor->imonitor)
    {
        FILE *fp = open_file("PoLiMEr_powercap-tags");
        if (fp == NULL)
            return 1;

#ifndef _HEADER_OFF
        fprintf(fp, "Tag ID\tZone\tTimestamp\tPower Cap Long (W)\tPower Cap Short (W)\tTime Window Long (s)\tTime Window Short (s)\tTime since start (s)\tPCAP FLAG\tNumber of active poli tags\tEmon tag list\n");
#endif

        int tag_num;
        for (tag_num = 0; tag_num < system_info->num_pcap_tags; tag_num++)
        {
            struct pcap_tag *tag = &system_info->pcap_tag_list[tag_num];
            double start_offset = tag->wtime - system_info->initial_mpi_wtime;

            char time_str_buffer[20];
            get_timestamp(start_offset, time_str_buffer, sizeof(time_str_buffer));

            fprintf(fp, "%d\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t", tag->id,
                tag->zone, time_str_buffer, tag->watts_long, tag->watts_short,
                tag->seconds_long, tag->seconds_short, start_offset,
                tag->pcap_flag, tag->num_active_poli_tags);

            int i;
            for (i = 0; i < tag->num_active_poli_tags; i++)
            {
                struct poli_tag etag = tag->active_poli_tags[i];
                if (i < tag->num_active_poli_tags - 1)
                    fprintf(fp, "\"%s,", etag.tag_name);
                else
                    fprintf(fp, "%s\"", etag.tag_name);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    return 0;
}

static int polling_info_to_file (void)
{
    if (monitor->imonitor)
    {
        if (poller->time_counter <= 0)
            return 0;

#ifndef _TIMER_OFF

        FILE *fp = open_file("PoLiMEr");
        if (fp == NULL)
            return 1;

        int zone;
#ifndef _HEADER_OFF
        fprintf(fp, "Count\tTimestamp\tTime since start (s)\t");
        if (!system_info->sysmsr->error_state)
        {
            fprintf(fp, "RAPL pkg E (J)\tRAPL pp0 E (J)\tRAPL pp1 E (J)\tRAPL platform E (J)\tRAPL dram E (J)\t");
            fprintf(fp, "RAPL pkg E since start (J)\tRAPL pp0 E since start (J)\tRAPL pp1 E since start (J)\tRAPL platform E since start (J)\tRAPL dram E since start (J)\t");
            fprintf(fp, "RAPL pkg P (W)\tRAPL pp0 P (W)\tRAPL pp1 P (W)\tRAPL platform P (W)\tRAPL dram P (W)");
        }
#ifdef _CRAY
        fprintf(fp, "\tCray node E (J)\tCray cpu E (J)\tCray memory E (J)\t");
        fprintf(fp, "Cray node E since start (J)\tCray cpu E since start (J)\tCray memory E since start (J)\t");
        fprintf(fp, "Cray node P (W)\tCray cpu P (W)\tCray memory P (W)\t");
        fprintf(fp, "Cray node P calc (W)\tCray cpu P calc (W)\tCray memory P calc (W)\t");
        fprintf(fp, "Cpufreq frequency (MHz)\tCray frequency (MHz)\t");
#else
#ifdef _BGQ
        write_bgq_header(&fp);
        write_bgq_ediff_header(&fp);
#endif
        fprintf(fp, "\tCpufreq frequency (MHz)\t");
#endif
        if (!system_info->sysmsr->error_state)
        {
            for (zone = 0; zone < system_info->sysmsr->num_zones - 1; zone++)
            {
                fprintf(fp, "%s power cap long (W)\t", zone_names[zone]);
                fprintf(fp, "%s power cap short (W)\t", zone_names[zone]);
            }
            fprintf(fp, "%s power cap long (W)\t", zone_names[system_info->sysmsr->num_zones - 1]);
            fprintf(fp, "%s power cap short (W)\n", zone_names[system_info->sysmsr->num_zones - 1]);
        }
#endif
        int counter;
        for (counter = 0; counter < poller->time_counter; counter++)
        {
            if (system_info->num_poli_tags > 0)
            {
                struct poli_tag *this_end_poli_tag = get_poli_tag_for_end_time_counter(counter);
                if (this_end_poli_tag != 0)
                    fprintf(fp, "--- EMON TAG END: %s\n", this_end_poli_tag->tag_name);
            }
            if (system_info->num_pcap_tags > 0)
            {
                struct pcap_tag *this_pcap = get_pcap_for_time_counter(counter);
                if (this_pcap != 0)
                    fprintf(fp, "*** SET POWER CAP TAG %d TO: %s, %lf\n",
                        this_pcap->id, this_pcap->zone, this_pcap->watts_long);
            }
            if (system_info->num_poli_tags > 0)
            {
                struct poli_tag *this_start_poli_tag = get_poli_tag_for_start_time_counter(counter);
                if (this_start_poli_tag != 0)
                    fprintf(fp, "--- EMON TAG START: %s\n", this_start_poli_tag->tag_name);
            }

            struct system_poll_info *info = &system_info->system_poll_list[counter];

            double time_from_start = info->wtime - system_info->initial_mpi_wtime;

            char time_str_buffer[20];
            get_timestamp(time_from_start, time_str_buffer, sizeof(time_str_buffer));

            fprintf(fp, "%d\t%s\t%lf\t", info->counter, time_str_buffer, time_from_start);

            if (!system_info->sysmsr->error_state)
            {
                struct rapl_energy *energy_j = &(info->current_energy.rapl_energy);
                struct rapl_energy *watts = &(info->computed_power.rapl_energy);

                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t", energy_j->package, energy_j->pp0, energy_j->pp1, energy_j->platform, energy_j->dram);
                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t", (energy_j->package - system_info->initial_energy.rapl_energy.package), (energy_j->pp0 - system_info->initial_energy.rapl_energy.pp0), (energy_j->pp1 - system_info->initial_energy.rapl_energy.pp1), (energy_j->platform - system_info->initial_energy.rapl_energy.platform), (energy_j->dram - system_info->initial_energy.rapl_energy.dram));
                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t", watts->package, watts->pp0, watts->pp1, watts->platform, watts->dram);
            }
#ifdef _CRAY
            struct cray_measurement *cmeasurement = &(info->current_energy.cray_meas);
            struct cray_measurement *cpower = &(info->computed_power.cray_meas);

            fprintf(fp, "%lf\t%lf\t%lf\t", cmeasurement->node_energy, cmeasurement->cpu_energy, cmeasurement->memory_energy);
            fprintf(fp, "%lf\t%lf\t%lf\t", (cmeasurement->node_energy - system_info->initial_energy.cray_meas.node_energy), (cmeasurement->cpu_energy - system_info->initial_energy.cray_meas.cpu_energy), (cmeasurement->memory_energy - system_info->initial_energy.cray_meas.memory_energy));
            fprintf(fp, "%lf\t%lf\t%lf\t", cmeasurement->node_power, cmeasurement->cpu_power, cmeasurement->memory_power);
            fprintf(fp, "%lf\t%lf\t%lf\t", cpower->node_measured_power, cpower->cpu_measured_power, cpower->memory_measured_power);
            fprintf(fp, "%lf\t%lf\t", info->freq.freq, info->freq.cray_freq);
#else
#ifdef _BGQ
            struct bgq_measurement *bgq_meas = &(info->current_energy.bgq_meas);
            write_bgq_output(&fp, bgq_meas);
            write_bgq_ediff(&fp, bgq_meas, &(system_info->initial_energy.bgq_meas));
#endif
            fprintf(fp, "%lf\t", info->freq.freq);
#endif
            for (zone = 0; zone < system_info->sysmsr->num_zones - 1; zone++)
            {
                fprintf(fp, "%lf\t", info->pcap_info_list[zone].watts_long);
                fprintf(fp, "%lf\t", info->pcap_info_list[zone].watts_short);
            }
            fprintf(fp, "%lf\t", info->pcap_info_list[system_info->sysmsr->num_zones - 1].watts_long);
            fprintf(fp, "%lf\n", info->pcap_info_list[system_info->sysmsr->num_zones - 1].watts_short);
        }
        fclose(fp);
#else //_TIMER_OFF is set
        return 0;
#endif
    }

    return 0;
}

static void finalize_power_interfaces (struct system_info_t * system_info)
{
    finalize_msrs(system_info);
#ifdef _CRAY
    finalize_cray_pm_counters(system_info);
#endif
    return;
}

static int finalize_tags (void)
{
    if (system_info->num_open_tags + system_info->num_closed_tags == 0)
        return 0;

    poli_log(WARNING, monitor, "There are unfinished tags! Attempting to close them...");

    int tag_num;

    for (tag_num = 0; tag_num < system_info->num_poli_tags; tag_num++)
    {
        struct poli_tag *this_poli_tag = &system_info->poli_tag_list[tag_num];
        if ( this_poli_tag->closed == 0)
        {
            poli_log(WARNING, monitor, "Tag %s is not finished. It will be closed automatically.", this_poli_tag->tag_name);
            end_existing_poli_tag(this_poli_tag);
        }
    }

    return 0;
}

int poli_finalize(void)
{
    poli_sync();

    poli_log(TRACE, monitor, "Entering %s", __FUNCTION__);

    if (monitor->imonitor)
    {
        poli_log(TRACE, monitor, "Finalizing: Getting application power, energy and time summary");

        poli_log(TRACE, monitor, "Checking if any poli tags are unfinished");

        struct poli_tag *app_summary = &system_info->poli_tag_list[0]; //need to close application summary tag which is the first one
        end_existing_poli_tag(app_summary);
        /* Check if any poli tags are unfinished */
        finalize_tags();

        poli_log(TRACE, monitor, "Resetting the system");

        /* Reset system power caps */
        if (poli_reset_system() != 0)
            if (!system_info->sysmsr->error_state)
                poli_log(ERROR, monitor, "Couldn't reset system!");

#ifndef _TIMER_OFF
        poli_log(TRACE, monitor, "Stopping timer");
        stop_timer();
#endif
        poli_log(TRACE, monitor, "Pushing results to file");
        file_handler();

        poli_log(TRACE, monitor,   "Closing frequency file");
        if (system_info->cur_freq_file)
            close(system_info->cur_freq_file);

        poli_log(TRACE, monitor,   "Cleaning up structures");

        /* Cleanup */
        if (system_info->poli_tag_list)
        {
            free(system_info->poli_tag_list);
            system_info->poli_tag_list = 0;
        }
        if (system_info->pcap_tag_list)
        {
            free(system_info->pcap_tag_list);
            system_info->pcap_tag_list = 0;
        }
        if (system_info->current_pcap_list) //this must be done after system reset
        {
            free(system_info->current_pcap_list);
            system_info->current_pcap_list = 0;
        }
#ifndef _TIMER_OFF
        if (system_info->system_poll_list)
        {
            free(system_info->system_poll_list);
            system_info->system_poll_list = 0;
        }
#endif
#ifdef _BENCH
        if (system_info->system_poll_list_em)
        {
            free(system_info->system_poll_list_em);
            system_info->system_poll_list_em = 0;
        }
#endif

        finalize_power_interfaces(system_info);

        if (poller)
            free(poller);

        if (system_info)
            free(system_info);

    }

#ifndef _NOMPI
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Comm_free(&monitor->mynode_comm);
#endif

    if (monitor)
        free(monitor);

    poli_log(TRACE, monitor,   "Finishing %s", __FUNCTION__);

    return 0;
}

/*                   END OF FINALIZING/CLEANUP                                */
