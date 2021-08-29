#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "PoLiMEr.h"
#include "PoLiLog.h"
#include "cray_pm-handler.h"

const char *path = "/sys/cray/pm_counters/";

const char *pm_filenames[NUM_COUNTERS] = {"energy", "power", "cpu_energy", "cpu_power",
    "memory_energy", "memory_power", "power_cap", "raw_scan_hz", "freshness", "generation",
    "version", "startup"};

static double evaluate_boundaries(double end, double start, int compute_power, double energy, double total_time);

int init_cray_pm_counters (struct system_info_t * system_info)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    system_info->syscray = malloc(sizeof(struct system_cray_info));
    system_info->syscray->counters = 0;
    system_info->syscray->num_counters = NUM_COUNTERS;

    system_info->syscray->counters = calloc(MAX_NUM_COUNTERS, sizeof(struct pm_counter));

    int counter;
    for (counter = 0; counter < NUM_COUNTERS; counter++)
    {
        char filename[100];
        memset(filename, '\0', 100);
        strcpy(filename, path);
        strcat(filename, pm_filenames[counter]);
        system_info->syscray->counters[counter].pm_filename = filename;

        if (strcmp(pm_filenames[counter], "energy") == 0)
            system_info->syscray->counters[counter].type = ENERGY;
        if (strcmp(pm_filenames[counter], "power") == 0)
            system_info->syscray->counters[counter].type = POWER;
        if (strcmp(pm_filenames[counter], "cpu_energy") == 0)
            system_info->syscray->counters[counter].type = CPU_ENERGY;
        if (strcmp(pm_filenames[counter], "cpu_power") == 0)
            system_info->syscray->counters[counter].type = CPU_POWER;
        if (strcmp(pm_filenames[counter], "memory_energy") == 0)
            system_info->syscray->counters[counter].type = MEMORY_ENERGY;
        if (strcmp(pm_filenames[counter], "memory_power") == 0)
            system_info->syscray->counters[counter].type = MEMORY_POWER;
        if (strcmp(pm_filenames[counter], "power_cap") == 0)
            system_info->syscray->counters[counter].type = POWER_CAP;
        if (strcmp(pm_filenames[counter], "raw_scan_hz") == 0)
            system_info->syscray->counters[counter].type = RAW_SCAN_HZ;
        if (strcmp(pm_filenames[counter], "freshness") == 0)
            system_info->syscray->counters[counter].type = FRESHNESS;
        if (strcmp(pm_filenames[counter], "generation") == 0)
            system_info->syscray->counters[counter].type = GENERATION;
        if (strcmp(pm_filenames[counter], "version") == 0)
            system_info->syscray->counters[counter].type = VERSION;
        if (strcmp(pm_filenames[counter], "startup") == 0)
            system_info->syscray->counters[counter].type = STARTUP;

        system_info->syscray->counters[counter].pm_file = cray_open_pm_file(system_info->syscray->counters[counter].pm_filename);
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

int finalize_cray_pm_counters (struct system_info_t * system_info)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    int counter;
    for (counter = 0; counter < NUM_COUNTERS; counter++)
    {
        struct pm_counter *pm_counter = &system_info->syscray->counters[counter];
        if (pm_counter->pm_file)
            close(pm_counter->pm_file);
    }
    free(system_info->syscray->counters);
    if (system_info->syscray)
        free(system_info->syscray);

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return 0;
}

int cray_open_pm_file (char *counter_name)
{
    poli_log(TRACE, NULL, "Entering %s", __FUNCTION__);

    int pm_file;

    pm_file = open(counter_name, O_RDONLY);
    if (pm_file < 0)
    {
        poli_log(ERROR, NULL, "Something went wrong with opening the cray counter %s", counter_name);
        return -1;
    }

    poli_log(TRACE, NULL, "Finishing %s", __FUNCTION__);

    return pm_file;
}

double cray_read_pm_counter(int pm_file)
{
    double result = 0.0;
    char buff[200];
    memset(buff, '\0', sizeof(buff));
    if (pm_file)
    {
        if(lseek(pm_file, 0, SEEK_SET) < 0)
            poli_log(ERROR, NULL, "Attempt to set file descriptor to beginning of file %d failed!", pm_file);
        if (read(pm_file, buff, sizeof(buff)))
        {
            char *token;
            token = strtok(buff, " \t\n");
            sscanf(token, "%lf", &result);
        }
    }
    return result;
}

int get_cray_measurement (struct cray_measurement *cm, struct system_info_t * system_info)
{
    cm->node_energy = -1;
    cm->node_power = -1;
    cm->cpu_energy = -1;
    cm->cpu_power = -1;
    cm->memory_energy = -1;
    cm->memory_power = -1;

    int i;
    for (i = 0; i < system_info->syscray->num_counters; i++)
    {
        double measurement = -1.0;
        //get measurement only from counters that report power and energy measurements
        if ((system_info->syscray->counters[i].type == ENERGY) || (system_info->syscray->counters[i].type == POWER) ||
            (system_info->syscray->counters[i].type == CPU_ENERGY) || (system_info->syscray->counters[i].type == CPU_POWER) ||
            (system_info->syscray->counters[i].type == MEMORY_POWER) || (system_info->syscray->counters[i].type == MEMORY_ENERGY))
        {
            measurement = cray_read_pm_counter(system_info->syscray->counters[i].pm_file);

            // find out what counter we just read and what measurement to update
            switch (system_info->syscray->counters[i].type)
            {
                case ENERGY:
                    cm->node_energy = measurement;
                    break;
                case POWER:
                    cm->node_power = measurement;
                    break;
                case CPU_POWER:
                    cm->cpu_power = measurement;
                    break;
                case CPU_ENERGY:
                    cm->cpu_energy = measurement;
                    break;
                case MEMORY_ENERGY:
                    cm->memory_energy = measurement;
                    break;
                case MEMORY_POWER:
                    cm->memory_power = measurement;
                    break;
                default:
                    poli_log(ERROR, NULL, "The counter %d did not report a valid measurement.", system_info->syscray->counters[i].type);
                    break;
            }
        }
    }

    if (cm->node_energy == -1 && cm->node_power == -1 && cm->cpu_energy == -1 &&
        cm->cpu_power == -1 && cm->memory_energy == -1 && cm->memory_power == -1)
        poli_log(ERROR, NULL, "%s: wasn't able to get any measurements from %d counters.\n", __FUNCTION__, system_info->syscray->num_counters);

    return 0;
}

int compute_cray_total_measurements (struct cray_measurement *cm, struct cray_measurement *end, struct cray_measurement *start, double total_time)
{
    double endnode_energy, endnode_power, endcpu_energy, endcpu_power, endmemory_energy, endmemory_power = 0.0;
    if (end)
    {
        endnode_energy = end->node_energy;
        endnode_power = end->node_power;
        endcpu_energy = end->cpu_energy;
        endcpu_power = end->cpu_power;
        endmemory_energy = end->memory_energy;
        endmemory_power = end->memory_power;
    }
    cm->node_energy = evaluate_boundaries(endnode_energy, start->node_energy, 0, -1, total_time);
    cm->node_power = evaluate_boundaries(endnode_power, start->node_power, 0, -1, total_time);
    cm->node_measured_power = evaluate_boundaries(-1, 0, 1, cm->node_energy, total_time);
    cm->cpu_energy = evaluate_boundaries(endcpu_energy, start->cpu_energy, 0, -1, total_time);
    cm->cpu_power = evaluate_boundaries(endcpu_power, start->cpu_power, 0, -1, total_time);
    cm->cpu_measured_power = evaluate_boundaries(-1, 0, 1, cm->cpu_energy, total_time);
    cm->memory_energy = evaluate_boundaries(endmemory_energy, start->memory_energy, 0, -1, total_time);
    cm->memory_power = evaluate_boundaries(endmemory_power, start->memory_power, 0, -1, total_time);
    cm->memory_measured_power = evaluate_boundaries(-1, 0, 1, cm->memory_energy, total_time);
    return 0;
}

static double evaluate_boundaries(double end, double start, int compute_power, double energy, double total_time)
{
    if (end >= start)
        return (end - start);
    if (compute_power && (energy > 0))
        return (energy / total_time);
    return 0.0;
}
