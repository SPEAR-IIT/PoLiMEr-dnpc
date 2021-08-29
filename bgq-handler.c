#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>

#include <spi/include/kernel/location.h>
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/memory.h>
#include <firmware/include/personality.h>
#include <hwi/include/bqc/nd_500_dcr.h>
#include <mpix.h>
#include <hwi/include/common/uci.h>

// Needed for timebase
//#include <hwi/include/bqc/A2_inlines.h>

#define EMON_DEFINE_GLOBALS
#include <spi/include/emon/emon.h>

#include "PoLiMEr.h"
#include "PoLiLog.h"
#include "bgq-handler.h"


// Info about the co-ordinates of the process
unsigned int row, col, midplane, nodeboard, computecard, core;

static double evaluate_boundaries(double end, double start);

int get_comm_split_color_bgq (struct monitor_t * monitor)
{
    int status;
    uint32_t procIDOnNode;

    MPIX_Hardware_t hw;
    MPIX_Hardware(&hw);
    Personality_t personality;
    status = Kernel_GetPersonality(&personality, sizeof(personality));
    int coord[hw.torus_dimension];

    procIDOnNode = Kernel_ProcessorID();

    bg_decodeComputeCardCoreOnNodeBoardUCI(personality.Kernel_Config.UCI, &row, &col,
        &midplane, &nodeboard, &computecard, &core);

    coord[0] = personality.Network_Config.cnBridge_A;
    coord[1] = personality.Network_Config.cnBridge_B;
    coord[2] = personality.Network_Config.cnBridge_C;
    coord[3] = personality.Network_Config.cnBridge_D;
    coord[4] = personality.Network_Config.cnBridge_E;

    //printf("BG_persInfo_init, my coords{%u,%u,%u,%u,%u} rankInPset %u,sizeOfPset %u,idOfPset %u\n",hw.Coords[0],hw.Coords[1],hw.Coords[2],hw.Coords[3],hw.Coords[4],hw.rankInPset,hw.sizeOfPset,hw.idOfPset);

    char temp[20];

    snprintf(temp, 20, "%d%d%d%d%d", hw.Coords[0],hw.Coords[1],hw.Coords[2],hw.Coords[3],hw.Coords[4]);

    monitor->color = atoi(temp);

    snprintf(monitor->my_host, sizeof(monitor->my_host),
      "RC%X%X-M%d-N%02d-CC%02d-A%dof%d-B%dof%d-C%dof%d-D%dof%d-E%dof%d",
      row, col, midplane, nodeboard, computecard,
      personality.Network_Config.Acoord, personality.Network_Config.Anodes,
      personality.Network_Config.Bcoord, personality.Network_Config.Bnodes,
      personality.Network_Config.Ccoord, personality.Network_Config.Cnodes,
      personality.Network_Config.Dcoord, personality.Network_Config.Dnodes,
      personality.Network_Config.Ecoord, personality.Network_Config.Enodes
);

    // Assign compute card 0 and procIDOnNode 0 to do monitoring
    if ((0 == computecard) && (0 == procIDOnNode))
        monitor->imonitor = 1;

    status = 0;
    if (monitor->imonitor)
        status = EMON_SetupPowerMeasurement();

    return status;
}

void init_bgq_measurement (struct bgq_measurement *bm)
{
    int i;
    for (i = 0; i < 14; i++)
    {
        bm->volts[i] = 0.0;
        bm->amps[i] = 0.0;
    }

    bm->card_power = 0.0;
    bm->cpu = 0.0;
    bm->dram = 0.0;
    bm->optics = 0.0;
    bm->pci = 0.0;
    bm->network = 0.0;
    bm->link_chip = 0.0;
    bm->sram = 0.0;

    bm->card_en = 0.0;
    bm->cpu_en = 0.0;
    bm->dram_en = 0.0;
    bm->optics_en = 0.0;
    bm->pci_en = 0.0;
    bm->network_en = 0.0;
    bm->link_chip_en = 0.0;
    bm->sram_en = 0.0;
}


int get_bgq_measurement (struct bgq_measurement *bm, struct system_info_t * system_info)
{
    double pw = EMON_GetPower_impl(bm->volts, bm->amps);
    if (pw == -1)
    {
        poli_log(ERROR, NULL, "EMON_GetPower_impl() failed!");
        return 1;
    }

    bm->card_power = pw;

    bm->k_const = (double) domain_info[0].k_const;
    bm->cpu += bm->volts[0] * bm->amps[0] * bm->k_const + bm->volts[1] * bm->amps[1] * bm->k_const;

    bm->k_const = domain_info[1].k_const;
    bm->dram += bm->volts[2] * bm->amps[2] * bm->k_const + bm->volts[3] * bm->amps[3] * bm->k_const;

    bm->k_const = domain_info[2].k_const;
    bm->optics += bm->volts[4] * bm->amps[4] * bm->k_const + bm->volts[5] * bm->amps[5] * bm->k_const;

    bm->k_const = domain_info[3].k_const;
    bm->pci += bm->volts[6] * bm->amps[6] * bm->k_const + bm->volts[7] * bm->amps[7] * bm->k_const;

    bm->k_const = domain_info[4].k_const;
    bm->network += bm->volts[8] * bm->amps[8] * bm->k_const + bm->volts[9] * bm->amps[9] * bm->k_const;

    bm->k_const = domain_info[5].k_const;
    bm->link_chip += bm->volts[10] * bm->amps[10] * bm->k_const + bm->volts[11] * bm->amps[11] * bm->k_const;

    bm->k_const = domain_info[6].k_const;
    bm->sram += bm->volts[12] * bm->amps[12] * bm->k_const + bm->volts[13] * bm->amps[13] * bm->k_const;

    return 0;
}


int compute_bgq_total_measurements (struct bgq_measurement *bm,
    struct bgq_measurement *end, struct bgq_measurement *start, double total_time)
{
    double endcard_power, endcpu, enddram, endoptics, endpci, endnetwork, endlink_chip, endsram = 0.0;
    if (end)
    {
        endcard_power = end->card_power;
        endcpu = end->cpu;
        enddram = end->dram;
        endoptics = end->optics;
        endpci = end->pci;
        endnetwork = end->network;
        endlink_chip = end->link_chip;
        endsram = end->sram;
    }
    bm->card_power = evaluate_boundaries(endcard_power, start->card_power);
    bm->cpu = evaluate_boundaries(endcpu, start->cpu);
    bm->dram = evaluate_boundaries(enddram, start->dram);
    bm->optics = evaluate_boundaries(endoptics, start->optics);
    bm->pci = evaluate_boundaries(endpci, start->pci);
    bm->network = evaluate_boundaries(endnetwork, start->network);
    bm->link_chip = evaluate_boundaries(endlink_chip, start->link_chip);
    bm->sram = evaluate_boundaries(endsram, start->sram);

    bm->card_en = bm->card_power * total_time;
    bm->cpu_en = bm->cpu * total_time;
    bm->dram_en = bm->dram * total_time;
    bm->optics_en = bm->optics * total_time;
    bm->pci_en = bm->pci * total_time;
    bm->network_en = bm->network_en * total_time;
    bm->link_chip_en = bm->link_chip * total_time;
    bm->sram_en = bm->sram * total_time;

    return 0;
}


static double evaluate_boundaries(double end, double start)
{
    if (end >= start)
        return (end - start);
    return start;
}


void write_bgq_header(FILE **fp)
{
    fprintf(*fp, "Row\tCol\tMidplane\tNodeboard\t");
    fprintf(*fp, "BGQ Node Card E (J)\tBGQ cpu E (J)\tBGQ dram E (J)\tBGQ optics E (J)\tBGQ pci E (J)\tBGQ network E (J)\tBGQ link chip E (J)\tBGQ sram E (J)\t");
    fprintf(*fp, "BGQ Node Card P (W)\tBGQ cpu P (W)\tBGQ dram P (W)\tBGQ optics P (W)\tBGQ pci P (W)\tBGQ network P (W)\tBGQ link chip P (W)\tBGQ sram P (W)\t");
}

void write_bgq_output(FILE **fp, struct bgq_measurement *bm)
{
    fprintf(*fp, "%X\t%X\t%d\t%02d\t", row, col, midplane, nodeboard);

    fprintf(*fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", bm->card_en, bm->cpu_en,
        bm->dram_en, bm->optics_en, bm->pci_en, bm->network_en, bm->link_chip_en, bm->sram_en);

    fprintf(*fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", bm->card_power, bm->cpu,
        bm->dram, bm->optics, bm->pci, bm->network, bm->link_chip, bm->sram);

}

void write_bgq_ediff_header(FILE **fp)
{
    fprintf(*fp, "BGQ Node Card E since start (J)\tBGQ cpu E since start (J)\tBGQ dram E since start (J)\tBGQ optics E since start (J)\tBGQ pci E since start (J)\tBGQ network E since start (J)\tBGQ link chip E since start (J)\tBGQ sram E since start (J)\t");
}

void write_bgq_ediff(FILE **fp, struct bgq_measurement *end, struct bgq_measurement *start)
{
    fprintf(*fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", end->card_power - start->card_power,
        end->cpu - start->cpu, end->dram - start->dram, end->optics - start->optics, end->pci - start->pci,
        end->network - start->network, end->link_chip - start->link_chip, end->sram - start->sram);
}

