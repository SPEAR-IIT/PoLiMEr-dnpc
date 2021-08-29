/* Largerly adopted from MonEQ https://doi.org/10.1016/j.parco.2016.05.015 */

#ifndef __BGQ_HANDLER_H
#define __BGQ_HANDLER_H

#ifdef __cplusplus
//extern "C"
//{
#endif

struct system_info_t;
struct monitor_t;


struct bgq_measurement
{
    double card_power;
    double volts[14];
    double amps[14];
    double cpu;
    double dram;
    double link_chip;
    double network;
    double optics;
    double pci;
    double sram;

    double k_const;

    double card_en;
    double cpu_en;
    double dram_en;
    double link_chip_en;
    double network_en;
    double optics_en;
    double pci_en;
    double sram_en;
};

int get_comm_split_color_bgq (struct monitor_t * monitor);

void init_bgq_measurement (struct bgq_measurement *bm);

int get_bgq_measurement (struct bgq_measurement *bm, struct system_info_t * system_info);
int compute_bgq_total_measurements (struct bgq_measurement *bm, struct bgq_measurement *end, struct bgq_measurement *start, double total_time);



void write_bgq_header(FILE **fp);
void write_bgq_output(FILE **fp, struct bgq_measurement *bm);
void write_bgq_ediff_header(FILE **fp);
void write_bgq_ediff(FILE **fp, struct bgq_measurement *end, struct bgq_measurement *start);

#ifdef __cplusplus
//}
#endif


#endif
