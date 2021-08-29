**NOTE:** This is a modified version of a repository that already existed.
- https://xgitlab.cels.anl.gov/open-source/PoLiMEr 
- Ivana Marincic, Venkatram Vishwanath, and Henry Hoffmann. 2017. PoLiMEr: An Energy Monitoring and Power Limiting Interface for HPC Applications. In <i>Proceedings of the 5th International Workshop on Energy Efficient Supercomputing</i> (<i>E2SC'17</i>). Association for Computing Machinery, New York, NY, USA, Article 7, 1–8. DOI:https://doi.org/10.1145/3149412.3149419

Everything below is by the original author.

# PoLiMEr - **Po**wer **Li**miting & **M**onitoring of **E**negy

Thank you for your interest in this project! PoLiMEr is still under development, so bear with me if the documentation is not in sync with the project.

#### Audience
PoLiMEr is suitable for all HPC users, experienced and beginners, even users working on small-scale projects.

# Prerequisites (Intel systems)

To use PoLiMEr, either the msr-safe module is needed, or read/write permissions to the MSRs (usually requires root access).

If you are on Theta, the msr-safe module is enabled by default and you can proceed to the next section.

Otherwise, if you are not sure, check the output of `ls /dev/cpu` (on compute nodes, if applicable) and verify that there is a `msr_whitelist` file. If yes, you should be able to proceed to the next section. If there isn't such a file, contact you system administrator to install the msr-safe module on your system.

If you want to set up the msr-safe module yourself, please refer to https://github.com/LLNL/msr-safe

# Getting Started

Clone and build the repo, then navigate to the Testing section down below.

## Clone the repo and build

Clone and set up directories:
```
git clone https://xgitlab.cels.anl.gov/open-source/PoLiMEr.git
cd PoLiMEr
mkdir bin
mkdir lib
```

If you are on Theta:
```
make CRAY=yes COBALT=yes
```
On a non-Argonne XC40 simply omit the Cobalt flag: `make CRAY=yes`.

If you are on a Blue Gene/Q system (if not using Argonne systems omit `COBALT=yes`):
```
make BGQ=yes COBALT=yes
```

If you are on JLSE:
```
make JLSE=yes
```

This will place the PoLiMEr static and shared library under `PoLiMEr/lib`. Configure the `LIBDIR` variable to change the path of the library.

#### Build Options

If you would like debug messages turned on use:
```
make DEBUG=yes
```
For more detailed debug messages:
```
make TRACE=yes
```

Other flags:

* `NOMPI=yes` builds PoLiMEr for profiling non-MPI applications. MPI is turned on by default.
* `NOOMP=yes` builds PoLiMEr for profiling non-OpenMP applications. OpenMP is turned on by default.
* `COBALT=yes` if your system uses the Cobalt scheduler.
* `PMI=yes` if the PMI library is present on your system (true for XC40). Not necessary if using the `CRAY` option.
* `HEADER_OFF=yes` to turn off printing of headers in output files
* `TIMER_OFF=yes` to turn off polling feature
* `BENCH=yes` to time specific PoLiMEr functions (used to measure PoLiMEr overhead)

# Testing

The `test` directoy included with PoLiMEr comes with a hello world application that demonstrates how to link with PoLiMEr and that you can run to test if PoLiMEr works properly on your system.

```
cd test/hello
```

For Theta run:
```
make CRAY=yes
qsub -A <your_project> simple_run.sh
```

On a different XC40, please pass the executable `./mpi_power_demo <num_threads>` to the scheduler of your system. You are free to choose the number of nodes, ranks and threads (e.g. 2 nodes, 8 ranks/node, 2 threads/rank).

For Blue Gene/Q run:
```
make BGQ=yes
qsub -A <your_project> simple_run_bgq.sh
```

Note that it is recommended to request at least 32 nodes for power measurements in BGQ systems.

For JLSE run: `make JLSE=yes`

# Linking

##### Static Linking

In the Makefile of the application set the include path: `-I${PATH_TO_POLIMER}/PoLiMEr/include`

Then specify the library, after the object file: `${PATH_TO_POLIMER}/PoLiMEr/lib/libpolimer.a`

##### Linking with shared library
In the Makefile of the application set the include path: `-I${PATH_TO_POLIMER}/PoLiMEr/include`

Then specify the library, after the object file: `-L${PATH_TO_POLIMER}/PoLiMEr/lib -lpolimer`

or `-L${PATH_TO_POLIMER}/PoLiMEr/lib ${PATH_TO_POLIMER}/PoLiMEr/lib/libpolimer.so`



# Usage

To use PoLiMEr at the bare minimum, add the following in your application:
```
#include "PoLiMEr.h"
...
MPI_Init(&argc, &argv);
...
poli_init(); //must be called *AFTER* MPI_Init if profiling an MPI application

...

poli_finalize(); //must be called *BEFORE* MPI_Finalize

MPI_Finalize();
```

If polling is not switched off (if `TIMER_OFF` is not set), this will produce a file `PoLiMEr_<node>_<jobid>.txt` which contains power and energy measurements polled at a specified interval during the application runtime.

Another file: `PoLiMEr_energy-tags_<node>_<jobid>.sh` contains total aggregate power, energy and time of the application. This file is always generated with the tag `application_summary`.

The third file `PoLiMEr_powercap-tags_<node>_<jobid>.sh` contains information about when and what power caps were set. This file is always generated marking that the system has been reset when PoLiMEr was finalized.

### Tagging

#### Basic tagging:

```
poli_start_tag("my_tag_name");

/*do some work here */

poli_end_tag("my_tag_name");
```

This will start energy, power and time measurements for that specific code block. In the energy tags file there will be a new entry called `my_tag_name`. In the poller file, the start and end of the tag will be indicated.

#### Profiling loops

Both of these are allowed.

```
1.
poli_start_tag("outer_loop");
for (int i = 0; i < 100; i++)
    /*do some work*/
poli_end_tag("outer_loop");

2.
for (int j = 0; j < 100; j++)
{
    poli_start_tag("inside_loop");
    /*do some work*/
    poli_end_tag("inside_loop");
}
```

The first example will generate a single tag called `outer_loop`.

The second example will generate 100 tags called `inside_loop` and place them all in the energy tag file. To get total energy, power and timing, the 100 tags can be aggregated, summing over the energy and power values, or taking the maximum time. Another option is to obtain the average energy, power and time.

#### Other combinations

Nested tags are allowed:

```
poli_start_tag("tag1");
...
poli_start_tag("tag2");
...
poli_start_tag("tag3");
...
poli_end_tag("tag3");
...
poli_end_tag("tag2");
...
poli_end_tag("tag1");
```

<aside class="warning">
Interleaved tags are **NOT** allowed! You are responsible for closing the tags in the correct order when using *distinct* nested tags. This choice was made for the sake of performance.
</aside>

Example of interleaved tags:
```
// Do not use interleaved tags
// PoLiMEr will not close these tags in the correct order.
poli_start_tag("tag1");
poli_start_tag("tag2");
poli_end_tag("tag1");
poli_start_tag("tag2");
```

In this example, when calling `poli_end_tag("tag1")`, PoLiMEr will actually close tag2 first, and later tag1. While PoLiMer will not throw an error and will handle the discrepancy in the best way it cam, the resulting energy, power and time values will not reflect your initial intent. It's best to avoid this scenario completely, and either separate these tags out, use nested tags, or simply repeat a tag.

In case you forget to close a tag, PoLiMEr will close it at the end, and notify you about which tag was not closed.

In case you forget to open a tag, PoLiMEr will close the last tag that was opened, or issue a warning if there were no open tags at all.

### Power Limiting/Capping

If you just want a general power cap applied without having to think about it use:
```
double watts = 120.0;
poli_set_power_cap(watts);
```
By default, power cap will be imposed on the package level (on XC40 the package corresponds to the processor die, ie a single CPU)

You can also specify extra power cap parameters:

```
poli_set_power_cap_with_params (char *zone_name, double watts_long, double watts_short, double seconds_long, double seconds_short);
```

The possible zone names are:
* `PACKAGE`
* `CORE`
* `UNCORE`
* `PLATFORM`
* `DRAM`

**Note that XC40 only supports `PACKAGE, CORE, DRAM`.**

** FURTHER NOTE: Most systems, including XC40 only support setting power cap of PACKAGE. Even though the other zones are active, enforcing power caps on them may not be supported.**

To read in the current value of a power cap use the default PACKAGE function:

```
double watts;
poli_get_power_cap(&watts);
```

Or to retrieve a specific parameter for a specific zone use:
```
poli_get_power_cap_for_param (char *zone_name, char *param, double *result);
```

The possible paramters are:
* `watts_long`
* `watts_short`
* `seconds_long`
* `seconds_short`
* `enabled_long`
* `enabled_short`
* `clamped_long`
* `clamped_short`


More instructions will be added later. For now, see `PoLiMEr.h` for the list of user-accessible functions.
