# Power System Dynamic Analysis of the IEEE 9-Bus System

Final project for **EEE 306: Power System I Laboratory** (BUET, July 2023). Group project, Section C2, Group 04.

MATLAB and Simulink analysis of the IEEE 9-bus system: load flow, time-domain dynamic simulation of the synchronous generators, and overcurrent relay behaviour under LG, LL, LLG and three-phase faults.

## What's inside

| Path | Description |
|---|---|
| `Main_File.m` | Main script: runs the load flow, builds the Y-bus and runs the time-domain simulation |
| `IEEE9Bus.m` | IEEE 9-bus case data (MATPOWER format) |
| `DAEs1.m`, `DAEs2.m`, `DAEs3.m` | System equations for the pre-fault, fault and post-fault stages (solved with `ode15s`) |
| `results/` | Load flow results |
| `report/` | Full project report |

## Requirements

MATLAB with Simulink and [MATPOWER](https://matpower.org/) on the MATLAB path. Run `Main_File.m`.
