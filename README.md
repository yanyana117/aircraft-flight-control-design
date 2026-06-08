# Aircraft Flight Control Design

MATLAB flight-control design project for an F-16A lateral/directional linear
model with actuator dynamics and digital controller simulation.

This repository is a cleaned portfolio version of an aircraft control project.
It keeps the controller design scripts and helper routines while removing
course handouts, grading sheets, temporary files, and unrelated homework.

## Highlights

- Modeled F-16A lateral/directional dynamics with aileron and rudder inputs.
- Added first-order actuator dynamics with position/rate command outputs.
- Designed digital controllers using LQR/SDR-style cost weighting.
- Implemented NZSP tracking, PI-SDR, PI-NZSP, and command/rate weighting cases.
- Simulated closed-loop state, control, command, and actuator-rate histories.

## Repository Structure

```text
controllers/
  sdr_regulator.m
  nzsp_tracking_control.m
  pi_sdr_integral_control.m
  pi_nzsp_tracking_control.m
  pif_nzsp_command_rate_weighting.m
utils/
  lqrdjv.m
  QPMCALC.M
docs/
  project_summary.md
```

## Requirements

- MATLAB
- Control System Toolbox

Add the helper folder before running a controller script:

```matlab
addpath("utils")
run("controllers/nzsp_tracking_control.m")
```

## Resume Summary

Designed and simulated digital flight-control laws for an F-16A
lateral/directional model in MATLAB, including actuator dynamics, LQR/SDR gain
selection, closed-loop modal analysis, and time-history evaluation of states,
commands, and control-surface rates.
