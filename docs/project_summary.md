# Project Summary

This project designs digital flight controllers for a linear F-16A
lateral/directional model. The model includes sideslip, roll rate, yaw rate,
roll angle, yaw angle, aileron actuator dynamics, and rudder actuator dynamics.

## Methods

- Continuous state-space aircraft model.
- First-order actuator augmentation.
- Zero-order-hold discretization.
- LQR/SDR gain design with continuous-to-discrete cost conversion.
- NZSP tracking and integral-control variants.
- Closed-loop eigenvalue, damping, and time-history analysis.

## Controller Scripts

- `sdr_regulator.m`: baseline SDR regulator design.
- `nzsp_tracking_control.m`: yaw-angle tracking with NZSP formulation.
- `pi_sdr_integral_control.m`: integral-state SDR controller.
- `pi_nzsp_tracking_control.m`: PI tracking formulation.
- `pif_nzsp_command_rate_weighting.m`: command/rate weighted tracking case.

## Interview Talking Points

- Why actuator dynamics are included before controller design.
- How weighting matrices affect LQR/SDR controller behavior.
- How discrete sample time changes the controller implementation.
- How closed-loop poles relate to damping ratio and settling behavior.
- How command and rate limits affect practical control-system design.
