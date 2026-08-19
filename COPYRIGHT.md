# Copyright

## Copyright Notice

Copyright (c) 2026 JAMSTEC, Hiroyasu Hasumi, and other copyright holders.

See the copyright holders listed below for details.

---

## Copyright Holders

The following organization and individuals hold copyright in all or part of COCO.

### Individual Copyright Holders

| Name              | Affiliation at the Time Copyright Arose                           | Current Affiliation                                              | Representative Contribution                                                  |
| ----------------- | ----------------------------------------------------------------- | ---------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| Hiroyasu Hasumi   | Center for Climate System Research (CCSR), The University of Tokyo       | Atmosphere and Ocean Research Institute (AORI), The University of Tokyo | Overall architecture and core framework                                      |
| Yasuhiro Yamanaka | Center for Climate System Research (CCSR), The University of Tokyo       | Hokkaido University                                              | Overall architecture and core framework inherited from the predecessor model |
| Hideyuki Nakano   | Center for Climate System Research (CCSR), The University of Tokyo       | Meteorological Research Institute, Japan Meteorological Agency (MRI/JMA)      | Bottom boundary layer scheme |
| Akira Oka         | Center for Climate System Research (CCSR), The University of Tokyo       | Atmosphere and Ocean Research Institute (AORI), The University of Tokyo | Noh and Kim surface mixed-layer scheme                                       |
| Yoshiki Komuro    | Center for Climate System Research (CCSR), The University of Tokyo       | Japan Agency for Marine-Earth Science and Technology (JAMSTEC)   | Support for CORE and other boundary-condition datasets                       |
| Hiroaki Tatebe    | Center for Climate System Research (CCSR), The University of Tokyo       | Japan Agency for Marine-Earth Science and Technology (JAMSTEC)   | Second-Order Moment tracer advection scheme                                  |
| Takao Kawasaki    | Atmosphere and Ocean Research Institute (AORI)*, The University of Tokyo | Japan Agency for Marine-Earth Science and Technology (JAMSTEC)   | Tidal mixing parameterization                                                |

* The Atmosphere and Ocean Research Institute (AORI) is the successor organization to the Center for Climate System Research (CCSR).

### Institutional Copyright Holder

| Copyright Holder                                               | Developers                                                                                             | Representative Contributions                                                                                                                                                                                                                                                                                                                                            |
| -------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Japan Agency for Marine-Earth Science and Technology (JAMSTEC) | Tatsuo Suzuki, Hiroaki Tatebe, Yoshiki Komuro, Masao Kuroki, Takao Kawasaki, Fuyuki Saito, Koji Ogochi | Fortran 90/95 conversion and modularization of the code base; tripolar grid support; time-staggered stepping scheme; sea-ice parameterizations; GLS surface mixed-layer scheme; mixed-layer eddy parameterization; support for JRA and other boundary-condition datasets; density-coordinate output; multi-format I/O framework; OpenACC support; GTOOL3 library integration; performance optimization. |

---

## Contributors

The following individuals have made significant technical contributions to COCO but are not listed above as copyright holders.

| Name               | Affiliation                                                         | Representative Contribution                                        |
| ------------------ | ------------------------------------------------------------------- | ------------------------------------------------------------------ |
| Atsushi Numaguchi  |                                              | Development of core infrastructure routines |
| Hiroyuki Tsujino   | Meteorological Research Institute, Japan Meteorological Agency (MRI/JMA)     | Early debugging, testing, and development of COCO                  |
| Shogo Urakawa      | Meteorological Research Institute, Japan Meteorological Agency (MRI/JMA)     | Extensive debugging and verification of the model             |
| Takateru Yamagishi | Research Organization for Information Science and Technology (RIST) | GPU support, including the implementation of OpenACC                        |

---

## Notes

* The representative contributions summarize the principal technical contributions of each copyright holder or contributor. They are intended to acknowledge major contributions and do not define or limit the scope of copyright.
* Affiliations are shown at the time of the principal contribution and, where applicable, the current affiliation.
* This document is maintained for attribution purposes and does not alter the ownership of copyright or the terms of the BSD 3-Clause License.
