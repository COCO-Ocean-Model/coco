# COCO

COCO is the ocean general circulation model. It\'s also the oceanic part of the coupled general circulation model MIROC and NICOCO.

## Overview

COCO is an open-source software package developed by contributors, including researchers affiliated with the Atmosphere and Ocean Research Institute (AORI), the University of Tokyo and the Japan Agency for Marine-Earth Science and Technology (JAMSTEC).

The software is intended for scientific research, operational applications, education, and commercial use.

For information about the management and maintainers of this repository, see `MAINTAINERS.md`.

## Open-source Release

**COCO v6.0.0 is the first version publicly released as open-source software.** The major version number reflects substantial updates accumulated since the previous closed-source release (v5.0.0), including changes that are not fully backward compatible.


## History

COCO traces its origins to a prototype ocean model developed by Nobuo Suginohara in the 1970s. The prototype and its successors were continuously used and developed within the ocean modeling group led by Professor Suginohara at the Center for Climate System Research (CCSR), the University of Tokyo, providing the foundation for subsequent generations of ocean models.
Building upon this legacy, Yasuhiro Yamanaka developed a substantially renewed modeling framework, introducing the numerical discretization and coding strategy that later served as an important reference in the design of COCO.

Based on these concepts, Hiroyasu Hasumi developed COCO (CCSR Ocean Component Model) as a full-fledged ocean general circulation model and authored its core architecture and implementation. Since then, COCO has been continuously developed and maintained through the contributions of members of the ocean modeling group at CCSR (now the Atmosphere and Ocean Research Institute, the University of Tokyo), together with collaborators at the Japan Agency for Marine-Earth Science and Technology (JAMSTEC). These collaborative efforts have expanded the model's physical capabilities, computational performance, and applicability to a wide range of oceanographic and climate research.

## Features

COCO is designed as a flexible and efficient ocean general circulation model suitable for a wide range of applications, from process studies to global climate simulations.

* Solves the hydrostatic, Boussinesq primitive equations on the sphere.
* Explicit free-surface formulation.
* Generalized curvilinear horizontal coordinate system, including support for tripolar grids.
* Hybrid vertical coordinate combining z-level and sigma coordinates, enabling improved representation of the upper ocean while avoiding surface layer outcropping.
* Barotropic/baroclinic mode splitting for computational efficiency.
* Modular physical parameterizations, including options for vertical and lateral mixing, sea surface forcing, and bottom boundary processes.
* Multi-category sea ice model with elastic-viscous-plastic (EVP) rheology.
* One-layer sea ice thermodynamics.
* Designed for applications ranging from regional process studies to global ocean and coupled climate simulations.
* Used as the ocean component of the MIROC coupled climate model.

## How to get

```bash
git clone https://github.com//coco-pub.git
cd coco-pub
```

## Usage

Detailed documentation is currently available in Japanese. See the documentation included in this repository for build instructions, configuration, and usage examples.

## License

This software is released under the BSD 3-Clause License.

See `LICENSE` for details.

## Citation

If this software contributes to scientific publications, technical reports, operational products, commercial services, or other publicly distributed works, the authors kindly request that the software and/or associated publications be cited.

Citation information is provided in `CITATION.md`

## Disclaimer

This software is provided "as is", without warranty of any kind, express or implied.

The authors, contributors, copyright holders, and their affiliated organizations shall not be liable for any claim, damages, or other liability arising from the use of this software.

## Endorsement

Use of this software does not imply endorsement by the authors, contributors, or their affiliated organizations.

Neither the names of the copyright holders, contributors, nor the names of their affiliated institutions (including JAMSTEC) may be used to endorse or promote products derived from this software without specific prior written permission.

## Copyright

Copyright in this software remains with the respective copyright holders.

Where copyright belongs to an organization under applicable laws, regulations, or institutional policies governing works made in the course of employment, that organization is recognized as the copyright holder.

The following copyright notice is used for the initial open-source release of COCO:

> Copyright (c) 2026 JAMSTEC, Hiroyasu Hasumi, and other copyright holders.
>
> See COPYRIGHT.md for details.

See `COPYRIGHT.md` for the list of copyright holders, contributors, and their representative contributions.
