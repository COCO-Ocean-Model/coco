# Citation Guide

Thank you for using **COCO**.

If you use COCO in scientific publications, presentations, technical reports, operational products, commercial services, or other publicly distributed works, please cite the software as described below.

## General Citation

To ensure the reproducibility and traceability of scientific results, **we strongly recommend citing the GitHub Release corresponding to the version of COCO used in your study**. A GitHub Release provides a permanent snapshot of the exact source code, documentation, and release history associated with your results.

Whenever possible, cite a tagged release rather than the development branch.

**Recommended software citation**

> COCO Development Team (Year). *COCO*, Version X.Y.Z. GitHub Release. https://github.com/ORG/COCO/releases/tag/vX.Y.Z

In addition to citing the software itself, users are encouraged to cite a publication describing the model that is appropriate for their application.

The original description of COCO is:

> Hasumi, H. (2006). *CCSR Ocean Component Model (COCO) Version 4.0*. CCSR Report No. 25, Center for Climate System Research, The University of Tokyo.

Although this report describes an earlier implementation of COCO, it remains the primary reference dedicated to the model itself. Since then, COCO has evolved substantially, including the transition from the original Fortran 77 implementation to the current Fortran 90/95 code base. At present, no dedicated model description paper is available for the modern implementation.

For studies using COCO as the ocean component of **MIROC6**, users may alternatively cite:

> Tatebe, H., et al. (2019). *Description and basic evaluation of simulated mean state, internal variability, and climate sensitivity in MIROC6*. Geoscientific Model Development, **12**, 2727–2765. https://doi.org/10.5194/gmd-12-2727-2019

---

## Component-specific Citations

Several components of COCO implement numerical methods and physical parameterizations that have been described in scientific publications. If your work relies substantially on one or more of these components, please cite the corresponding reference(s) in addition to the software itself.

| Component / Feature                           | Recommended Citation         |
| --------------------------------------------- | ---------------------------- |
| COCO software                                 | GitHub Release (recommended) |
| COCO model description                        | Hasumi (2006)                |
| MIROC6 ocean component                        | Tatebe et al. (2019)         |
| Tidal mixing parameterization                 | Kawasaki et al. (2021)       |
| Tracer advection (Second-Order Moment scheme) | Tatebe and Hasumi (2010)     |
| Bottom boundary layer scheme                  | Nakano and Suginohara (2002) |

---

## Why cite a GitHub Release?

Unlike traditional publications, a GitHub Release identifies the exact version of the source code used in a study. This improves the reproducibility, transparency, and traceability of scientific results by providing a permanent record of the software implementation.

---

## Future DOI Support

If persistent identifiers (e.g., DOIs through Zenodo or similar services) become available for COCO releases in the future, we recommend citing the DOI associated with the corresponding release instead of the GitHub URL.

---

## Questions

If you are uncertain which references are most appropriate for your work, we recommend citing:

1. the GitHub Release corresponding to the version used; and
2. either **Hasumi (2006)** or **Tatebe et al. (2019)**, depending on your application;

together with the publications describing the model components or parameterizations that are central to your study.
