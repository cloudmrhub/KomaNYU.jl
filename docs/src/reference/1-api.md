# API Overview

```@raw html
<p align="center"><img width="100%" src="../../assets/koma-schema-subdirs.svg"/></p>
```

KomaNYU is divided into the following sub-packages:
- [KomaNYUBase](2-koma-base.md): Custom types and functions
- [KomaNYUCore](3-koma-core.md): Simulation functions
- [KomaNYUFiles](4-koma-files.md): File I/O functions
- [KomaNYUPlots](5-koma-plots.md): Plotting functions
- [KomaNYU](6-koma-mri.md): User Interface

The idea of separating the package into sub-packages is to make it easier to maintain and to allow users to use only the parts they need. Two common use-cases can be:
- **GUI users**: They will use the `KomaNYU` package to interact with the GUI. Internally this includes al the others.
- **Advanced users**: They will use the sub-packages directly to build their own scripts or notebooks, for simulation only `KomaNYUCore` is required.