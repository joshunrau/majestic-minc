# Majestic MINC

This repository is a fork of [`minc-toolkit`](https://github.com/BIC-MNI/minc-toolkit-v2) v1.9.19. In addition, it also includes `conda-forge` recipes for a number of packages requiring native code that was previously unavailable for `osx-arm64` targets.

## Building Conda Package

### Create Environment
```shell
conda create -c minc-forge --name minc conda-build anaconda-client minc-toolkit-v2 r-base=4.1.3
```

### Compile Packages

```shell
conda-build -c minc-forge --R 4.1.3 conda/r-rminc
```