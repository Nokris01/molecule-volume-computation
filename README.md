# Molecular Volume Maps

A C++17 command-line tool that calculates Gaussian atom-density volumes from
PDB and mmCIF coordinate files. It can write legacy VTK maps, OpenDX maps, and
per-element MVOL density channels. The project builds directly with `g++`;
CMake is not required.

## Features

- Reads `.pdb` and coordinate-bearing `.cif` files.
- Calculates atom-number density on a configurable 3D grid.
- Writes legacy VTK output, with optional OpenDX (`--dx`) and MVOL (`--mvol`)
  output.
- Uses OpenMP to parallelize voxel calculation.
- Includes GLM in `glm/`, so GLM does not need to be installed separately.
- Can write additional per-voxel atom-number density data, including per-element channels (--mvol).

## Requirements

- A C++17-capable GCC toolchain (`g++`).
- OpenMP support in that toolchain.
- The GLM headers included in this repository.

## Build

Run the command for the target platform from the project root. For a release
binary, build for the target operating system and CPU architecture; binaries
are not generally portable between platforms.

### Linux x86-64

```sh
g++ -std=c++17 -O3 -fopenmp -I. main.cpp -o volume_map
```

### Linux ARM64 (AArch64)

On an ARM64 Linux machine, use the same command:

```sh
g++ -std=c++17 -O3 -fopenmp -I. main.cpp -o volume_map
```

This produces a Linux ARM64 executable when `g++` is the ARM64 compiler. For
release builds, avoid `-march=native`: that can make the binary depend on CPU
features specific to the machine that built it. To cross-compile from another
Linux architecture, use an AArch64 cross-compiler and ensure its OpenMP runtime
is available on the target system, for example:

```sh
aarch64-linux-gnu-g++ -std=c++17 -O3 -fopenmp -I. main.cpp -o volume_map
```

### Windows x64

Use a MinGW-w64 GCC toolchain with OpenMP support, such as the MSYS2 UCRT64
environment. In its UCRT64 terminal, install GCC if needed and compile:

```sh
pacman -S mingw-w64-ucrt-x86_64-gcc
g++ -std=c++17 -O3 -fopenmp -I. main.cpp -o volume_map.exe
```

Run the executable from PowerShell or Command Prompt:

```powershell
.\volume_map.exe
```

The Windows build may require GCC runtime DLLs, including the OpenMP runtime,
to be installed or included alongside the executable. Test the release on a
clean Windows machine; do not assume that an `.exe` built on Linux will run on
Windows.

### Linux AppImage

The repository includes an `AppDir/` staging directory for the Linux x86-64
command-line app. Build the executable into it from the project root:

```sh
mkdir -p AppDir/usr/bin
g++ -std=c++17 -O3 -fopenmp -I. main.cpp -o AppDir/usr/bin/volume_map
```

Install `linuxdeploy` and its AppImage output plugin, then deploy the executable
and create the AppImage:

```sh
./linuxdeploy-x86_64.AppImage \
  --appdir AppDir \
  --executable AppDir/usr/bin/volume_map \
  --desktop-file AppDir/usr/share/applications/volume-map.desktop \
  --icon-file AppDir/usr/share/icons/hicolor/scalable/apps/volume-map.svg \
  --output appimage
```

The generated AppImage is for Linux x86-64 only. Build it on an older
glibc-based distribution (or compatible build environment) to maximize
compatibility with older Linux distributions. `linuxdeploy` bundles the
application's non-system shared-library dependencies; glibc is provided by the
host system and must not be bundled. Test the output on target distributions.

## Run

By default, the program reads each `.pdb` and `.cif` file in `batch/` and writes
one legacy VTK structured-points volume (`.vtk`) per input to `results/`.
Override the directories and density parameters with:

```sh
./volume_map INPUT_DIR OUTPUT_DIR \
  --spacing 0.4 --sigma 0.8 --cutoff 4
```

On Windows, replace `./volume_map` with `.\volume_map.exe`.

`--spacing` sets voxel size in Å. Smaller spacing gives finer spatial
resolution and substantially larger/slower files; dimensions are derived from
the structure bounds and the kernel padding. For example, reducing spacing
from `0.5` to `0.4` Å gives about 1.95 times as many voxels for the same
volume. Increase spacing for smaller maps.

`--sigma` controls Gaussian smoothing, not grid resolution. A smaller sigma
makes narrower peaks and increases their maximum values; a larger sigma makes
the density smoother and lowers the peaks. `--cutoff` is the number of sigma
used to truncate each kernel and pad the grid, not an intensity control. Keep
it around `4` unless you specifically want to truncate more of the Gaussian
tails.

The Gaussian calculation is parallelized over exclusive 16×16×16 voxel
bricks. Each voxel is owned by one worker, so overlapping atom kernels no
longer contend on atomic updates. Gaussian factors along each axis are
precomputed per atom, avoiding an exponential evaluation for every voxel.
Grid values use x-contiguous storage to improve calculation and VTK/MVOL
output locality; DX serialization uses the order expected by the renderer's
loader. Larger bricks were benchmarked but did not improve this workload.
OpenMP uses the available thread count by default; override it with
`OMP_NUM_THREADS`, for example:

```sh
OMP_NUM_THREADS=8 ./volume_map --dx
```

The program reports the OpenMP thread count at startup.

For visualization-only intensity adjustment, use `--density-scale FACTOR`.
This multiplies all density values written to VTK, DX, and MVOL without
changing grid geometry or smoothing. It defaults to `1` and must be positive.
For example, `--density-scale 10` makes values ten times larger. This is a
display scale, not a change in the underlying atom-density model; the DX
header records the factor and labels its density values as scaled.

Pass `--dx` to also convert each VTK map generated in that run to an OpenDX
`.dx` map after the batch finishes. The `.vtk` files are retained:

```sh
./volume_map --dx
```

Pass `--mvol` to write an additional versioned `.mvol` file per structure:

```sh
./volume_map --mvol
```

Each `.mvol` contains a total-density channel and a separate additive density
channel for every chemical element present in that structure (for example,
`C`, `H`, `N`, `O`, and `S`). At each voxel, a type channel is the sum of the
Gaussian contributions from all selected atoms of that element. The channels
are atom-number densities weighted by occupancy, not normalized probabilities
between zero and one. The sum of element channels equals the total density up
to float32 rounding.

### MVOL binary format, version 1.0

All integers and floating-point values are little-endian. Floating-point
values use IEEE-754 binary32 or binary64. There is no compression. The file is
written as a fixed header, channel-name table, then one dense float32 array
per channel.

| Byte offset | Type | Field |
|---:|---|---|
| 0 | 4 bytes | Magic ASCII `MVOL` |
| 4 | u16 | Major version (`1`) |
| 6 | u16 | Minor version (`0`) |
| 8 | u32 | Header size in bytes, including the channel-name table |
| 12 | u32 × 3 | Grid dimensions `nx, ny, nz` |
| 24 | float64 × 3 | Position of the first voxel center in Å (`x, y, z`) |
| 48 | float64 × 3 | Grid spacing in Å (`dx, dy, dz`) |
| 72 | u32 | Density semantics (`1` = Gaussian atom-number density) |
| 76 | u32 | Voxel data type (`1` = IEEE-754 float32) |
| 80 | u32 | Voxel layout (`1` = channel-major, x-index fastest) |
| 84 | u32 | Number of channels, including total density |
| 88 | u64 | Number of voxels per channel (`nx * ny * nz`) |
| 96 | variable | Channel-name table |

The channel-name table has one entry per channel, in data order. Each entry
is a one-byte unsigned name length followed by that many ASCII bytes. The
first entry is `density`; subsequent entries are distinct element symbols,
sorted lexicographically (for example, `C`, `H`, `N`, `O`, `S`). Thus
`header_size` is `96` plus the size of all length-prefixed names, and the
voxel payload begins at that offset.

The payload contains `channel_count` consecutive arrays, each with
`voxel_count` little-endian float32 values. Within each array, x varies
fastest, then y, then z. Coordinates are at voxel centers:
`origin + (x * dx, y * dy, z * dz)`. The `density` channel and element
channels use the same fixed sigma, cutoff, and occupancy weighting as the VTK
output. A reader should validate the version, header size, dimensions,
channel-name lengths, and expected file length before accessing voxel data.

Coordinates and spacing are in angstroms. Each atom contributes a normalized,
isotropic Gaussian with the specified fixed sigma; atom occupancy is its
weight. The grid covers all selected atom centers plus `cutoff * sigma`
padding. VTK data is binary, big-endian 32-bit floating-point density and
sampled at voxel centers.

The optional OpenDX output uses the same grid origin, spacing and scalar
density values as its corresponding VTK file. It uses the text-based
gridpositions/gridconnections format and named field footer. DX values are
written with the last grid index (`z`) varying fastest, matching the OpenDX
array ordering and the renderer's `loadDx` index mapping. This differs from
the VTK payload order. Density values are Gaussian atom-number density in
Å⁻³; they are not rescaled for display unless `--density-scale` is specified.
Choose the renderer's contour/transfer-function range for the map's reported
density range.

PDB/mmCIF coordinate files do not define a separate world-space `origin`;
their atom coordinates are already expressed in Cartesian Å. The DX `origin`
is the coordinate of the first grid point, placed at the lowest atom
coordinate minus the Gaussian padding. It is therefore usually not zero and
is not expected to equal an origin field from the structure. Both the atom
coordinates and the map grid remain in the same Cartesian frame, so load them
without independently translating or recentering either one.

Only the first model is used. For alternative atom locations, the location
with the highest occupancy is selected (ties prefer a blank location, then
the lexicographically earlier location). ATOM and HETATM records, including
hydrogen atoms and waters, are included. B-factors are not used to set sigma,
and assembly operations are not applied; an already expanded assembly file
uses the coordinates in that file.

Inputs with no coordinate-bearing atom-site records are reported as failures.
In the supplied batch, `9TNI-sf.cif` is a structure-factor file rather than a
coordinate model, so it has no volume map to calculate. Compressed `.cif.gz`
files are not selected; the batch contains an uncompressed assembly CIF too.
The program continues after per-file failures, reports the batch summary, and
returns a nonzero exit code if any input could not be converted.
