# Double-Slit Diffraction (Rayleigh-Sommerfeld, Fresnel and Fraunhofer)

MATLAB simulation of the diffraction of a monochromatic plane wave passing
through a **rectangular double slit**. The intensity is computed by numerically
integrating the Rayleigh-Sommerfeld diffraction integral and compared with the
**Fresnel** and **Fraunhofer** approximations.

## Physical model

Diffracted field at a point `(x, y, z)` (Rayleigh-Sommerfeld integral):

    U(x,y,z) = −1/(2π) ∬ [ z · e^(i·k·n·r) / r² ] · ( i·k·n − 1/r ) dx' dy'

    r = √( (x − x')² + (y − y')² + z² )

The integral is evaluated over the two rectangular apertures, and the
intensity is `I = |U|²`. In the plots it is normalized by the intensity at the
center, `I₀ = I(x = 0)`.

**Approximations compared:**

- **Fresnel:** paraxial expansion of `r` up to second order.
- **Fraunhofer:** far field, which reduces the problem to a Fourier transform
  of the aperture.

## Geometry and parameters

| Parameter | Default | Meaning |
|---|---|---|
| `lambda` | 650 nm | Wavelength |
| `L` | 0.015 mm | Width of each slit |
| `h` | 0.02 mm | Height of each slit |
| `d` | 0.06 mm | Distance between slit centers |
| `z` | 1 m | Distance to the screen |
| `n` | 1 | Refractive index |
| `AbsTol` | 1e-4 | Absolute tolerance of `integral2` |

## What it does

The script is organized in six blocks:

1. **Method comparison:** intensity profile along `x` (with `y = 0`) using
   Rayleigh-Sommerfeld, Fresnel and Fraunhofer, in absolute and normalized
   form (`I/I₀` as a function of `d·x/(z·λ)`).
2. **Varying `z`:** effect of bringing the screen closer to the slits.
3. **Varying `d`:** effect of the slit separation.
4. **Varying `L`:** effect of the slit width.
5. **Varying `λ`:** effect of the wavelength.
6. **3D representation:** `I/I₀` intensity map on the `(x, y)` plane using
   Rayleigh-Sommerfeld (400 × 400 grid).

The computation is **parallelized** with `parfor`, distributing the
observation points among the available cores.

## Usage

1. Open MATLAB and move to the project folder.
2. Run:

```matlab
practica1codigo
```

3. Edit the parameters at the top of the file (`p` and `config` structures) to
   explore other cases.

A pool of parallel workers (`parpool`) is started on the first run. With
`config.MAXIMIZAR_IMAGEN = 1`, figures are maximized before being saved.

> **Note:** the computation is expensive, especially the 3D map (160,000
> points, two double integrals each). For a quick test, reduce the resolution
> of `xpoints3D` / `ypoints3D` or the number of points in `xpoints`.

## Generated files

Figures are saved to the working folder as `.png` files:

| File | Contents |
|---|---|
| `1evaluacioninicial.png` | Rayleigh-Sommerfeld intensity |
| `2intensidadnormalizada.png` | Normalized intensity |
| `3comparacionmetodos_x.png` | Method comparison vs. `x` |
| `3comparacionmetodos_norm.png` | Method comparison vs. `d·x/(z·λ)` |
| `evaluacion_z_1..3.png` | Sweep in `z` |
| `evaluacion_d_1..3.png` | Sweep in `d` |
| `evaluacion_L_1..3.png` | Sweep in `L` |
| `evaluacion_lambda_1..3.png` | Sweep in `λ` |
| `3intensidadtridimensional.png` | 3D intensity map |

<img width="3000" height="1505" alt="3comparacionmetodos_x" src="https://github.com/user-attachments/assets/ce1894ba-b377-4f8d-99bf-0cf0725093ed" />


## Requirements

- MATLAB R2019b or later (uses the `WindowState` property).
- **Parallel Computing Toolbox** (`parpool`, `parfor`).
