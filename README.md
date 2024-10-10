# Rend-a-Pixel (Raytracer)

**Rend-a-Pixel** is a raytracing rendering engine developed on top of the **Lightwave** Framework as the final project for the [Computer Graphics course at Saarland University](https://graphics.cg.uni-saarland.de/) lectured by [Prof. Dr.-Ing. Philipp Slusallek](https://graphics.cg.uni-saarland.de/people/slusallek.html) during the Winter Semester 2023/2024. Some of the implemented features are showcased below:

## Area Lights 

<table>
<tr>
  <td align="center">No Area Lights</td>
  <td align="center">Uniform Sphere Sampling</td>
</tr>
<tr>
    <td><img src="images/area_lights/no_area_lights_ref.jpeg"</td>
    <td><img src="images/area_lights/area_lights_sphere_uniform_sampling.jpeg"</td>
</tr>
  <tr>
  <td align="center">Cosine-Weighted Sampling</td>
  <td align="center">Subtended-Cone Sampling</td>
</tr>
<tr>
    <td><img src="images/area_lights/area_lights_sphere_cosine_weighted_sampling.jpeg"</td>
    <td><img src="images/area_lights/area_lights_sphere_subtended_cone_sampling.jpeg"</td>
</tr>
</table>

## Multiple Importance Sampling

<table>
<tr>
  <td align="center">BSDF Sampling</td>
  <td align="center">NEE</td>
  <td align="center">MIS</td>
</tr>
<tr>
  <td><img src="images/mis_path_tracer/veach_bsdf.jpeg"</td>
  <td><img src="images/mis_path_tracer/veach_nee.jpeg"</td>
  <td><img src="images/mis_path_tracer/veach_mis.jpeg"</td>
</tr>
</table>

## Shading Normals

<table>
<tr>
  <td align="center">No Normal Mapping</td>
  <td align="center">Normal Mapping</td>
</tr>
<tr>
    <td><img src="images/shading_normals/no_normal_mapping.jpeg"</td>
    <td><img src="images/shading_normals/normal_mapping.jpeg"</td>
</tr>
</table>

## A Thinlens Camera Model

<table>
<tr>
  <td align="center">Perspective</td>
  <td align="center">With Thinlens</td>
</tr>
<tr>
    <td><img src="images/thinlens_camera_model/no_thinlens_ducks.jpg"</td>
    <td><img src="images/thinlens_camera_model/thinlens_ducks.jpg"</td>
</tr>
</table>

## Alpha Masking

<table>
<tr>
  <td align="center">Without Alpha Masking</td>
  <td align="center">With Alpha Masking</td>
</tr>
<tr>
    <td><img src="images/alpha_masking/no_alpha_masking.jpeg"</td>
    <td><img src="images/alpha_masking/alpha_masking.jpeg"</td>
</tr>
</table>

## Custom Bokeh Shapes

<table>
<tr>
  <td align="center">Sphere-looking Lights</td>
  <td align="center">3-bladed Bokeh Shapes</td>
</tr>
<tr>
    <td><img src="images/custom_bokeh_shapes/simple_bokeh_example_without_bokeh.jpeg"</td>
    <td><img src="images/custom_bokeh_shapes/simple_bokeh_example_test.jpeg"</td>
</tr>
</table>

## Image Denoising

<div align="center">
  <table width="100%">
    <tr>
      <td align="center">Noisy</td>
      <td align="center">Denoised</td>
    </tr>
    <tr>
        <td><img src="images/denoising/noisy_pathtracing_lights.jpeg"</td>
        <td><img src="images/denoising/denoised_pathtracing_lights.jpeg"</td>
    </tr>
  </table>
</div>

## Halton Sampler

<table>
<tr>
  <td align="center">Independent Sampling</td>
  <td align="center">Normal Halton Sampling</td>
</tr>
<tr>
    <td><img src="images/halton_sampler/bunny_constant_independent.jpeg"</td>
    <td><img src="images/halton_sampler/bunny_constant_halton_no_permutation.jpeg"</td>
</tr>
  <tr>
  <td align="center">Digit-permutated Halton Sampling</td>
  <td align="center">Owen-scrambled Halton Sampling</td>
</tr>
<tr>
    <td><img src="images/halton_sampler/bunny_constant_halton_digit_permutations.jpeg"</td>
    <td><img src="images/halton_sampler/bunny_constant_halton_owen_scramble.jpeg"</td>
</tr>
</table>

### Halton Sampler (In detail)

<table>
<tr>
  <td align="center">Independent Sampling</td>
  <td align="center">Normal Halton Sampling</td>
</tr>
<tr>
    <td><img src="images/halton_sampler/bunny_constant_detail_independent.png"</td>
    <td><img src="images/halton_sampler/bunny_constant_detail_halton.png"</td>
</tr>
  <tr>
  <td align="center">Digit-permutated Halton Sampling</td>
  <td align="center">Owen-scrambled Halton Sampling</td>
</tr>
<tr>
    <td><img src="images/halton_sampler/bunny_constant_detail_permutedigits.png"</td>
    <td><img src="images/halton_sampler/bunny_constant_detail_owen.png"</td>
</tr>
</table>

## Summary of Features
- [x] Camera Models
  - [x] Basic Perspective Camera
  - [x] Thinlens Camera
- [x] Basic Primitives
  - [x] Sphere
  - [x] Rectangles
  - [x] Triangle/Generic Meshes
- [x] Integrators
  - [x] Albedo
  - [x] Normals
  - [x] Direct Lighting
  - [x] Path Tracing 
- [x] BSDFs & Lighting Models:
  - [x] Materials: 
    - [x] Diffuse
    - [x] Conductor
    - [x] Rough Conductor
    - [x] Dielectric
    - [x] Principled
  - [x] Lambertian Emission
- [x] Textures:
  - [x] Checkerboard Texture
  - [x] Image Texture
- [x] Lights:
  - [x] Environment Map
    - [ ] Improved Environment Sampling
  - [x] Area Lights
    - [x] Uniform Sphere Sampling
    - [x] Cosine-Weighted Sampling
    - [x] Subtended-Cone Sampling  
  - [x] Point Light
  - [x] Directional Light
- [x] Sampling:
  - [x] BSDF Sampling 
  - [x] Next Event Estimation (NEE)
  - [x] Multiple Importance Sampling (MIS) 
- [x] Image denoising using [Intel&reg; Open Image Denoise](https://www.openimagedenoise.org/)
- [x] Acceleration Structures:
  - [x] SAH Bounding Volume Hierarchy
- [x] Shading Normals
- [x] Alpha Masking
- [x] Custom Bokeh Shapes

## Copyright & Credits
&copy; The Lightwave Framework was written by [Alexander Rath](https://graphics.cg.uni-saarland.de/people/rath.html), with contributions from [Ömercan Yazici](https://graphics.cg.uni-saarland.de/people/yazici.html) and [Philippe Weier](https://graphics.cg.uni-saarland.de/people/weier.html). Their support was invaluable in the coding of these features. The scenes showcasing the features were provided by their team, and should be used under permission. Many textures and models were taken from [Poly Haven](https://polyhaven.com.)'s extensive library.
                    
## References

<ol style="padding-left: 80px;">
  <li>Tomas Akenine-Mller, Eric Haines, and Naty Hoffman. <i>Real-Time Rendering</i> (4th ed.) USA: A. K. Peters, Ltd., 2018. ISBN: 0134997832. 
  <li>John F. Hughes et al. <i>Computer Graphics - Principles and Practice</i>, (3rd ed.) Addison-Wesley, 2014. ISBN: 978-0-321-39952-6.
  <li>Matt Pharr, Wenzel Jakob, and Greg Humphreys. <i>Physically Based Rendering: From Theory to Implementation</i>, (3rd ed.) Morgan Kaufmann Publishers Inc., 2016, ISBN: 978-0128006450.
  <li>Matt Pharr, Wenzel Jakob, and Greg Humphreys. <i>Physically Based Rendering: From Theory to Implementation</i>, (4th ed.) MIT Press, 2023, ISBN: 978-0262048026.</li>
</ol>
