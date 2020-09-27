
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// cameras/perspective_dist.cpp*
#include "cameras/perspective_dist.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "light.h"
#include "stats.h"

namespace pbrt {

// PerspectiveDistCamera Method Definitions
PerspectiveDistCamera::PerspectiveDistCamera(const AnimatedTransform &CameraToWorld,
                                     const Bounds2f &screenWindow,
                                     Float shutterOpen, Float shutterClose,
                                     Float lensRadius, Float focalDistance,
                                     Float fov, Film *film,
                                     const Medium *medium)
    : PerspectiveDistCamera(CameraToWorld, PerspectiveDist(fov, 1e-2f, 1000.f),
                       screenWindow, shutterOpen, shutterClose, lensRadius,
                       focalDistance, film, medium) {
    // Compute differential changes in origin for perspective camera rays
    dxCamera =
        (RasterToCamera(Point3f(1, 0, 0)) - RasterToCamera(Point3f(0, 0, 0)));
    dyCamera =
        (RasterToCamera(Point3f(0, 1, 0)) - RasterToCamera(Point3f(0, 0, 0)));

    // Compute image plane bounds at $z=1$ for _PerspectiveDistCamera_
    Point2i res = film->fullResolution;
    Point3f pMin = RasterToCamera(Point3f(0, 0, 0));
    Point3f pMax = RasterToCamera(Point3f(res.x, res.y, 0));
    pMin /= pMin.z;
    pMax /= pMax.z;
    A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));
}
// return distortion factor
Float applyDistortion(Float r_c) {
  //  Float kc[] = {0.5960, 0.3106, 0.1182}; // tmp turn this into params
    Float kc[] = {0.5960, 0.3106, 0.1182}; // tmp
    Float r2 = r_c*r_c;
    return 1 + r2*(kc[0] + r2*kc[1] + r2*r2*kc[2]);
}

Float applyDistortionAngle(Float theta){
    Float kc[] = {0.3812e-4, 0.0, 0.0};
    Float theta2 = theta*theta;
    return theta*(kc[0]+theta2*(kc[1]+theta2*kc[2]));
}
// return inverse distortion factor
Float invertDistortion(Float r_d) {
    Float kc[] = {-0.1688, 1.1354, -0.0165};
//    {-0.1688, 1.1354, -0.0165}; // tmp turn this into params
    int it = 0;
    Float r = r_d;
    while (true) {
        Float r2 = r*r,
              f  = r*(1 + r2*(kc[0] + r2*kc[1] + r2*r2*kc[2])) - r_d,
              df = 1 + r2*(3*kc[0] + 5*kc[1]*r2 + 7*kc[2]*r2*r2);

        r -= f / df;

        if (std::abs(f) < 1e-6 || ++it > 4)
            break;
    }
//    printf( "%6.4lf \n", r_d/r );
    return r_d/r;
//    return r/r_d;
}

Float PerspectiveDistCamera::GenerateRay(const CameraSample &sample,
                                     Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);
    // Apply distortion -----tmp zhenyi
    Float x_d = pCamera.x/pCamera.z;
    Float y_d = pCamera.y/pCamera.z;
    Float r_d = sqrt(x_d*x_d + y_d*y_d);

    Float correction_factor = applyDistortion(r_d);
    if (r_d == 0) {
        correction_factor=1;
    }
    Vector3f dir = Normalize(Vector3f(pCamera.x*correction_factor, pCamera.y*correction_factor, pCamera.z));
//    Float r = sqrt(pCamera.x*pCamera.x+ pCamera.y*pCamera.y + pCamera.z*pCamera.z);
//    Float psi = atan2(pCamera.y, pCamera.x);
//    Float theta = acos(pCamera.z/r);
//    Float imageHieght = applyDistortionAngle(theta);
//    Vector3f dir = Normalize(Vector3f(imageHieght*cos(psi), imageHieght*sin(psi), pCamera.z));
    //---------------------------------
    *ray = Ray(Point3f(0, 0, 0), dir);
//    *ray = Ray(Point3f(0, 0, 0), Normalize(Vector3f(pCamera)));
    //  tmp: apply distortion to ray direction here

    //    }
    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
}

Float PerspectiveDistCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);
//    Vector3f dir = Normalize(Vector3f(pCamera.x, pCamera.y, pCamera.z));
    //--------------------------------------------------------------------------
    // add distortion tmp --zhenyi
    // test without using a params input
    // the sampled pixel is the distorted pixel, so we need to inverse the distortion
    // to find the correct ray direction
    Float x_d = pCamera.x/pCamera.z;
    Float y_d = pCamera.y/pCamera.z;
    float r_d = sqrt(x_d*x_d + y_d*y_d);
    Float correction_factor = applyDistortion(r_d);
    if (r_d == 0) {
        correction_factor=1;
    }
//    printf( "%6.4lf \n", correction_factor );
    Vector3f dir = Normalize(Vector3f(pCamera.x*correction_factor, pCamera.y*correction_factor, pCamera.z));
//    Float r = sqrt(pCamera.x*pCamera.x+ pCamera.y*pCamera.y + pCamera.z*pCamera.z);
//    Float psi = atan2(pCamera.y, pCamera.x);
//    Float theta = acos(pCamera.z/r);
//    Float imageHieght = applyDistortionAngle(theta);
//    Vector3f dir = Normalize(Vector3f(theta*cos(psi), theta*sin(psi), pCamera.z));
    //---------------------------------------------------------------------------
    Float w = ray->wavelength; // Save the wavelength information
    *ray = RayDifferential(Point3f(0, 0, 0), dir);
    ray->wavelength = w; // Reassign
    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    // Compute offset rays for _PerspectiveDistCamera_ ray differentials
    if (lensRadius > 0) {
        // Compute _PerspectiveDistCamera_ ray differentials accounting for lens

        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Vector3f dx = Normalize(Vector3f(pCamera + dxCamera));
        Float ft = focalDistance / dx.z;
        Point3f pFocus = Point3f(0, 0, 0) + (ft * dx);
        ray->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

        Vector3f dy = Normalize(Vector3f(pCamera + dyCamera));
        ft = focalDistance / dy.z;
        pFocus = Point3f(0, 0, 0) + (ft * dy);
        ray->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
    } else {
        ray->rxOrigin = ray->ryOrigin = ray->o;
        ray->rxDirection = Normalize(Vector3f(pCamera) + dxCamera);
        ray->ryDirection = Normalize(Vector3f(pCamera) + dyCamera);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    ray->hasDifferentials = true;
    return 1;
}

Spectrum PerspectiveDistCamera::We(const Ray &ray, Point2f *pRaster2) const {
    // Interpolate camera matrix and check if $\w{}$ is forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    // apply distortion --zhenyi
//    dir = applyDistortion(<#Float r_c#>)
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) return 0;

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return raster position if requested
    if (pRaster2) *pRaster2 = Point2f(pRaster.x, pRaster.y);

    // Return zero importance for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y)
        return 0;

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;

    // Return importance for point on image plane
    Float cos2Theta = cosTheta * cosTheta;
    return Spectrum(1 / (A * lensArea * cos2Theta * cos2Theta));
}

void PerspectiveDistCamera::Pdf_We(const Ray &ray, Float *pdfPos,
                               Float *pdfDir) const {
    // Interpolate camera matrix and fail if $\w{}$ is not forward-facing
    Transform c2w;
    CameraToWorld.Interpolate(ray.time, &c2w);
    Float cosTheta = Dot(ray.d, c2w(Vector3f(0, 0, 1)));
    if (cosTheta <= 0) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Map ray $(\p{}, \w{})$ onto the raster grid
    Point3f pFocus = ray((lensRadius > 0 ? focalDistance : 1) / cosTheta);
    Point3f pRaster = Inverse(RasterToCamera)(Inverse(c2w)(pFocus));

    // Return zero probability for out of bounds points
    Bounds2i sampleBounds = film->GetSampleBounds();
    if (pRaster.x < sampleBounds.pMin.x || pRaster.x >= sampleBounds.pMax.x ||
        pRaster.y < sampleBounds.pMin.y || pRaster.y >= sampleBounds.pMax.y) {
        *pdfPos = *pdfDir = 0;
        return;
    }

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdfPos = 1 / lensArea;
    *pdfDir = 1 / (A * cosTheta * cosTheta * cosTheta);
}

Spectrum PerspectiveDistCamera::Sample_Wi(const Interaction &ref, const Point2f &u,
                                      Vector3f *wi, Float *pdf,
                                      Point2f *pRaster,
                                      VisibilityTester *vis) const {
    // Uniformly sample a lens interaction _lensIntr_
    Point2f pLens = lensRadius * ConcentricSampleDisk(u);
    Point3f pLensWorld = CameraToWorld(ref.time, Point3f(pLens.x, pLens.y, 0));
    Interaction lensIntr(pLensWorld, ref.time, medium);
    lensIntr.n = Normal3f(CameraToWorld(ref.time, Vector3f(0, 0, 1)));

    // Populate arguments and compute the importance value
    *vis = VisibilityTester(ref, lensIntr);
    *wi = lensIntr.p - ref.p;
    Float dist = wi->Length();
    *wi /= dist;

    // Compute PDF for importance arriving at _ref_

    // Compute lens area of perspective camera
    Float lensArea = lensRadius != 0 ? (Pi * lensRadius * lensRadius) : 1;
    *pdf = (dist * dist) / (AbsDot(lensIntr.n, *wi) * lensArea);
    return We(lensIntr.SpawnRay(-*wi), pRaster);
}

bool PerspectiveDistCamera::CanSample_Wi() const {
    return true;
}


PerspectiveDistCamera *CreatePerspectiveDistCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film->fullResolution.x) / Float(film->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    Float fov = params.FindOneFloat("fov", 90.);
    Float halffov = params.FindOneFloat("halffov", -1.f);
    if (halffov > 0.f)
        // hack for structure synth, which exports half of the full fov
        fov = 2.f * halffov;
    return new PerspectiveDistCamera(cam2world, screen, shutteropen, shutterclose,
                                 lensradius, focaldistance, fov, film, medium);
}

}  // namespace pbrt
