
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

// Added by Zhenyi for Lidar Simulation
// lights/laser.cpp*
#include "lights/laser.h"
#include "paramset.h"
#include "sampling.h"
#include "reflection.h"
#include "stats.h"

namespace pbrt {

// LaserLight Method Definitions
LaserLight::LaserLight(Transform &LightToWorld,
                     const MediumInterface &mediumInterface, const Spectrum &I,
                     Float totalWidth, Float falloffStart, Transform LaserToWorld,
                       Transform WorldToLaser)
    : Light((int)LightFlags::DeltaPosition, LightToWorld, mediumInterface),
      pLight(LaserToWorld(Point3f(0, 0, 0))),
      I(I),
      cosTotalWidth(std::cos(Radians(totalWidth))),
      cosFalloffStart(std::cos(Radians(falloffStart))) {}

Spectrum LaserLight::Sample_Li(const Interaction &ref, const Point2f &u,
                              Vector3f *wi, Float *pdf,
                              VisibilityTester *vis) const {
    ProfilePhase _(Prof::LightSample);
    *wi = Normalize(pLight - ref.p);
    *pdf = 1.f;
    *vis =
        VisibilityTester(ref, Interaction(pLight, ref.time, mediumInterface));
    return I * Falloff(-*wi) / DistanceSquared(pLight, ref.p);
}

Float LaserLight::Falloff(const Vector3f &w) const {
    Vector3f wl = Normalize(WorldToLaser(w));
    Float cosTheta = wl.z;
    if (cosTheta < cosTotalWidth) return 0;
    if (cosTheta >= cosFalloffStart) return 1;
    // Compute falloff inside Laserlight cone
    Float delta =
        (cosTheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);
    return (delta * delta) * (delta * delta);
}

Spectrum LaserLight::Power() const {
    return I * 2 * Pi * (1 - .5f * (cosFalloffStart + cosTotalWidth));
}

Float LaserLight::Pdf_Li(const Interaction &, const Vector3f &) const {
    return 0.f;
}

Spectrum LaserLight::Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                              Ray *ray, Normal3f *nLight, Float *pdfPos,
                              Float *pdfDir) const {
    ProfilePhase _(Prof::LightSample);
    Vector3f w = UniformSampleCone(u1, cosTotalWidth);
    *ray = Ray(pLight, LaserToWorld(w), Infinity, time, mediumInterface.inside);
    *nLight = (Normal3f)ray->d;
    *pdfPos = 1;
    *pdfDir = UniformConePdf(cosTotalWidth);
    return I * Falloff(ray->d);
}

void LaserLight::Pdf_Le(const Ray &ray, const Normal3f &, Float *pdfPos,
                       Float *pdfDir) const {
    ProfilePhase _(Prof::LightPdf);
    *pdfPos = 0;
    *pdfDir = (CosTheta(WorldToLaser(ray.d)) >= cosTotalWidth)
                  ? UniformConePdf(cosTotalWidth)
                  : 0;
}
    
void LaserLight::SetLaserToWorld(Point3f &newDir, Point3f &newFrom) {
    // for a laser, newFrom defines the origin, newDir is the dir;
    Vector3f dir = Normalize(Vector3f(newDir));
    Vector3f du, dv;
    CoordinateSystem(dir, &du, &dv);
    Transform dirToZ =
        Transform(Matrix4x4(du.x, du.y, du.z, 0., dv.x, dv.y, dv.z, 0., dir.x,
                            dir.y, dir.z, 0., 0, 0, 0, 1.));
    Transform currentTranform = Transform();
    LaserToWorld =
        currentTranform * Translate(Vector3f(newFrom.x, newFrom.y, newFrom.z))* Inverse(dirToZ);
    WorldToLaser = Inverse(LaserToWorld);
    pLight = LaserToWorld(Point3f(0, 0, 0));
}
    
//std::unique_ptr<Light> LaserLight::Clone(int Seed) {
//    return std::unique_ptr<Light>(new LaserLight(*this));
//}

    
std::shared_ptr<LaserLight> CreateLaserLight(const Transform &l2w,
                                           const Medium *medium,
                                           const ParamSet &paramSet) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Float coneangle = paramSet.FindOneFloat("coneangle", 30.);
    Float conedelta = paramSet.FindOneFloat("conedeltaangle", 0);
    // Compute Laserlight world to laser transformation
    Transform LightToWorld = l2w ;
    Transform LaserToWorld = LightToWorld;
    Transform WorldToLaser = Inverse(LaserToWorld);
    return std::make_shared<LaserLight>(LightToWorld, medium, I * sc, coneangle,
                                       coneangle - conedelta, LaserToWorld, WorldToLaser);

}

}  // namespace pbrt
