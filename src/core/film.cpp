
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


// core/film.cpp*
#include "film.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"

// We need this for writing out an appropriate .dat file for the spectral image:
#include <sstream>
#include <iostream>
#include <fstream>

namespace pbrt {
    
STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);
    
// Film Method Definitions
Film::Film(const Point2i &resolution, const Bounds2f &cropWindow,
               std::unique_ptr<Filter> filt, Float diagonal,
               const std::string &filename, Float scale, bool sf, std::string dt, Float maxSampleLuminance)
    : fullResolution(resolution),
      diagonal(diagonal * .001),
      filter(std::move(filt)),
      filename(filename),
      scale(scale),
      maxSampleLuminance(maxSampleLuminance),
      spectralFlag(sf),
      datatype(dt){
    // Compute film image bounds
    croppedPixelBounds =
        Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                         std::ceil(fullResolution.y * cropWindow.pMin.y)),
                 Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                         std::ceil(fullResolution.y * cropWindow.pMax.y)));
    LOG(INFO) << "Created film with full resolution " << resolution <<
        ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
        croppedPixelBounds;

    // Allocate film image storage
    pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
    filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

    // Precompute filter weight table
    int offset = 0;
    for (int y = 0; y < filterTableWidth; ++y) {
        for (int x = 0; x < filterTableWidth; ++x, ++offset) {
            Point2f p;
            p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
            p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
            filterTable[offset] = filter->Evaluate(p);
        }
    }
}
    
Bounds2i Film::GetSampleBounds() const {
    Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
                                Vector2f(0.5f, 0.5f) - filter->radius),
                            Ceil(Point2f(croppedPixelBounds.pMax) -
                                Vector2f(0.5f, 0.5f) + filter->radius));
    return (Bounds2i)floatBounds;
}
    
Bounds2f Film::GetPhysicalExtent() const {
    Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}
    
std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
    Point2i(1, 1);
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
    return std::unique_ptr<FilmTile>(new FilmTile(
                                                    tilePixelBounds, filter->radius, filterTable, filterTableWidth,
                                                    maxSampleLuminance));
}
    
void Film::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int c = 0; c < 3; ++c)
            pixel.splatXYZ[c] = pixel.xyz[c] = 0;
        pixel.filterWeightSum = 0;
        pixel.L = 0;
    }
}
    
void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->pixelBounds;
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
        for (int i = 0; i < 3; ++i) mergePixel.xyz[i] += xyz[i];
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
        // If we are using a spectral film, save also save values directly
        if(spectralFlag){
            mergePixel.L += tilePixel.contribSum;
        }
            
    }
}
    
void Film::SetImage(const Spectrum *img) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        img[i].ToXYZ(p.xyz);
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}
    
void Film::AddSplat(const Point2f &p, Spectrum v) {
    ProfilePhase pp(Prof::SplatFilm);
        
    if (v.HasNaNs()) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with NaN values "
                                    "at (%f, %f)", p.x, p.y);
        return;
    } else if (v.y() < 0.) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with negative "
                                    "luminance %f at (%f, %f)", v.y(), p.x, p.y);
        return;
    } else if (std::isinf(v.y())) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with infinite "
                                    "luminance at (%f, %f)", p.x, p.y);
        return;
    }
        
    if (!InsideExclusive((Point2i)p, croppedPixelBounds)) return;
    if (v.y() > maxSampleLuminance)
        v *= maxSampleLuminance / v.y();
    Float xyz[3];
    v.ToXYZ(xyz);
    Pixel &pixel = GetPixel((Point2i)p);
    for (int i = 0; i < 3; ++i) pixel.splatXYZ[i].Add(xyz[i]);
}
    
void Film::WriteImage(Float splatScale) {
        
    if(!spectralFlag){
        // This is the standard way of writing out an image into RGB format.
            
        // Convert image to RGB and compute final pixel values
        LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
        std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
        int offset = 0;
        for (Point2i p : croppedPixelBounds) {
            // Convert pixel XYZ color to RGB
            Pixel &pixel = GetPixel(p);
            XYZToRGB(pixel.xyz, &rgb[3 * offset]);
                
            // Normalize pixel with weight sum
            Float filterWeightSum = pixel.filterWeightSum;
            if (filterWeightSum != 0) {
                Float invWt = (Float)1 / filterWeightSum;
                rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
                rgb[3 * offset + 1] =
                std::max((Float)0, rgb[3 * offset + 1] * invWt);
                rgb[3 * offset + 2] =
                std::max((Float)0, rgb[3 * offset + 2] * invWt);
            }
                
            // Add splat value at pixel
            Float splatRGB[3];
            Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                pixel.splatXYZ[2]};
            XYZToRGB(splatXYZ, splatRGB);
            rgb[3 * offset] += splatScale * splatRGB[0];
            rgb[3 * offset + 1] += splatScale * splatRGB[1];
            rgb[3 * offset + 2] += splatScale * splatRGB[2];
                
            // Scale pixel value by _scale_
            rgb[3 * offset] *= scale;
            rgb[3 * offset + 1] *= scale;
            rgb[3 * offset + 2] *= scale;
            ++offset;
        }
            
        // Write RGB image
        LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
        pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
    } else {
        // Otherwise, we write it out as a multispectral image with no filter weighting.
                     
        // spectralData holds all the values in the multispectral image
        
        std::unique_ptr<Float[]> spectralData(new Float[nSpectralSamples * croppedPixelBounds.Area()]);
        
        // Choose write out number of samples based on datatype. --Zhenyi
        int nDataSamples = 31;
        if (datatype.compare("depth")==0 || datatype.compare("mesh")==0
            || datatype.compare("material")==0) {
            nDataSamples = 1;}
        else if (datatype.compare("coordinates")==0) {
            nDataSamples = 3;}
        // lidar output: [x, y, z, reflectance, irradiance, instance ID]
        else if (datatype.compare("pointcloud")==0) {
            nDataSamples = 6;}

        int offset = 0;
        for (Point2i p : croppedPixelBounds) {
                
            // Get spectrum directly
            Pixel &pixel = GetPixel(p);
            Spectrum currSpectrum = pixel.L;
            //Spectrum currSpectrum = Spectrum::FromXYZ(pixel.xyz);

            
            // Loop through the current spectrum and put each value into spectralData
            for(int i = 0; i < nDataSamples; i++){
                Float filterWeightSum = pixel.filterWeightSum;
                Float invWt = (Float)1 / filterWeightSum;
                currSpectrum[i] = currSpectrum[i] * invWt;
                spectralData[offset*nDataSamples + i] = currSpectrum[i];
            }
                
            // We're not going to weight the output data with the filter. In other words, the more rays, the higher the output value will be. We can scale to an appropriate luminance later in ISET.
                
            // Do we even need splatting? What exactly is it?
            // Add splat value at every spectral sample
            Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                pixel.splatXYZ[2]};
            Spectrum splatSpectrum = Spectrum::FromXYZ(splatXYZ);
                
            for(int i = 0; i < nDataSamples; i++){
                spectralData[offset*nDataSamples + i] += splatScale * splatSpectrum[i];
                spectralData[offset*nDataSamples + i] *= scale;
            }
                
            ++offset;
        }
            
        // TODO: Turn this into function analogous to pbrt::WriteImage
            
        // Write RGB image
        LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
            
        std::ofstream myfile;
        int extPos = filename.find_last_of(".");
        std::string datFilename = filename.substr(0,extPos) + ".dat"; // Filename is now going to be xxx.dat
        myfile.open(datFilename.c_str());
        // Print out dimensions of the image
        myfile << croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x << " " << croppedPixelBounds.pMax.y - croppedPixelBounds.pMin.y << " " << nDataSamples << "\n";
            
        // Print out focal length and field of view information
        // TODO: Is it possible for us to do this here? We don't have access to focal length in the film class...
        // myfile << focalLength << " " << fStop << " " << fieldOfView << "\n";
            
        // Note: To be compatible with pbrt-v2-spectral we need to write image out column by column, e.g.
        // 1 2 3
        // 4 5 6
        // 7 8 9
        // Needs to be written out as [1 4 7 2 5 8 3 6 9]
        // TODO: How do we do this? Careful indexing I imagine...
        // For now let's just change piReadDAT so it can read it correctly. We therefore print out a v3 flag so piReadDat knows what to do.
        myfile << "v3 \n";
            
        myfile.close();
            
        //Open file for binary writing
        FILE * spectralDataBin;
        spectralDataBin = fopen(datFilename.c_str(), "a");
            
        //Write binary image
        
        
        for (int i = 0; i < nDataSamples; i++)
        {
            for (int j = 0; j < croppedPixelBounds.Area(); j++)
            {
                double r = (double)spectralData[nDataSamples * j + i];
                fwrite((void*)(&r), sizeof(r), 1, spectralDataBin);
            }
        }
            
        // piReadDAT expects the data to be serialized wavelength by wavelength. In other words, you would go through all the rows and columns of wavelength index = 1 first, then move on to all the pixels for wavelength index = 2 next, etc.

        fclose(spectralDataBin);
            
    }
            
}
    
    
Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    // Intentionally use FindOneString() rather than FindOneFilename() here
    // so that the rendered image is left in the working directory, rather
    // than the directory the scene file lives in.
    std::string filename = params.FindOneString("filename", "");
    if (PbrtOptions.imageFile != "") {
        if (filename != "") {
            Warning(
                    "Output filename supplied on command line, \"%s\", ignored "
                    "due to filename provided in scene description file, \"%s\".",
                    PbrtOptions.imageFile.c_str(), filename.c_str());
        } else
            filename = PbrtOptions.imageFile;
    }
    if (filename == "") filename = "pbrt.exr";
        
    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop(Point2f(0, 0), Point2f(1, 1));
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);
        
    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                    Infinity);
    
    bool spectralFlag = params.FindOneBool("spectralFlag", true); //delete me
    
    std::string datatype = params.FindOneString("datatype", "spectral");
    
    return new Film(Point2i(xres, yres), crop, std::move(filter), diagonal,
                    filename, scale, spectralFlag, datatype, maxSampleLuminance);
}
    
}  // namespace pbrt
