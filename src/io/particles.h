#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <map>
#include <vector>


class particles
{


public:

    size_t nParticles, nSteps, nVertices, nVertexAllLayers;
    std::vector<double> xPos, yPos, zPos;
    std::vector<double> xVertex, yVertex, zVertex,
                    velocityVx, velocityVy, velocityVz,
                    zTopVertex, zLevelParticle;
    std::vector<int> trace_sizes_global;

    particles();

    void load_particle_data(const std::string& filename);
    // void load_particle_data2(const std::string& filename);
    

};

#endif
