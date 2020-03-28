#ifndef MPASO_H
#define MPASO_H

#include <string>
#include <map>
#include <vector>

using namespace std;

class mpaso
{


public:

    size_t nCells, nEdges, nVertices, nVertLevels;
    std::vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
    std::vector<double> xyzCell; // xCell, yCell, zCell;
    std::vector<int> indexToVertexID, indexToCellID;
    std::vector<int> verticesOnEdge, cellsOnVertex;
    std::vector<double> velocityX, velocityY, velocityZ;
    std::vector<double> velocityXv, velocityYv, velocityZv;
    std::vector<double> zTop;
    std::vector<double> zTopVertex, zTopVertexNorm;
    std::map<int, int> vertexIndex, cellIndex;

    std::vector<double> xParticle, yParticle, zParticle; // for LIGHT
    size_t Time, nParticles; // for LIGHT

    std::vector<int> graph_info_part0, graph_info_part1, graph_info_part2,
                    graph_info_part3;

    std::vector<int> partitions;

    std::vector<double> uVertexVelocity, vVertexVelocity, wVertexVelocity;
    std::vector<double> vertVelocityTop;

    mpaso();

    void loadMeshFromNetCDF_CANGA(const std::string& filename, size_t time_id=0);
    void write_webvis_files(string &opfile_centers, string &opfile_vertices, string &opfile_edge_ids, string &opfile_vertvel);
    void load_partition_data(const std::string& filename);


    // void read_graph_info_part(const std::string& filename);
    // void load_light_particle_trace(const std::string& filename);
};

#endif
