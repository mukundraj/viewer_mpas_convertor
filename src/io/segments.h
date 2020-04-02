#ifndef SEGMENTS_H
#define SEGMENTS_H

#include <string>
#include <map>
#include <vector>



class segments
{


public:

    size_t nSteps, nSegments;
    std::vector<double> px, py, pz;
    std::vector<int> step;
    std::vector<int> segsizes;
    std::vector<int> seedid;
    std::vector<std::vector<double>> formatted_pxyz; // reformatted px,py,pz (merged segments)


    segments();

    void load_segments_data(const std::string& filename);
    void write_webvis_files(const std::string& opfile_traces);


};

#endif
