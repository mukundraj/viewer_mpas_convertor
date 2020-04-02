#include "particles.h"
#include <netcdf.h>
#include "def.h"
#include <iostream>
#include <fstream>
#include <string>
#include "segments.h"
#include <algorithm>
#include "misc.h"
#include <unordered_map>


using namespace std;

segments::segments(){

}

void segments::write_webvis_files(const std::string& opfile_traces){

    ofstream myfile (opfile_traces.c_str(), ios::out|ios::trunc);
    string sep = "";
    if (myfile.is_open())
    {
      for(size_t i = 0; i < formatted_pxyz.size(); i += 1){
        if (formatted_pxyz[i].size()>0){
            for (size_t j=0; j < formatted_pxyz[i].size(); j+=3){
              myfile <<sep<<formatted_pxyz[i][j] << " " << formatted_pxyz[i][j+1]<<" "<< formatted_pxyz[i][j+2];
              sep = " ";
            }
            myfile <<"\n";
        }


      }
      myfile.close();
    }

}
void segments::load_segments_data(const std::string& filename){

    int ncid;
    int dimid_x, dimid_y, dimid_z, dimid_segsize;


    int varid_segsizes, varid_px, varid_py, varid_pz, varid_seedid, varid_step;

    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "x", &dimid_x) );

    // NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "y", &dimid_y) );

    // NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "z", &dimid_z) );

    // NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "segsize", &dimid_segsize) );
    

    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_x, &nSteps) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_segsize, &nSegments) );

    NC_SAFE_CALL( nc_inq_varid(ncid, "px", &varid_px) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "py", &varid_py) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "pz", &varid_pz) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "segsizes", &varid_segsizes) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "seedid", &varid_seedid) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "step", &varid_step) );

    px.resize(nSteps);
    py.resize(nSteps);
    pz.resize(nSteps);
    segsizes.resize(nSegments);
    seedid.resize(nSegments);
    step.resize(nSegments);

    size_t start[1] = {0}, count[1] = {nSteps};
    NC_SAFE_CALL(nc_get_vara_double(ncid, varid_px, start, count, &px[0]));
    NC_SAFE_CALL(nc_get_vara_double(ncid, varid_py, start, count, &py[0]));
    NC_SAFE_CALL(nc_get_vara_double(ncid, varid_pz, start, count, &pz[0]));

    start[0] = 0; count[0] = nSegments;
    NC_SAFE_CALL(nc_get_vara_int(ncid, varid_segsizes, start, count, &segsizes[0]));
    NC_SAFE_CALL(nc_get_vara_int(ncid, varid_seedid, start, count, &seedid[0]));
    NC_SAFE_CALL(nc_get_vara_int(ncid, varid_step, start, count, &step[0]));

    for (size_t i =0 ;i<segsizes.size(); i++){
        fprintf(stderr, "seg %d\n", segsizes[i]);
    }


//    for (size_t i =0 ;i<px.size(); i++){
//        fprintf(stderr, "px %f\n", px[i]);
//    }


    NC_SAFE_CALL( nc_close(ncid));

    int numseeds = 1 + *std::max_element(seedid.begin(), seedid.end());

    formatted_pxyz.resize(numseeds);

    

    std::unordered_map<int, int> trace_sizes;
    for (size_t i=0; i<nSegments; i++){

        int cur_seed_id = seedid[i];
        if (trace_sizes.find(cur_seed_id)==trace_sizes.end())
            trace_sizes[cur_seed_id] = segsizes[i];
        else
             trace_sizes[cur_seed_id] += segsizes[i];

    }

    int idx = 0;
    for (size_t i=0; i<nSegments; i++){
        
        int cur_seg_size = segsizes[i];
        int cur_seed_id = seedid[i];
        int start_step = step[i];
        // dprint("ss %d %d", cur_seg_size, cur_seed_id);
        if (cur_seed_id>-1){
            
            formatted_pxyz[cur_seed_id].resize(3*trace_sizes[cur_seed_id]);
            
            for (int j=3*start_step; j<3*(start_step+cur_seg_size); j+=3){
                
                formatted_pxyz[cur_seed_id][j] = px[idx];
                formatted_pxyz[cur_seed_id][j+1] = py[idx];
                formatted_pxyz[cur_seed_id][j+2] = pz[idx];

                idx++;

                // formatted_pxyz[cur_seed_id].push_back(px[j]);
                // formatted_pxyz[cur_seed_id].push_back(py[j]);
                // formatted_pxyz[cur_seed_id].push_back(pz[j]);
                // dprint(" %f %f %f, %d, %d", px[j], py[j], pz[j], idx, cur_seg_size);
            }
            // std::reverse(formatted_pxyz[cur_seed_id].begin(), formatted_pxyz[cur_seed_id].end()); // to fix the disconneced segments in viewer
        }else{
            idx += 1;
        }
        // idx += cur_seg_size;

    }



    /*
    int ncid;

    int dimid_nStep, dimid_nParticle, dimid_nVertices, dimid_nVertexAllLayers;
    int varid_xPos, varid_yPos, varid_zPos, varid_tsg,
            varid_xVertex, varid_yVertex, varid_zVertex, varid_zTopVertex,
            varid_velocityVx, varid_velocityVy, varid_velocityVz; // trace_sizes_global



    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nStep", &dimid_nStep) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nParticle", &dimid_nParticle) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_nVertices) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertexAllLayers", &dimid_nVertexAllLayers) );

    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nParticle, &nParticles) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nStep, &nSteps) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertices, &nVertices) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertexAllLayers, &nVertexAllLayers) );

    fprintf(stderr, "particles.cpp load_particle_data %ld %ld\n", nParticles, nSteps);

    NC_SAFE_CALL( nc_inq_varid(ncid, "xPos", &varid_xPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "yPos", &varid_yPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zPos", &varid_zPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "trace_sizes_global", &varid_tsg) );

    NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );

    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVx", &varid_velocityVx) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVy", &varid_velocityVy) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVz", &varid_velocityVz) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zTopVertex", &varid_zTopVertex) );



    const size_t start_pos[2] = {0, 0};
    const size_t count_pos[2] = {nSteps, nParticles};

    xPos.resize(nSteps * nParticles);
    yPos.resize(nSteps * nParticles);
    zPos.resize(nSteps * nParticles);
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xPos, start_pos, count_pos, &xPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yPos, start_pos, count_pos, &yPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zPos, start_pos, count_pos, &zPos[0]) );

    const size_t start_tsg[1] = {0};
    const size_t count_tsg[1] = {nParticles};
    trace_sizes_global.resize(nParticles);

    NC_SAFE_CALL( nc_get_vara_int(ncid, varid_tsg, start_tsg, count_tsg, &trace_sizes_global[0]) );

    const size_t start_vertices[1] = {0};
    const size_t count_vertices[1] = {nVertices};
    xVertex.resize(nVertices);
    yVertex.resize(nVertices);
    zVertex.resize(nVertices);

    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, count_vertices, &xVertex[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, count_vertices, &yVertex[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, count_vertices, &zVertex[0]) );




    const size_t start_velocityV[1] = {0};
    const size_t count_velocityV[1] = {nVertexAllLayers};
    velocityVx.resize(nVertexAllLayers);
    velocityVy.resize(nVertexAllLayers);
    velocityVz.resize(nVertexAllLayers);
    zTopVertex.resize(nVertexAllLayers);

//    fprintf(stderr, "\n nvertexVelocity %ld %d", nVertexVelocity, varid_velocityVx);

    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVx, start_velocityV, count_velocityV, &velocityVx[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVy, start_velocityV, count_velocityV, &velocityVy[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVz, start_velocityV, count_velocityV, &velocityVz[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zTopVertex, start_velocityV, count_velocityV, &zTopVertex[0]) );

//    for (size_t i=0;i<xPos.size();i++){
//        fprintf(stderr, "%f ", yPos[i] );
//    }
//    fprintf(stderr, "trace_sizes_global ");
//    for (size_t i=0;i<nParticles;i++){
//        fprintf(stderr, "%d ", trace_sizes_global[i] );
//    }



    NC_SAFE_CALL( nc_close(ncid));

//    Eigen::MatrixXd M(3, nVertices);
//            for (int i=0;i<nVertices;i++){
//                M(0,i) = xVertex[i];
//                M(1,i) =yVertex[i];
//                M(2,i) = zVertex[i];
//            }

//    Nabo::NNSearchD* nns = Nabo::NNSearchD::createKDTreeLinearHeap(M);

//    // get nearest vertex neighbor ids using mpas_c
//        int K = 6;
//        Eigen::VectorXi nearest_idx(K);
//        Eigen::VectorXd dists2(K);

//        Eigen::VectorXd q(3);
//        // q<<X[0], X[1], X[2];
////        q << -5682112.802109 , -619878.622472 , 2814587.637284;
////        q << -5676505.755330 , -501327.815854 , 2849300.230092;
////        q << -5676505.755330 , -501327.815854 , 2849300.230092;
//          q << -2089982.455199 , 3032836.244409 , -5198695.665453;

//        nns->knn(q, nearest_idx, dists2, K, 0, Nabo::NNSearchF::SORT_RESULTS);

////        std::cout<<"nearest_idx\n"<<nearest_idx<<"\n\n";
////       int nearest_idx_array[K];
////            for (int i=0;i<K;i++){
////                nearest_idx_array[i] = nearest_idx[i];
////                fprintf(stderr, "nearest %d %f %f %f %f\n", nearest_idx_array[i], xVertex[i],
////                    yVertex[i], zVertex[i], dists2[i]);

////            }
//            this->nearest_idx = nearest_idx;

*/

}
