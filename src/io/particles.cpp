#include "particles.h"
#include <netcdf.h>
#include "def.h"
#include <iostream>
#include <fstream>
#include <string>

particles::particles(){

}



void particles::load_particle_data(const std::string& filename){


    int ncid;

    int dimid_nStep, dimid_nParticle, dimid_nVertices, dimid_nVertexAllLayers;
    int varid_xPos, varid_yPos, varid_zPos, varid_tsg,
            varid_xVertex, varid_yVertex, varid_zVertex, varid_zTopVertex,
            varid_velocityVx, varid_velocityVy, varid_velocityVz; // trace_sizes_global
    int varid_zLevelParticle;



    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "Time", &dimid_nStep) );
    NC_SAFE_CALL( nc_inq_dimid(ncid, "nParticles", &dimid_nParticle) );
//    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_nVertices) );
//    NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertexAllLayers", &dimid_nVertexAllLayers) );

    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nParticle, &nParticles) );
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nStep, &nSteps) );
//    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertices, &nVertices) );
//    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_nVertexAllLayers, &nVertexAllLayers) );

    fprintf(stderr, "particles.cpp load_particle_data %ld %ld\n", nParticles, nSteps);

    NC_SAFE_CALL( nc_inq_varid(ncid, "xParticle", &varid_xPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "yParticle", &varid_yPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zParticle", &varid_zPos) );
    NC_SAFE_CALL( nc_inq_varid(ncid, "zLevelParticle", &varid_zLevelParticle) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "trace_sizes_global", &varid_tsg) );

//    NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );

//    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVx", &varid_velocityVx) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVy", &varid_velocityVy) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "velocityVz", &varid_velocityVz) );
//    NC_SAFE_CALL( nc_inq_varid(ncid, "zTopVertex", &varid_zTopVertex) );



    const size_t start_pos[2] = {0, 0};
    const size_t count_pos[2] = {nSteps, nParticles};

    xPos.resize(nSteps * nParticles);
    yPos.resize(nSteps * nParticles);
    zPos.resize(nSteps * nParticles);
    zLevelParticle.resize(nSteps * nParticles);
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xPos, start_pos, count_pos, &xPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yPos, start_pos, count_pos, &yPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zPos, start_pos, count_pos, &zPos[0]) );
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zLevelParticle, start_pos, count_pos, &zLevelParticle[0]) );

    int val=100;
    fprintf(stderr, "\nlen %ld (%.15f %.15f %.15f) (%.15f %.15f %.15f)\n", xPos.size(), xPos[0+val], yPos[0+val], zPos[0+val], xPos[235160*5+val], yPos[235160*5+val], zPos[235160*5+val]);

    // for (size_t i=0; i<235160*10; i+=235160)
    //         fprintf(stderr, "zLev %f\n", zLevelParticle[i]);

//    const size_t start_tsg[1] = {0};
//    const size_t count_tsg[1] = {nParticles};
//    trace_sizes_global.resize(nParticles);

//    NC_SAFE_CALL( nc_get_vara_int(ncid, varid_tsg, start_tsg, count_tsg, &trace_sizes_global[0]) );

//    const size_t start_vertices[1] = {0};
//    const size_t count_vertices[1] = {nVertices};
//    xVertex.resize(nVertices);
//    yVertex.resize(nVertices);
//    zVertex.resize(nVertices);

//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, count_vertices, &xVertex[0]) );
//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, count_vertices, &yVertex[0]) );
//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, count_vertices, &zVertex[0]) );




//    const size_t start_velocityV[1] = {0};
//    const size_t count_velocityV[1] = {nVertexAllLayers};
//    velocityVx.resize(nVertexAllLayers);
//    velocityVy.resize(nVertexAllLayers);
//    velocityVz.resize(nVertexAllLayers);
//    zTopVertex.resize(nVertexAllLayers);

////    fprintf(stderr, "\n nvertexVelocity %ld %d", nVertexVelocity, varid_velocityVx);

//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVx, start_velocityV, count_velocityV, &velocityVx[0]) );
//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVy, start_velocityV, count_velocityV, &velocityVy[0]) );
//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityVz, start_velocityV, count_velocityV, &velocityVz[0]) );
//    NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zTopVertex, start_velocityV, count_velocityV, &zTopVertex[0]) );

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

}
