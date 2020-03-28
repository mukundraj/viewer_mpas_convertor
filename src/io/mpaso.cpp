#include "mpaso.h"
#include <netcdf.h>
#include "def.h"
// #include "interpolators.h"
#include <iostream>
#include <fstream>
#include <string>
#include "misc.h"


mpaso::mpaso(){
}

void mpaso::write_webvis_files(string &opfile_centers, string &opfile_vertices, string &opfile_edge_ids, string &opfile_vertvel){

    ofstream myfile (opfile_centers.c_str(), ios::out|ios::trunc);
    if (myfile.is_open())
    {
      for(size_t i = 0; i < xyzCell.size(); i += 3){
          myfile << xyzCell[i] << " " << xyzCell[i+1]<<" "<< xyzCell[i+2]<<" "<<partitions[i/3]<<"\n";
      }
      myfile.close();
    }


    myfile.open(opfile_vertices.c_str(), ios::out|ios::trunc);
    if (myfile.is_open()){

      for (size_t i=0; i<xVertex.size(); i+=1){
          myfile << xVertex[i]<<" "<<yVertex[i]<<" "<<zVertex[i]<<"\n";
      }

      myfile.close();
    }


    myfile.open(opfile_edge_ids.c_str(), ios::out|ios::trunc);
    if (myfile.is_open()){

      for (size_t i=0; i<verticesOnEdge.size(); i+=2){
          myfile << verticesOnEdge[i]<<" "<<verticesOnEdge[i+1]<<"\n";
      }

      myfile.close();
    }


    // write out vertex velocity files

    // for (int i=0; i<nVertLevels; i++){
    //   string filename = opfile_vertvel + to_string(i) + ".csv";

    //   myfile.open(filename.c_str(), ios::out|ios::trunc);
    //   if (myfile.is_open()){
    //     for (int j=0; j<nVertices; j++){
    //       myfile << uVertexVelocity[j*nVertLevels + i]<< " "<<vVertexVelocity[j*nVertLevels + i]<<" "<<wVertexVelocity[j*nVertLevels + i]<<"\n";

    //     }
    //     myfile.close();
    //   }

    // }
    

}

void mpaso::load_partition_data(const std::string& filename){

  std::vector<std::vector<int>> in_data = read_csv(filename);

  for (size_t i=0; i<in_data.size(); i++)
    partitions.push_back(in_data[i][0]);

  dprint("psize %ld", in_data.size());

}

void mpaso::loadMeshFromNetCDF_CANGA(const std::string& filename, size_t time_id){


  int ncid;
  int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels;
  int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
      varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell,
      varid_verticesOnEdge, varid_cellsOnVertex,
      varid_indexToVertexID, varid_indexToCellID,
      varid_velocityX, varid_velocityY, varid_velocityZ, 
      varid_uVertexVelocity, varid_vVertexVelocity, varid_wVertexVelocity;
  int varid_zTop;




  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  NC_SAFE_CALL( nc_inq_dimid(ncid, "nCells", &dimid_cells) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nEdges", &dimid_edges) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertices", &dimid_vertices) );
  NC_SAFE_CALL( nc_inq_dimid(ncid, "nVertLevels", &dimid_vertLevels) );
   NC_SAFE_CALL( nc_inq_dimid(ncid, "nEdges", &dimid_edges) );


  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_cells, &nCells) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_edges, &nEdges) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertices, &nVertices) );
  NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_vertLevels, &nVertLevels) );

  NC_SAFE_CALL( nc_inq_varid(ncid, "indexToVertexID", &varid_indexToVertexID) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "indexToCellID", &varid_indexToCellID) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "latCell", &varid_latCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonCell", &varid_lonCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xCell", &varid_xCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yCell", &varid_yCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zCell", &varid_zCell) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "latVertex", &varid_latVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "lonVertex", &varid_lonVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "xVertex", &varid_xVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "yVertex", &varid_yVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zVertex", &varid_zVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "cellsOnVertex", &varid_cellsOnVertex) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityX", &varid_velocityX) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityY", &varid_velocityY) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "velocityZ", &varid_velocityZ) );
  NC_SAFE_CALL( nc_inq_varid(ncid, "zTop", &varid_zTop));

  NC_SAFE_CALL( nc_inq_varid(ncid, "uVertexVelocity", &varid_uVertexVelocity));
  NC_SAFE_CALL( nc_inq_varid(ncid, "vVertexVelocity", &varid_vVertexVelocity));
  NC_SAFE_CALL( nc_inq_varid(ncid, "wVertexVelocity", &varid_wVertexVelocity));



  const size_t start_cells[1] = {0}, size_cells[1] = {nCells};

  indexToCellID.resize(nCells);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToCellID, start_cells, size_cells, &indexToCellID[0]) );
  for (int i=0; i<nCells; i++) { // shoudn't i vary from 1 to nCells here? then cellIndex[indexToCellID[i-1]] = i-1
    cellIndex[indexToCellID[i]] = i;
    // fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
    //std::cout<<i<<" "<<indexToCellID[i]<<" "<<cellIndex[i]<<"\n";
  }

  std::vector<double> coord_cells;
  coord_cells.resize(nCells);
  xyzCell.resize(nCells*3);
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3] = coord_cells[i];
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3+1] = coord_cells[i];
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &coord_cells[0]) );
  for (int i=0; i<nCells; i++) xyzCell[i*3+2] = coord_cells[i];

  const size_t start_vertices[1] = {0}, size_vertices[1] = {nVertices};
  latVertex.resize(nVertices);
  lonVertex.resize(nVertices);
  xVertex.resize(nVertices);
  yVertex.resize(nVertices);
  zVertex.resize(nVertices);
  indexToVertexID.resize(nVertices);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]) );


  const size_t start_time_vertex_level[3] = {time_id, 0, 0}, size_time_vertex_level[3] = {1, nVertices, nVertLevels};

  uVertexVelocity.resize(nVertices*nVertLevels);
  vVertexVelocity.resize(nVertices*nVertLevels);
  wVertexVelocity.resize(nVertices*nVertLevels);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_uVertexVelocity, start_time_vertex_level, size_time_vertex_level, &uVertexVelocity[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_vVertexVelocity, start_time_vertex_level, size_time_vertex_level, &vVertexVelocity[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_wVertexVelocity, start_time_vertex_level, size_time_vertex_level, &wVertexVelocity[0]) );



  for (int i=0; i<nVertices; i++) {
    vertexIndex[indexToVertexID[i]] = i;
    // fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
//    if (i<20){
//        fprintf(stderr, "\n%d %f %f %f", i, xVertex[i], yVertex[i], zVertex[i]);
//    }

  }

  const size_t start_edges2[2] = {0, 0}, size_edges2[2] = {nEdges, 2};
  verticesOnEdge.resize(nEdges*2);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]) );

  // for (int i=0; i<nEdges; i++)
  //   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);

  const size_t start_vertex_cell[2] = {0, 0}, size_vertex_cell[2] = {nVertices, 3};
  cellsOnVertex.resize(nVertices*3);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid_cellsOnVertex, start_vertex_cell, size_vertex_cell, &cellsOnVertex[0]) );

  // getting velocity at one time at cell centers
  const size_t start_time_cell_level[3] = {time_id, 0, 0}, size_time_cell_level[3] = {1, nCells, nVertLevels};
  velocityX.resize(nCells*nVertLevels);
  velocityY.resize(nCells*nVertLevels);
  velocityZ.resize(nCells*nVertLevels);



  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityX, start_time_cell_level, size_time_cell_level, &velocityX[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityY, start_time_cell_level, size_time_cell_level, &velocityY[0]) );
  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_velocityZ, start_time_cell_level, size_time_cell_level, &velocityZ[0]) );

  zTop.resize(nCells*nVertLevels);
  zTopVertex.resize(nVertices*nVertLevels);
  zTopVertexNorm.resize(nVertices*nVertLevels);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid_zTop, start_time_cell_level, size_time_cell_level, &zTop[0]) );

  // derive velocity on cell verticies (on only one altitude level here)
  velocityXv.resize(nVertices*nVertLevels);
  velocityYv.resize(nVertices*nVertLevels);
  velocityZv.resize(nVertices*nVertLevels);




  NC_SAFE_CALL( nc_close(ncid));



}

