#include <string>
#include <iostream>
#include <vector>
#include "io/mpaso.h"
#include "io/particles.h"
#include "io/segments.h"

using namespace std;

int main(int argc, char **argv){

	string ipfile = "/Users/mukundraj/Desktop/work/datasets/mpas/MPAS-O_V6.0_EC60to30/output.nc";
	string opfile_centers = "/Users/mukundraj/Desktop/work/projects/viewer_mpas_convertor/output/centers.csv";
	string opfile_vertices = "/Users/mukundraj/Desktop/work/projects/viewer_mpas_convertor/output/vertices.csv";
	string opfile_edge_ids = "/Users/mukundraj/Desktop/work/projects/viewer_mpas_convertor/output/edgeids.csv";
	string opfile_vertvel = "/Users/mukundraj/Desktop/work/projects/viewer_mpas_convertor/output/vertvel";

	string ip_particles = "/Users/mukundraj/Desktop/work/datasets/mpas/MPAS-O_V6.0_EC60to30/particles.nc";
	string ip_partinfo = "/Users/mukundraj/Desktop/work/datasets/mpas/MPAS-O_V6.0_EC60to30/graph.info.part.4";
	

	mpaso mpas1;

	particles par1;
	par1.load_particle_data(ip_particles);

	mpas1.loadMeshFromNetCDF_CANGA(ipfile, 0);
	mpas1.load_partition_data(ip_partinfo);
	cout<<mpas1.xyzCell.size()<<" "<<mpas1.nCells<<endl;



	mpas1.write_webvis_files(opfile_centers, opfile_vertices, opfile_edge_ids, opfile_vertvel);

	string ip_segments = "/Users/mukundraj/Desktop/work/datasets/mpas/MPAS-O_V6.0_EC60to30/segments.nc";
	segments s1;
	s1.load_segments_data(ip_segments);


	string op_traces = "/Users/mukundraj/Desktop/work/projects/viewer_mpas_convertor/output/traces.csv";
	s1.write_webvis_files(op_traces);



	return 0;
}