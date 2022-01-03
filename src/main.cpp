// Create simple 3d highway environment using PCL

#include "highway.h"
//#include <fstream>

int main(int argc, char* argv[])
{
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);

	// set camera position and angle
	viewer->initCameraParameters();
	float x_pos = 0;
	viewer->setCameraPosition ( x_pos-26, 0, 15.0, x_pos+25, 0, 0, 0, 0, 1);

	Highway highway(viewer);

	//initHighway(viewer);

	int frame_per_sec = 30;
	int sec_interval = 10;
	int frame_count = 0;
	int time_us = 0;

	double egoVelocity = 25;

	while (frame_count < (frame_per_sec*sec_interval))
	{
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

		//stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		highway.stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		viewer->spinOnce(1000/frame_per_sec);
		frame_count++;
		time_us = 1000000*frame_count/frame_per_sec;
		
	}
	/*
    std::vector<double> r_nis = std::move(highway.traffic[0].ukf.radar_nis);
    std::vector<double> l_nis = std::move(highway.traffic[0].ukf.lidar_nis);

    std::ofstream fs;
    fs.open("lidar_nis.csv");
    for (int i=0; i<l_nis.size(); ++i)
    {
        fs << l_nis[i] << "\n";
    }
    fs.close();

    fs.open("radar_nis.csv");
    for (int i=0; i<r_nis.size(); ++i)
    {
        fs << r_nis[i] << "\n";
    }
    fs.close();
    */
    return 0;
}
