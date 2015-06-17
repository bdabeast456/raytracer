all: raytracer.cpp ray_util.cpp lodepng.cpp
	  g++ -g -Wall -o raytracer raytracer.cpp ray_util.cpp lodepng.cpp

