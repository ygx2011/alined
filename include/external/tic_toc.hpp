/*
* Timer tool inspired by MATLAB's tic-toc function.
*
* Copyright (C) 2016 Andrea Luca Lampart <lamparta at student dot ethz dot ch> (ETH Zurich)
* For more information see <https://github.com/andrealampart/elised_plus_plus>
*/

#pragma once

#include <iostream>
#include <chrono>


//******************************************************************************************
//	TIMER:
//	Use this function as an equivalent to the matlab function:
//	Place tic() at the beginning of the code
//	Place toc() at the end of the code
//	Time in between tic() and toc() is measured in milli- or nanoseconds
//  @param option       "nsec" for nanoseconds, "msec" for milliseconds
//  @param annotation   "Annotation:" custom annotation
//  @param color        "red, green, yellow, blue" are supported colors
//******************************************************************************************

//set start time
static std::chrono::high_resolution_clock::time_point start_time;
static std::string TIMER_ANNOTATION;
static std::string TIMER_COLOR;
static bool TIMER_OPTION;

inline void tic(std::string option = "msec" ,std::string annotation = "", std::string color = "green"){
  TIMER_ANNOTATION = annotation;
  if(color == "red"){
    TIMER_COLOR = "\033[1;31m";
  }
  else if(color == "green"){
    TIMER_COLOR = "\033[1;32m";
  }
  else if(color == "yellow"){
    TIMER_COLOR = "\033[1;33m";
  }
  else if(color == "blue"){
    TIMER_COLOR = "\033[1;34m";
  }
  if(option == "msec"){
    TIMER_OPTION = true;
  }
  else if(option == "nsec"){
    TIMER_OPTION = false;
  }

	start_time = std::chrono::high_resolution_clock::now();

}



inline uint64_t toc(std::string muted = ""){
	std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now()-start_time;
	//Use colored console output
  if(TIMER_OPTION)
  {
    if(muted != "muted")
    std::cout << TIMER_COLOR << TIMER_ANNOTATION << " Elapsed time is " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() <<" milliseconds. \033[0m\n" << std::endl;

    return (uint64_t)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
  }
  else{
    if(muted != "muted")
    std::cout << TIMER_COLOR << TIMER_ANNOTATION << " Elapsed time is " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() <<" nanoseconds. \033[0m\n" << std::endl;

    return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
  }
}


