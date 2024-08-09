#ifndef profiler_H
#define profiler_H

#include <chrono>
#include <string>
#include <iostream>
using namespace std;

//!A class for simple timing within the codebase
/*!
Initialize the profiler with some string to be used in printing, then call it
with start() and end().  Asking the profiler for timing() will return the
average time taken between start and end blocks,and print() will print some
information to screen for convenience.  chrono is used to do the timing
*/
class profiler
    {
    public:
        profiler(){functionCalls = 0; timeTaken = 0;};
        profiler(string profilerName) : name(profilerName) {functionCalls = 0; timeTaken = 0;};

        void start()
            {
            startTime = chrono::high_resolution_clock::now();
            };
        void end()
            {
            endTime = chrono::high_resolution_clock::now();
            chrono::duration<double> difference = endTime-startTime;
            timeTaken += difference.count();
            functionCalls +=1;
            };

        double timing()
            {
            if(functionCalls>0)
                return timeTaken/functionCalls;
            else
                return 0;
            };

        void print()
            {
            cout << "profiler \"" << name << "\" took an average of " << timing() << " per call over " << functionCalls << " calls...total time = "<<timing()*functionCalls << endl;
            }

        void setName(string _name){name=_name;};
        chrono::time_point<chrono::high_resolution_clock>  startTime;
        chrono::time_point<chrono::high_resolution_clock>  endTime;
        int functionCalls;
        double timeTaken;
        string name;
    };
#endif
