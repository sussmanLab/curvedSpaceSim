#ifndef dataType_H
#define dataType_H

/*! \file dataType.h
defines other custom simple data types used in the code 
*/

class int3
    {
    public:
        int x;
        int y;
        int z;

        int3& operator=(const int3& a)
            {
            this->x = a.x;
            this->y = a.y;
            this->z = a.z;
            return *this;
            };
    };

class double3
    {
    public:
        double x;
        double y;
        double z;

        double3& operator=(const double3& a)
            {
            this->x = a.x;
            this->y = a.y;
            this->z = a.z;
            return *this;
            };
    };

class double4
    {
    public:
        double x;
        double y;
        double z;
        double w;

        double4& operator=(const double4& a)
            {
            this->x = a.x;
            this->y = a.y;
            this->z = a.z;
            this->w = a.w;
            return *this;
            };
    };

#endif
