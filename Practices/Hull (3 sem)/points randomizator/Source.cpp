#ifndef PointPairFile
#define PointPairFile


#include "Source1.cpp"

struct PointPair {
public:
		Point first;
		Point second;

	PointPair(const Point firstPoint, const Point secondPoint) {
		first = firstPoint;
		second = secondPoint;
	}
};

#endif