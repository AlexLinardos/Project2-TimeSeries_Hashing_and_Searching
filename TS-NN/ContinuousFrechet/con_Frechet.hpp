#include <math.h>
#include <vector>
#include <algorithm>
#include "../utils.hpp"

#include "./Fred/frechet.hpp"

namespace cF
{
    Curve vecToCurve(Item &item)
    {
        Points points(2);
        for (int i = 0; i < item.xij.size(); i++)
        {
            // Coordinates coordinates;
            // coordinates.push_back((coordinate_t)item.xij[i]);
            // coordinates.push_back((coordinate_t)i);
            Point point(2);
            point.set(0, (coordinate_t)item.xij[i]);
            point.set(1, (coordinate_t)i);
            points.add(point);
        }

        // Curve * curve = new Curve(points, item.id);
        Curve curve(points, item.id);
        return curve;
    }

    Curve vecToCurve(curves::Curve2d &curve)
    {
        Points points(2);
        for (int i = 0; i < curve.data.size(); i++)
        {
            // Coordinates coordinates;
            // coordinates.push_back((coordinate_t)item.xij[i]);
            // coordinates.push_back((coordinate_t)i);
            Point point(2);
            point.set(0, (coordinate_t)curve.data[i].x);
            point.set(1, (coordinate_t)curve.data[i].y);
            points.add(point);
        }

        // Curve * fred_curve = new Curve(points, item.id);
        Curve fred_curve(points, curve.id);
        return fred_curve;
    }

    double distance(Item &item1,  Item &item2)
    {
        return Frechet::Continuous::distance(vecToCurve(item1) , vecToCurve(item2)).value;
    }
    
    double c_distance(curves::Curve2d &item1,  curves::Curve2d &item2)
    {
        return Frechet::Continuous::distance(vecToCurve(item1) , vecToCurve(item2)).value;
    }
    
    curves::Curve2d * filter_curve(curves::Curve2d &curve, double e) // for any consecutive points a, b, c, if |a − b| ≤ ε and |b − c| ≤ ε then remove b
    {
        curves::Curve2d * filtered_curve = new curves::Curve2d(curve.id);
        filtered_curve->data
        int i=0;
        while(1)
        {

            i++;
        }
        return 
    }


}