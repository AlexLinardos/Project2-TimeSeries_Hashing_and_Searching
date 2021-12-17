#include <math.h>
#include <vector>
#include <algorithm>
#include "../utils.hpp"
#include "../curves.hpp"

#include "./Fred/frechet.hpp"

namespace cF
{
    // Curve vecToCurve(Item &item)
    // {
    //     Points points(2);
    //     for (int i = 0; i < item.xij.size(); i++)
    //     {
    //         // Coordinates coordinates;
    //         // coordinates.push_back((coordinate_t)item.xij[i]);
    //         // coordinates.push_back((coordinate_t)i);
    //         Point point(2);
    //         point.set(0, (coordinate_t)item.xij[i]);
    //         point.set(1, (coordinate_t)i);
    //         points.add(point);
    //     }

    //     // Curve * curve = new Curve(points, item.id);
    //     Curve curve(points, item.id);
    //     return curve;
    // }

    Curve fredCurve(curves::Curve2d &curve)
    {
        Points points(1);
        for (int i = 0; i < curve.data.size(); i++)
        {
            Point point(1);
            point.set(0, (coordinate_t)curve.data[i].x);

            Point point2(1);
            point2.set(0, (coordinate_t)curve.data[i].y);

            points.add(point);
            points.add(point2);
        }

        // Curve * fred_curve = new Curve(points, item.id);
        Curve fred_curve(points, curve.id);
        return fred_curve;
    }

    // double distance(Item &item1,  Item &item2)
    // {
    //     return Frechet::Continuous::distance(vecToCurve(item1), vecToCurve(item2)).value;
    // }
    
    double c_distance(curves::Curve2d &item1,  curves::Curve2d &item2)
    {
        return Frechet::Continuous::distance(fredCurve(item1) , fredCurve(item2)).value;
    }
}