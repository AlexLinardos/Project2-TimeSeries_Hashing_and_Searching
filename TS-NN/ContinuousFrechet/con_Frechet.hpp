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

    double distance(Item &item1, Item &item2)
    {
        return Frechet::Continuous::distance(vecToCurve(item1), vecToCurve(item2)).value;
    }
}