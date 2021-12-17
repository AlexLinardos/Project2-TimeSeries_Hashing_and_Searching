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

    // minima maxima y removing from π˜ = v1, . . . , v` any vi s.t. vi ∈ [min{vi−1, vi+1}, max{vi−1, vi+1}] 
    // στο min/max αφαιρεις ενα στοιχειο αν ειναι μεγαλυτερο απο το προηγουμενο του και μικροτερο απο το επομενο του ή αν ειναι ισο με ενα απο τα 2
    void minima_maxima(vector<double> &p)
    {
        int og_size = p.size();
        int i=0;
        int cuts=0;
        int total_cuts=0;
        while(i+2 < p.size())
        {
            cuts=0;
            if ((p[i+1] >= min(p[i], p[i+2])) || (p[i+1] >= max(p[i], p[i+2])))
            {
                p.erase(p.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
            if(cuts==0)
                i++;
        }
        while(cuts!=0)
        {
            cuts=0;
            i--;
            if ((p[i+1] >= min(p[i], p[i+2])) || (p[i+1] >= max(p[i], p[i+2])))
            {
                p.erase(p.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "minima_maxima cut " << total_cuts << " points out of "<< og_size << endl;
        return;
    }

    
    curves::Curve2d * filter_curve(curves::Curve2d &curve, double e) // for any consecutive points a, b, c, if |a − b| ≤ ε and |b − c| ≤ ε then remove b
    {
        curves::Curve2d * filtered_curve = new curves::Curve2d(curve.id);
        filtered_curve->data=curve.data;
        int i=0;
        int cuts=0;
        int total_cuts=0;
        while(i+2 < filtered_curve->data.size())
        {
            cuts=0;
            if((point2d_L2(filtered_curve->data[i], filtered_curve->data[i+1]) <= e) && (point2d_L2(filtered_curve->data[i+1], filtered_curve->data[i+2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
            if(cuts==0)
                i++;
        }
        while(cuts!=0)
        {
            cuts=0;
            i--;
            if((point2d_L2(filtered_curve->data[i], filtered_curve->data[i+1]) <= e) && (point2d_L2(filtered_curve->data[i+1], filtered_curve->data[i+2]) <= e))
            {
                (filtered_curve->data).erase(filtered_curve->data.begin()+i+1); // erase element i+1
                cuts++;
                total_cuts++;
            }
        }

        cout << "filtering cut " << total_cuts << " points out of "<< curve.data.size() << endl;
        return filtered_curve;
    }

    vector<curves::Curve2d> * filter_curves(vector<curves::Curve2d> &curve_dataset, double e)
    {
        vector<curves::Curve2d> * filtered_curve_dataset = new vector<curves::Curve2d>;
        for(int i=0; i < curve_dataset.size(); i++)
        {
            cout << "curve#" << i << ": ";
            filtered_curve_dataset->push_back(*filter_curve(curve_dataset[i], e));
        }
        return filtered_curve_dataset;
    }

    // Item * to1dcurves(vector<curves::Curve2d> &curve_dataset)
    // {
    //     Item * curves = new Iteam 
    // }


}