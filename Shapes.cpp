#include <iostream>
#include <vector>
#include <math.h>

class CPoint
{
public:
    explicit CPoint(double x = 0, double y = 0)
        : x_(x), y_(y) {}
    CPoint(const CPoint &other)
        : x_(other.x_), y_(other.y_) {}
    CPoint &operator=(const CPoint &other)
    {
        if (&other == this)
            return *this;
        x_ = other.x_;
        y_ = other.y_;
    }
    ~CPoint() {}
    double X() const
    {
        return x_;
    }
    double Y() const
    {
        return y_;
    }
    double radius() const
    {
        return std::sqrt(x_ * x_ + y_ * y_);
    }

private:
    double x_;
    double y_;
};

class CPolyline
{
public:
    explicit CPolyline(unsigned number_points, std::vector<CPoint> points)
        : number_points_(number_points), points_(points) {}
    explicit CPolyline(std::vector<CPoint> points)
        : number_points_(points.size()), points_(points) {}
    CPolyline(const CPolyline &other)
        : number_points_(other.number_points_), points_(other.points_) {}
    CPolyline &operator=(const CPolyline &other)
    {
        if (&other == this)
            return *this;
        number_points_ = other.number_points_;
        points_ = other.points_;
    }
    ~CPolyline() {}
    unsigned Number_points()
    {
        return number_points_;
    }

protected:
    unsigned number_points_;
    std::vector<CPoint> points_;
};

int main()
{
    CPoint p1(1, 2);
    std::cout << p1.radius() << std::endl;
}