#include <iostream>
#include <cmath>

#define PI 3.14159265

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
        {
            return *this;
        }
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
    double Radius() const
    {
        return std::sqrt(x_ * x_ + y_ * y_);
    }

private:
    double x_;
    double y_;
};

double Distance(CPoint p1, CPoint p2)
{
    return std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
}

class CPolyline
{
public:
    explicit CPolyline()
        : number_points_(0), points_(new CPoint[number_points_]) {}
    explicit CPolyline(size_t number_points, CPoint *points)
        : number_points_(number_points), points_(new CPoint[number_points])
    {
        memcpy(points_, points, sizeof(CPoint) * number_points_);
    }
    CPolyline(const CPolyline &other)
        : number_points_(other.number_points_), points_(new CPoint[other.number_points_])
    {
        memcpy(points_, other.points_, sizeof(CPoint) * number_points_);
    }
    CPolyline &operator=(const CPolyline &other)
    {
        if (&other == this)
        {
            return *this;
        }
        number_points_ = other.number_points_;
        delete[] points_;
        points_ = new CPoint[number_points_];
        memcpy(points_, other.points_, sizeof(CPoint) * number_points_);
        return *this;
    }
    virtual ~CPolyline()
    {
        delete[] points_;
    }
    virtual size_t NumberPoints() const
    {
        return number_points_;
    }
    virtual CPoint Point(size_t i) const
    {
        return points_[i];
    }
    virtual double Length() const
    {
        double length = 0;
        for (size_t i = 1; i < number_points_; ++i)
        {
            length += Distance(points_[i - 1], points_[i]);
        }
        return length;
    }

private:
    size_t number_points_;
    CPoint *points_;
};

class CClosedPolyline : public CPolyline
{
public:
    explicit CClosedPolyline() : CPolyline() {}
    explicit CClosedPolyline(size_t number_points, CPoint *points)
        : CPolyline(number_points, points) {}
    CClosedPolyline(const CClosedPolyline &other)
        : CPolyline(other) {}
    CClosedPolyline &operator=(const CClosedPolyline &other)
    {
        CPolyline::operator=(other);
        return *this;
    }
    ~CClosedPolyline() {}
    double Length() const
    {
        size_t number_points = this->NumberPoints();
        if (number_points < 2)
            return 0;
        double length_ = 0;
        for (size_t i = 1; i < number_points; ++i)
        {
            length_ += Distance(this->Point(i), this->Point(i - 1));
        }
        length_ += Distance(this->Point(0), this->Point(number_points - 1));
        return length_;
    }
};

class CPolygon
{
public:
    explicit CPolygon() : outline_() {}
    explicit CPolygon(size_t number_points, CPoint *points)
        : outline_(number_points, points) {} // checking for self-intersections
    CPolygon(const CPolygon &other)
        : outline_(other.outline_) {}
    CPolygon &operator=(const CPolygon &other)
    {
        if (this == &other)
        {
            return *this;
        }
        outline_ = other.outline_;
        return *this;
    }
    ~CPolygon() {}
    virtual double Perimeter() const
    {
        return outline_.Length();
    }
    virtual double Square() const
    {
        double square = 0;
        for (size_t i = 0; i < outline_.NumberPoints() - 1; ++i)
        {
            CPoint p = outline_.Point(i);
            CPoint p_next = outline_.Point(i + 1);
            square += (p.X() + p_next.X()) * (p.Y() - p_next.Y());
        }
        CPoint p_end = outline_.Point(outline_.NumberPoints() - 1);
        CPoint p_start = outline_.Point(0);
        square += (p_end.X() + p_start.X()) * (p_end.Y() - p_start.Y());
        square = std::fabs(square * 0.5);
        return square;
    }
    size_t NumberPoints() const
    {
        return outline_.NumberPoints();
    }

private:
    CClosedPolyline outline_;
};

class CRegularPolygon : public CPolygon
{
    explicit CRegularPolygon() : CPolygon(), side_(0) {}
    explicit CRegularPolygon(size_t number_points, CPoint *points)
        : CPolygon(number_points, points), side_(0)
    {
        if (number_points != 0)
        {
            bool is_regular = true;
            side_ = Distance(points[number_points - 1], points[0]);
            for (size_t i = 1; i < number_points; ++i)
            {
                double dist = Distance(points[i - 1], points[i]);
                if (dist != side_)
                {
                    is_regular = false;
                    break;
                }
            }
            if (!is_regular)
                std::cerr << "This polygon is not regular!";
        }
    }
    CRegularPolygon(const CRegularPolygon &other)
        : CPolygon(other), side_(other.side_) {}
    CRegularPolygon &operator=(const CRegularPolygon &other)
    {
        if (&other == this)
        {
            return *this;
        }
        CPolygon::operator=(other);
        side_ = other.side_;
        return *this;
    }
    ~CRegularPolygon() {}
    double Perimeter() const
    {
        double perimeter = this->NumberPoints() * side_;
        return perimeter;
    }
    double Square() const
    {
        double square = (this->NumberPoints() * side_ * side_) / (4 * std::tan(PI / this->NumberPoints()));
        return square;
    }

private:
    double side_;
};

class CTriangle : public CPolygon
{
    explicit CTriangle() : CPolygon() {}
    explicit CTriangle(size_t number_points, CPoint *points)
        : CPolygon(number_points, points)
    {
        if (number_points != 3)
        {
            std::cerr << "It is not a triangle!";
        }
    }
    explicit CTriangle(CPoint point1, CPoint point2, CPoint point3)
    {
        CPoint *points = new CPoint[3];
        points[0] = point1;
        points[1] = point2;
        points[2] = point3;
        CPolygon(3, points);
    }
    CTriangle(const CTriangle &other)
        : CPolygon(other) {}
    CTriangle &operator=(const CTriangle &other)
    {
        CPolygon::operator=(other);
        return *this;
    }
    ~CTriangle() {}
};

class CTrapezoid : public CPolygon
{
    explicit CTrapezoid() : CPolygon() {}
    explicit CTrapezoid(size_t number_points, CPoint *points)
        : CPolygon(number_points, points)
    {
        bool is_parallel = false, not_parallel = false;
        double k1, k2, k3, k4;
        k1 = (points[0].Y() - points[1].Y()) / (points[0].X() - points[1].X());
        k2 = (points[1].Y() - points[2].Y()) / (points[1].X() - points[2].X());
        k3 = (points[2].Y() - points[3].Y()) / (points[2].X() - points[3].X());
        k4 = (points[0].Y() - points[3].Y()) / (points[0].X() - points[3].X());
        if (std::fabs(k1 - k3) < 0.000001)
        {
            is_parallel = true;
        }
        else
        {
            not_parallel = true;
        }
        if (std::fabs(k2 - k4) < 0.000001)
        {
            is_parallel = true;
        }
        else
        {
            not_parallel = true;
        }
        if (number_points != 4 || !is_parallel || !not_parallel)
        {
            std::cerr << "It is not a triangle!";
        }
    }
    explicit CTrapezoid(CPoint point1, CPoint point2, CPoint point3, CPoint point4)
    {
        CPoint *points = new CPoint[4];
        points[0] = point1;
        points[1] = point2;
        points[2] = point3;
        points[3] = point4;
        CTrapezoid(4, points);
    }
    CTrapezoid(const CTrapezoid &other)
        : CPolygon(other) {}
    CTrapezoid &operator=(const CTrapezoid &other)
    {
        CPolygon::operator=(other);
        return *this;
    }
    ~CTrapezoid() {}
};

int main()
{
    double a = 1.6565695, b = 2.56;
    std::cout << a / b << std::endl;
    CPoint p1(1, 2);
    std::cout << p1.Radius() << std::endl;
    return 0;
}