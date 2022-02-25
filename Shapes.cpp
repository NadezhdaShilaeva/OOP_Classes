#include <iostream>
#include <cmath>
#include <string>

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
        return *this;
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

private:
    double x_;
    double y_;
};

double Distance(CPoint p1, CPoint p2)
{
    return std::sqrt(std::pow(p1.X() - p2.X(), 2) + std::pow(p1.Y() - p2.Y(), 2));
}

double Area(CPoint p1, CPoint p2, CPoint p3)
{
    return (p2.X() - p1.X()) * (p3.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (p3.X() - p1.X());
}

bool IntersectLine(double p1, double p2, double p3, double p4)
{
    if (p1 > p2)
    {
        std::swap(p1, p2);
    }
    if (p3 > p4)
    {
        std::swap(p3, p4);
    }
    return std::max(p1, p3) <= std::min(p2, p4);
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
        double length = 0;
        for (size_t i = 1; i < number_points; ++i)
        {
            length += Distance(this->Point(i), this->Point(i - 1));
        }
        length += Distance(this->Point(0), this->Point(number_points - 1));
        return length;
    }
};

class CPolygon
{
public:
    explicit CPolygon() : outline_() {}
    explicit CPolygon(size_t number_points, CPoint *points)
        : outline_(number_points, points)
    {
        if (number_points < 3)
        {
            std::cerr << "It is not a polygon!" << std::endl;
        }
        bool self_intersection = false;
        for (size_t i = 3; i < number_points; ++i)
        {
            CPoint p1 = points[i - 1];
            CPoint p2 = points[i];
            for (size_t j = 1; j < i - 1; ++j)
            {
                CPoint p3 = points[j - 1];
                CPoint p4 = points[j];
                if (IntersectLine(p1.X(), p2.X(), p3.X(), p4.X()) && IntersectLine(p1.Y(), p2.Y(), p3.Y(), p4.Y()) 
                && Area(p1, p2, p3) * Area(p1, p2, p4) <= 0 && Area(p3, p4, p1) * Area(p3, p4, p2) <= 0)
                {
                    self_intersection = true;
                    break;
                }
            }
        }
        CPoint p1 = points[number_points - 1];
        CPoint p2 = points[0];
        for (size_t j = 2; j < number_points - 1; ++j)
        {
            CPoint p3 = points[j - 1];
            CPoint p4 = points[j];
            if (IntersectLine(p1.X(), p2.X(), p3.X(), p4.X()) && IntersectLine(p1.Y(), p2.Y(), p3.Y(), p4.Y()) 
            && Area(p1, p2, p3) * Area(p1, p2, p4) <= 0 && Area(p3, p4, p1) * Area(p3, p4, p2) <= 0)
            {
                self_intersection = true;
                break;
            }
        }
        if (self_intersection)
        {
            std::cerr << "It is not a polygon!" << std::endl;
        }
    }
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
    virtual ~CPolygon() {}
    virtual double Perimeter() const
    {
        return outline_.Length();
    }
    virtual double Square() const
    {
        if (NumberPoints() < 3)
        {
            return 0;
        }
        double square = 0;
        for (size_t i = 0; i < NumberPoints() - 1; ++i)
        {
            CPoint p = outline_.Point(i);
            CPoint p_next = outline_.Point(i + 1);
            square += (p.X() + p_next.X()) * (p.Y() - p_next.Y());
        }
        CPoint p_end = outline_.Point(NumberPoints() - 1);
        CPoint p_start = outline_.Point(0);
        square += (p_end.X() + p_start.X()) * (p_end.Y() - p_start.Y());
        square = std::fabs(square * 0.5);
        return square;
    }
    virtual const std::string Name() const
    {
        return "Polygon";
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
public:
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
            {
                std::cerr << "This polygon is not regular!" << std::endl;
            }
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
        double perimeter = NumberPoints() * side_;
        return perimeter;
    }
    double Square() const
    {
        if (NumberPoints() < 3)
        {
            return 0;
        }
        double square = (NumberPoints() * side_ * side_) / (4 * std::tan(PI / NumberPoints()));
        return square;
    }
    const std::string Name() const
    {
        return "Regular polygon";
    }

private:
    double side_;
};

class CTriangle : public CPolygon
{
public:
    explicit CTriangle() : CPolygon() {}
    explicit CTriangle(size_t number_points, CPoint *points)
        : CPolygon(number_points, points)
    {
        if (number_points != 3)
        {
            std::cerr << "It is not a triangle!" << std::endl;
        }
    }
    CTriangle(const CTriangle &other)
        : CPolygon(other) {}
    CTriangle &operator=(const CTriangle &other)
    {
        CPolygon::operator=(other);
        return *this;
    }
    ~CTriangle() {}
    const std::string Name() const
    {
        return "Triangle";
    }
};

class CTrapezoid : public CPolygon
{
public:
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
            std::cerr << "It is not a triangle!" << std::endl;
        }
    }
    CTrapezoid(const CTrapezoid &other)
        : CPolygon(other) {}
    CTrapezoid &operator=(const CTrapezoid &other)
    {
        CPolygon::operator=(other);
        return *this;
    }
    ~CTrapezoid() {}
    const std::string Name() const
    {
        return "Trapezoid";
    }
};

int main()
{
    size_t n;
    std::cin >> n;
    CPolygon **shapes = new CPolygon *[n];
    for (size_t i = 0; i < n; ++i)
    {
        size_t num_points;
        std::cin >> num_points;
        CPoint *points = new CPoint[num_points];
        for (size_t j = 0; j < num_points; ++j)
        {
            double x, y;
            std::cin >> x >> y;
            points[j] = CPoint{x, y};
        }
        if (num_points == 3)
        {
            shapes[i] = new CTriangle{num_points, points};
        }
        else
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
            if (num_points == 4 && is_parallel && not_parallel)
            {
                shapes[i] = new CTrapezoid{num_points, points};
            }
            else
            {
                bool is_regular = true;
                double side = Distance(points[num_points - 1], points[0]);
                for (size_t i = 1; i < num_points; ++i)
                {
                    double dist = Distance(points[i - 1], points[i]);
                    if (dist != side)
                    {
                        is_regular = false;
                        break;
                    }
                }
                if (is_regular)
                {
                    shapes[i] = new CRegularPolygon{num_points, points};
                }
                else
                {
                    shapes[i] = new CPolygon{num_points, points};
                }
            }
        }
    }
    for (size_t i; i < n; ++i)
    {
        std::cout << shapes[i]->Name() << " " << shapes[i]->Perimeter() << " " << shapes[i]->Square() << std::endl;
    }
    return 0;
}

// 4
// 4
// 2 3 4 1 -2 -1 -1 2
// 3
// 2 1 7 8 4 5
// 4
// 0 3 3 0 0 -3 -3 0
// 5
// 3 4 5 11 12 8 9 5 5 6