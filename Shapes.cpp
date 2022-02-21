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
    explicit CPolyline()
        : number_points_(0), points_(new CPoint[number_points_]) {}
    explicit CPolyline(size_t number_points, CPoint *points)
        : number_points_(number_points), points_(new CPoint[number_points])
    {
        if (sizeof(points) != number_points * sizeof(CPoint))
        {
            std::cerr << "Error! Different number of points!" << std::endl;
        }
        else
        {
            memcpy(points_, points, sizeof(CPoint) * number_points_);
        }
    }
    CPolyline(const CPolyline &other)
        : number_points_(other.number_points_), points_(new CPoint[other.number_points_])
    {
        memcpy(points_, other.points_, sizeof(CPoint) * number_points_);
    }
    CPolyline &operator=(const CPolyline &other)
    {
        if (&other == this)
            return *this;
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
    virtual size_t GetNumberPoints() const
    {
        return number_points_;
    }
    virtual CPoint GetPoint(size_t i) const
    {
        return points_[i];
    }
    virtual double Length() const
    {
        double length_ = 0;
        for (int i = 1; i < number_points_; ++i)
        {
            length_ += sqrt(pow((points_[i - 1].X() - points_[i].X()), 2) + (pow((points_[i - 1].Y() - points_[i].Y()), 2)));
        }
        return length_;
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
        if (this->GetNumberPoints() < 2)
            return 0;
        double length_ = 0;
        for (int i = 1; i < this->GetNumberPoints(); ++i)
        {
            length_ += std::sqrt(std::pow(this->GetPoint(i - 1).X() - this->GetPoint(i).X(), 2) + std::pow(this->GetPoint(i - 1).Y() - this->GetPoint(i).Y(), 2));
        }
        length_ += std::sqrt(std::pow(this->GetPoint(this->GetNumberPoints() - 1).X() - this->GetPoint(0).X(), 2) + std::pow(this->GetPoint(this->GetNumberPoints() - 1).Y() - this->GetPoint(0).Y(), 2));
        return length_;
    }
};

class CPolygon : public CClosedPolyline
{
public:
    explicit CPolygon() : CClosedPolyline() {}
    explicit CPolygon(size_t number_points, CPoint *points)
        : CClosedPolyline(number_points, points) {} // checking for self-intersections
    CPolygon(const CPolygon &other)
        : CClosedPolyline(other) {}
    CPolygon &operator=(const CPolygon &other)
    {
        CPolyline::operator=(other);
        return *this;
    }
    ~CPolygon() {}
    virtual double Perimeter() const
    {
        return this->Length();
    }
    virtual double Square() const
    {
        double square = 0;
        for (int i = 0; i < this->GetNumberPoints() - 1; ++i)
        {
            CPoint p = this->GetPoint(i);
            CPoint p_next = this->GetPoint(i + 1);
            square += p.X() * p_next.Y() - p.Y() * p_next.X();
        }
        CPoint p_end = this->GetPoint(this->GetNumberPoints() - 1);
        CPoint p_start = this->GetPoint(0);
        square += p_end.X() * p_start.Y() - p_end.Y() * p_start.X();
        square = std::fabs(square * 0.5);
        return square;
    }
};

class CRegularPolygon : CPolygon
{
    explicit CRegularPolygon() : CPolygon(), side_(0) {}
    explicit CRegularPolygon(size_t number_points, CPoint *points)
        : CPolygon(number_points, points), side_(0)
    {
        if (number_points != 0)
        {
            bool is_regular = true;
            side_ = std::sqrt(std::pow(points[number_points - 1].X() - points[0].X(), 2) + std::pow(points[number_points - 1].Y() - points[0].Y(), 2));
            for (int i = 1; i < number_points; ++i)
            {
                double dist = std::sqrt(std::pow(points[i - 1].X() - points[i].X(), 2) + std::pow(points[i - 1].Y() - points[i].Y(), 2));
                if (dist != side_)
                {
                    is_regular = false;
                    break;
                }
            }
            if (!is_regular)
                std::cerr << "Your polygon is not regular!";
        }
    }
    CRegularPolygon(const CRegularPolygon &other)
        : CPolygon(other), side_(other.side_) {}
    CRegularPolygon &operator=(const CRegularPolygon &other)
    {
        CPolyline::operator=(other);
        if (&other == this)
            return *this;
        this->side_ = other.side_;
        return *this;
    }
    ~CRegularPolygon() {}
    double Perimeter() const
    {
        double perimeter = this->GetNumberPoints() * side_;
        return perimeter;
    }
    double Square() const
    {
        double square = (this->GetNumberPoints() * side_ * side_) / (4 * std::tan(PI / this->GetNumberPoints()));
        return square;
    }

private:
    double side_;
};

class CTriangle : CPolygon
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
    CTriangle(const CTriangle &other)
        : CPolygon(other) {}
    CTriangle &operator=(const CTriangle &other)
    {
        CPolyline::operator=(other);
        return *this;
    }
    ~CTriangle() {}
};

class CTrapezoid : CPolygon
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
    CTrapezoid(const CTrapezoid &other)
        : CPolygon(other) {}
    CTrapezoid &operator=(const CTrapezoid &other)
    {
        CPolyline::operator=(other);
        return *this;
    }
    ~CTrapezoid() {}
};

int main()
{
    double a = 1.6565695, b = 2.56;
    std::cout << a/b << std::endl;
    CPoint p1(1, 2);
    std::cout << p1.radius() << std::endl;
    return 0;
}