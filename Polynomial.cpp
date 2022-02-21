#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

class CPolynomial
{
public:
    explicit CPolynomial()
        : degree_(0), coeff_(new double[degree_ + 1])
    {
        coeff_[0] = 0;
    }
    explicit CPolynomial(size_t degree, double *coeff)
    {
        degree_ = 0;
        for (int i = 0; i <= degree; ++i)
        {
            if (coeff[i] != 0)
            {
                degree_ = i;
            }
        }
        coeff_ = new double[degree_ + 1];
        memcpy(coeff_, coeff, sizeof(double) * (degree_ + 1));
    }
    CPolynomial(const CPolynomial &other)
        : degree_(other.degree_), coeff_(new double[degree_ + 1])
    {
        memcpy(coeff_, other.coeff_, sizeof(double) * (degree_ + 1));
    }
    ~CPolynomial()
    {
        delete[] coeff_;
    }
    CPolynomial &operator=(const CPolynomial &other)
    {
        if (&other == this)
        {
            return *this;
        }
        degree_ = other.degree_;
        delete[] coeff_;
        coeff_ = new double[degree_];
        memcpy(coeff_, other.coeff_, sizeof(double) * (degree_ + 1));
        return *this;
    }
    bool operator!=(const CPolynomial &right) const
    {
        if (degree_ != right.degree_)
        {
            return true;
        }
        for (int i = 0; i <= degree_; ++i)
        {
            if (std::fabs(coeff_[i] - right.coeff_[i]) > 0.000001)
            {
                return true;
            }
        }
        return false;
    }
    bool operator==(const CPolynomial &right) const
    {
        return !(*this == right);
    }
    CPolynomial operator+(const CPolynomial &right) const
    {
        size_t degree = std::max(degree_, right.degree_);
        double *coeff = new double[degree + 1];
        for (int i = 0; i <= degree; ++i)
        {
            if (i <= std::min(degree_, right.degree_))
            {
                coeff[i] = coeff_[i] + right.coeff_[i];
            }
            else
            {
                coeff[i] = degree_ > right.degree_ ? coeff_[i] : right.coeff_[i];
            }
        }
        CPolynomial result{degree, coeff};
        return result;
    }
    CPolynomial operator-(const CPolynomial &right) const
    {
        size_t degree = std::max(degree_, right.degree_);
        double *coeff = new double[degree + 1];
        for (int i = 0; i <= degree; ++i)
        {
            if (i <= std::min(degree_, right.degree_))
            {
                coeff[i] = coeff_[i] - right.coeff_[i];
            }
            else
            {
                coeff[i] = degree_ > right.degree_ ? coeff_[i] : right.coeff_[i];
            }
        }
        CPolynomial result{degree, coeff};
        return result;
    }
    CPolynomial operator-() const
    {
        CPolynomial result{*this};
        for (int i = 0; i <= degree_; ++i)
        {
            result.coeff_[i] = -coeff_[i];
        }
        return result;
    }
    CPolynomial &operator+=(const CPolynomial &right)
    {
        CPolynomial result = *this + right;
        *this = result;
        return *this;
    }
    CPolynomial &operator-=(const CPolynomial &right)
    {
        CPolynomial result = *this - right;
        *this = result;
        return *this;
    }
    CPolynomial operator*(double right) const
    {
        size_t degree = degree_;
        double *coeff = new double[degree + 1];
        for (int i = 0; i <= degree; ++i)
        {
            coeff[i] = coeff_[i] * right;
        }
        CPolynomial result{degree, coeff};
        return result;
    }
    friend CPolynomial operator*(double left, const CPolynomial &right)
    {
        return right * left;
    }
    CPolynomial operator*(const CPolynomial &right) const
    {
        size_t degree = degree_ + right.degree_;
        double *coeff = new double[degree + 1];
        for (int i = 0; i <= degree; ++i)
        {
            coeff[i] = 0;
        }
        for (int i = 0; i <= degree_; ++i)
        {
            for (int j = 0; j <= right.degree_; ++j)
            {
                coeff[i + j] += coeff_[i] * right.coeff_[j];
            }
        }
        CPolynomial result{degree, coeff};
        return result;
    }
    CPolynomial operator/(double right) const
    {
        if (right == 0)
        {
            CPolynomial result;
            return result;
        }
        CPolynomial result;
        result.degree_ = degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] /= right;
        }
        return result;
    }
    CPolynomial &operator*=(double right)
    {
        CPolynomial result = *this * right;
        *this = result;
        return *this;
    }
    CPolynomial &operator*=(const CPolynomial &right)
    {
        CPolynomial result = *this * right;
        *this = result;
        return *this;
    }
    CPolynomial &operator/=(double right)
    {
        CPolynomial result = *this / right;
        *this = result;
        return *this;
    }
    const double &operator[](size_t i) const
    {
        return coeff_[i];
    }
    friend std::ostream &operator<<(std::ostream &stream, const CPolynomial &polynomial);
    friend CPolynomial &operator>>(std::istream &stream, CPolynomial &polynomial);

private:
    size_t degree_;
    double *coeff_;
};

std::ostream &operator<<(std::ostream &stream, const CPolynomial &polynomial)
{
    if (polynomial.degree_ == 0)
    {
        stream << polynomial.coeff_[0];
    }
    else
    {
        stream << polynomial.coeff_[polynomial.degree_] << "x^" << polynomial.degree_;
        for (int i = polynomial.degree_ - 1; i > 0; --i)
        {
            if (polynomial.coeff_[i] > 0)
            {
                stream << ' + ' << polynomial.coeff_[i] << "x^" << i;
            }
            if (polynomial.coeff_[i] < 0)
            {
                stream << ' - ' << -polynomial.coeff_[i] << "x^" << i;
            }
        }
        if (polynomial.coeff_[0] > 0)
        {
            stream << " + " << polynomial.coeff_[0];
        }
        if (polynomial.coeff_[0] < 0)
        {
            stream << " - " << -polynomial.coeff_[0];
        }
    }
    return stream;
}
CPolynomial &operator>>(std::istream &stream, CPolynomial &polynomial)
{
    double coeff;
    stream >> coeff;
    char symbol = stream.get();
    if(symbol == '\n')
    {
        polynomial.degree_ = 0;
        delete[] polynomial.coeff_;
        polynomial.coeff_ = new double[1];
        polynomial.coeff_[0] = coeff;
    }
    while(symbol!='\n')
    {
        while(symbol)
        symbol = stream.get();
    }
    return polynomial;
}

int main()
{
    return 0;
}