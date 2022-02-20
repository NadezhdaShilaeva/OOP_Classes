#include <iostream>
#include <algorithm>
#include <cmath>

class CPolynomial
{
public:
    explicit CPolynomial()
        : degree_(0), coeff_(new double[degree_ + 1])
    {
        coeff_[0] = 0;
    }
    explicit CPolynomial(unsigned degree, double *coeff)
        : degree_(degree), coeff_(new double[degree_ + 1]) ///////////
    {
        if (sizeof(coeff) / sizeof(double) != degree + 1)
        {
            std::cout << "Error! Different degrees!" << std::endl;
        }
        else
        {
            memcpy(coeff_, coeff, sizeof(double) * (degree_ + 1));
        }
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
    friend bool operator!=(const CPolynomial &left, const CPolynomial &right)
    {
        if (left.degree_ != right.degree_)
        {
            return true;
        }
        for (int i = 0; i <= left.degree_; ++i)
        {
            if (left.coeff_[i] != right.coeff_[i])
            {
                return true;
            }
        }
        return false;
    }
    friend bool operator==(const CPolynomial &left, const CPolynomial &right)
    {
        return !(left == right);
    }
    friend CPolynomial operator+(const CPolynomial &left, const CPolynomial &right)
    {
        CPolynomial result;
        result.degree_ = std::max(right.degree_, left.degree_);
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            if (i <= std::min(left.degree_, right.degree_))
            {
                result.coeff_[i] = left.coeff_[i] + right.coeff_[i];
            }
            else
            {
                result.coeff_[i] = left.degree_ > right.degree_ ? left.coeff_[i] : right.coeff_[i];
            }
        }
        return result;
    }
    friend CPolynomial operator-(const CPolynomial &left, const CPolynomial &right)
    {
        CPolynomial result;
        result.degree_ = std::max(right.degree_, left.degree_);
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            if (i <= std::min(left.degree_, right.degree_))
            {
                result.coeff_[i] = left.coeff_[i] - right.coeff_[i];
            }
            else
            {
                result.coeff_[i] = left.degree_ > right.degree_ ? left.coeff_[i] : right.coeff_[i];
            }
        }
        return result;
    }
    friend CPolynomial operator-(const CPolynomial &other)
    {
        CPolynomial result;
        result.degree_ = other.degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] = -other.coeff_[i];
        }
        return result;
    }
    CPolynomial &operator+=(const CPolynomial &other) ////////
    {
        degree_ = std::max(right.degree_, left.degree_);
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            if (i <= std::min(left.degree_, right.degree_))
            {
                result.coeff_[i] = left.coeff_[i] + right.coeff_[i];
            }
            else
            {
                result.coeff_[i] = left.degree_ > right.degree_ ? left.coeff_[i] : right.coeff_[i];
            }
        }
        return *this;
    }
    friend CPolynomial operator*(const CPolynomial &left, double num)
    {
        CPolynomial result;
        result.degree_ = left.degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] = left.coeff_[i] * num;
        }
        return result;
    }
    friend CPolynomial operator*(double num, const CPolynomial &right)
    {
        CPolynomial result;
        result.degree_ = right.degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] = right.coeff_[i] * num;
        }
        return result;
    }
    friend CPolynomial operator*(const CPolynomial &left, const CPolynomial &right) ////////
    {
        CPolynomial result;
        result.degree_ = left.degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] *= num;
        }
        return result;
    }
    friend CPolynomial operator/(const CPolynomial &left, double num)
    {
        CPolynomial result;
        result.degree_ = left.degree_;
        result.coeff_ = new double[result.degree_];
        for (int i = 0; i <= result.degree_; ++i)
        {
            result.coeff_[i] /= num;
        }
        return result;
    }
    CPolynomial &operator*=(double num)
    {
        for (int i = 0; i <= degree_; ++i)
        {
            coeff_[i] *= num;
        }
        return *this;
    }
    CPolynomial &operator *=(const CPolynomial& other)
    {
        //////////
    }
    CPolynomial &operator/=(double num)
    {
        for (int i = 0; i <= degree_; ++i)
        {
            coeff_[i] /= num;
        }
        return *this;
    }
    const double &operator[](unsigned i) const
    {
        return coeff_[i];
    }

private:
    unsigned degree_;
    double *coeff_;
};

int main()
{
    return 0;
}