#ifndef LINEAR_REGRESSION_H
#define LINEAR_REGRESSION_H

#include <cassert>
//--
#include "Eigen/Eigen"
#include "Eigen/Dense"
//--
#define ei_assert assert
#define EIGEN2_SUPPORT
#include "Eigen/LeastSquares"

/**
 *  Helper functions for fitting a line to a set of points in
 *  R^2.
 */
namespace LinearRegression{
    
    template<typename T>
    struct Line{
        T m,b; //Slope and y intercept, respectively

        /**
         *  Get the y-value of the line at a specific point.
         *  @param x is the x value.  Returns f(x) = mx + b
         */
        inline T operator()(const T x)const{ return m * x + b; }

        /// \brief Parameterized constructor.
        Line(T myM = (T)1.0, T myB = (T)0.0) :
            m(myM), b(myB) {}
        
        /// \brief Overloaded stream insertion operator
        friend std::ostream& operator<<(std::ostream& os, const Line& line){
            os << "y = " << line.m << "x + " << line.b;
            return os;
        }

    };

    
    /**
     *  Fit a line to a set of points.  Optinally return the residual.
     *
     *  @param N is the number of points.  Must be >= 2.
     *  @param dataPoints is a pointer to a vector of points in R^2.
     *  @param residual returns the residual error.  If NULL this is ignored.  Defaults to NULL.
     */
    template<typename T>
    struct Line<T> linearRegression(
        int N,
        const Eigen::Matrix<T,Eigen::Dynamic,2>* dataPoints,
        T* residual = NULL)
    {
        assert(N >= 2);
        assert(dataPoints != NULL);
        assert(dataPoints->rows() >= N);

        Line<T> retLine; //Declare line to return
        
        //Compute parameters of the line
        //TODO: Delete unnedded vars
        T xBar, yBar, xSum, ySum, xySum, xSqSum;
        xBar = yBar = xSum = ySum = xySum = xSqSum = static_cast<T>(0.0);
        const T invNf = static_cast<T>(1.0) / static_cast<float>(N);
        const T Nf    = static_cast<float>(N);
        for(int i = 0; i < N; i++){
            const T x = (*dataPoints)(i,0);
            const T y = (*dataPoints)(i,1);
            xSum   += x;
            ySum   += y;
            xySum  += x*y;
            xSqSum += x*x;
        }
        xBar = xSum * invNf;
        yBar = ySum * invNf;

        const T m = (Nf * xySum - (xSum * ySum)) / (Nf * xSqSum - (xSum * xSum));
        const T b = invNf * (ySum - m*xSum);

        retLine.m = m;
        retLine.b = b;

        //potentially compute ressidual
        if(residual != NULL){
            *residual = (T)(0.0);
            for(int i = 0; i < N; i++){
                const T lineVal = retLine( (*dataPoints)(i,0) );
                const T goalVal = (*dataPoints)(i,1);
                *residual +=
                    (goalVal - lineVal) * (goalVal - lineVal);
            }
        }

        return retLine;
    }
}


#endif //LINEAR_REGRESSION_H
