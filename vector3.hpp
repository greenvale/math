/*
Vector3 class
William Denny, 29th August 2022
*/

namespace mathlib
{

    class Vector3
    {
    private:
        double m_x;
        double m_y;
        double m_z;
    
    public:
        Vector3();
        Vector3(const double& x, const double& y, const double& z);
        
        double getX();
        double getY();
        double getZ();
        
        void setX(const double& x);
        void setY(const double& y);
        void setZ(const double& z);
       
    }

}
