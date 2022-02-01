using System;

namespace sim
{
    public class SimplePend
    {
        int n = 2;
        private double len = 1.1;
        private double g = 9.81;

        private double[] x;
        private double[] f;

        //--------------------------------------------------------------------
        //SimplePend Constructor
        //--------------------------------------------------------------------
        public SimplePend()
        {
            x = new double[n];
            f = new double[n];

            x[0] = 1.0;
            x[1] = 0.0;
        }

        //--------------------------------------------------------------------
        //rhsFunc: Calculates right hand side of the pendulum's ODE.
        //--------------------------------------------------------------------
        public void rhsFunc(double[] st, double[] ff)
        {
            ff[0] = st[1];
            ff[1] = -(g/len)*Math.Sin(st[0]);
        }

        //--------------------------------------------------------------------
        //step: Performs ONE integration step via Euler's Method.
        //--------------------------------------------------------------------
        public void step(double dt)
        {
            rhsFunc(x, f);
            int i;
            for(i=0; i<n; ++i)
            {
                x[i] = x[i] + f[i]*dt;
            }
        }

        //--------------------------------------------------------------------
        //Getters and Setters Methods
        //--------------------------------------------------------------------
        public double L
        {
            get {return(len);}

            set
            {
                if(value > 0.0)
                {
                    len = value;
                }
            }
        }

        public double G
        {
            get {return(g);}

            set
            {
                if(value > 0.0)
                {
                    g = value;
                }
            }
        }

        public double theta
        {
            get {return x[0];}

            set {x[0] = value;}
        }

        public double thetaDot
        {
            get {return x[1];}

            set {x[1] = value;}
        }
    }

}