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
        //rk4: Performs ONE integration step via RK4 Method.
        //--------------------------------------------------------------------
        public void rk4(double dt)
        {
            double[,] sl = new double[4,2] {{0.0 , 0.0} , {0.0 , 0.0} , {0.0 , 0.0} , {0.0 , 0.0}};
            double[] xi = new double[n];

            int i;
            //Solving for kA, kB, kC, and kD
            rhsFunc(x, f);
            for(i=0; i<n; ++i)
            {
                sl[0,i] = f[i];
                xi[i] = x[i] + sl[0,i]*0.5*dt;
            }

            rhsFunc(xi, f);
            for(i=0; i<n; ++i)
            {
                sl[1,i] = f[i];
                xi[i] = x[i] + sl[1,i]*0.5*dt;
            }

            rhsFunc(xi, f);
            for(i=0; i<n; ++i)
            {
                sl[2,i] = f[i];
                xi[i] = x[i] + sl[2,i]*dt;
            }

            rhsFunc(xi, f);
            for(i=0; i<n; ++i)
            {
                sl[3,i] = f[i];
                xi[i] = x[i] + sl[3,i]*0.5*dt;
            }

            for(i=0; i<n; ++i)
            {
                x[i] = x[i] + (sl[0,i] + 2*sl[1,i] + 2*sl[2,i] + sl[3,i])*dt/6;
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