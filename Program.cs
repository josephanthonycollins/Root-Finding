using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;


//Name: Joseph Collins
//ID: 98718584
//Module: AM6004
//Assignment 1

namespace ConsoleApplication1
{
    //we will use delegates to implement the functions
    public delegate double funcref(double val);

    class Program
    {
        static void Main(string[] args)
        {
            Solver s = new Solver();

            funcref f1 = new funcref(val => Math.Exp(val) + Math.Pow(2, -val) + 2 * Math.Sin(val) - 6);
            funcref f2 = new funcref(val => 2 * val * Math.Cos(2 * val) - (val - 2) * (val - 2));
            funcref f3 = new funcref(val => Math.Exp(val) - 3 * Math.Pow(val, 2));
            funcref f = new funcref(val => Math.Exp(val) + Math.Pow(2, -val) + 2 * Math.Sin(val) - 6);

            double linterval = 0;
            double uinterval = 0;

            //Here we allow the user to specify which one of the three equations they would like to find the roots of
            string str = "";
            int choice = 0;
            Console.Write("Three functions.\n\nf1(x) = e^x+2^-x+2Sin(x)-6\nf2(x) = 2xCos(2x) - (x-2)^2\nf3(x) = e^x -3x^2\n\nPlease enter an integer between 1 and 3:");
            str = Console.ReadLine();
            bool b = int.TryParse(str, out choice);
            while (!b || ((choice < 1) || (choice > 3)))
            {
                Console.Write("\nInvalid Input. Please renter an integer between 1 and 3: ");
                str = Console.ReadLine();
                b = int.TryParse(str, out choice);
            }

            //user chooses which function to examine
            switch (choice)
            {
                case 1:
                    Console.Write("\nWe will work with function f1(x). Enter number of intervals:");
                    f = f1;
                    linterval = -4;
                    uinterval = 4;
                    break;
                case 2:
                    Console.Write("\nWe will work with function f2(x).Enter number of intervals:");
                    f = f2;
                    linterval = 0;
                    uinterval = 4;
                    break;
                case 3:
                    Console.Write("\nWe will work with function f3(x).Enter number of intervals:");
                    linterval = -2;
                    uinterval = 5;
                    f = f3;
                    break;
            }

            //user chooses the number of intervals
            int intervals = 0;
            str = Console.ReadLine();
            b = int.TryParse(str, out intervals);
            while (!b || ((intervals < 1) || (intervals > 1000)))
            {
                Console.Write("\nInvalid Input. Please enter an integer between 1 and 1000: ");
                str = Console.ReadLine();
                b = int.TryParse(str, out intervals);
            }

            //user chooses a tolerance level
            Console.Write("\nNow enter the tolerance level:");
            double tolerance = 0;
            str = Console.ReadLine();
            b = double.TryParse(str, out tolerance);
            while (!b || (tolerance <= 0))
            {
                Console.Write("\nInvalid Input. Please enter a double >=0: ");
                str = Console.ReadLine();
                b = double.TryParse(str, out tolerance);
            }

            //bracketing  
            Console.Write("\nBracketing Output\n");
            Console.Write("*****************\n");
            s.Bracketing(uinterval, linterval, intervals, f);
            s.print_bracketing_roots();

            //Bisection
            Console.Write("\nBisection Method\n");
            Console.Write("*****************\n");
            s.Bisection(tolerance, f);
            s.print_roots();

            //Newton Raphson
            Console.Write("\nNewton Raphson Method\n");
            Console.Write("*****************\n");
            s.NewtonRaphson(tolerance, f);
            s.print_roots_iterations();

            int milliseconds = 3000;
            Thread.Sleep(milliseconds);

            //Richardson extrapolation with error control
            //user chooses a tolerance level
            Console.Write("\nRichardson Extrapolation with error Control\n");
            Console.Write("*********************************************\n");
            Console.Write("Provide a double at which you wish to evaluate the derivative:");
            double x = 0;
            str = Console.ReadLine();
            b = double.TryParse(str, out x);
            while (!b)
            {
                Console.Write("\nInvalid Input. Please enter a double: ");
                str = Console.ReadLine();
                b = double.TryParse(str, out x);
            }
            s.richardson_extrapolation_errorcontrol(f, x, 0.5, 40, tolerance);
            s.print_richardson_output();

            milliseconds = 5000;
            Thread.Sleep(milliseconds);

            //Advanced test
            Console.Write("\n\n**********************************************************************");
            Console.Write("\n**********************OPTICAL WAGEGUIDE*******************************\n");
            Console.Write("\nBracketing Output\n");
            Console.Write("*****************\n");

            Waveguide w = new Waveguide();
            funcref func = new funcref(w.Calc);

            double lambda = 1.3, nc = 3.44, ns = 3.34;
            double k0 = 2.0 * Math.PI / lambda;
            linterval = k0 * ns;
            uinterval = k0 * nc;
            //bracketing  
            intervals = 20;
            s.Bracketing(uinterval, linterval, intervals, func);
            s.print_bracketing_roots();
            milliseconds = 3000;
            Thread.Sleep(milliseconds);
            //Bisection
            Console.Write("\nBiscection Output\n");
            Console.Write("*****************\n");
            s.Bisection(tolerance, func);
            s.print_roots();
            milliseconds = 3000;
            Thread.Sleep(milliseconds);

            //Newton Raphson Richardson Error Control
            Console.Write("\nNewton Raphson Error Control Output\n");
            Console.Write("*************************************\n");
            s.NewtonRaphson_Richardson_ErrorControl(tolerance, func);
            s.print_roots_iterations();
            milliseconds = 3000;
            Thread.Sleep(milliseconds);
            Console.Write("\n\n*************************************\n");
            Console.Write("End of Program\n");
            Console.Write("*************************************\n");
            Console.ReadLine();
        }
    }

    //this class is used to create a function which passed as a delegate
    class Waveguide
    {

        public double Calc(double val)
        {
            double lambda = 1.3, W = 3.0, nc = 3.44, ns = 3.34, ncl = 1.0;
            double k0 = 2.0 * Math.PI / lambda;
            double linterval = k0 * ns;
            double uinterval = k0 * nc;
            double qb = Math.Sqrt(val * val - k0 * k0 * ncl * ncl);
            double pb = Math.Sqrt(val * val - k0 * k0 * ns * ns);
            double hb = Math.Sqrt(k0 * k0 * nc * nc - val * val);
            double mult1 = 1 - (pb * qb) / (hb * hb);
            double mult2 = ((qb) / (hb)) + ((pb) / (hb));
            return mult1 * Math.Sin(W * hb) - mult2 * Math.Cos(W * hb);
        }

        //end of Waveguide class
    }

    //this class will contain the bracketing, bisection, newtron raphson, richardson methods
    //it will also contain attributes which list the roots, the errrors, the intervals etc
    class Solver
    {
        //this list will contain the right intervals, will be utilised when we want to search for roots
        List<double> right = new List<double>();

        public List<double> Right
        {
            get { return right; }
            set { right = value; }
        }

        //this list will contain the left intervals, will be utilised when we want to search for roots
        List<double> left = new List<double>();

        public List<double> Left
        {
            get { return left; }
            set { left = value; }
        }

        //this list will contain the roots
        List<double> roots = new List<double>();

        public List<double> Roots
        {
            get { return roots; }
            set { roots = value; }
        }

        //this list will contain the number of iterations required required for the NewtonRaphson to converge
        List<double> nriterations = new List<double>();

        public List<double> NRiterations
        {
            get { return nriterations; }
            set { nriterations = value; }
        }

        //this list will contain the output associated with running the richardson_extrapolation_errorcontrol() method
        List<double> richardsonoutput = new List<double>();

        public List<double> Richardsonoutput
        {
            get { return richardsonoutput; }
            set { richardsonoutput = value; }
        }

        //bracketing() method
        public void Bracketing(double u, double l, double intervals, funcref thefunc)
        {
            double upper = 0;
            double lower = 0;
            double dx = 0;
            //firstly we need to ensure that the Right and Left lists are empty
            if (Right.Count() > 0)
            {
                Right.Clear();
                Left.Clear();
            }

            dx = (u - l) / intervals;

            lower = l;
            upper = lower + dx;
            while (upper <= u)
            {
                //test the bracket
                if (thefunc(lower) * thefunc(upper) < 0)
                {
                    //there is a root in this interval so add the values to the relevants list
                    Left.Add(lower);
                    Right.Add(upper);
                }
                else
                {
                    //do nothing
                }

                //move to the next bracket
                lower = upper;
                upper = upper + dx;
            }

            //end of bracketing() method
        }

        //method to print to screen output of bracketing() method
        public void print_bracketing_roots()
        {
            if (Right.Count() > 0)
            {
                int z = 0;
                Console.WriteLine("\nThere are {0} roots.\n", Right.Count());
                foreach (double d in Right)
                {
                    Console.Write("Root in interval ({0,3:F6}, {1,3:F6})\n", Left[z], d);
                    z++;
                }
            }
            else
                Console.WriteLine("No roots were found.\n");
            //end of print_bracketing_roots() method
        }

        //bisection() method
        public void Bisection(double tol, funcref thefunc)
        {
            //clear the roots list
            if (Roots.Count() > 0)
                Roots.Clear();

            //ensure that there are intervals upon which to apply the bisection method)
            if (Right.Count() == 0)
            {
                Console.WriteLine("No intervals are given, hence we cannot run the bisection method.\n");
            }
            else
            {
                int z = 0;
                foreach (double d in Right)
                {
                    double upper = d;
                    double lower = Left[z];

                    int counter = 0;

                    //check upper and lower intervals are closer than the tolerance level
                    while ((upper - lower > tol) & (counter < 1000))
                    {
                        double m = (upper + lower) / 2;
                        double s = thefunc(lower) * thefunc(m);
                        if (s < 0)
                        {
                            upper = m;
                        }
                        else
                        {
                            lower = m;
                        }

                        //prevent an infinite loop
                        counter++;
                        if (counter > 1000)
                            Console.WriteLine("The Bisection method has done more than 1000 iterations and is still not within the specified tolerance.We will now assume that the root is given by the last iteration\n");
                    }

                    double val = (upper + lower) / 2;
                    Roots.Add(val);
                    z++;
                }

            }
            //end of Bisection() method    
        }

        //method to print to screen output roots
        public void print_roots()
        {
            if (Roots.Count() > 0)
            {
                int z = 0;
                Console.WriteLine("There are {0} roots.\n", Roots.Count());
                foreach (double d in Roots)
                {
                    Console.Write("Root {0,3:F6} in interval ({1,3:F6}, {2,3:F6})\n", d, Left[z], Right[z]);
                    z++;
                }
            }
            else
                Console.Write("No roots were found.\n");
            //end of print_roots() method
        }

        //method to print to screen Roots and NRiterations
        public void print_roots_iterations()
        {
            if (Roots.Count() > 0)
            {
                int z = 0;
                Console.WriteLine("There are {0} roots.\n", Roots.Count());
                foreach (double d in Roots)
                {
                    Console.Write("Root {0,3:F6} in interval ({1,3:F6}, {2,3:F6}), {3} iterations\n", d, Left[z], Right[z], NRiterations[z]);
                    z++;
                }
            }
            else
                Console.WriteLine("No roots were found.\n");
            //end of print_roots_iterations() method
        }

        //recursive method to implement Richardsons extrapolation
        public double richardson_extrapolation(funcref thefunc, double x, double h, double j)
        {
            double dval1 = 0;
            double dval2 = 0;
            if (j == 1)
                return (1 / (2 * h)) * (thefunc(x + h) - thefunc(x - h));
            else
            {
                dval1 = richardson_extrapolation(thefunc, x, h / 2.0, j - 1);
                dval2 = richardson_extrapolation(thefunc, x, h, j - 1);
                return dval1 + (1 / (Math.Pow(4, (double)j - 1) - 1)) * (dval1 - dval2);
            }

            //end of richardson_extrapolation() method
        }

        //method to implement richardson extrapolation with error control
        public void richardson_extrapolation_errorcontrol(funcref thefunc, double x, double h, int maxiter, double toler)
        {
            double deriv = 0;
            double derivold = 0;
            double err = 0;
            double errold = 1000000000;
            int itercount = 0;
            int count = 0;

            //firstly ensure that the "Richardsonoutput" list is empty
            if (Richardsonoutput.Count() > 0)
                Richardsonoutput.Clear();

            //next run the algorithm and populate the l8ist
            for (count = 1; count <= maxiter; count++)
            {
                //get the derivative estimate at iteration "count"
                deriv = richardson_extrapolation(thefunc, x, h, count);
                err = Math.Abs(deriv - derivold);
                //is the error less than a specified tolerance
                if (err < toler)
                {
                    itercount = count;
                    Richardsonoutput.Add(deriv);
                    Richardsonoutput.Add(err);
                    Richardsonoutput.Add((double)count);
                    break;
                }
                //is the latest error greatest than the previous error
                if (err > errold)
                {
                    deriv = derivold;
                    err = errold;
                    itercount = count - 1;
                    Richardsonoutput.Add(deriv);
                    Richardsonoutput.Add(err);
                    Richardsonoutput.Add((double)count);
                }
                derivold = deriv;
                errold = err;
                Richardsonoutput.Add(deriv);
                Richardsonoutput.Add(err);
                Richardsonoutput.Add((double)count);
            }

        }

        //method to print to screen output associated with richardson_extrapolation_errorcontrol() method
        public void print_richardson_output()
        {
            if (Richardsonoutput.Count() > 0)
            {
                int z = Richardsonoutput.Count();
                int counter = z / 3;
                int i = 0;
                for (i = 0; i < counter; i++)
                {
                    Console.Write("\n(deriv,err,iter) = ({0,3:F6},{1,3:F6},{2})", Richardsonoutput[0 + i * 3], Richardsonoutput[1 + i * 3], Richardsonoutput[2 + i * 3]);
                }
            }
            else
                Console.WriteLine("No output associated with Richardson extrapoloation with error control, please ensure the method() is called.\n");
            //end of method
        }


        //NewtonRaphson() method
        public void NewtonRaphson(double tol, funcref thefunc)
        {
            //clear the roots list
            if (Roots.Count() > 0)
                Roots.Clear();
            if (NRiterations.Count() > 0)
                NRiterations.Clear();

            //ensure that there are intervals upon which to apply the Newton Raphson method)
            if (Right.Count() == 0)
            {
                Console.WriteLine("No intervals are given, hence we cannot run the Newton Raphson method.\n");
            }
            else
            {
                int z = 0;

                //running Newton Raphson on each interval
                foreach (double d in Right)
                {
                    double upper = d;
                    double lower = Left[z];
                    double rinit = (upper + lower) / 2;
                    double rold = rinit;
                    double numer = 0, denom = 0, dx = 0;
                    double rnew = 0;

                    int counter = 0;

                    //we will use "counter" to ensure Newton Raphson does not loop infinitely
                    while (counter < 1000)
                    {
                        numer = thefunc(rold);
                        /***********************************************************make to call Richardsonextrapolation in the correct manner*****/
                        denom = richardson_extrapolation(thefunc, rold, 0.5, 10);
                        if (Double.IsNaN(denom))
                        {
                            Console.WriteLine("Stopping Newton Raphson in interval ({0},{1}) as f'[x] at {2} is not defined", Left[z], d, rold);
                            break;
                        }
                        if (denom == 0.0)
                        {
                            Console.WriteLine("Stopping Newton Raphson in interval ({0},{1}) as f'[x] at {2} equals 0", Left[z], d, rold);
                            break;
                        }
                        dx = numer / denom;
                        rnew = rold - dx;
                        //if the change "dx" is less than the tolerance, then we have found an estimate of the root
                        if (Math.Abs(dx) < Math.Abs(tol))
                        {
                            Roots.Add(rnew);
                            NRiterations.Add(counter + 1);
                            break;
                        }

                        rold = rnew;

                        //prevent an infinite loop
                        counter++;
                        if (counter > 1000)
                        {
                            Console.WriteLine("The Newton Raphson method has done more than 1000 iterations and is still not within the specified tolerance.We will now assume that the root is given by the last iteration\n");
                            Roots.Add(rnew);
                            NRiterations.Add(counter);
                        }
                    }

                    double val = rnew;
                    z++;
                    //end of Newton Raphson on interval
                }

            }
            //end of NewtonRaphson() method    
        }

        //NewtonRaphson_Richardson_ErrorControl() method
        public void NewtonRaphson_Richardson_ErrorControl(double tol, funcref thefunc)
        {
            //clear the roots list
            if (Roots.Count() > 0)
                Roots.Clear();
            if (NRiterations.Count() > 0)
                NRiterations.Clear();

            //ensure that there are intervals upon which to apply the Newton Raphson method)
            if (Right.Count() == 0)
            {
                Console.WriteLine("No intervals are given, hence we cannot run the Newton Raphson method.\n");
            }
            else
            {
                int z = 0;

                //running Newton Raphson on each interval
                foreach (double d in Right)
                {
                    double upper = d;
                    double lower = Left[z];
                    double rinit = (upper + lower) / 2;
                    double rold = rinit;
                    double numer = 0, denom = 0, dx = 0;
                    double rnew = 0;

                    int counter = 0;

                    //we will use "counter" to ensure Newton Raphson does not loop infinitely
                    while (counter < 1000)
                    {
                        numer = thefunc(rold);
                        richardson_extrapolation_errorcontrol(thefunc, rold, 0.025, 20, 0.0001);
                        int k = Richardsonoutput.Count();
                        denom = Richardsonoutput[k - 3];
                        if (Double.IsNaN(denom))
                        {
                            Console.WriteLine("Stopping Newton Raphson in interval ({0},{1}) as f'[x] at {2} is not defined", Left[z], d, rold);
                            break;
                        }
                        if (denom == 0.0)
                        {
                            Console.WriteLine("Stopping Newton Raphson in interval ({0},{1}) as f'[x] at {2} equals 0", Left[z], d, rold);
                            break;
                        }
                        dx = numer / denom;
                        rnew = rold - dx;
                        //if the change "dx" is less than the tolerance, then we have found an estimate of the root
                        if (Math.Abs(dx) < Math.Abs(tol))
                        {
                            Roots.Add(rnew);
                            NRiterations.Add(counter + 1);
                            break;
                        }

                        rold = rnew;

                        //prevent an infinite loop
                        counter++;
                        if (counter > 1000)
                        {
                            Console.WriteLine("The Newton Raphson method has done more than 1000 iterations and is still not within the specified tolerance.We will now assume that the root is given by the last iteration\n");
                            Roots.Add(rnew);
                            NRiterations.Add(counter);
                        }
                    }

                    double val = rnew;
                    z++;
                    //end of Newton Raphson on interval
                }

            }
            //end of NewtonRaphson() method    
        }


        //end of Solver class
    }
    //end of namespace
}
