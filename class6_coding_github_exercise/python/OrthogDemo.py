# Define a function which verifies whether two sine waves are orthogonal or not

def OrthogDemo(n, m):    
    
    # Import required packages
    import numpy as np
    from scipy.integrate import quad
    
    # Define a linear time array of 100 values
    t = np.arange(0, 1, 0.01)
    
    # Define waves 1 and 2 based on the inputs of the function
    wave1 = np.sin(2*np.pi*n*t)
    wave2 = np.sin(2*np.pi*m*t)
    
    # Multiply the waves together 
    product = wave1*wave2
    
    # Define a function to calculate the integral 
    def product_fcn(t):
        return np.sin(2*np.pi*n*t)*np.sin(2*np.pi*m*t)

    # Calculate the integral and error of the product 
    integral,err = quad(product_fcn,0,1)
    
    # Define an output class out and store the data in that
    class output:
        pass
    
    out = output()
    
    out.time = t
    out.wave1 = wave1
    out.wave2 = wave2
    out.product = product
    out.integral = integral
    
    return out

