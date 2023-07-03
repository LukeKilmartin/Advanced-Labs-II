from scipy import odr
import numpy as np
#Used to change python arrays into a latex table format
def longtable_2f(array_of_arrays):
    transposed_array_of_arrays=array_of_arrays.T
    num_of_arrays,length_of_each_array=np.shape(transposed_array_of_arrays)
    length=0
    num=0
    while num < num_of_arrays-1:
        length=0
        while length < length_of_each_array-1:
            print(f"{transposed_array_of_arrays[num,length]:.2f}",end="")
            print(" & ",end="")
            length+=1
        print(f"{transposed_array_of_arrays[num,length]:.2f}",end="")
        print(" \\\\")
        num+=1
    length=0
    while length < length_of_each_array-1:
        print(f"{transposed_array_of_arrays[num,length]:.2f}",end="")
        print(" & ",end="")
        length+=1
    print(f"{transposed_array_of_arrays[num,length]:.2f}",end="")
    num+=1

#Linear Chi Square
def lin_chi_sqr(x, y, parameter_array, y_errs): 
    #Defining my own (generalised) function for calculating the (linear fit) 𝜒2  sum so I can use it in future.
    #See background information in cells at start if reading this in future.    
    
    #imports
    from scipy.stats import chi2  
    import numpy as np
    from scipy.optimize import curve_fit

    #p0 is an array of starting values for the best-fit parameter which are close to the optimum values. 
    #For linear fit, this means [slope, y-intercept].
    def func(x, m, c):
        return m*x+c
    popt, pcov = curve_fit(func, x, y, p0=parameter_array, sigma=y_errs, absolute_sigma=True)
    slope = popt[0]
    y_int = popt[1]
    y_fit = func(x, slope, y_int)
    return np.sum(((y-y_fit)/y_errs)**2)
#Linear Reduced Chi Square
def lin_reduced_chi_sqr(x, y, parameter_array, y_errs):
    #Defining a function for calculating the reduced chi-squared value.
    #𝜈 = 𝑁−𝑛_𝑐
    #𝜈 is known as the number of degrees of freedom
    #N = data points
    #𝑛_𝑐  is the number of constraints derived from the data (free parameters in the fit)
    #Note that np.size(parameter_array) == n_c
    #Reduced chi-squared value = 𝜒2/𝜈
    import numpy as np
    d_o_f = np.size(y)-np.size(parameter_array)
    return lin_chi_sqr(x, y, parameter_array, y_errs)/d_o_f
#Quadratic Chi Square
def quad_chi_sqr(x, y, parameter_array, y_errs): 
    #imports
    from scipy.stats import chi2  
    import numpy as np
    from scipy.optimize import curve_fit

    #p0 is an array of starting values for the best-fit parameter which are close to the optimum values. 
    #For linear fit, this means [slope, y-intercept].
    #There are 3 parameters for a quadratic in the form a+bx+cx^2.
    def quad_func(x, a2, slope, y_int):
        return y_int+slope*x+a2*x**2
    popt, pcov = curve_fit(quad_func, x, y, p0=parameter_array, sigma=y_errs, absolute_sigma=True)
    a2 = popt[0]
    slope = popt[1]
    y_int = popt[2]
    y_fit = quad_func(x, a2, slope, y_int)
    return np.sum(((y-y_fit)/y_errs)**2)
#Quadratic Reduced Chi Square
def quad_reduced_chi_sqr(x, y, parameter_array, y_errs):
    #Defining a function for calculating the reduced chi-squared value.
    #𝜈 = 𝑁−𝑛_𝑐
    #𝜈 is known as the number of degrees of freedom
    #N = data points
    #𝑛_𝑐  is the number of constraints derived from the data (free parameters in the fit)
    #Note that np.size(parameter_array) == n_c
    #Reduced chi-squared value = 𝜒2/𝜈
    import numpy as np
    d_o_f = np.size(y)-np.size(parameter_array)
    return quad_chi_sqr(x, y, parameter_array, y_errs)/d_o_f
#Linear curve_fit
def linear_best_fit(xarray,yarray,yerrs,init_estimate_m,init_estimate_c,title="Linear Fit",xlabel="xlabel",ylabel="ylabel",precision=3):
    '''Uses scipy.optimize.curvefit() to fit a straight line to the data.
    Must use np arrays. 'precision' is the number of decimal places. ''' #Docstring
    
    #imports
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from scipy.stats import chi2

    #defining linear function
    def func(x, m, c):
        return m*x+c
    
    initial_estimate = np.array([init_estimate_m, init_estimate_c])
    
    x = xarray
    y = yarray
    
    #curve_fit() to find our line of best fit
    popt, pcov = curve_fit(func, x, y, p0=initial_estimate, sigma=yerrs, absolute_sigma=True)
    #print(*popt) #[m, c]
    
    #Best fit y-values
    y_fit = func(x, *popt)
    
    #Plotting the errors on the best fit plots
    def gr(x,m,c): #mx+c
        dfdm=x
        dfdc=1
        return np.array([[dfdm],[dfdc]])
    
    error=[np.sqrt(float(gr(i,*popt).T @ pcov @ gr(i,*popt))) for i in x]
    #print(error)
    upper=y_fit+error
    lower=y_fit-error
    
    #Plotting 
    plt.figure(figsize=(6,4), dpi=1000)
    
    #Curve of best fit plot
    plt.plot(x, y_fit,label="Line of best fit")
    plt.fill_between(x,upper,lower,facecolor='blue',alpha=0.25,zorder=11)
    plt.errorbar(x,y,yerr=yerrs,fmt='.',label="Data (with error bars)",c="r")
    plt.xlabel(f"xlabel")
    plt.ylabel(f"ylabel")
    plt.legend()
    plt.title(f"title")
    plt.show()
    
    print(f"The best-fit slope is {popt[0]:.{precision}f}, with an error of {np.sqrt(pcov[0,0]):.{precision}f}.\nThe y-intercept is {popt[1]:.{precision}f}, with an error of {np.sqrt(pcov[1,1]):.{precision}f}.")
    
    #Chi squared
    chi2_value = lin_chi_sqr(x, y, initial_estimate,yerrs)
    #Degrees of freedom
    d_o_f = np.size(y)-np.size(initial_estimate)
    #Reduced chi squared
    reduced_chi2_value = lin_reduced_chi_sqr(x, y, initial_estimate,yerrs)
    #Chi squared P-value
    P=chi2.sf(chi2_value, d_o_f)
 
    print(f"The chi squared value is {chi2_value:.5f}.\nThe reduced chi squared value is {reduced_chi2_value:.5f}.\nThe chi squared P-value is {P:.5f}.")
#Quadratic curve_fit
def quad_best_fit(xarray,yarray,yerrs,init_estimate_a2,init_estimate_m,init_estimate_c,title="Quadratic Fit",xlabel="xlabel",ylabel="ylabel",precision=3):
    '''Uses scipy.optimize.curvefit() to fit a quadratic curve to the data.
    Must use np arrays. 'precision' is the number of decimal places. ''' #Docstring
    
    #imports
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from scipy.stats import chi2

    #defining quadratic function
    def quad_func(x, a2, slope, y_int):
        return y_int+slope*x+a2*x**2
    
    initial_estimate = np.array([init_estimate_a2, init_estimate_m, init_estimate_c])
    
    x = xarray
    y = yarray
    
    #curve_fit() to find our line of best fit
    popt, pcov = curve_fit(quad_func, x, y, p0=initial_estimate, sigma=yerrs, absolute_sigma=True)
    #print(*popt) #a2, m, c
    
    #Best fit y-values
    y_fit = quad_func(x, *popt)
    
    #Plotting the errors on the best fit plots
    def gr(x,a2,m,c): #mx+c
        dfda=x**2
        dfdm=x
        dfdc=1
        return np.array([[dfda],[dfdm],[dfdc]])
    
    error=[np.sqrt(float(gr(i,*popt).T @ pcov @ gr(i,*popt))) for i in x]
    #print(error)
    upper=y_fit+error
    lower=y_fit-error
    
    #Plotting 
    plt.figure(figsize=(6,4), dpi=1000)
    #Plotted best fit
    plotrange = np.linspace(min(x),max(x),1000)
    plotyfit = quad_func(plotrange, *popt)
    
    #Curve of best fit plot
    plt.plot(plotrange, plotyfit,label="Quadratic curve of best fit")
    plt.fill_between(x,upper,lower,facecolor='blue',alpha=0.25,zorder=11)
    plt.errorbar(x,y,yerr=yerrs,fmt='.',label="Data (with error bars)",c="r")
    plt.xlabel(f"xlabel")
    plt.ylabel(f"ylabel")
    plt.legend()
    plt.title(f"title")
    plt.show()
    
    print(f"The best-fit curve is {popt[2]:.{precision}f}±{np.sqrt(pcov[2,2]):.{precision}f}+{popt[1]:.{precision}f}±{np.sqrt(pcov[1,1]):.{precision}f}*x+{popt[0]:.{precision}f}±{np.sqrt(pcov[0,0]):.{precision}f}*x^2")
    
    #Chi squared
    chi2_value = quad_chi_sqr(x, y, initial_estimate,yerrs)
    #Degrees of freedom
    d_o_f = np.size(y)-np.size(initial_estimate)
    #Reduced chi squared
    reduced_chi2_value = quad_reduced_chi_sqr(x, y, initial_estimate,yerrs)
    #Chi squared P-value
    P=chi2.sf(chi2_value, d_o_f)
 
    print(f"The chi squared value is {chi2_value:.5f}.\nThe reduced chi squared value is {reduced_chi2_value:.5f}.\nThe chi squared P-value is {P:.5f}.")

#WODR section (Weighted Orthogonal Distance Regression)
#wodr is the general call, and then executes either linear (default), quadratic or 
#exponential depending on given argument
def linear(xarray,yarray,xerrs,yerrs,init_estimate,title,xlabel,ylabel,precision):
    #import
    from scipy import odr
    import numpy as np
    import matplotlib.pyplot as plt
    
    if np.size(init_estimate)!=2:
        raise ValueError("Invalid initial guess array size. Expected np.size=2. m*x+c -> [m,c]")
    else:
        def func(init_estimate, x):
            m,c=init_estimate
            return m*x+c
        
    # Model object
    model = odr.Model(func)
    
    # test data and error
    x = xarray
    y = yarray
    noise_x = xerrs
    noise_y = yerrs
    
    # Create a RealData object
    data = odr.RealData(x, y, sx=noise_x, sy=noise_y)
    
    initial_estimate = init_estimate
    
    # Set up ODR with the model and data.
    odr1 = odr.ODR(data, model, beta0=initial_estimate)
    
    # Run the regression.
    out = odr1.run()
    
    #print fit parameters and 1-sigma estimates
    #(see scipy.odr.Output)
    popt = out.beta     #Estimated parameter values
    pcov = out.cov_beta #Covariance matrix of the estimated parameters
    
    #Plotting the lines of best fit with their errors (using the Matrix Approach), now that we have the popt array for each of them:

    #Best fit y-values for each line
    y_fit = func(popt,x)
    
    #Plotting the errors on the best fit plots
    def gr(x,m,c): #m*x+c
        dfdm=x
        dfdc=1
        return np.array([[dfdm],[dfdc]])
    
    error=[np.sqrt(float(gr(i,*popt).T @ pcov @ gr(i,*popt))) for i in x]
    upper=y_fit+error
    lower=y_fit-error
    
    #Plotting 
    plt.figure(figsize=(6,4), dpi=1000)
    
    #Curve of best fit plot
    plt.plot(x, y_fit,label="Line of best fit")
    plt.fill_between(x,upper,lower,facecolor='blue',alpha=0.25,zorder=11,label="Confidence interval")
    plt.errorbar(x,y,xerr=noise_x, yerr=noise_y, barsabove=True, fmt="r.", label="Error bars")
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    plt.legend()
    plt.title(f"{title}")
    plt.show()

    error_m,error_c=np.sqrt(np.diag(pcov))
    print(f"The line of best fit is given by ({popt[0]:.{precision}f}±{error_m:.{precision}f})*x + {popt[1]:.{precision}f}±{error_c:.{precision}f}")

def quadratic(xarray,yarray,xerrs,yerrs,init_estimate,title,xlabel,ylabel,precision):
    #import
    from scipy import odr
    import numpy as np
    import matplotlib.pyplot as plt
    
    if np.size(init_estimate)!=3:
        raise ValueError("Invalid initial guess array size. Expected np.size=3. ax^2+mx+c -> [a,m,c].")
    else:
        def quad_func(init_estimate, x): #ax^2+mx+c
            a,m,c=init_estimate
            return a*x**2+m*x+c
    # Model object
    model = odr.Model(quad_func)
    
    # test data and error
    x = xarray
    y = yarray
    noise_x = xerrs
    noise_y = yerrs
    
    # Create a RealData object
    data = odr.RealData(x, y, sx=noise_x, sy=noise_y)
    
    initial_estimate = init_estimate
    
    # Set up ODR with the model and data.
    odr1 = odr.ODR(data, model, beta0=initial_estimate)
    
    # Run the regression.
    out = odr1.run()
    
    #print fit parameters and 1-sigma estimates
    #(see scipy.odr.Output)
    popt = out.beta     #Estimated parameter values
    pcov = out.cov_beta #Covariance matrix of the estimated parameters
    
    #Plotting the lines of best fit with their errors (using the Matrix Approach), now that we have the popt array for each of them:

    #Best fit y-values for each line
    y_fit = quad_func(popt,x)
    
    #Plotting the errors on the best fit plots
    def gr(x,a,m,c): #ax^2+mx+c
        dfda=x**2
        dfdm=x
        dfdc=1
        return np.array([[dfda],[dfdm],[dfdc]])
    
    error=[np.sqrt(float(gr(i,*popt).T @ pcov @ gr(i,*popt))) for i in x]
    upper=y_fit+error
    lower=y_fit-error
    
    #Plotting 
    plt.figure(figsize=(6,4), dpi=1000)
    #Plotted best fit
    plotrange = np.linspace(min(x),max(x),1000)
    plotyfit = quad_func(popt,plotrange)
    
    #Curve of best fit plot
    plt.plot(plotrange, plotyfit,label="Quadratic curve of best fit")
    plt.fill_between(x,upper,lower,facecolor='blue',alpha=0.25,zorder=11,label="Confidence interval")
    plt.errorbar(x,y,xerr=noise_x, yerr=noise_y, barsabove=True, fmt="r.", label="Error bars")
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    plt.legend()
    plt.title(f"{title}")
    plt.show()

    error_a,error_m,error_c=np.sqrt(np.diag(pcov))
    print(f"The curve of best fit is given by ({popt[0]:.{precision}f}±{error_a:.{precision}f})*x^2 + ({popt[1]:.{precision}f}±{error_m:.{precision}f})*x + {popt[2]:.{precision}f}±{error_c:.{precision}f}")

def exponential(xarray,yarray,xerrs,yerrs,init_estimate,title,xlabel,ylabel,precision):
    #import
    from scipy import odr
    import numpy as np
    import matplotlib.pyplot as plt
    
    if np.size(init_estimate)!=3:
        raise ValueError("Invalid initial guess array size. Expected np.size=3. a*np.exp(b*x)+c -> [a,b,c]")
    else:
        def exp_func(init_estimate,x): #a*np.exp(b*x)+c
            a,b,c=init_estimate
            return a*np.exp(b*x)+c
    # Model object
    model = odr.Model(exp_func)
    
    # test data and error
    x = xarray
    y = yarray
    noise_x = xerrs
    noise_y = yerrs
    
    # Create a RealData object
    data = odr.RealData(x, y, sx=noise_x, sy=noise_y)
    
    initial_estimate = init_estimate
    
    # Set up ODR with the model and data.
    odr1 = odr.ODR(data, model, beta0=initial_estimate)
    
    # Run the regression.
    out = odr1.run()
    
    #print fit parameters and 1-sigma estimates
    #(see scipy.odr.Output)
    popt = out.beta     #Estimated parameter values
    pcov = out.cov_beta #Covariance matrix of the estimated parameters
    
    #Plotting the lines of best fit with their errors (using the Matrix Approach), now that we have the popt array for each of them:

    #Best fit y-values for each line
    y_fit = exp_func(popt,x)
    
    #Plotting the errors on the best fit plots
    def gr(x,a,b,c): #ax^2+mx+c
        dfda=np.exp(b*x)
        dfdb=a*np.exp(b*x)*x
        dfdc=1
        return np.array([[dfda],[dfdb],[dfdc]])
    
    error=[np.sqrt(float(gr(i,*popt).T @ pcov @ gr(i,*popt))) for i in x]
    upper=y_fit+error
    lower=y_fit-error
    
    #Plotting 
    plt.figure(figsize=(6,4), dpi=1000)
    #Plotted best fit
    plotrange = np.linspace(min(x),max(x),1000)
    plotyfit = exp_func(popt,plotrange)
    
    #Curve of best fit plot
    plt.plot(plotrange, plotyfit,label="Exponential curve of best fit")
    plt.fill_between(x,upper,lower,facecolor='blue',alpha=0.25,zorder=11,label="Confidence interval")
    plt.errorbar(x,y,xerr=noise_x, yerr=noise_y, barsabove=True, fmt="r.", label="Error bars")
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    plt.legend()
    plt.title(f"{title}")
    plt.show()

    error_a,error_b,error_c=np.sqrt(np.diag(pcov))
    print(f"The exponential curve of best fit is given by ({popt[0]:.{precision}f}±{error_a:.{precision}f}) * e^[({popt[1]:.{precision}f}±{error_b:.{precision}f})*x] + {popt[2]:.{precision}f}±{error_c:.{precision}f}")

def wodr(xarray,yarray,xerrs,yerrs,init_estimate,function='linear',title="Title",xlabel="xlabel",ylabel="ylabel",precision=3):
    '''Using the Matrix Approach with the Weighted Orthogonal Distance Regression method of finding lines of best fit.'''

    function_types = ['linear', 'quadratic', 'exponential']
    if function in function_types:
        eval(function+"(xarray,yarray,xerrs,yerrs,init_estimate,title,xlabel,ylabel,precision)")
    else:
        raise ValueError("Invalid function type. Expected one of: %s" % function_types)