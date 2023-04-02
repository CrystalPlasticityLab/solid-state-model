# Tensor
Base class `container` may contain of arrays with any dimension and rank (even rank = 0 for scalar arrays).
Based on `container` class `object` provides basic tensor-vector and scalar arrays calculus in N dimensional spatial.

Class `object` provide +/-/scal product operators between each other. So you should not keep in mind at how basis are component of current tensor/vector, 
internal functionality does operation correctly.

Using smart pointers for manage basis and component containers allows to have many objects (tensor/vector) at the same basis object. 
So all objects linked to the one basis will change at the same time when basis will be changed.
Basis object will live till the last tensor/vector object will be destroyed.

There is possibility to handle array of scalar or even just scalar variable like tensor object wut without basis poiter.
It is an advantage for abstraction and construction of any measures (see the next paragraph).

# Measure
The main purposes the next features are imtroduced for fast implementation of any mathematical model based on [state variable approach](https://en.wikipedia.org/wiki/State_variable). The main points:
- state of a system may be described by finite set of tensor or/and scalar state variables (**SV**)
- each state variable may have evolutionary equation (usually - differential equation)
- each state variable may have dependecies on other strate variables
- 
Class `StateMeasure` based on `object` class is providing unified logic for working with state variable:
- it has `rate` - value of rate of **SV**
- it has `value` - value **SV**
- also it contains previous value and rate of **SV** for ability to implement of more convinience numerical schemas 
- `rate_equation()` - evolution equation in rate form (the rate of state variable dependency on set of parameters and set of **SV**'s)
- `finit_equation()` - evolution equation in finite form (describes dependency of finite value on set of parameters and set of **SV**'s)
- `calc_rate()` - by default it is the first order schema to calculate rate, may be overriden
- `integrate_value(T dt)` - by default the first order (Euler) schema to integrate value, may be overriden
- 
# Numerical schema
Base abstract class `AbstractSchema` describes main functionality:
- `init()` - assumed that it is called before calc rate/finite and integration procedure
- `step()` - the main stage of numerical schema, `rate` or `value` on the end of current step are calculated
- `finalize()` - actions after if necessery

There are a few types of numerical schema:
- `RATE_ASSIGN`: *dX(n+1) := dx_new, X(n+1) = X(n) + dX(n+1)dt*
- `RATE_CALCULATE`: *dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)dt*
- `FINITE_ASSIGN`: *X(n+1) := x_new, dX(n+1) = (X(n+1)-X(n))/dt*
- `FINITE_CALCULATE`: *X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt*
- 
Class `DefaultSchema` implements (abstract methods of `AbstractSchema`) a plain first order numerical schema over `StateMeasure`.

# State
Class `state` has array of `measure` that describe a state of an object.
State has only one `basis` object so any tensor measures are linked to the same basis, at the same time indifferent scalar measures are not linked to any basis.
